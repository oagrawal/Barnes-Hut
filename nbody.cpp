#include <mpi.h>
#include <iostream>
#include <iomanip> // For std::setprecision
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <getopt.h>
#include <windows.h>

// Constants for the simulation
const double G = 0.0001;    // Gravitational constant
const double RLIMIT = 0.03; // Minimum distance to avoid infinite forces
const double DOMAIN_MIN = 0.0;
const double DOMAIN_MAX = 4.0;
const int WINDOW_SIZE = 800; // Window size in pixels

// Structure definitions
struct Body
{
    int index;     // -1 for internal, correct for external nodes
    double px, py; // Position
    double mass;   // -1 if out of bounds
    double vx, vy; // Velocity
    double fx, fy; // Force
};

struct Node
{
    struct Body *b; // Pointer to body (can represent internal / external nodes)
    struct Node *NW, *NE, *SW, *SE;
    double minx, maxx, miny, maxy;
};

// Global variables for visualization
HWND hwnd = NULL;
HDC hdc = NULL;
HBRUSH hBrush = NULL;
HPEN hPen = NULL;

// Window procedure
LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    switch (uMsg)
    {
    case WM_CLOSE:
        DestroyWindow(hwnd);
        break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hwnd, uMsg, wParam, lParam);
    }
    return 0;
}

// Initialize visualization window
bool initVisualization()
{
    WNDCLASSEXA wc = {0};
    wc.cbSize = sizeof(WNDCLASSEXA);
    wc.lpfnWndProc = WindowProc;
    wc.hInstance = GetModuleHandle(NULL);
    wc.lpszClassName = "NBodyVisualization";

    if (!RegisterClassExA(&wc))
    {
        return false;
    }

    // Create the window
    hwnd = CreateWindowExA(
        0,
        "NBodyVisualization",
        "N-Body Simulation",
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT,
        WINDOW_SIZE, WINDOW_SIZE,
        NULL,
        NULL,
        GetModuleHandle(NULL),
        NULL);

    if (!hwnd)
    {
        return false;
    }

    // Show the window
    ShowWindow(hwnd, SW_SHOW);
    UpdateWindow(hwnd);

    // Get the device context
    hdc = GetDC(hwnd);

    // Create drawing objects
    hBrush = CreateSolidBrush(RGB(0, 0, 0));           // Black brush for bodies
    hPen = CreatePen(PS_SOLID, 1, RGB(128, 128, 128)); // Gray pen for grid

    return true;
}

// Clean up visualization resources
void cleanupVisualization()
{
    if (hPen)
        DeleteObject(hPen);
    if (hBrush)
        DeleteObject(hBrush);
    if (hdc)
        ReleaseDC(hwnd, hdc);
    if (hwnd)
        DestroyWindow(hwnd);
}

void visualizeBodies(const std::vector<Body> &bodies, int rank)
{
    // Only rank 0 should visualize
    if (rank != 0)
        return;

    // Process Windows messages
    MSG msg;
    while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

    // Clear the window
    RECT rect;
    GetClientRect(hwnd, &rect);
    FillRect(hdc, &rect, (HBRUSH)GetStockObject(WHITE_BRUSH));

    // Draw grid lines
    SelectObject(hdc, hPen);
    for (int i = 0; i <= 10; i++)
    {
        int pos = i * (WINDOW_SIZE / 10);
        // Vertical lines
        MoveToEx(hdc, pos, 0, NULL);
        LineTo(hdc, pos, WINDOW_SIZE);
        // Horizontal lines
        MoveToEx(hdc, 0, pos, NULL);
        LineTo(hdc, WINDOW_SIZE, pos);
    }

    // Find max mass for scaling
    double maxMass = 0.0;
    for (const Body &b : bodies)
    {
        if (b.mass > 0 && b.mass != -1)
        {
            maxMass = std::max(maxMass, b.mass);
        }
    }

    // Draw bodies
    SelectObject(hdc, hBrush);
    for (const Body &b : bodies)
    {
        if (b.mass == -1)
            continue; // Skip lost particles

        // Map position from [0,4] to [0,WINDOW_SIZE]
        int x = static_cast<int>((b.px / DOMAIN_MAX) * WINDOW_SIZE);
        int y = static_cast<int>((1.0 - b.py / DOMAIN_MAX) * WINDOW_SIZE);

        // Calculate radius based on mass (min 3, max 15 pixels)
        const int minRadius = 3;
        const int maxRadius = 15;
        int radius = minRadius;
        if (maxMass > 0)
        {
            radius = minRadius + static_cast<int>((b.mass / maxMass) * (maxRadius - minRadius));
        }

        // Draw body as a circle
        Ellipse(hdc,
                x - radius, y - radius,
                x + radius, y + radius);
    }
}

std::vector<Body> readBodiesFromFile(const std::string &filename);
Node *constructTreeHelper(Node *root, const Body &b, double minx, double maxx, double miny, double maxy);
Node *constructTree(std::vector<Body> &bodies);

bool parseArguments(int argc, char *argv[], std::string &inputFile, std::string &outputFile,
                    int &steps, double &theta, double &dt, bool &visualization)
{
    // Set default values
    inputFile = "input/nb-10.txt";
    outputFile = "output.txt";
    steps = 100;
    theta = 0.5; // Common default for Barnes-Hut
    dt = 0.005;
    visualization = false;
    int opt;
    while ((opt = getopt(argc, argv, "i:o:s:t:d:V")) != -1)
    {
        switch (opt)
        {
        case 'i':
            inputFile = optarg;
            break;
        case 'o':
            outputFile = optarg;
            break;
        case 's':
            steps = std::stoi(optarg);
            break;
        case 't':
            theta = std::stod(optarg);
            break;
        case 'd':
            dt = std::stod(optarg);
            break;
        case 'V':
            visualization = true;
            break;
        default:
            return false;
        }
    }

    // Validate values
    if (steps <= 0 || dt <= 0.0 || theta < 0.0)
    {
        std::cerr << "Error: Invalid parameter values. Steps, dt, and theta must be positive." << std::endl;
        return false;
    }

    return true;
}

Node *insertCorrect(Node *root, const Body &b)
{
    if (!root)
        return nullptr;

    double midx = (root->maxx + root->minx) / 2;
    double midy = (root->maxy + root->miny) / 2;

    if (b.px <= midx && b.py > midy)
    {
        root->NW = constructTreeHelper(root->NW, b, root->minx, midx, midy, root->maxy);
    }
    else if (b.px > midx && b.py > midy)
    {
        root->NE = constructTreeHelper(root->NE, b, midx, root->maxx, midy, root->maxy);
    }
    else if (b.px <= midx && b.py <= midy)
    {
        root->SW = constructTreeHelper(root->SW, b, root->minx, midx, root->miny, midy);
    }
    else if (b.px > midx && b.py <= midy)
    {
        root->SE = constructTreeHelper(root->SE, b, midx, root->maxx, root->miny, midy);
    }
    return root;
}

Node *constructTreeHelper(Node *root, const Body &b, double minx, double maxx, double miny, double maxy)
{
    // Skip lost particles
    if (b.mass == -1)
    {
        return root;
    }

    if (root == nullptr)
    {
        // Allocate new node with 'new'
        root = new Node();
        root->NW = nullptr;
        root->NE = nullptr;
        root->SW = nullptr;
        root->SE = nullptr;
        root->minx = minx;
        root->maxx = maxx;
        root->miny = miny;
        root->maxy = maxy;
        root->b = new Body(b);
        return root;
    }

    // If this is an empty node, just put the body here
    if (root->b->mass == 0)
    {
        *root->b = b;
    }
    // If this is an external node (leaf with a single body)
    else if (root->NW == nullptr && root->NE == nullptr &&
             root->SW == nullptr && root->SE == nullptr)
    {
        // Store existing body
        Body b2 = *root->b;

        // Mark as internal node
        root->b->index = -1;

        // Insert both bodies
        root = insertCorrect(root, b);
        root = insertCorrect(root, b2);

        // Update center of mass and total mass
        double totalMass = root->b->mass + b.mass;
        if (totalMass > 0)
        { // Prevent division by zero
            // Use a more stable formula for weighted average to avoid overflow
            double w1 = root->b->mass / totalMass;
            double w2 = b.mass / totalMass;
            root->b->px = w1 * root->b->px + w2 * b.px;
            root->b->py = w1 * root->b->py + w2 * b.py;

            // root->b->px = ((root->b->px*root->b->mass) + (b.px*b.mass)) / (b.mass + root->b->mass);
            // root->b->py = ((root->b->py*root->b->mass) + (b.py*b.mass)) / (b.mass + root->b->mass);

            root->b->mass = totalMass;
        }
        else
        {
            // If both masses are zero, just use the position of the new body
            root->b->px = b.px;
            root->b->py = b.py;
            root->b->mass = 0;
        }
    }
    // If this is an internal node
    else
    {
        // Update center of mass and total mass
        double totalMass = root->b->mass + b.mass;
        if (totalMass > 0)
        { // Prevent division by zero
            // Use a more stable formula for weighted average to avoid overflow
            double w1 = root->b->mass / totalMass;
            double w2 = b.mass / totalMass;
            root->b->px = w1 * root->b->px + w2 * b.px;
            root->b->py = w1 * root->b->py + w2 * b.py;
            // root->b->px = ((root->b->px*root->b->mass) + (b.px*b.mass)) / (b.mass + root->b->mass);
            // root->b->py = ((root->b->py*root->b->mass) + (b.py*b.mass)) / (b.mass + root->b->mass);

            root->b->mass = totalMass;
        }
        else
        {
            // If both masses are zero, just use the position of the new body
            root->b->px = b.px;
            root->b->py = b.py;
            root->b->mass = 0;
        }

        // Insert the new body
        root = insertCorrect(root, b);
    }
    return root;
}

Node *constructTree(std::vector<Body> &bodies)
{
    Node *root = nullptr;

    for (size_t i = 0; i < bodies.size(); i++)
    {
        root = constructTreeHelper(root, bodies[i], DOMAIN_MIN, DOMAIN_MAX, DOMAIN_MIN, DOMAIN_MAX);
    }
    return root;
}

void destroyTree(Node *root)
{
    if (root == nullptr)
        return;

    // Recursively destroy all children
    destroyTree(root->NW);
    destroyTree(root->NE);
    destroyTree(root->SW);
    destroyTree(root->SE);

    // Free the body pointer
    delete root->b;

    // Free this node
    delete root;
}

void printTree(Node *root, int level = 0, const std::string &prefix = "")
{
    if (root == nullptr)
    {
        std::cout << prefix << "└── [empty]" << std::endl;
        return;
    }

    // Print current node
    std::cout << prefix;
    if (level == 0)
    {
        std::cout << "Root: ";
    }
    else
    {
        std::cout << "└── ";
    }

    if (root->b->index == -1)
    {
        std::cout << "Internal Node (mass=" << root->b->mass
                  << ", pos=(" << root->b->px << "," << root->b->py << ")"
                  << ", bounds=[" << root->minx << "," << root->maxx << "]["
                  << root->miny << "," << root->maxy << "])" << std::endl;
    }
    else
    {
        std::cout << "Body " << root->b->index
                  << " (mass=" << root->b->mass
                  << ", pos=(" << root->b->px << "," << root->b->py << ")"
                  << ", vel=(" << root->b->vx << "," << root->b->vy << "))" << std::endl;
    }

    // Print children with increased indentation
    std::string newPrefix = prefix + "    ";
    printTree(root->NW, level + 1, newPrefix + "NW: ");
    printTree(root->NE, level + 1, newPrefix + "NE: ");
    printTree(root->SW, level + 1, newPrefix + "SW: ");
    printTree(root->SE, level + 1, newPrefix + "SE: ");
}

void calcForce(Body *b, Node *root, double theta)
{
    if (root == nullptr)
    {
        return;
    }

    // Skip calculation if this is the same body or if either body is lost
    if ((root->b->index == b->index && root->b->index != -1) ||
        b->mass == -1 || root->b->mass == -1)
    {
        return;
    }

    // Calculate distance vector (preserving direction)
    double dx = root->b->px - b->px;
    double dy = root->b->py - b->py;
    double d = sqrt(dx * dx + dy * dy);

    // Check if node is external (leaf node)
    if (root->b->index != -1)
    {
        // External node (and not the body itself)
        if (d < RLIMIT)
        {
            d = RLIMIT; // Prevent division by zero/very small numbers
        }

        // Calculate gravitational force using the new formula
        // F = G*M0*M1*d / d^3
        // Projected forces: Fx = G*M0*M1*dx / d^3, Fy = G*M0*M1*dy / d^3
        double d3 = d * d * d;
        b->fx += (G * b->mass * root->b->mass * dx) / d3;
        b->fy += (G * b->mass * root->b->mass * dy) / d3;

        // double d2 = d * d;
        // double d3 = d2 * d;
        // // Optional: Consider factoring out common terms
        // double gm = G * b->mass * root->b->mass;
        // b->fx += gm * dx / d3;
        // b->fy += gm * dy / d3;

        return;
    }
    else
    {
        // Internal node - check if it's far enough away
        double s = root->maxx - root->minx; // Size of region

        if ((s / d) < theta)
        {
            // Far enough away, treat as single body
            if (d < RLIMIT)
            {
                d = RLIMIT;
            }

            // Calculate gravitational force using the new formula
            // F = G*M0*M1*d / d^3
            // Projected forces: Fx = G*M0*M1*dx / d^3, Fy = G*M0*M1*dy / d^3
            double d3 = d * d * d;
            b->fx += (G * b->mass * root->b->mass * dx) / d3;
            b->fy += (G * b->mass * root->b->mass * dy) / d3;

            // double d2 = d * d;
            // double d3 = d2 * d;
            // // Optional: Consider factoring out common terms
            // double gm = G * b->mass * root->b->mass;
            // b->fx += gm * dx / d3;
            // b->fy += gm * dy / d3;

            return;
        }
        else
        {
            // Too close, recursively check children
            calcForce(b, root->NW, theta);
            calcForce(b, root->NE, theta);
            calcForce(b, root->SW, theta);
            calcForce(b, root->SE, theta);
        }
    }
}

void calculateForcesParallel(Node *root, std::vector<Body> &bodies, double theta, int rank, int size)
{
    // Reset forces for all bodies
    for (auto &body : bodies)
    {
        body.fx = 0.0;
        body.fy = 0.0;
    }

    // Each rank calculates forces for its assigned subset of bodies (round-robin)
    for (size_t i = rank; i < bodies.size(); i += size)
    {
        // Skip lost particles
        if (bodies[i].mass == -1)
        {
            continue;
        }
        calcForce(&bodies[i], root, theta);
    }

    // Create arrays for force aggregation
    std::vector<double> local_fx(bodies.size(), 0.0);
    std::vector<double> local_fy(bodies.size(), 0.0);
    std::vector<double> global_fx(bodies.size(), 0.0);
    std::vector<double> global_fy(bodies.size(), 0.0);

    // Copy local results to the arrays
    for (size_t i = rank; i < bodies.size(); i += size)
    {
        local_fx[i] = bodies[i].fx;
        local_fy[i] = bodies[i].fy;
    }

    // Combine results from all processes
    MPI_Allreduce(local_fx.data(), global_fx.data(), bodies.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local_fy.data(), global_fy.data(), bodies.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // Update bodies with the combined forces
    for (size_t i = 0; i < bodies.size(); i++)
    {
        bodies[i].fx = global_fx[i];
        bodies[i].fy = global_fy[i];
    }
}

void updateBodies(std::vector<Body> &bodies, double dt)
{
    double ax, ay;
    for (auto &body : bodies)
    {
        if (body.mass == 0 || body.mass == -1)
        {
            // Skip bodies with zero mass or lost particles
            continue;
        }
        ax = body.fx / body.mass;
        ay = body.fy / body.mass;
        body.vx = body.vx + ax * dt;
        body.vy = body.vy + ay * dt;
        body.px = body.px + body.vx * dt + 0.5 * ax * (dt * dt);
        body.py = body.py + body.vy * dt + 0.5 * ay * (dt * dt);

        // Check if particle is outside the domain
        if (body.px < DOMAIN_MIN || body.px > DOMAIN_MAX ||
            body.py < DOMAIN_MIN || body.py > DOMAIN_MAX)
        {
            // Mark as lost particle
            body.mass = -1;
        }
    }
}

void writeBodiestoFile(const std::vector<Body> &bodies, const std::string &outputFile, int rank)
{
    // Only rank 0 should write to the output file to avoid conflicts
    if (rank != 0)
    {
        return;
    }

    // Open the output file
    std::ofstream outFile(outputFile);
    if (!outFile.is_open())
    {
        std::cerr << "Error: Could not open output file " << outputFile << std::endl;
        return;
    }

    // Count non-lost bodies
    int validBodies = 0;
    for (const auto &body : bodies)
    {
        // if (body.mass != -1) {
        //     validBodies++;
        // }
        validBodies++;
    }

    // Write the number of bodies
    outFile << validBodies << std::endl;

    // Write each body's data
    for (const auto &body : bodies)
    {
        // Skip lost particles
        // if (body.mass == -1) {
        //     continue;
        // }

        // Write data with scientific notation and specific precision
        outFile << body.index << "  "
                << std::scientific << std::setprecision(6) << body.px << "  "
                << std::scientific << std::setprecision(6) << body.py << "  "
                << std::scientific << std::setprecision(6) << body.mass << "  "
                << std::scientific << std::setprecision(6) << body.vx << "  "
                << std::scientific << std::setprecision(6) << body.vy << std::endl;

        // outFile << body.index << "  "
        //         << body.px << "  "
        //         << body.py << "  "
        //         << body.mass << "  "
        //         << body.vx << "  "
        //         << body.vy << std::endl;
    }

    outFile.close();
    std::cout << "Wrote " << validBodies << " bodies to " << outputFile << std::endl;
}

int main(int argc, char *argv[])
{
    // Initialize MPI environment
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Parse command line arguments
    std::string inputFile, outputFile;
    int steps;
    double theta, dt;
    bool visualization;

    if (!parseArguments(argc, argv, inputFile, outputFile, steps, theta, dt, visualization))
    {
        if (rank == 0)
        {
            std::cerr << "Usage: " << argv[0] << " -i <input_file> -o <output_file> "
                      << "-s <steps> -t <theta> -d <dt> [-V]" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    // Print parameters (only from rank 0)
    if (rank == 0)
    {
        std::cout << "Parameters:" << std::endl;
        std::cout << "  Input file: " << inputFile << std::endl;
        std::cout << "  Output file: " << outputFile << std::endl;
        std::cout << "  Steps: " << steps << std::endl;
        std::cout << "  Theta: " << theta << std::endl;
        std::cout << "  Timestep: " << dt << std::endl;
        std::cout << "  Visualization: " << (visualization ? "On" : "Off") << std::endl;
    }

    // Initialize visualization if requested
    if (visualization && rank == 0)
    {
        if (!initVisualization())
        {
            std::cerr << "Failed to initialize visualization" << std::endl;
            MPI_Finalize();
            return 1;
        }
    }

    // Read input file
    std::vector<Body> bodies;
    bodies = readBodiesFromFile(inputFile);
    if (rank == 0)
    {
        std::cout << "Read " << bodies.size() << " bodies from " << inputFile << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Main simulation loop
    for (int step = 0; step < steps; step++)
    {
        Node *root = constructTree(bodies);
        calculateForcesParallel(root, bodies, theta, rank, size);
        updateBodies(bodies, dt);

        if (visualization)
        {
            visualizeBodies(bodies, rank);
            Sleep(5); // Reduced delay for faster animation
        }

        destroyTree(root);
    }

    // Write final state to output file
    writeBodiestoFile(bodies, outputFile, rank);

    // Clean up visualization
    if (visualization && rank == 0)
    {
        cleanupVisualization();
    }

    MPI_Finalize();
    return 0;
}

std::vector<Body> readBodiesFromFile(const std::string &filename)
{
    std::vector<Body> bodies;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return bodies;
    }

    // Read number of bodies
    int numBodies;
    file >> numBodies;

    // Reserve space to avoid reallocations
    bodies.reserve(numBodies);

    // Read each body
    for (int i = 0; i < numBodies; i++)
    {
        Body body;
        file >> body.index >> body.px >> body.py >> body.mass >> body.vx >> body.vy;

        // Initialize forces to zero
        body.fx = 0.0;
        body.fy = 0.0;

        bodies.push_back(body);
    }

    file.close();
    return bodies;
}
