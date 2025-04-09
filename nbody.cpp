#include <mpi.h>
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <getopt.h>
#include <chrono>

const double G = 0.0001;  
const double RLIMIT = 0.03;   
const double DOMAIN_MIN = 0.0;
const double DOMAIN_MAX = 4.0;

struct Body {
    int index;          // -1 for internal, correct for external nodes
    double px, py;      // Position
    double mass;        // -1 if out of bounds
    double vx, vy;      // Velocity
    double fx, fy;      // Force
};

struct Node {
    struct Body* b;      // Pointer to body (can represent internal / external nodes)
    struct Node *NW, *NE, *SW, *SE;
    double minx, maxx, miny, maxy;
};


std::vector<Body> readBodiesFromFile(const std::string& filename);
Node* constructTreeHelper(Node* root, const Body& b, double minx, double maxx, double miny, double maxy);
Node* constructTree(std::vector<Body>& bodies);

bool parseArguments(int argc, char* argv[], std::string& inputFile, std::string& outputFile, 
                   int& steps, double& theta, double& dt, bool& visualization, bool& sequential) {

    inputFile = "input/nb-10.txt";
    outputFile = "output.txt";
    steps = 100;
    theta = 0.5; 
    dt = 0.005;
    visualization = false;
    sequential = false;
    int opt;
    while ((opt = getopt(argc, argv, "i:o:s:t:d:VS")) != -1) {
        switch (opt) {
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
            case 'S':
                sequential = true;
                break;
            default:
                return false;
        }
    }
    
    if (steps <= 0 || dt <= 0.0 || theta < 0.0) {
        std::cerr << "Error: Invalid parameter values. Steps, dt, and theta must be positive." << std::endl;
        return false;
    }
    
    return true;
}


Node* insertCorrect(Node* root, const Body& b) {
    if(!root) return nullptr;
    
    double midx = (root->maxx + root->minx) / 2;
    double midy = (root->maxy + root->miny) / 2;
    
    if(b.px <= midx && b.py > midy) {
        root->NW = constructTreeHelper(root->NW, b, root->minx, midx, midy, root->maxy);
    } else if(b.px > midx && b.py > midy) {
        root->NE = constructTreeHelper(root->NE, b, midx, root->maxx, midy, root->maxy);
    } else if(b.px <= midx && b.py <= midy) {
        root->SW = constructTreeHelper(root->SW, b, root->minx, midx, root->miny, midy);
    } else if(b.px > midx && b.py <= midy) {
        root->SE = constructTreeHelper(root->SE, b, midx, root->maxx, root->miny, midy);
    }
    return root;
}

Node* constructTreeHelper(Node* root, const Body& b, double minx, double maxx, double miny, double maxy) {

    if (b.mass == -1) {
        return root;
    }
    
    if(root == nullptr) {
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
    
    if(root->b->mass == 0) {
        *root->b = b;
    } 
    else if(root->NW == nullptr && root->NE == nullptr && 
            root->SW == nullptr && root->SE == nullptr) {
        Body b2 = *root->b;
        root->b->index = -1;
        root = insertCorrect(root, b);
        root = insertCorrect(root, b2);
        
        double totalMass = root->b->mass + b.mass;
        if (totalMass > 0) { 
            double w1 = root->b->mass / totalMass;
            double w2 = b.mass / totalMass;
            root->b->px = w1 * root->b->px + w2 * b.px;
            root->b->py = w1 * root->b->py + w2 * b.py;

            root->b->mass = totalMass;
        } else {
            root->b->px = b.px;
            root->b->py = b.py;
            root->b->mass = 0;
        }
    } 
    else {
        double totalMass = root->b->mass + b.mass;
        if (totalMass > 0) {  
            double w1 = root->b->mass / totalMass;
            double w2 = b.mass / totalMass;
            root->b->px = w1 * root->b->px + w2 * b.px;
            root->b->py = w1 * root->b->py + w2 * b.py;

            root->b->mass = totalMass;
        } else {
            root->b->px = b.px;
            root->b->py = b.py;
            root->b->mass = 0;
        }
        
        root = insertCorrect(root, b);
    }
    return root;
}

Node* constructTree(std::vector<Body>& bodies) {
    Node* root = nullptr;
    
    for(size_t i = 0; i < bodies.size(); i++) {
        root = constructTreeHelper(root, bodies[i], DOMAIN_MIN, DOMAIN_MAX, DOMAIN_MIN, DOMAIN_MAX);
    }
    return root;
}

void destroyTree(Node* root) {
    if(root == nullptr) return;
    
    destroyTree(root->NW);
    destroyTree(root->NE);
    destroyTree(root->SW);
    destroyTree(root->SE);
    
    delete root->b;
    
    delete root;
}

void printTree(Node* root, int level = 0, const std::string& prefix = "") {
    if (root == nullptr) {
        std::cout << prefix << "└── [empty]" << std::endl;
        return;
    }
    
    std::cout << prefix;
    if (level == 0) {
        std::cout << "Root: ";
    } else {
        std::cout << "└── ";
    }
    
    if (root->b->index == -1) {
        std::cout << "Internal Node (mass=" << root->b->mass 
                  << ", pos=(" << root->b->px << "," << root->b->py << ")"
                  << ", bounds=[" << root->minx << "," << root->maxx << "][" 
                  << root->miny << "," << root->maxy << "])" << std::endl;
    } else {
        std::cout << "Body " << root->b->index 
                  << " (mass=" << root->b->mass 
                  << ", pos=(" << root->b->px << "," << root->b->py << ")"
                  << ", vel=(" << root->b->vx << "," << root->b->vy << "))" << std::endl;
    }
    
    std::string newPrefix = prefix + "    ";
    printTree(root->NW, level + 1, newPrefix + "NW: ");
    printTree(root->NE, level + 1, newPrefix + "NE: ");
    printTree(root->SW, level + 1, newPrefix + "SW: ");
    printTree(root->SE, level + 1, newPrefix + "SE: ");
}

void calcForce(Body* b, Node* root, double theta) {
    if (root == nullptr) {
        return;
    }
    
    if ((root->b->index == b->index && root->b->index != -1) || 
        b->mass == -1 || root->b->mass == -1) {
        return;
    }
    
    double dx = root->b->px - b->px;
    double dy = root->b->py - b->py;
    double d = sqrt(dx*dx + dy*dy);
    
    if (root->b->index != -1) {
        if (d < RLIMIT) {
            d = RLIMIT;  
        }
        double d3 = d * d * d;
        b->fx += (G * b->mass * root->b->mass * dx) / d3;
        b->fy += (G * b->mass * root->b->mass * dy) / d3;
        
        return;
    } else {
        double s = root->maxx - root->minx;  
        
        if ((s / d) < theta) {
            if (d < RLIMIT) {
                d = RLIMIT;
            }
            double d3 = d * d * d;
            b->fx += (G * b->mass * root->b->mass * dx) / d3;
            b->fy += (G * b->mass * root->b->mass * dy) / d3;

            return;
        } else {
            calcForce(b, root->NW, theta);
            calcForce(b, root->NE, theta);
            calcForce(b, root->SW, theta);
            calcForce(b, root->SE, theta);
        }
    }
}

void calculateForcesParallel(Node* root, std::vector<Body>& bodies, double theta, int rank, int size) {
    for (auto& body : bodies) {
        body.fx = 0.0;
        body.fy = 0.0;
    }
    
    for (size_t i = rank; i < bodies.size(); i += size) {
        if (bodies[i].mass == -1) {
            continue;
        }
        calcForce(&bodies[i], root, theta);
    }
    
    std::vector<double> local_fx(bodies.size(), 0.0);
    std::vector<double> local_fy(bodies.size(), 0.0);
    std::vector<double> global_fx(bodies.size(), 0.0);
    std::vector<double> global_fy(bodies.size(), 0.0);
    
    for (size_t i = rank; i < bodies.size(); i += size) {
        local_fx[i] = bodies[i].fx;
        local_fy[i] = bodies[i].fy;
    }
    
    MPI_Allreduce(local_fx.data(), global_fx.data(), bodies.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local_fy.data(), global_fy.data(), bodies.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    for (size_t i = 0; i < bodies.size(); i++) {
        bodies[i].fx = global_fx[i];
        bodies[i].fy = global_fy[i];
    }
}

void calculateForcesSeq(Node* root, std::vector<Body>& bodies, double theta) {
    for (auto& body : bodies) {
        body.fx = 0.0;
        body.fy = 0.0;
    }
    
    for (size_t i = 0; i < bodies.size(); i += 1) {
        if (bodies[i].mass == -1) {
            continue;
        }
        calcForce(&bodies[i], root, theta);
    }
    
}


void updateBodies(std::vector<Body>& bodies, double dt){
    double ax, ay;
    for (auto& body : bodies) {
        if (body.mass == 0 || body.mass == -1) {
            continue;
        }
        ax = body.fx / body.mass;
        ay = body.fy / body.mass;  
        body.vx = body.vx + ax*dt;
        body.vy = body.vy + ay*dt;
        body.px = body.px + body.vx*dt + 0.5*ax*(dt*dt);
        body.py = body.py + body.vy*dt + 0.5*ay*(dt*dt);
        
        if (body.px < DOMAIN_MIN || body.px > DOMAIN_MAX || 
            body.py < DOMAIN_MIN || body.py > DOMAIN_MAX) {
            body.mass = -1;
        }
    }
}

void writeBodiestoFile(const std::vector<Body>& bodies, const std::string& outputFile, int rank) {
    if (rank != 0) {
        return;
    }
    
    std::ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open output file " << outputFile << std::endl;
        return;
    }
    
    int validBodies = 0;
    for (const auto& body : bodies) {
        validBodies++;
    }
    
    outFile << validBodies << std::endl;
    
    for (const auto& body : bodies) {
        outFile << body.index << "  "
                << std::scientific << std::setprecision(6) << body.px << "  "
                << std::scientific << std::setprecision(6) << body.py << "  "
                << std::scientific << std::setprecision(6) << body.mass << "  "
                << std::scientific << std::setprecision(6) << body.vx << "  "
                << std::scientific << std::setprecision(6) << body.vy << std::endl;

    }
    
    outFile.close();
}


std::vector<Body> readBodiesFromFile(const std::string& filename) {
    std::vector<Body> bodies;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return bodies;
    }
    
    int numBodies;
    file >> numBodies;
    
    bodies.reserve(numBodies);
    
    for (int i = 0; i < numBodies; i++) {
        Body body;
        file >> body.index >> body.px >> body.py >> body.mass >> body.vx >> body.vy;
        
        body.fx = 0.0;
        body.fy = 0.0;
        
        bodies.push_back(body);
    }
    
    file.close();
    return bodies;
}


int parallel_imp(int argc, char* argv[], std::string& inputFile, std::string& outputFile, 
                   int& steps, double& theta, double& dt, bool& visualization, bool& sequential){

    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<Body> bodies;
    bodies = readBodiesFromFile(inputFile);
    MPI_Barrier(MPI_COMM_WORLD);

    Node* root = constructTree(bodies);
    
    MPI_Barrier(MPI_COMM_WORLD);

    double start_time = MPI_Wtime();


    for (int step = 0; step < steps; step++) {
        calculateForcesParallel(root, bodies, theta, rank, size);

        updateBodies(bodies, dt);

        MPI_Barrier(MPI_COMM_WORLD);
        
        destroyTree(root);
        root = constructTree(bodies);
    }
    
    double end_time = MPI_Wtime();

    
    writeBodiestoFile(bodies, outputFile, rank);

    destroyTree(root);
    if (rank == 0) {
        std::cout << std::fixed << std::setprecision(6) << (end_time - start_time) << std::endl;
    }


    MPI_Finalize();
    return 0;

}


int sequential_imp(int argc, char* argv[], std::string& inputFile, std::string& outputFile, 
                   int& steps, double& theta, double& dt, bool& visualization, bool& sequential){


    // Read input file
    std::vector<Body> bodies;
    bodies = readBodiesFromFile(inputFile);
    // if (rank == 0) {
    //     std::cout << "Read " << bodies.size() << " bodies from " << inputFile << std::endl;
    // }

    // Initialize the Barnes-Hut tree
    Node* root = constructTree(bodies);
    
    // Print the tree structure for debugging (only from rank 0)
    // if (rank == 0) {
    //     std::cout << "\nBarnes-Hut Tree Structure (BEFORE):" << std::endl;
    //     printTree(root);
    //     std::cout << std::endl;
    // }

    // Use standard C++ timing instead of MPI_Wtime
    auto start_time = std::chrono::high_resolution_clock::now();

    // Main simulation loop
    for (int step = 0; step < steps; step++) {
        // TODO: Calculate forces on each body
        // TODO: Update positions and velocities
            // if (rank == 0) {
            //   std::cout << "step: " << step << std::endl;  
            // }

        // Calculate forces in parallel
        calculateForcesSeq(root, bodies, theta);
        // Update positions and velocities
        updateBodies(bodies, dt);

        
        destroyTree(root);
        root = constructTree(bodies);
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;
    
    // if (rank == 0) {
    //     std::cout << "\nBarnes-Hut Tree Structure (AFTER):" << std::endl;
    //     printTree(root);
    //     std::cout << std::endl;
    // }

    // Write final state to output file - use rank 0 for sequential version
    writeBodiestoFile(bodies, outputFile, 0);

    // Clean up
    destroyTree(root);
    std::cout << std::fixed << std::setprecision(6) << elapsed_seconds.count() << std::endl;

    return 0;
}

int main(int argc, char* argv[]) {
    std::string inputFile, outputFile;
    int steps;
    double theta, dt;
    bool visualization;
    bool sequential;
    
    if (!parseArguments(argc, argv, inputFile, outputFile, steps, theta, dt, visualization, sequential)) {
        std::cerr << "Usage: " << argv[0] << " -i <input_file> -o <output_file> "
                    << "-s <steps> -t <theta> -d <dt> [-V] [-S]" << std::endl;
        return 1;
    }
    
    if (!sequential){
        return parallel_imp(argc, argv, inputFile, outputFile, steps, theta, dt, visualization, sequential);
    }else{
        return sequential_imp(argc, argv, inputFile, outputFile, steps, theta, dt, visualization, sequential);
    }
    
}

