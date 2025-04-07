#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <getopt.h>

// Constants for the simulation
const double G = 0.0001;  // Gravitational constant
const double RLIMIT = 0.03;    // Minimum distance to avoid infinite forces
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
    struct Body b;      // can represent internal / external nodes
    struct Node *NW, *NE, *SW, *SE;
    double minx, maxx, miny, maxy;
};


// Function prototypes
std::vector<Body> readBodiesFromFile(const std::string& filename);
Node* constructTreeHelper(Node* root, const Body& b, double minx, double maxx, double miny, double maxy);

// Parse command line arguments
bool parseArguments(int argc, char* argv[], std::string& inputFile, std::string& outputFile, 
                   int& steps, double& theta, double& dt, bool& visualization) {
    // Set default values
    inputFile = "input/nb-10.txt";
    outputFile = "output.txt";
    steps = 100;
    theta = 0.5;  // Common default for Barnes-Hut
    dt = 0.005;
    visualization = false;
    
    int opt;
    while ((opt = getopt(argc, argv, "i:o:s:t:d:V")) != -1) {
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
            default:
                return false;
        }
    }
    
    // Validate values
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

/*

Constructing the Barnes-Hut tree. To construct the Barnes-Hut tree, 
insert the bodies one after another. To insert a body b into the tree 
rooted at node x, use the following recursive procedure:

If node x does not contain a body, put the new body b here.

If node x is an internal node, update the center-of-mass and total mass of x. 
Recursively insert the body b in the appropriate quadrant.

If node x is an external node, say containing a body named c, then there are 
two bodies b and c in the same region. Subdivide the region further by creating 
four children. Then, recursively insert both b and c into the appropriate 
quadrant(s). Since b and c may still end up in the same quadrant, there may 
be several subdivisions during a single insertion. Finally, update the 
center-of-mass and total mass of x.

*/
Node* constructTreeHelper(Node* root, const Body& b, double minx, double maxx, double miny, double maxy) {
    if(root == nullptr) {
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
        root->b = b;
        return root;
    }
    
    // If this is an empty node, just put the body here
    if(root->b.mass == 0) {
        root->b = b;
    } 
    // If this is an external node (leaf with a single body)
    else if(root->NW == nullptr && root->NE == nullptr && 
            root->SW == nullptr && root->SE == nullptr) {
        // Store existing body
        Body b2 = root->b;
        
        // Mark as internal node
        root->b.index = -1;
        
        // Insert both bodies
        root = insertCorrect(root, b);
        root = insertCorrect(root, b2);
        
        // Update center of mass and total mass
        double totalMass = root->b.mass + b.mass;
        if (totalMass > 0) {  // Prevent division by zero
            // Use a more stable formula for weighted average to avoid overflow
            double w1 = root->b.mass / totalMass;
            double w2 = b.mass / totalMass;
            root->b.px = w1 * root->b.px + w2 * b.px;
            root->b.py = w1 * root->b.py + w2 * b.py;
            root->b.mass = totalMass;
        } else {
            // If both masses are zero, just use the position of the new body
            root->b.px = b.px;
            root->b.py = b.py;
            root->b.mass = 0;
        }
    } 
    // If this is an internal node
    else {
        // Update center of mass and total mass
        double totalMass = root->b.mass + b.mass;
        if (totalMass > 0) {  // Prevent division by zero
            // Use a more stable formula for weighted average to avoid overflow
            double w1 = root->b.mass / totalMass;
            double w2 = b.mass / totalMass;
            root->b.px = w1 * root->b.px + w2 * b.px;
            root->b.py = w1 * root->b.py + w2 * b.py;
            root->b.mass = totalMass;
        } else {
            // If both masses are zero, just use the position of the new body
            root->b.px = b.px;
            root->b.py = b.py;
            root->b.mass = 0;
        }
        
        // Insert the new body
        root = insertCorrect(root, b);
    }
    return root;
}

Node* constructTree(const std::vector<Body>& bodies) {
    Node* root = nullptr;
    
    for(size_t i = 0; i < bodies.size(); i++) {
        root = constructTreeHelper(root, bodies[i], DOMAIN_MIN, DOMAIN_MAX, DOMAIN_MIN, DOMAIN_MAX);
    }
    return root;
}

void destroyTree(Node* root) {
    if(root == nullptr) return;
    
    // Recursively destroy all children
    destroyTree(root->NW);
    destroyTree(root->NE);
    destroyTree(root->SW);
    destroyTree(root->SE);
    
    // Free this node
    delete root;
}

// Function to print the tree structure for debugging
void printTree(Node* root, int level = 0, const std::string& prefix = "") {
    if (root == nullptr) {
        std::cout << prefix << "└── [empty]" << std::endl;
        return;
    }
    
    // Print current node
    std::cout << prefix;
    if (level == 0) {
        std::cout << "Root: ";
    } else {
        std::cout << "└── ";
    }
    
    if (root->b.index == -1) {
        std::cout << "Internal Node (mass=" << root->b.mass 
                  << ", pos=(" << root->b.px << "," << root->b.py << ")"
                  << ", bounds=[" << root->minx << "," << root->maxx << "][" 
                  << root->miny << "," << root->maxy << "])" << std::endl;
    } else {
        std::cout << "Body " << root->b.index 
                  << " (mass=" << root->b.mass 
                  << ", pos=(" << root->b.px << "," << root->b.py << ")"
                  << ", vel=(" << root->b.vx << "," << root->b.vy << "))" << std::endl;
    }
    
    // Print children with increased indentation
    std::string newPrefix = prefix + "    ";
    printTree(root->NW, level + 1, newPrefix + "NW: ");
    printTree(root->NE, level + 1, newPrefix + "NE: ");
    printTree(root->SW, level + 1, newPrefix + "SW: ");
    printTree(root->SE, level + 1, newPrefix + "SE: ");
}

int main(int argc, char* argv[]) {
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
    
    if (!parseArguments(argc, argv, inputFile, outputFile, steps, theta, dt, visualization)) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " -i <input_file> -o <output_file> "
                      << "-s <steps> -t <theta> -d <dt> [-V]" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }
    
    // Print parameters (only from rank 0)
    if (rank == 0) {
        std::cout << "Parameters:" << std::endl;
        std::cout << "  Input file: " << inputFile << std::endl;
        std::cout << "  Output file: " << outputFile << std::endl;
        std::cout << "  Steps: " << steps << std::endl;
        std::cout << "  Theta: " << theta << std::endl;
        std::cout << "  Timestep: " << dt << std::endl;
        std::cout << "  Visualization: " << (visualization ? "On" : "Off") << std::endl;
    }
    
    // Read input file
    std::vector<Body> bodies;
    if (rank == 0) {
        bodies = readBodiesFromFile(inputFile);
        std::cout << "Read " << bodies.size() << " bodies from " << inputFile << std::endl;
    }

    Node* root;
    // Print the tree structure for debugging (only from rank 0)
    if (rank == 0) {
        // Initialize the Barnes-Hut tree
        root = constructTree(bodies);
        std::cout << "\nBarnes-Hut Tree Structure:" << std::endl;
        printTree(root);
        std::cout << std::endl;
    }

    for (int i = 0; i < steps; i++){
        // main simulation loop
        // For each new step, you would destroy the old tree and create a new one
        // after the bodies have moved
        
        // If you need to rebuild the tree each step:
        // destroyTree(root);
        // root = constructTree(bodies);
    }

    // Clean up by freeing all allocated memory
    destroyTree(root);
    
    MPI_Finalize();
    return 0;
}


// Read bodies from input file
std::vector<Body> readBodiesFromFile(const std::string& filename) {
    std::vector<Body> bodies;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return bodies;
    }
    
    // Read number of bodies
    int numBodies;
    file >> numBodies;
    
    // Reserve space to avoid reallocations
    bodies.reserve(numBodies);
    
    // Read each body
    for (int i = 0; i < numBodies; i++) {
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

