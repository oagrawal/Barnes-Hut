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
    int index;
    double x, y;       // Position
    double mass;
    double vx, vy;     // Velocity
    double fx, fy;     // Force
};


// Function prototypes
std::vector<Body> readBodiesFromFile(const std::string& filename);
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
    
    // Start timing
    // double startTime = MPI_Wtime();
    
    // Read input file (only rank 0 needs to do this initially)
    std::vector<Body> bodies;
    if (rank == 0) {
        bodies = readBodiesFromFile(inputFile);
        if (rank == 0) {
            std::cout << "Read " << bodies.size() << " bodies from " << inputFile << std::endl;
        }
    }


    

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
        file >> body.index >> body.x >> body.y >> body.mass >> body.vx >> body.vy;
        
        // Initialize forces to zero
        body.fx = 0.0;
        body.fy = 0.0;
        
        bodies.push_back(body);
    }
    
    file.close();
    return bodies;
}

