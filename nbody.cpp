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

// Structure to represent a body in the simulation
struct Body {
    int index;
    double x, y;       // Position
    double mass;
    double vx, vy;     // Velocity
    double fx, fy;     // Force
};

// Function prototypes
std::vector<Body> readBodiesFromFile(const std::string& filename);
void writeBodiestoFile(const std::string& filename, const std::vector<Body>& bodies);
void computeForces(std::vector<Body>& bodies, int startIdx, int endIdx, double theta);
void updateBodies(std::vector<Body>& bodies, double dt, int startIdx, int endIdx);

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
    double startTime = MPI_Wtime();
    
    // Read input file (only rank 0 needs to do this initially)
    std::vector<Body> bodies;
    if (rank == 0) {
        bodies = readBodiesFromFile(inputFile);
        if (rank == 0) {
            std::cout << "Read " << bodies.size() << " bodies from " << inputFile << std::endl;
        }
    }
    


    // Main simulation loop
    for (int step = 0; step < steps; step++) {
    }
    
    // Write output file (only rank 0 needs to do this)
    if (rank == 0) {
        writeBodiestoFile(outputFile, bodies);
        std::cout << "Wrote " << bodies.size() << " bodies to " << outputFile << std::endl;
    }
    
    // End timing and print elapsed time
    double endTime = MPI_Wtime();
    if (rank == 0) {
        std::cout << "Simulation completed in " << endTime - startTime << " seconds" << std::endl;
    }
    
    // Free MPI datatype
    // MPI_Type_free(&bodyType);
    
    // Finalize MPI
    MPI_Finalize();
    return 0;
}

// Read bodies from input file
std::vector<Body> readBodiesFromFile(const std::string& filename) {
    std::vector<Body> bodies;
    // Implementation would go here
    return bodies;
}

// Write bodies to output file
void writeBodiestoFile(const std::string& filename, const std::vector<Body>& bodies) {
    return;
}

// Calculate forces between bodies (naive implementation - replace with Barnes-Hut)
void computeForces(std::vector<Body>& bodies, int startIdx, int endIdx, double theta) {
    return;
}

// Update positions and velocities using Leapfrog-Verlet integration
void updateBodies(std::vector<Body>& bodies, double dt, int startIdx, int endIdx) {
    return;
}
