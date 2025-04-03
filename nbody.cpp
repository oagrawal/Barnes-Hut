#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

// Constants for the simulation
const double G = 6.67430e-11;  // Gravitational constant
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
void computeForces(std::vector<Body>& bodies, int startIdx, int endIdx);
void updateBodies(std::vector<Body>& bodies, double dt, int startIdx, int endIdx);

int main(int argc, char* argv[]) {
    // Initialize MPI environment
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Parse command line arguments
    if (argc != 5) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <num_iterations> <timestep> <input_file> <output_file>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }
    
    // Parse command line arguments according to specification
    double iterations = std::stod(argv[1]);
    double dt = std::stod(argv[2]);
    std::string inputFile = argv[3];
    std::string outputFile = argv[4];
    
    // Convert iterations to int for the loop
    int steps = static_cast<int>(iterations);
    
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
    
    // Broadcast the number of bodies to all processes
    int numBodies = 0;
    if (rank == 0) {
        numBodies = bodies.size();
    }
    MPI_Bcast(&numBodies, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Create MPI datatype for Body struct
    MPI_Datatype bodyType;
    int blocklengths[7] = {1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype types[7] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[7];
    
    offsets[0] = offsetof(Body, index);
    offsets[1] = offsetof(Body, x);
    offsets[2] = offsetof(Body, y);
    offsets[3] = offsetof(Body, mass);
    offsets[4] = offsetof(Body, vx);
    offsets[5] = offsetof(Body, vy);
    offsets[6] = offsetof(Body, fx);
    
    MPI_Type_create_struct(7, blocklengths, offsets, types, &bodyType);
    MPI_Type_commit(&bodyType);
    
    // If not rank 0, resize the bodies vector
    if (rank != 0) {
        bodies.resize(numBodies);
    }
    
    // Broadcast all bodies to all processes
    MPI_Bcast(bodies.data(), numBodies, bodyType, 0, MPI_COMM_WORLD);
    
    // Determine body range for this process
    int bodiesPerProcess = numBodies / size;
    int remainder = numBodies % size;
    int startIdx = rank * bodiesPerProcess + std::min(rank, remainder);
    int endIdx = startIdx + bodiesPerProcess + (rank < remainder ? 1 : 0);
    
    // Main simulation loop
    for (int step = 0; step < steps; step++) {
        // Reset forces
        for (int i = startIdx; i < endIdx; i++) {
            bodies[i].fx = bodies[i].fy = 0.0;
        }
        
        // Compute forces - this is where you would implement Barnes-Hut algorithm
        computeForces(bodies, startIdx, endIdx);
        
        // Collect all force calculations from all processes
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                      bodies.data(), numBodies, bodyType, MPI_COMM_WORLD);
        
        // Update positions and velocities
        updateBodies(bodies, dt, startIdx, endIdx);
        
        // Sync updated positions and velocities
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                     bodies.data(), numBodies, bodyType, MPI_COMM_WORLD);
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
    MPI_Type_free(&bodyType);
    
    // Finalize MPI
    MPI_Finalize();
    return 0;
}

// Read bodies from input file
std::vector<Body> readBodiesFromFile(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error opening input file: " << filename << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    int n;
    inFile >> n;
    
    std::vector<Body> bodies(n);
    for (int i = 0; i < n; i++) {
        inFile >> bodies[i].index >> bodies[i].x >> bodies[i].y 
               >> bodies[i].mass >> bodies[i].vx >> bodies[i].vy;
        bodies[i].fx = bodies[i].fy = 0.0;
    }
    
    inFile.close();
    return bodies;
}

// Write bodies to output file
void writeBodiestoFile(const std::string& filename, const std::vector<Body>& bodies) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error opening output file: " << filename << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    outFile << bodies.size() << std::endl;
    for (const auto& body : bodies) {
        outFile << body.index << "\t" 
                << std::scientific << body.x << "\t" 
                << std::scientific << body.y << "\t" 
                << std::scientific << body.mass << "\t" 
                << std::scientific << body.vx << "\t" 
                << std::scientific << body.vy << std::endl;
    }
    
    outFile.close();
}

// Calculate forces between bodies (naive implementation - replace with Barnes-Hut)
void computeForces(std::vector<Body>& bodies, int startIdx, int endIdx) {
    // TODO: Replace this with Barnes-Hut algorithm
    // This is a naive O(n²) direct calculation
    for (int i = startIdx; i < endIdx; i++) {
        for (int j = 0; j < bodies.size(); j++) {
            if (i == j) continue;
            
            double dx = bodies[j].x - bodies[i].x;
            double dy = bodies[j].y - bodies[i].y;
            double dist = sqrt(dx*dx + dy*dy);
            
            // Apply minimum distance constraint
            if (dist < RLIMIT) dist = RLIMIT;
            
            // Calculate gravitational force
            double force = G * bodies[i].mass * bodies[j].mass / (dist * dist);
            
            // Compute force components
            bodies[i].fx += force * dx / dist;
            bodies[i].fy += force * dy / dist;
        }
    }
}

// Update positions and velocities using Leapfrog-Verlet integration
void updateBodies(std::vector<Body>& bodies, double dt, int startIdx, int endIdx) {
    for (int i = startIdx; i < endIdx; i++) {
        // Calculate acceleration
        double ax = bodies[i].fx / bodies[i].mass;
        double ay = bodies[i].fy / bodies[i].mass;
        
        // Update position
        bodies[i].x += bodies[i].vx * dt + 0.5 * ax * dt * dt;
        bodies[i].y += bodies[i].vy * dt + 0.5 * ay * dt * dt;
        
        // Update velocity
        bodies[i].vx += ax * dt;
        bodies[i].vy += ay * dt;
        
        // Check if body is outside domain
        if (bodies[i].x < DOMAIN_MIN || bodies[i].x > DOMAIN_MAX || 
            bodies[i].y < DOMAIN_MIN || bodies[i].y > DOMAIN_MAX) {
            bodies[i].mass = -1;  // Mark as lost
        }
    }
}
