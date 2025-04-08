# Compiler and flags
CXX = mpicxx
CXXFLAGS = -std=c++11 -Wall -O3 -g
LDFLAGS = 

# Source files and target
SRC = nbody.cpp
HEADERS = 
OBJ = $(SRC:.cpp=.o)
TARGET = nbody

# Default target
all: $(TARGET)

# Linking
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compiling
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJ) $(TARGET) *~

# Run example (not called by default)
run: $(TARGET)
	mpirun -n 4 ./$(TARGET) -i input/nb-10.txt -o output/output.txt -s 100 -t 0.005 -d 0.01

.PHONY: all clean run
