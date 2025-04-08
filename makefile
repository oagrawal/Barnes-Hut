# Compiler and flags
CXX = g++
MSMPI_INC = "C:\Program Files (x86)\Microsoft SDKs\MPI\Include"
MSMPI_LIB64 = "C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64"
CXXFLAGS = -O2 -std=c++11 -I$(MSMPI_INC)
LDFLAGS = -L$(MSMPI_LIB64) -lmsmpi -lgdi32

# Source files and target
SRC = nbody.cpp
HEADERS = 
OBJ = $(SRC:.cpp=.o)
TARGET = nbody.exe

# Default target
all: $(TARGET)

# Linking
$(TARGET): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

# Compiling
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	del /F /Q $(OBJ) $(TARGET) *~

# Run example (not called by default)
run: $(TARGET)
	mpiexec -n 4 .\$(TARGET) -i input\nb-10.txt -o output\output.txt -s 100 -t 0.005 -d 0.01 -V

.PHONY: all clean run
