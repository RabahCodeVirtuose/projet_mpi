# ==============================
#       Makefile MPI propre
# ==============================

CXX = mpic++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -pedantic
LIBS = -lcgraph

SRC = main_mpi.cpp \
      ForGraphMPI.cpp \
      Distribution.cpp \
      ParallelFWBlocks.cpp\
      Utils.cpp\

OBJ = $(SRC:.cpp=.o)

TARGET = main_mpi

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGET)
