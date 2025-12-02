# Makefile pour construire le programme MPI qui génère le graphe .dot
# à partir du fichier FASTA (BuildDotMPI.cpp)

# Compilateur MPI C++
CXX     = mpic++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pedantic

# Nom de l'exécutable
TARGET  = build_dot

# Sources
SRCS    = BuildMatrixMPI.cpp
OBJS    = $(SRCS:.cpp=.o)

# Règle principale
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Nettoyage
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
