# ==============================================================================
# Makefile pour Needleman-Wunsch OpenMP
# ==============================================================================

CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -march=native
OMPFLAGS = -fopenmp

SRC = build_matrix_needleman.cpp
TARGET = build_matrix_needleman

DATASET = ../../DATA/dataset_2000seq.fa

.PHONY: all clean run benchmark

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $< -o $@

run: $(TARGET)
	OMP_NUM_THREADS=12 ./$(TARGET) $(DATASET)

benchmark: $(TARGET)
	@echo "Threads | Temps (s)"
	@echo "--------|----------"
	@for T in 1 2 4 8 12; do \
		echo -n "  $$T     | "; \
		OMP_NUM_THREADS=$$T ./$(TARGET) $(DATASET) 2>&1 | grep "Temps:" | awk '{print $$2}'; \
	done

clean:
	rm -f $(TARGET) ../../DATA/Resulat_sequence_by_premier_algo.dot