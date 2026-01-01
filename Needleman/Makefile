CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pedantic
OMPFLAGS = -fopenmp

all: test_seq test_omp

test_seq: test_seq.o needleman.o needleman_debug.o
	$(CXX) $(CXXFLAGS) -o $@ $^

test_omp: test_omp.o needleman_omp.o needleman_debug.o
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $^

needleman.o: needleman.cpp needleman_common.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

needleman_omp.o: needleman_omp.cpp needleman_common.hpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -c $< -o $@

needleman_debug.o: needleman_debug.cpp needleman_common.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

test_seq.o: test_seq.cpp needleman_common.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

test_omp.o: test_omp.cpp needleman_common.hpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -c $< -o $@

clean:
	rm -f *.o test_seq test_omp

.PHONY: all clean
