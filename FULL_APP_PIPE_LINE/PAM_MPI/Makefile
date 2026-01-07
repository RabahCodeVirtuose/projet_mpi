CXX      = mpic++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pedantic -fopenmp
LDFLAGS  = -fopenmp

SRCS  = main.cpp PAM.cpp Utils.cpp
OBJS  = $(SRCS:.cpp=.o)
TARGET = pam_mpi

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
