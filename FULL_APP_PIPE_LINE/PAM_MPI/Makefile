CXX      = mpic++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pedantic

SRCS  = main.cpp PAM.cpp Utils.cpp
OBJS  = $(SRCS:.cpp=.o)
TARGET = pam_mpi

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f $(OBJS) $(TARGET)
