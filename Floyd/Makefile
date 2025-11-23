CC := g++
CFLAGS := -std=c++11 -Wall -Wextra -pedantic -O3
LDFLAGS := -lcgraph

SRC := $(wildcard *.cpp)
OBJ := $(SRC:.cpp=.o)

.PHONY: all clean

all: main

main: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o main *~
