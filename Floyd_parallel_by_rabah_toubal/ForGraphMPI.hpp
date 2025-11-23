#ifndef FORGRAPHMPI_HPP
#define FORGRAPHMPI_HPP

#include <iostream>
#include <map>
#include <graphviz/cgraph.h>

int* lectureGrapheMPI(char* f, int* nb_nodes, std::map<std::string,int>* my_nodes);

#endif
