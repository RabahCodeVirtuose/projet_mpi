#ifndef FORGRAPHMPI_HPP
#define FORGRAPHMPI_HPP

#include <iostream>
#include <map>
#include <string>
#include <graphviz/cgraph.h>

// On évite using namespace std; dans un header pour ne pas polluer partout
// On utilise std::map et std::string explicitement.

int* lectureGrapheMPI(char* f, int* nb_nodes, std::map<std::string,int>* my_nodes);

#endif
