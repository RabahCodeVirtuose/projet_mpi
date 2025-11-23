#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

std::vector<string> Lecture(string file_name)
{
    std::vector<string> V;
    ifstream f;
    f.open(file_name, ifstream::in);
    string line;
    string seq = "";
    getline(f, line);
    getline(f, line);
    while (!f.eof()) {
        if (line[0] != '>')
            seq += line;
        else {
            V.push_back(seq);
            seq = "";
        }
        getline(f, line);
    }
    V.push_back(seq);
    f.close();
    return V;
}
