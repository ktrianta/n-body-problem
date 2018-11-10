#ifndef IO_HPP
#define IO_HPP

#include <fstream>
#include "types.hpp"


typedef sim::data_type dtype;

std::ofstream openFileToWrite(std::string filename, std::string dirname = ".");
int readDataFromFile(std::string fname, int N, dtype *m, dtype (*r)[3], dtype (*u)[3]);
void writeDataToFile(int N, dtype (*a)[3], std::ofstream& fs);
void writeDataToFile(int N, dtype (*a)[3], dtype (*b)[3], std::ofstream& fs);

#endif  // IO_HPP
