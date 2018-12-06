#ifndef IO_HPP
#define IO_HPP

#include <fstream>
#include "types.hpp"


typedef sim::data_type dtype;

void openFileToWrite(std::ofstream& file, std::string filename, std::string dirname = ".");
int readDataFromFile(std::string fname, int N, dtype *m, dtype (*rx),dtype (*ry),dtype (*rz), dtype (*u)[3]);
int readDataFromFile(std::string fname, int N, dtype *m, dtype (*r)[3], dtype (*u)[3]);
int readDataFromFile(std::string fname, int N, dtype (*a)[7]);
void writeDataToFile(int N, dtype (*a)[3], std::ofstream& fs);
void writeDataToFile(int N, dtype (*ax),dtype (*ay),dtype (*az), dtype (*b)[3], std::ofstream& fs);
void writeDataToFile(int N, dtype (*a)[3], dtype (*b)[3], std::ofstream& fs);
void writeDataToFile(int N, dtype (*a)[7], std::ofstream& fs);

#endif  // IO_HPP
