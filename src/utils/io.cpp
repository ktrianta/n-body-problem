#include "io.hpp"


std::string baseName(std::string const& path) {
    return path.substr(path.find_last_of("/") + 1);
}

void openFileToWrite(std::ofstream& fs, std::string filename, std::string dirname) {
    std::string out_filename = dirname + "/" + baseName(filename) + ".out";
    fs.open(out_filename);
}



int readDataFromFile(std::string filename, int N, dtype *m, dtype (*rx), dtype (*ry),dtype (*rz), dtype (*u)[3]) {
    std::ifstream fs;

    fs.open(filename);
    if (!fs.is_open()) {
        return -1;
    }

    for (size_t i = 0; i < N; i++) {
        fs >> m[i]
           >> rx[i] >> ry[i] >> rz[i]
           >> u[i][0] >> u[i][1] >> u[i][2];
    }

    fs.close();
    return 0;
}


int readDataFromFile(std::string filename, int N, dtype *m, dtype (*r)[3], dtype (*u)[3]) {
    std::ifstream fs;

    fs.open(filename);
    if (!fs.is_open()) {
        return -1;
    }

    for (size_t i = 0; i < N; i++) {
        fs >> m[i]
           >> r[i][0] >> r[i][1] >> r[i][2]
           >> u[i][0] >> u[i][1] >> u[i][2];
    }

    fs.close();
    return 0;
}

int readDataFromFile(std::string filename, int N, dtype (*a)[7]) {
    std::ifstream fs;

    fs.open(filename);
    if (!fs.is_open()) {
        return -1;
    }

    for (size_t i = 0; i < N; i++) {
        fs >> a[i][0]
           >> a[i][1] >> a[i][2] >> a[i][3]
           >> a[i][4] >> a[i][5] >> a[i][6];
    }

    fs.close();
    return 0;
}

void writeDataToFile(int N, dtype (*a)[3], std::ofstream& fs) {
    for (size_t i = 0; i < N; i++) {
        fs << a[i][0] << "    "
           << a[i][1] << "    "
           << a[i][2] << "\n";
    }
}

void writeDataToFile(int N, dtype (*ax),dtype (*ay),dtype (*az), dtype (*b)[3], std::ofstream& fs) {
    for (size_t i = 0; i < N; i++) {
        fs << ax[i] << "    "
           << ay[i] << "    "
           << az[i] << "    "
           << b[i][0] << "    "
           << b[i][1] << "    "
           << b[i][2] << "\n";
    }
}


void writeDataToFile(int N, dtype (*a)[3], dtype (*b)[3], std::ofstream& fs) {
    for (size_t i = 0; i < N; i++) {
        fs << a[i][0] << "    "
           << a[i][1] << "    "
           << a[i][2] << "    "
           << b[i][0] << "    "
           << b[i][1] << "    "
           << b[i][2] << "\n";
    }
}


void writeDataToFile(int N, dtype (*a)[7], std::ofstream& fs) {
    for (size_t i = 0; i < N; i++) {
        fs << a[i][1] << "    "
           << a[i][2] << "    "
           << a[i][3] << "    "
           << a[i][4] << "    "
           << a[i][5] << "    "
           << a[i][6] << "\n";
    }
}
