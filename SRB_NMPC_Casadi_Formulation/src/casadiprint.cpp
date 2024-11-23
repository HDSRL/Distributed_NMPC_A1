#include <sstream>
#include <iostream>
#include <fstream>
#include "nmpc_params.h"
#include <casadi/casadi.hpp>

// Utility function for converting SX to MX and printing
std::string casadiPrintSXorMX(const casadi::SX& sx) {
    casadi::MX mx = casadi::MX(sx);  // Convert SX to MX
    std::ostringstream stream;
    mx.disp(stream);  // Display MX into stream
    return stream.str();
}

// // Overloaded functions for writing to file
// void writeMatrixToFile(const casadi::SX& matrix, const std::string& filename) {
//     std::ofstream file(filename);
//     if (file.is_open()) {
//         file << matrix << std::endl;
//         file.close();
//     } else {
//         std::cerr << "Unable to open file" << std::endl;
//     }
// }

// void writeMatrixToFile(const casadi::DM& matrix, const std::string& filename) {
//     std::ofstream file(filename);
//     if (file.is_open()) {
//         file << matrix << std::endl;
//         file.close();
//     } else {
//         std::cerr << "Unable to open file" << std::endl;
//     }
// }
