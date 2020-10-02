#pragma once

#include <fstream>
#include <string>

#include "matrix.h"
#include "amg_level.h"

void print_coarsening(const Matrix& matrix, int Nr, int Nz, std::string fn);
void print_matrix(const Matrix& matrix, std::string fn);
void print_vector(const std::vector<double>& vec, int Nr, int Nz, std::string fn);
void print_hierarchy(AMGhierarchy& hier, int Nr, int Nz, std::string fn);

