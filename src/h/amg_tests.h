#pragma once

#include <vector>
#include "matrix.h"
#include "amg_solve_components.h"
#include "amg_level.h"

std::vector<std::vector<double>> matrix_vector_to_dense(const Matrix& m1, const Matrix& m2, int m1size1, int m1size2, int m2size1,int m2size2);
void continued_restriction(AMGhierarchy& hier, std::vector<double>& vec, int imin, int imax);
void continued_interpolation(AMGhierarchy& hier, std::vector<double>& vec, int imax, int imin);
void test_smoothing(AMGhierarchy& hier, double maxerror, int maxsteps, double omega);
