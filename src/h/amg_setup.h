#pragma once 

#include "matrix.h"

//beta = 0.35 good choice via GMG_2_AMG
void coarsen_matrix(Matrix& matrix, double thresh=0.25, double beta=0.35);

Matrix create_interpolation_matrix(Matrix& matrix, double thresh);

//for Poisson-like problems thresh=0.25 is good guess (ruge/stueben, "amg introduction")
Matrix create_coarse_matrix(const Matrix& matrix, const Matrix& int_matrix);
