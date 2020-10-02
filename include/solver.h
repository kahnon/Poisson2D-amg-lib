#pragma once
#include <vector>

void calculate_potential( std::vector<double>& phi, std::vector<double>& rhs, double& field_time ) ;

void init_solver(const std::vector<double>& Eps, int NZ, int NR);

void destroy_solver();
