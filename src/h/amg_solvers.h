#pragma once

#include <vector>
#include <functional>
#include "matrix.h"
#include "amg_level.h"


void sor(const Matrix& matrix, const std::vector<double>& rhs, std::vector<double>& phi, int error_interval=10, double maxerror=1e-12, int maxstep=2000, double omega=1.7);

void jacobi(const Matrix& matrix, const std::vector<double>& rhs, std::vector<double>& phi, int error_interval=10, double maxerror=1e-12, int maxstep=5000, double omega=1);

void cg(const Matrix& matrix, const std::vector<double>& rhs, std::vector<double>& phi, int maxsteps, double acc=1e-12);

void solve_v_cycle(AMGhierarchy& Hier, const std::vector<double>& rho, std::vector<double>& phi, int maxsteps, 
		   int smooth_steps=5, bool gen_flag=0, double maxerror=1e-12);

//TODO pass parameters of amg
void amg_preconditioned_cg(AMGhierarchy& hier, const std::vector<double>& rhs, std::vector<double>& phi, 
			   std::function<void(AMGhierarchy&,const std::vector<double>&,std::vector<double>&,int,int,bool,double)> precondition, 
			   bool gen_flag, int maxsteps, double acc=1e-12);
