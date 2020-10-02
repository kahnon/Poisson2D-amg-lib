#pragma once

#include <vector>
#include <string>
#include <math.h>
#include <functional>
#include "matrix.h"

//these base iterations should only be used for smoothing or as solver kernel
//no error checks necessary here

//sor iteration with system matrix "matrix", rhs "rhs"
//and guess phi which should be 0 for coarse levels which afterwards stores the solution
//phi_old is passed only to assure the same function pointer as jacobi may be used
//XXX runtime critical
//only base iteration for inner domain, not for boundary points!
inline void sor_base(int i, const Matrix& matrix, const std::vector<double>& rhs, const std::vector<double>& phi_old, std::vector<double>& phi, double omega=1){
  const line& ln = matrix[i];
  //matrix vector multiplication
  double phi0 = rhs[i] - ln.val*phi[ln.ind];
  for(auto & cn : ln.conns) phi0 -= std::get<double>(cn) * phi[std::get<int>(cn)];
  if(ln.val!=0) phi0/=ln.val;
  else{
    std::cout<<"diagonal value is 0 at "<<i<<std::endl;
    exit(-1);
  }
  phi[i] += omega * phi0;
}


//with boundary points
inline void sor_base_with_bounds(int i, const Matrix& matrix, const std::vector<double>& rhs, const std::vector<double>& phi_old, std::vector<double>& phi, double omega=1){
  const line& ln = matrix[i];
  //apply constant bc if necessary
  if (ln.val == 1){ 
    phi[i] = rhs[i];
  }else{
    sor_base(i,matrix,rhs,phi_old,phi,omega);
  }
}


//jacobi iteration with system matrix "matrix", rhs "rhs", damping parameter omega
//and guess phi which should be 0 for coarse levels which afterwards stores the solution
//XXX runtime critical
//TODO put together with sor_base
inline void jacobi_base(int i, const Matrix& matrix, const std::vector<double>& rhs, const std::vector<double>& phi_old, std::vector<double>& phi, double omega=0.6666){
  const line& ln = matrix[i];
  //matrix vector multiplication
  double phi0 = rhs[i] - ln.val*phi[ln.ind];
  for(auto & cn : ln.conns) phi0 -= std::get<double>(cn) * phi_old[std::get<int>(cn)];
  phi0/=ln.val;
  phi[i] += omega * phi0;
}


//with boundary points
inline void jacobi_base_with_bounds(int i, const Matrix& matrix, const std::vector<double>& rhs, const std::vector<double>& phi_old, std::vector<double>& phi, double omega=0.6666){
  const line& ln = matrix[i];
  //apply constant bc if necessary
  if (ln.val == 1){ 
    phi[i] = rhs[i];
  }else{
    jacobi_base(i,matrix,rhs,phi_old,phi,omega);
  }
}


//TODO use template magic to choose between different smoothing strategies
//at compile time
//right now cf smoothing with sor
inline void smoothing(AMGlvl& lvl, int smoothing_steps=1, double omega=1){
  std::vector<double> phi_old;
  for(auto k=0;k<smoothing_steps;++k){
#if 1
    //two sweeps over fine points possibly improves convergence
    for(auto& i:lvl.boundary_list) lvl.phi[i] = lvl.rhs[i];
    for(auto& i:lvl.coarse_list) sor_base(i, lvl.matrix, lvl.rhs, phi_old, lvl.phi);
    for(auto& i:lvl.fine_list) sor_base(i, lvl.matrix, lvl.rhs, phi_old, lvl.phi);
    //for(auto& i:lvl.fine_list) sor_base(i, lvl.matrix, lvl.rhs, phi_old, lvl.phi);
#else
    for(auto k=0;k<lvl.boundary_list.size();++k){ 
      auto& i=lvl.boundary_list[k];
      lvl.phi[i] = lvl.rhs[i];
    }
    for(auto k=0;k<lvl.coarse_list.size();++k){
      auto& i=lvl.coarse_list[k];
      sor_base(i, lvl.matrix, lvl.rhs, phi_old, lvl.phi);
    }
    for(auto k=0;k<lvl.fine_list.size();++k){
      auto& i=lvl.fine_list[k];
      sor_base(i, lvl.matrix, lvl.rhs, phi_old, lvl.phi);
    }
    for(auto k=0;k<lvl.fine_list.size();++k){
      auto& i=lvl.fine_list[k];
      sor_base(i, lvl.matrix, lvl.rhs, phi_old, lvl.phi);
    }
#endif
  }
}
