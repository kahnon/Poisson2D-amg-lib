#pragma once

#include <vector>
#include <tuple>

#include "matrix.h"
#include "amg_level.h"


double go_down(int i, AMGhierarchy& Hier, double maxerror, int maxsteps);
void go_up(int i, AMGhierarchy& Hier, double maxerror, int maxsteps);

//inlining might offer performance benefit
//interpolation from coarser to finer level
//assumes interpolation matrix is passed
//XXX runtime critical
inline std::vector<double> interpolation(const Matrix& matrix, const std::vector<double>& vec){
  std::vector<double> int_vec(matrix.size(),0);

  for(auto i=0u;i<matrix.size();++i){
    double& element = int_vec[i];
    for(auto& cn:matrix[i].conns){
      int ind=std::get<int>(cn);
      double val=std::get<double>(cn);
      element += val*vec[ind];
    }
  }
  return int_vec;
}

//restriction from finer to coarser level
//assumes restriction matrix is transpose of interpolation matrix
//XXX runtime critical
inline std::vector<double> restriction(const AMGlvl& lvl, const std::vector<double>& vec){
  int coarse_size = lvl.matrix.size();
  const Matrix& matrix=lvl.interpolation;
  std::vector<double> res_vec(coarse_size,0);

  for(auto i=0u;i<matrix.size();++i){
    auto& ln = matrix[i];
    for(auto& cn:ln.conns){
      int ind=std::get<int>(cn);
      double val=std::get<double>(cn);
      double& element = res_vec[ind];
      element += vec[i] * val;
    }
  }
  return res_vec;
}


