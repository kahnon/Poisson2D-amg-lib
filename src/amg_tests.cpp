#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>

#include "matrix.h"
#include "amg_solve_components.h"
#include "amg_level.h"
#include "smoothing.h"



//TODO should be removed when solver is finished
std::vector<std::vector<double>> matrix_vector_to_dense(const Matrix& m1, const Matrix& m2, int m1size1, int m1size2, int m2size1,int m2size2){
  //first convert m1 and m2 to dense matrices then multiply

  auto find_matrix_val = [](const Matrix& matrix, int i, int j){
    if(i==j) return matrix[i].val;
    else{
      for(auto& cn:matrix[i].conns){
	if(std::get<int>(cn)==j) return std::get<double>(cn);
      }
      return (double) 0;
    }
  };

  //initialize matrix
  assert(m1size2 == m2size1 && "cant multiply dense matrices");
  std::vector<std::vector<double>> tmp(m1size1);
  std::for_each(tmp.begin(),tmp.end(),[&m2size2](std::vector<double>& v){v.resize(m2size2,0);});

  //calc m1*m2
  for(auto i=0u;i<tmp.size();++i){
    for(auto j=0u;j<tmp[0].size();++j){
      double sum=0;
      for(auto k=0;k<m1size2;++k){
	sum += find_matrix_val(m1,i,k)*find_matrix_val(m2,k,j);
      }
      tmp[i][j]=sum;
    }
  }


  std::vector<std::vector<double>> res(m2size2);
  std::for_each(res.begin(),res.end(),[&m2size2](std::vector<double>& v){v.resize(m2size2,0);});
  //calc m2^T*(m1*m2)
  for(auto i=0u;i<res.size();++i){
    for(auto j=0u;j<res[0].size();++j){
      double sum=0;
      for(auto k=0;k<m1size1;++k){
	sum+=0.5 * find_matrix_val(m2,k,i) * tmp[k][j];
      }
      res[i][j]=sum;
    }
  }
  return res;
}


//restrict from level imin up to level imax
void continued_restriction(AMGhierarchy& hier, std::vector<double>& vec, int imin, int imax){
  assert(vec.size() == hier[imin].matrix.size() && "sizes in continued restriction dont match");
  std::vector<std::vector<double>> phi_vec(1,vec);

  for(auto i=imin+1;i<=imax;++i){
    int size = hier[i].matrix.size();
    phi_vec.emplace_back(std::vector<double>(size,0));
  }

  for(auto i=imin;i<imax;++i){ 
    AMGlvl& nextlvl=hier[i+1];
    const std::vector<double>& current_phi = phi_vec[i];
    std::vector<double>& next_phi = phi_vec[i+1];
    next_phi = restriction(nextlvl,current_phi);
#if 0
    for(auto k=0;k<next_phi.size();++k){
      std::cout<<"k="<<k<<"   phi[k]="<<next_phi[k]<<std::endl;
    }
    std::cout<<std::endl<<std::endl;
#endif
  }
  vec = phi_vec.back();
}


//restrict from level imin up to level imax
void continued_interpolation(AMGhierarchy& hier, std::vector<double>& vec, int imax, int imin){
  assert(vec.size() == hier[imax].matrix.size() && "sizes in continued interpolation dont match");


  std::vector<std::vector<double>> phi_vec;
  for(auto i=imin;i<=imax-1;++i){
    int size = hier[i].matrix.size();
    phi_vec.emplace_back(std::vector<double>(size,0));
  }
  phi_vec.push_back(vec);

  for(auto i=imax;i>imin;--i){ 
    AMGlvl& currentlvl=hier[i];
    const std::vector<double>& current_phi = phi_vec[i];
    std::vector<double>& next_phi = phi_vec[i-1];
    next_phi = interpolation(currentlvl.interpolation,current_phi);
#if 0
    for(auto k=0;k<next_phi.size();++k){
      std::cout<<"k="<<k<<"   phi[k]="<<next_phi[k]<<std::endl;
    }
    std::cout<<std::endl<<std::endl;
#endif
  }
  vec = phi_vec.front();
}


//code duplication not a big deal here
void test_smoothing(AMGhierarchy& hier, double maxerror, int maxsteps, double omega){
  //testing gauss-seidel iteration
  std::cout<<"SOR smoothing:"<<std::endl;
  for(auto& lvl:hier){
    int size = lvl.matrix.size();
    std::vector<double> phi(size,1);
    std::vector<double> phi_old(phi);
    std::vector<double> rhs(size,1);
    std::cout<<"AMG level="<<lvl.level<<"   ";
    for(auto k=0;k<maxsteps;++k){
      for(auto i=0u;i<phi.size();++i){
	sor_base_with_bounds(i,lvl.matrix,rhs,phi_old,phi,omega);
      }
    }
    std::cout<<"error = "<<calc_max_norm_res(lvl.matrix,rhs,phi).second<<std::endl;
  }

  std::cout<<std::endl<<"Jacobi smoothing:"<<std::endl;
  for(auto& lvl:hier){
    int size = lvl.matrix.size();
    std::vector<double> phi(size,1);
    std::vector<double> phi_old(phi);
    std::vector<double> rhs(size,1);
    std::cout<<"AMG level="<<lvl.level<<"   ";
    for(auto k=0;k<maxsteps;++k){
      for(auto i=0u;i<phi.size();++i){
	jacobi_base_with_bounds(i,lvl.matrix,rhs,phi_old,phi,omega);
      }
      std::copy(phi.begin(),phi.end(),phi_old.begin());
    }
    std::cout<<"error = "<<calc_max_norm_res(lvl.matrix,rhs,phi).second<<std::endl;
  }
  std::cout<<std::endl;
}
