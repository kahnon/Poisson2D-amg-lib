#include <vector>
#include <utility>
#include <math.h>
#include <assert.h>
#include <iostream>

#include "matrix.h"


//set strong and weak connections within a matrix
void set_matrix_connections(Matrix& matrix, double thresh){
  //check if connections are strong or weak
  std::vector<double> max_conns(matrix.size(),0);
  for(auto i=0u;i<matrix.size();++i){
    line& ln=matrix[i];
    for(auto& cn:ln.conns){
      double val=std::get<double>(cn);
      //only account for negative connections
      if(val<max_conns[i]) max_conns[i]=val;
    }
  }

  //connection is strong if value of connection is sufficiently strong
  for(auto& ln : matrix){
    for(auto& cn : ln.conns){
      double val=std::get<double>(cn);
      double ind=std::get<int>(cn);
      if(val <= thresh*max_conns[ind]) std::get<ConnType>(cn) = STRONG;
    }
  }
}


//set strong and weak connections within a matrix
void set_matrix_transpose_connections(Matrix& matrix, double thresh){
  //check if connections are strong or weak
  std::vector<double> max_conns(matrix.size(),0);
  for(auto i=0u;i<matrix.size();++i){
    line& ln=matrix[i];
    for(auto& cn:ln.conns){
      double val=std::get<double>(cn);
      //only account for negative connections
      if(val<max_conns[i]) max_conns[i]=val;
    }
  }

  //get strong transpose connections
  //--> one point has strong influence on another one, if that other point's connection to this point is sufficiently large
  for(auto& ln : matrix){
    for(auto& cn : ln.conns){
      int ind=std::get<int>(cn);
      line& c_ln=matrix[ind];
      for(auto& conn : c_ln.conns){
	if(std::get<int>(conn) != ln.ind) continue;
	double val=std::get<double>(conn);
	if(val <= thresh*max_conns[ind]) std::get<ConnType>(cn) = STRONG;
      }
    }
  }
}


//calculate max norm of the residual
std::pair<int,double> calc_max_norm_res(const Matrix& matrix, const std::vector<double>& rhs, const std::vector<double>& phi){
  double error=0;
  int k=0;
  std::vector<double> res = matrix_vector(matrix,phi);
  for(auto i=0u;i<res.size();++i){
    double new_error=fabs(res[i] - rhs[i]);
    if( new_error > error){ 
      error=new_error;
      k=i;
    }
  }
  return std::make_pair(k,error);
}


//test matrix setup
void matrix_setup(Matrix& matrix, int Nr, int Nz, double thresh){
  int Ng=Nr*Nz;
  if(!matrix.empty()) matrix.clear();

  matrix.resize(Ng);

  for(auto r=0;r<Nr;++r){
    for(auto z=0;z<Nz;++z){
      int ind=r*Nz+z;
      line& ln = matrix[ind];
      ln.ind=ind;
      ln.status=UNDECIDED;

      if(z==0 || r==0 || r==Nr-1 || z==Nz-1){ 
	ln.val=1;
	ln.status=BOUNDARY;
      }
      else{
	ln.val=4;
	ln.conns.emplace_back(ind-Nz,-1,WEAK);
	ln.conns.emplace_back(ind-1,-1,WEAK);
	ln.conns.emplace_back(ind+1,-1,WEAK);
	ln.conns.emplace_back(ind+Nz,-1,WEAK);
      }
    }
  }
  set_matrix_connections(matrix,thresh);
}


void fill_rhs(std::vector<double>& dens, int Nr, int Nz){
  if(!dens.empty()) dens.clear();

  int Ng = Nr*Nz;
  dens.resize(Ng);

  for(auto r=0;r<Nr;++r){
    for(auto z=0;z<Nz;++z){
      int i=r*Nz+z;
      //dens[i] = 10 * (r+z)/(double) Ng;
      if(r==0 || z==0 || r==Nr-1 || z==Nz-1) dens[i]=0;
      else dens[i]=1;
    }
  }
}
