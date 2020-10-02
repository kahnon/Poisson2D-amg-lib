#pragma once

#include <vector>
#include <utility>
#include <tuple>
#include <assert.h>


//matrix definition in here as well as some useful matrix functions

//for matrix coarsening
enum Status{
  UNDECIDED=0, //not yet processed
  COARSE, //coarse point
  SFINE, //strongly connected to coarse point
  WFINE, //not strongly connected to coarse poit, but ensure strong connection to fine point
  BOUNDARY //boundary (no strong connections)
};

//for matrix connections
enum ConnType{
  WEAK=0,
  STRONG,
};

//TODO introduce base line class for matrices that do not need coarsening
//then create the following line class by inheritance
//saves some memory
struct line{
  int ind;
  double val;
  Status status;

  std::vector<std::tuple<int,double,ConnType>> conns;

  line():ind(-1),val(0),status(UNDECIDED){}

  line(int _ind, double _val, Status stat = UNDECIDED, std::vector<std::tuple<int,double,ConnType>> _conns = std::vector<std::tuple<int,double,ConnType>>())
    : ind(_ind),val(_val),status(stat){
    conns = _conns;
  }

  //in bytes
  long int memory_consumption(){
    double mem = sizeof(double) + sizeof(int) + sizeof(Status);
    for(auto i = 0u; i<conns.size(); ++i){
      mem += sizeof(int) + sizeof(double) + sizeof(ConnType);
    }
    return mem;
  }
};

typedef std::vector<line> Matrix;




//set strong and weak connections within a matrix
void set_matrix_connections(Matrix& matrix, double thresh);

//set weak and strong transpose connections within matrix
void set_matrix_transpose_connections(Matrix& matrix, double thresh);

//calculate max norm of the residual
std::pair<int,double> calc_max_norm_res(const Matrix& matrix, const std::vector<double>& rhs, const std::vector<double>& phi);

//set up test-system matrix
void matrix_setup(Matrix& matrix, int Nr, int Nz, double thresh);

//set up test-system rhs vector
void fill_rhs(std::vector<double>& dens, int Nr, int Nz);



//some inlined functions, XXX all runtime critical for cg methods
inline double scalar_product(const std::vector<double>& v1, const std::vector<double>& v2){
  assert(v1.size() == v2.size() && "vector dimensions dont match in scalar_product");
  double result=0;
  for(auto i=0u;i<v1.size();++i) result+=v1[i]*v2[i];
  return result;
}

//matrix vector multiplication for square matrices
inline std::vector<double> matrix_vector(const Matrix& matrix, const std::vector<double>& vec){
  std::vector<double> result(matrix.size(),0);

  for(auto& ln : matrix){
    //diagonal element
    result[ln.ind] += ln.val * vec[ln.ind];
    for(auto& cn:ln.conns){
      int i = std::get<int>(cn);
      double val = std::get<double>(cn);
      result[ln.ind]+=val*vec[i];
    }
  }
  return result;
}


inline double find_matrix_val(const Matrix& matrix, int i, int j){
  if(i==j) return matrix[i].val;
  else{
    for(auto& cn:matrix[i].conns){
      if(std::get<int>(cn)==j) return std::get<double>(cn);
    }
    return (double) 0;
  }
}

