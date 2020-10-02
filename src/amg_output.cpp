#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

#include "matrix.h"
#include "amg_level.h"
#include "amg_tests.h"

void print_coarsening(const Matrix& matrix, int Nr, int Nz, std::string fn){
  std::ofstream ofs(fn);
  for(auto r=0;r<Nr;++r){
    for(auto z=0;z<Nz;++z){
      int i=r*Nz+z;
      ofs<<std::setw(6)<<(int) matrix[i].status;
    }
    ofs<<std::endl;
  }
}

void print_matrix(const Matrix& matrix, std::string fn){
  std::ofstream ofs(fn);
  for(auto i=0u;i<matrix.size();++i){
    auto& ln = matrix[i];
    ofs<< ln.ind<<"   "<<std::setw(6)<< ln.val<<std::setw(6)<< (int) ln.status<<std::endl;
    for(auto& cn : ln.conns){
      ofs<<"-->"<<"    "<<std::get<int>(cn)<<"    "<<std::setw(6)<<std::get<double>(cn)<<std::setw(6)<<(int) std::get<ConnType>(cn)<<std::endl;
    }
    ofs<<std::endl;
  }
}


void print_vector(const std::vector<double>& vec, int Nr, int Nz, std::string fn){
  std::ofstream ofs(fn);
  for(auto r=0;r<Nr;++r){
    for(auto z=0;z<Nz;++z){
      int i=r*Nz+z;
      ofs<<std::setw(10)<<vec[i];
    }
    ofs<<std::endl;
  }
}


void print_hierarchy(AMGhierarchy& hier, int Nr, int Nz, std::string fn){
  size_t size = hier.front().size();
  std::vector<int> hmap(size,0);

  for(auto k=1u;k<hier.size();++k){
    Matrix& intp=hier[k].interpolation;
    size_t ctr=0;

    for(auto i=0u;i<size && ctr<intp.size() ;++i){
      if(hmap[i]==static_cast<int>(k)-1){
	line& ln = intp[ctr++];
	if((ln.conns.size() == 1) && std::get<double>(ln.conns.front())==1){
	  hmap[i]=k;
	}
      }
    }
  }
  std::ofstream ofs(fn);
  for(auto r=0;r<Nr;++r){
    for(auto z=0;z<Nz;++z){
      int i=r*Nz+z;
      ofs<<hmap[i]<<" ";
    }
    ofs<<std::endl;
  }
}


//to transfer matrix from PIC code to test
//assumes one line in file for each matrix line. diagonal entries come first
Matrix read_input_matrix(std::string fn){
  Matrix inmatrix;
  std::ifstream ifs(fn);

  if(ifs.good()){
    std::string line;
    while(std::getline(ifs,line)){
      double val;
      int ind;
      std::vector<std::tuple<int,double,ConnType>> conns;
      std::istringstream iss(line);
      //read first two values
      iss >> ind;
      iss >> val;

      int conn_ind=-1;
      while(iss >> conn_ind){
	double conn_val=0;
	iss >> conn_val;
	conns.emplace_back(conn_ind,conn_val,WEAK);
      }
      inmatrix.emplace_back(ind,val,UNDECIDED,conns);
    }
  }
  return inmatrix;
}

//to transfer rhs from PIC code to test
//assumes linear array (one line)
std::vector<double> read_input_vector(std::string fn){
  std::vector<double> invec;
  std::ifstream ifs(fn);

  if(ifs.good()){
    std::string line;
    while(std::getline(ifs,line)){
      std::istringstream iss(line);
      double val=0;
      while(iss >> val){
	invec.push_back(val);
      }
    }
  }
  return invec;
}
