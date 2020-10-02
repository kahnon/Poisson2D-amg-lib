#pragma once

#include <vector>
#include <tuple>
#include <iostream>
#include "matrix.h"

//TODO find more elegant way of constructing hierarchy other than copying (not rly critical for runtime)
struct AMGlvl{
  Matrix matrix;
  Matrix interpolation;
  int level;
  int max_level;
  std::vector<double> phi;
  std::vector<double> rhs;

  //lists to make iterations over only coarse/fine points and boundary points easier
  std::vector<int> coarse_list;
  std::vector<int> fine_list;
  std::vector<int> boundary_list;

  //add Matrix&& constructor so that initially calculated matrix doesnt have to be deleted
  AMGlvl(Matrix& _matrix, int lvl):level(lvl){
    matrix = _matrix;
    interpolation = Matrix();
    phi.resize(matrix.size(),0);
    rhs.resize(matrix.size(),0);
  }

  AMGlvl(Matrix& _matrix, Matrix& _interpolation, int lvl):AMGlvl(_matrix,lvl){
    interpolation = _interpolation;
  }

  AMGlvl():level(-1),max_level(-1){
  }

  int size(){
    return matrix.size();
  }

  //in bytes
  long int memory_consumption(){
    long int mem=0;
    for(auto& ln:matrix) mem+=ln.memory_consumption();
    for(auto& ln : interpolation) mem += ln.memory_consumption();
    mem += 2 * sizeof(int);
    mem += sizeof(phi) + sizeof(rhs);
    mem += phi.capacity()*sizeof(double);
    mem += rhs.capacity()*sizeof(double);
    mem += coarse_list.capacity()*sizeof(int);
    mem += fine_list.capacity()*sizeof(int);
    mem += boundary_list.capacity()*sizeof(int);
    return mem;
  }

  //for each point the sum of inteprolation weights should be 1
  //not sure if useful at all
  double interpolation_check(){
    if(!interpolation.empty()){
      double error=0;
      double sum=0;

      for(auto i=0u;i<interpolation.size();++i){
	auto& ln = interpolation[i];
	sum=0;
	for(auto& cn:ln.conns) sum+=std::get<double>(cn);
	if(sum!=0) sum=1-sum;
	//std::cout<<"i="<<i<<" sum="<<sum<<std::endl;
	if(sum>error) error=sum;
      }
      return sum;
    }else return 0;
  }

  void init_lists(){
    for(auto& ln:matrix){
      if(ln.status == UNDECIDED){
	std::cout<<"lists of AMGlvl "<<level<<" cannot be initialized. Matrix not yet coarsened"<<std::endl;
	break;
      }else if(ln.status == COARSE) coarse_list.push_back(ln.ind);
      else if(ln.status == SFINE || ln.status==WFINE) fine_list.push_back(ln.ind);
      else if(ln.status == BOUNDARY) boundary_list.push_back(ln.ind);
    }
  }
};


typedef std::vector<AMGlvl> AMGhierarchy;

long int memory_consumption(AMGhierarchy& hier);
long int number_of_points(AMGhierarchy& hier);

//TODO decide where to set thresh/beta as default values
AMGlvl create_coarser_level(AMGlvl& finelvl, double thresh=0.25, double beta=0.35);
AMGhierarchy create_amg_hierarchy(Matrix& _matrix, int max_lvl, int minsize, double thresh=0.25, double beta=0.35);
