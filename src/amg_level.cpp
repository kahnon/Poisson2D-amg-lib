#include <algorithm>
#include <iostream>

#include "matrix.h"
#include "amg_level.h"
#include "amg_setup.h"
#include "amg_benchmark.h"
#include "amg_define.h"
#include "amg_output.h"


//memory consumption in MB
long int memory_consumption(AMGhierarchy& hier){
  long int mem=sizeof(hier);
  for(auto& lvl:hier) mem+=lvl.memory_consumption();
  return mem/(1024*1024);
}


//total number of gridpoints of the amg_hierarchy
long int number_of_points(AMGhierarchy& hier){
  long int ctr=0;
  for(auto& lvl:hier) ctr += lvl.matrix.size();
  return ctr;
}


//first coarsen matrix of finelvl
//then create interpolation matrix and assign it to the coarser lvl
//then create new system matrix for the new lvl and assign it as well
AMGlvl create_coarser_level(AMGlvl& finelvl, double thresh, double beta){
  //profiling setup times
  int id = AMGbench::start();
  coarsen_matrix(finelvl.matrix, thresh, beta);
  coarsening_time+=AMGbench::stop(id);

  id = AMGbench::start();
  Matrix new_interpolation = create_interpolation_matrix(finelvl.matrix,thresh);
  interpolation_creation_time+=AMGbench::stop(id);

  id = AMGbench::start();
  Matrix coarse_matrix = create_coarse_matrix(finelvl.matrix,new_interpolation);
  coarse_creation_time+=AMGbench::stop(id);

  return AMGlvl(coarse_matrix,new_interpolation,finelvl.level+1);
}


AMGhierarchy create_amg_hierarchy(Matrix& _matrix, int max_lvl, int minsize, double thresh, double beta){
  int id = AMGbench::start();
  AMGhierarchy hierarchy;

  //initialize finest level
  hierarchy.emplace_back(AMGlvl(_matrix,0));

  //initialize benchmark creation variables
  init_benchmark_variables();

  for(auto i=1;i<max_lvl && hierarchy.back().size()>minsize; ++i){
    std::cout<<"creating AMGlvl "<<i<<std::endl;
    AMGlvl& last_lvl = hierarchy[i-1];
    hierarchy.emplace_back(create_coarser_level(last_lvl,thresh,beta));
    print_matrix(hierarchy.back().matrix,"matrix_"+std::to_string(i)+".dat");
    print_matrix(hierarchy.back().interpolation,"interpolation_"+std::to_string(i)+".dat");
  }
  std::cout<<std::endl;

  //set maxlevel for each hierarchy, and init point lists
  std::for_each(hierarchy.begin(),hierarchy.end(),[&hierarchy](AMGlvl& lvl){
      lvl.max_level=hierarchy.size();
      //do not init lists of coarsest grid
      if(lvl.level != lvl.max_level-1) lvl.init_lists();
    }
  );

  std::cout<<"*********************AMG SETUP TIMES*********************"<<std::endl;
  std::cout<<"   Setup Phase took "<<AMGbench::stop(id)/1000.<<"s."<<std::endl;
  std::cout<<"      Coarsening          took "<<coarsening_time<<"ms."<<std::endl;
  std::cout<<"      Interpolation setup took "<<interpolation_creation_time<<"ms."<<std::endl;
  std::cout<<"      Coarse matrix setup took "<<coarse_creation_time<<"ms."<<std::endl;
  std::cout<<"*********************************************************"<<std::endl<<std::endl;

  return hierarchy;
}


