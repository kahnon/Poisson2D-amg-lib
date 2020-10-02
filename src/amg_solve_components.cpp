#include <vector>
#include <iostream>
#include <math.h>

#include "matrix.h"
#include "amg_benchmark.h"
#include "amg_solve_components.h"
#include "smoothing.h"


//phi_vec and rhs_vec need to be initialized properly in amg cycle routine
//move "down" from finer to coarser level in the amg hierarchy
double go_down(int i, AMGhierarchy& Hier, double maxerror, int maxsteps){
  AMGlvl& fine_lvl = Hier[i];
  AMGlvl& coarse_lvl = Hier[i+1];
  double error=1;
  if(i==0) error=0;

  assert(i < fine_lvl.max_level-1 && "in AMG::go_down the index was chosen too high");

  //make it better readable
  std::vector<double>& fine_phi = fine_lvl.phi;
  const std::vector<double>& fine_rhs = fine_lvl.rhs;
  std::vector<double>& coarse_rhs = coarse_lvl.rhs;

  //maxstep-times smoothing iteration
  //TODO make this smoothing an inlined function and pass here
  int id = AMGbench::start();
  smoothing(fine_lvl,1);
  pre_smoothing_time+=AMGbench::stop(id);
  
  //calculate residual on current level
  id = AMGbench::start();
  std::vector<double> residual = matrix_vector(fine_lvl.matrix,fine_phi);
  for(auto i=0u;i<residual.size();++i) residual[i] = fine_rhs[i] - residual[i];
  residual_time+=AMGbench::stop(id);

  //on finest level calc max residual
  id = AMGbench::start();
  if(i==0){
    for(auto i=0u; i<residual.size();++i){
      double locerror = residual[i]>=0 ? residual[i] : -1*residual[i];
      if(locerror > error) error = locerror;
    }
  }
  error_time+=AMGbench::stop(id);
  
  //restrict residual to next level
  id = AMGbench::start();
  coarse_rhs = restriction(coarse_lvl,residual);
  restriction_time+=AMGbench::stop(id);
  return error;
}

//move "up" one level in the amg hierarchy
void go_up(int i, AMGhierarchy& Hier, double maxerror, int maxsteps){
  AMGlvl& coarse_lvl = Hier[i];
  AMGlvl& fine_lvl = Hier[i-1];

  assert(i > 0 && "in AMG::go_up the index was chosen to be 0");

  std::vector<double>& coarse_phi = coarse_lvl.phi;
  std::vector<double>& fine_phi = fine_lvl.phi;
  std::vector<double> fine_phi_old;
  
  //interpolation to finer lvl
  int id = AMGbench::start();
  std::vector<double> phi_correction = interpolation(coarse_lvl.interpolation,coarse_phi);
  interpolation_time += AMGbench::stop(id);

  //jacobi sweep after interpolation improves convergence
  //TODO pass number of jacobi sweeps as parameter
  id = AMGbench::start();
#if 0
  fine_phi_old = phi_correction;
  for(auto k=0;k<maxsteps;++k){
    for(auto& i: fine_lvl.fine_list) jacobi_base(i,fine_lvl.matrix,fine_phi_old,fine_phi_old,phi_correction,1);
    for(auto& i: fine_lvl.fine_list) fine_phi_old[i]=phi_correction[i];
  }
#endif
  jacobi_sweep_time+=AMGbench::stop(id);
  
  //add coarse grid correction to current iterate
  for(auto i=0u;i<fine_phi.size();++i) fine_phi[i]+=phi_correction[i];

  //if we are not on finest level, assign new first guess for next iteration
  std::fill(coarse_phi.begin(),coarse_phi.end(),0);

  //TODO make this smoothing an inlined function and pass here
  id = AMGbench::start();
  smoothing(fine_lvl,1);
  post_smoothing_time+=AMGbench::stop(id);
}


