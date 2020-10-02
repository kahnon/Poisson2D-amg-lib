#include <vector>
#include <limits>
#include <math.h>

#include "amg_solvers.h"
#include "smoothing.h"
#include "amg_solve_components.h"
#include "matrix.h"
#include "amg_benchmark.h"


//sor solver with system matrix "matrix", rhs "rhs", damping parameter omega
//and guess phi. Iteration until max_norm of residual is smaller than maxerror
//or maxstep is reached
void sor(const Matrix& matrix, const std::vector<double>& rhs, std::vector<double>& phi, int error_interval, double maxerror, int maxstep, double omega){
  assert(matrix.size() == rhs.size() && "matrix & vector size dont match in sor(...)");
  assert(matrix.size() == phi.size() && "matrix & vector size dont match in sor(...)");
  int steps=0;
  double error=1;
  std::vector<double> phi_old(0);
  do{
    //checkerboard pattern
    for(auto i=0u;i<rhs.size();i+=2) sor_base_with_bounds(i,matrix,rhs,phi_old,phi,omega);
    for(auto i=1u;i<rhs.size();i+=2) sor_base_with_bounds(i,matrix,rhs,phi_old,phi,omega);
    ++steps;
    if(steps%error_interval == 0) error = calc_max_norm_res(matrix,rhs,phi).second;
  }while(error > maxerror && steps<maxstep);
  std::cout<<"sor: error="<<error<<"  steps="<<steps<<std::endl;
}


//jacobi solver with system matrix "matrix", rhs "rhs", damping parameter omega
//and guess phi. Iteration until max_norm of residual is smaller than maxerror
//or maxstep is reached
void jacobi(const Matrix& matrix, const std::vector<double>& rhs, std::vector<double>& phi, int error_interval, double maxerror, int maxstep, double omega){
  assert(matrix.size() == rhs.size() && "matrix & vector size dont match in jacobi(...)");
  assert(matrix.size() == phi.size() && "matrix & vector size dont match in jacobi(...)");
  int steps=0;
  double error=1;
  std::vector<double> phi_old(phi);

  do{
    for(auto i=0u;i<rhs.size();++i) jacobi_base_with_bounds(i,matrix,rhs,phi_old,phi,1);
    ++steps;
    if(steps%error_interval == 0) error = calc_max_norm_res(matrix,rhs,phi).second;
    std::copy(phi.begin(),phi.end(),phi_old.begin());
  }while(error > maxerror && steps<maxstep);
  std::cout<<"jacobi: error="<<error<<"  steps="<<steps<<std::endl;
}


//non pre-conditioned cg method to solve the AMG system on the lowest level
//(should be) guaranteed to converge after dim(phi) steps
void cg(const Matrix& matrix, const std::vector<double>& rhs, std::vector<double>& phi, int maxsteps, double acc){
  assert(matrix.size() == rhs.size() && "matrix & vector size dont match in cg(...)");
  assert(matrix.size() == phi.size() && "matrix & vector size dont match in cg(...)");

  double error=0;
  int steps=0;

  std::vector<double> residual=matrix_vector(matrix,phi);
  for(auto i=0u;i<residual.size();++i) residual[i] = rhs[i] - residual[i];

  std::vector<double> direction=residual;
  std::vector<double> mvprod(direction.size(),0);
  std::vector<double> new_residual(residual);

  for(auto i=0;i<maxsteps;++i){
    error=0;
    ++steps;

    //save this mvproduct
    mvprod = matrix_vector(matrix,direction);

    //find phi in search direction and update gradient and residual
    double alpha = scalar_product(residual,residual);
    alpha /= scalar_product(direction,mvprod);
    for(auto k=0u;k<phi.size();++k) phi[k]+=alpha*direction[k];
    for(auto k=0u;k<residual.size();++k){ 
      new_residual[k] = residual[k] - alpha*mvprod[k];
      //find maximum norm of residual
      if(fabs(new_residual[k])>error) error=fabs(new_residual[k]);
    }
    if(error<acc) break;

    //update search direction
    double beta = scalar_product(new_residual,new_residual);
    beta /= scalar_product(residual,residual);
    for(auto k=0u;k<direction.size();++k) direction[k] = new_residual[k] + beta*direction[k];
    std::copy(new_residual.begin(),new_residual.end(),residual.begin());
  }
  //std::cout<<"   cg: error="<<error<<"  steps="<<steps<<std::endl;
}


//v-cycle: move straight down to coarsest level, 
//then back up and repeat until convergence
void solve_v_cycle(AMGhierarchy& Hier, const std::vector<double>& rho, std::vector<double>& phi, int maxsteps, 
		   int smooth_steps, bool gen_flag, double maxerror){
  init_benchmark_variables();
  int solve_id = AMGbench::start();

  int max_level=Hier[0].max_level;
  //int max_level=2;
  int steps=0;
  std::vector<int> stepsvec(max_level,smooth_steps);
  //if gen flag=true, use generalized cycle
  if(gen_flag){
    for(auto i=0;i<max_level;++i){
      stepsvec[i] = (i+1) * smooth_steps;
    }
  }

  AMGlvl& coarsest_lvl = Hier[max_level-1];

  //init first guess for each level
  Hier[0].phi = phi;
  Hier[0].rhs = rho;
  for(auto i=1;i<max_level;++i){
    std::fill(Hier[i].phi.begin(),Hier[i].phi.end(),0);
    std::fill(Hier[i].rhs.begin(),Hier[i].rhs.end(),0);
  }

  auto v_iteration = [&](){
    double error = 0;
    //go "down" from i to next coarser lvl all the way to coarsest level
    for(auto i=0;i<max_level-1;++i){ 
      int id = AMGbench::start();
      double cycerror = go_down(i,Hier,maxerror,stepsvec[i]);
      go_down_time+=AMGbench::stop(id);
      
      id=AMGbench::start();
      if(i==0){
	error=cycerror;
	if(error < maxerror) return error;
      }
      error_time+=AMGbench::stop(id);
    }

    //almost exact solution on coarsest grid
    int id = AMGbench::start();
    //cg(coarsest_lvl.matrix,coarsest_lvl.rhs,coarsest_lvl.phi,std::numeric_limits<int>::max());
    cg(coarsest_lvl.matrix,coarsest_lvl.rhs,coarsest_lvl.phi,500);
    //sor(coarsest_lvl.matrix, coarsest_lvl.rhs, coarsest_lvl.phi, 10, 1e-12, 500, 1);
    cg_time+=AMGbench::stop(id);

    //go back "up" again from i t0 next finer lvl
    id=AMGbench::start();
    for(auto i=max_level-1;i>0;--i) go_up(i,Hier,maxerror,stepsvec[i]);
    go_up_time+=AMGbench::stop(id);
    return error;
  };

  double error=1;
  do{
    error = v_iteration();
    ++steps;
    std::cout<<"v-cycle step "<<steps<<", error "<<error<<std::endl;
  }while(error>maxerror && steps<maxsteps);

  phi=Hier[0].phi;
  //std::cout<<"AMG V-cycle: cycles="<<steps<<"  "<<"error="<<error<<std::endl;
  
  long int pts = number_of_points(Hier);

  std::cout<<"*********************AMG V-CYCLE TIMES*********************"<<std::endl;
  std::cout<<"   Solve took "<<AMGbench::stop(solve_id)/1000.<<"s with "<<steps<<" cycles and error="<<error<<std::endl;
  std::cout<<"      Go Down      took "<<1e6*go_down_time/(steps*pts)<<" ps/(cyc*pt)."<<std::endl;
  std::cout<<"         Pre-Smoothing        took "<<1e6*pre_smoothing_time/(steps*pts)<<" ps/(cyc*pt)."<<std::endl;
  std::cout<<"         Residual Calculation took "<<1e6*residual_time/(steps*pts)<<" ps/(cyc*pt)."<<std::endl;
  std::cout<<"         Restriction          took "<<1e6*restriction_time/(steps*pts)<<" ps/(cyc*pt)."<<std::endl<<std::endl;
  std::cout<<"      CG solve     took "<<1e6*cg_time/(steps*pts)<<" ps/(cyc.*pt)"<<std::endl<<std::endl;
  std::cout<<"      Go Up        took "<<1e6*go_up_time/(steps*pts)<<" ps/(cyc*pt)."<<std::endl;
  std::cout<<"         Post-Smoothing       took "<<1e6*post_smoothing_time/(steps*pts)<<" ps/(cyc*pt)."<<std::endl;
  std::cout<<"         Jacobi Sweeps        took "<<1e6*jacobi_sweep_time/(steps*pts)<<" ps/(cyc*pt)."<<std::endl;
  std::cout<<"         Interpolation        took "<<1e6*interpolation_time/(steps*pts)<<" ps/(cyc*pt)."<<std::endl;
  std::cout<<"      Error Checks took "<<1e6*error_time/(steps*pts)<<" ps/(cyc*pt)."<<std::endl;
  std::cout<<"***********************************************************"<<std::endl<<std::endl;
}



//TODO pass parameters of amg preconditioning
void amg_preconditioned_cg(AMGhierarchy& hier, const std::vector<double>& rhs, std::vector<double>& phi, 
			   std::function<void(AMGhierarchy&,const std::vector<double>&,std::vector<double>&,int,int,bool,double)> precondition, 
			   bool gen_flag, int maxsteps, double acc){

  //system matrix the same for both systems
  Matrix& matrix = hier[0].matrix;
  assert(matrix.size() == rhs.size() && "matrix & vector size dont match in cg(...)");
  assert(matrix.size() == phi.size() && "matrix & vector size dont match in cg(...)");

  double error=0;
  int steps=0;

  //move out of here and cg() from above
  auto scalar_product = [](const std::vector<double>& v1, const std::vector<double>& v2){
    assert(v1.size() == v2.size() && "vector dimensions dont match in scalar_product");
    double result=0;
    for(auto i=0u;i<v1.size();++i) result+=v1[i]*v2[i];
    return result;
  };

  //needed to set residual as rhs for preconditioner
  std::vector<double> residual=matrix_vector(matrix,phi);
  for(auto i=0u;i<residual.size();++i) residual[i] = rhs[i] - residual[i];

  //preconditioned variable used here
  std::vector<double> prec_residual(residual);
  precondition(hier,residual,prec_residual,30,2,gen_flag,1e-12);
  //search direction in preconditioned residual
  std::vector<double> direction(prec_residual);
  std::vector<double> mvprod(direction.size(),0);

  double delta = scalar_product(residual,prec_residual);
  double new_delta = delta;


  for(auto i=0;i<maxsteps;++i){
    error=0;
    ++steps;

    //save this mvproduct
    mvprod = matrix_vector(matrix,direction);

    //find phi in search direction and update residual
    double alpha = delta;
    alpha /= scalar_product(direction,mvprod);

    for(auto k=0u;k<phi.size();++k) phi[k]+=alpha*direction[k];
    for(auto k=0u;k<residual.size();++k){ 
      residual[k] -= alpha*mvprod[k];
      if(fabs(residual[k]) > error) error = fabs(residual[k]);
    }
    if(error < acc) break;

    //use preconditioning
    precondition(hier,residual,prec_residual,30,2,gen_flag,1e-12);
    
    //update coefficients for new search direction
    new_delta = scalar_product(residual,prec_residual);
    double beta = new_delta / delta;
    delta = new_delta;

    //std::cout<<"step="<<steps<<"  error="<<error<<"  alpha="<<alpha<<"  beta="<<beta<<"  delta="<<delta<<std::endl;

    //update search direction
    for(auto k=0u;k<direction.size();++k) direction[k] = prec_residual[k] + beta*direction[k];
  }

  std::cout<<"amg_preconditioned_cg: error="<<error<<"  steps="<<steps<<std::endl;
}
