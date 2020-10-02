#define DEFINE_AMG
#include <iostream>
#include <tuple>
#include <vector>

#include "amg_benchmark.h"
#include "amg_level.h"
#include "amg_output.h"
#include "amg_solvers.h"
#include "amg_tests.h"
#include "benchmark.h"
#include "matrix.h"
#include "solver.h"


AMGhierarchy amg;

//Eps matrix must be provided by implementation
//it contains the boundary conditions of my domain:
//-100: phi = const. (also within the domain)
//-200: grad(phi) = const. (only at domain boundaries)
//1 = vacuum
//other values>0: relative permittivity
void init_solver(const std::vector<double>& Eps, int NZ, int NR){
  //setup routine not mpi_parallel

  //size parameters
  int NG = NZ*NR;

  assert(Eps.size() == static_cast<size_t>(NG) && "Eps matrix has wrong dimensions!");

  //time measurement (independent of amg profiling)
  double factorize_time = 0.;
  clock_t t_factorize0 = clock();

  // local crm 
  Matrix matrix;

  // count of non null entries
  for (int i = 0; i < NR; ++i)
    for (int j = 0; j < NZ; ++j) {
      int i_row = j + i * NZ;
      // Fixed potential on metal
      if (Eps[j + i * NZ] == -100.){
	//bingo
	matrix.emplace_back(i_row,1);

      // Fixed gradient
      } else if (Eps[j + i * NZ] == -200.){
	std::vector<std::tuple<int,double,ConnType>> conns;
	//west
	conns.push_back(std::make_tuple(i_row-1,-1,WEAK));
	//bingo
	matrix.emplace_back(i_row,1,UNDECIDED,conns);
	
      // Axis
      } else if (i == 0){
	double e_iP1_jM1 = Eps[(i + 1) * NZ + j - 1];
	if (e_iP1_jM1 <= 0.) e_iP1_jM1 = Eps[i * NZ + j];

	double e_iP1_jP1 = Eps[(i + 1) * NZ + j + 1];
	if (e_iP1_jP1 <= 0.) e_iP1_jP1 = Eps[i * NZ + j];

	double west = -1*e_iP1_jM1;  // a,b,c,d,e    *= dz*dz
	double east = -1*e_iP1_jP1;
	double north = -2. * (e_iP1_jM1 + e_iP1_jP1);
	double bingo = 3. * (e_iP1_jM1 + e_iP1_jP1);

	std::vector<std::tuple<int,double,ConnType>> conns;
	//west
	conns.push_back(std::make_tuple(i_row-1,west,WEAK));
	//east
	conns.push_back(std::make_tuple(i_row+1,east,WEAK));
	//north
	conns.push_back(std::make_tuple(i_row+NZ,north,WEAK));
	//bingo
	matrix.emplace_back(i_row,bingo,UNDECIDED,conns);

      // the rest: channell, dielectric, boundaries, etc
      }else{
	double e_iP1_jM1 = Eps[(i + 1) * NZ + j - 1];
	if (e_iP1_jM1 <= 0.) e_iP1_jM1 = Eps[i * NZ + j];
	if (e_iP1_jM1 == 0) {
	  if (Eps[(i + 1) * NZ + j] > 0)
	    e_iP1_jM1 = Eps[(i + 1) * NZ + j];
	  else
	    e_iP1_jM1 = Eps[i * NZ + j - 1];
	}

	double e_iP1_jP1 = Eps[(i + 1) * NZ + j + 1];
	if (e_iP1_jP1 <= 0.) e_iP1_jP1 = Eps[i * NZ + j];
	if (e_iP1_jP1 == 0) {
	  if (Eps[(i + 1) * NZ + j] > 0)
	    e_iP1_jP1 = Eps[(i + 1) * NZ + j];
	  else
	    e_iP1_jP1 = Eps[i * NZ + j + 1];
	}

	double e_iM1_jM1 = Eps[(i - 1) * NZ + j - 1];
	if (e_iM1_jM1 <= 0.) e_iM1_jM1 = Eps[i * NZ + j];
	if (e_iM1_jM1 == 0) {
	  if (Eps[(i - 1) * NZ + j] > 0)
	    e_iM1_jM1 = Eps[(i - 1) * NZ + j];
	  else
	    e_iM1_jM1 = Eps[i * NZ + j - 1];
	}

	double e_iM1_jP1 = Eps[(i - 1) * NZ + j + 1];
	if (e_iM1_jP1 <= 0.) e_iM1_jP1 = Eps[i * NZ + j];
	if (e_iM1_jP1 == 0) {
	  if (Eps[(i - 1) * NZ + j] > 0)
	    e_iM1_jP1 = Eps[(i - 1) * NZ + j];
	  else
	    e_iM1_jP1 = Eps[i * NZ + j + 1];
	}

	double north = -1*(0.25 * (e_iP1_jM1 + e_iP1_jP1) / i + 0.5 * (e_iP1_jM1 + e_iP1_jP1));
	double south = -1*(0.5 * (e_iM1_jM1 + e_iM1_jP1) - 0.25 * (e_iM1_jM1 + e_iM1_jP1) / i);
	double bingo = -1*(0.25 * (e_iM1_jM1 + e_iM1_jP1 - e_iP1_jM1 - e_iP1_jP1) / i -
		(e_iP1_jM1 + e_iP1_jP1 + e_iM1_jM1 + e_iM1_jP1));
	double west = -1*0.5 * (e_iP1_jM1 + e_iM1_jM1);
	double east = -1*0.5 * (e_iP1_jP1 + e_iM1_jP1);

	std::vector<std::tuple<int,double,ConnType>> conns;
	//west
	conns.push_back(std::make_tuple(i_row-1,west,WEAK));
	//east
	conns.push_back(std::make_tuple(i_row+1,east,WEAK));
	//north
	conns.push_back(std::make_tuple(i_row+NZ,north,WEAK));
	//south
	conns.push_back(std::make_tuple(i_row-NZ,south,WEAK));
	//bingo
	matrix.emplace_back(i_row,bingo,UNDECIDED,conns);
      }//j
    }//i
#if DEBUG_MATRIX_CREATION
  print_matrix(matrix,"amg_system_matrix.dat");
#endif

  //actual creation of amg hierarchy
  //magic numbers = amg initialization parameters
  //good set for parameters: create_amg_hierarchy(matrix,5,2048,0,0.35)
  amg = create_amg_hierarchy(matrix,16,64,0.25,0.35);
  print_matrix(matrix,"amg_system_matrix.dat");
  print_hierarchy(amg,NR,NZ,"amg_hierarchy.dat");

  factorize_time += ((double)(clock() - t_factorize0)) / CLOCKS_PER_SEC;

#if 0
  {
    std::vector<double> test_vec(amg[0].matrix.size(),1);
    continued_restriction(amg,test_vec,0,amg.size()-1);
    print_vector(test_vec,test_vec.size(),1,"vector_restricted.dat");
    continued_interpolation(amg,test_vec,amg.size()-1,0);
    print_vector(test_vec,test_vec.size(),1,"vector_interpolated.dat");
    print_vector(test_vec,NR,NZ,"vector_interpolated_formatted.dat");
    
  }
#endif
}

void destroy_solver(){
  if(!amg.empty()) amg.resize(0);
}

//phi contains the solution, rhs is rhs of Poisson equation
//must be supplied with correct dimensions by the user
//init_solver must be called before computation of the solution
void calculate_potential( std::vector<double>& phi, std::vector<double>& rhs, double& field_time ) 
{
  assert(phi.size() == rhs.size() && "Vector dimensions don't match!");

  auto start_time = Benchmark::start();
  solve_v_cycle(amg,rhs,phi,100,1,0,1e-6);
  field_time += Benchmark::stop(start_time);
}
