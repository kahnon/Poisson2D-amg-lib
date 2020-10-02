#pragma once

/*
 * =====================================================================================
 *
 *       Filename:  Benchmark.hpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  22/10/12 21:19:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Stefan Kemnitz (), kemnitz.stefan@googlemail.com
 *        Company:  
 *
 * =====================================================================================
 */

#include <iostream>
#include <sys/time.h>
#include <ctime>
#include <chrono>
#include <string>
#include <map>
#include "amg_define.h"

//profiling variables
XTRN_AMG double coarsening_time, interpolation_creation_time, coarse_creation_time;
XTRN_AMG double interpolation_time, restriction_time, residual_time, pre_smoothing_time, post_smoothing_time, jacobi_sweep_time; 
XTRN_AMG double go_up_time, go_down_time, cg_time, error_time;

inline void init_benchmark_variables(){
  coarsening_time=0; interpolation_creation_time=0; coarse_creation_time=0;
  interpolation_time=0; restriction_time=0; residual_time=0; pre_smoothing_time=0; post_smoothing_time=0, jacobi_sweep_time=0; 
  go_up_time=0; go_down_time=0; cg_time=0; error_time=0;
}

namespace AMGbench{

    XTRN_AMG std::map<int,decltype(std::chrono::high_resolution_clock::now())> begins;
    /// 
    /// @brief stores the actual time on call 
    /// @returns inline void
    ///
    inline int start() {
	auto begin = std::chrono::high_resolution_clock::now();
	static int ctr = 0;
	int id = ctr++;
	begins[id] = begin;
	return id;
    } // end of function start

    inline long stop( int id ) {
	auto begin = begins[id];
	auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin);
	begins.erase( id );
	return time.count() / 1000.0 / 1000.0;
    } // end of function stop

}
