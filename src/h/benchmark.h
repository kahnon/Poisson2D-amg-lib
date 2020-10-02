#pragma once
#include <chrono>

namespace Benchmark{
  typedef decltype(std::chrono::high_resolution_clock::now()) BenchStart;

  auto start() {
    return std::chrono::high_resolution_clock::now();
  }  

  // return the duration
  auto stop(BenchStart begin) {
    auto time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - begin);
    return time.count();  // returns time in seconds as a double
  }  
}
