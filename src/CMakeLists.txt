set( SOURCES_AMG  solver.cpp matrix.cpp amg_level.cpp amg_output.cpp amg_setup.cpp
		  amg_solve_components.cpp amg_solvers.cpp amg_tests.cpp )

add_library(solver_amg SHARED ${SOURCES_AMG})

target_include_directories(solver_amg PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_include_directories(solver_amg PUBLIC "./h")
