# cmake executable path: /tmp/cmake-3.30.8-linux-x86_64/bin/cmake

set( CONFIG_NAME "ubuntu-lsl" )

# Set compilers path
set(CMAKE_C_COMPILER "/usr/bin/gcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE PATH "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

# Set paths to mpi
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/bin/mpicxx" CACHE PATH "")
set(MPIEXEC "/usr/bin/mpirun" CACHE PATH "")

# Set paths to blas and lapack
set( BLAS_LIBRARIES "/usr/lib/x86_64-linux-gnu/libblas.so" CACHE PATH "" FORCE )
set( LAPACK_LIBRARIES "/usr/lib/x86_64-linux-gnu/liblapack.so" CACHE PATH "" FORCE )

# Cuda and openMP
set( ENABLE_CUDA OFF CACHE PATH "" FORCE )
set( ENABLE_OPENMP OFF CACHE PATH "" FORCE )

# TPLs
set( ENABLE_TRILINOS OFF CACHE PATH "" FORCE )
set( ENABLE_HYPRE ON CACHE BOOL "" FORCE )
set( ENABLE_CALIPER OFF CACHE BOOL "" FORCE )
set( ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE)
set( ENABLE_MATHPRESSO OFF CACHE BOOL "" FORCE )
set( ENABLE_SCOTCH ON CACHE BOOL "" FORCE )
set( ENABLE_WARNINGS_AS_ERRORS OFF CACHE BOOL "" FORCE )

if(NOT ( EXISTS "${GEOS_TPL_DIR}" AND IS_DIRECTORY "${GEOS_TPL_DIR}" ) )
   set(GEOS_TPL_DIR "${CMAKE_SOURCE_DIR}/../../thirdPartyLibs/install-${CONFIG_NAME}-release" CACHE PATH "" FORCE )
endif()

include(${CMAKE_CURRENT_LIST_DIR}/tpls.cmake)
