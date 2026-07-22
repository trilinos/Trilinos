# Rapid Optimization Library (ROL)

![Rapid Optimization Library](https://raw.githubusercontent.com/sandialabs/rol/refs/heads/develop/rol.png)

**ROL** (as in rock and _roll_) is a high-performance C++ library for numerical optimization.
ROL brings an extensive collection of state-of-the-art optimization algorithms to virtually
any application. Its programming interface supports any computational hardware, including
heterogeneous many-core systems with digital and analog accelerators. ROL has been used with
great success for optimal control, optimal design, inverse problems, image processing and
mesh optimization, in application areas including geophysics, structural dynamics, fluid
dynamics, electromagnetics, quantum computing, hypersonics and geospatial imaging.

For additional details, see [https://rol.sandia.gov](https://rol.sandia.gov).

Feature highlights:

1. Vector abstractions and matrix-free interface for universal applicability
2. Modern algorithms for unconstrained and constrained optimization
3. Easy-to-use methods for stochastic and risk-aware optimization
4. Fast and robust algorithms for nonsmooth optimization
5. Trust-region methods for inexact and adaptive computations
6. PDE-OPT application development kit for PDE-constrained optimization
7. Interfaces and algorithms for optimal experimental design


## Getting Started
ROL is a C++ library with a cmake build system.
There are minimal third-party library requirements consisting of BLAS and LAPACK
which also require a suitable Fortran compiler on the system.
On Linux-based systems with an apt package manager,
the dependencies may be installed with
```
apt install cmake gfortran libopenblas-dev liblapack3
```

Some typical configure scripts may be found in the `.github/workflows/` directory.
An in-source release build including ROL's examples and tests may be configured
from within the source directory with
```
cmake -D CMAKE_BUILD_TYPE:STRING=RELEASE \
      -D ENABLE_EXAMPLES:BOOL=ON \
      -D ENABLE_TESTS:BOOL=ON \
      -B build \
      -D CMAKE_INSTALL_PREFIX=build/install \
      .
```

After a successful configure, ROL is built by then changing
to the build directory and running `make`
```
cd build
make
```

One can then install it 
```
make install
```

## Building against ROL
Assuming we have a directory called 'roltest' which contains files we would like to run using ROL. 
In /dir/to/roltest, create a CMakeLists.txt with
the following: 

```
cmake_minimum_required(VERSION 3.10)  
project(test)
set(CMAKE_CXX_STANDARD 17) 

enable_language(C)
enable_language(CXX)

message(STATUS "ROL directory is: ${ROL_DIR}/include")
include_directories("${ROL_DIR}/include") #one down for all files
find_library(ROL_LIB NAMES rol librol PATHS ${ROL_DIR}/lib NO_DEFAULT_PATH)
if(ROL_LIB)
  message(STATUS "We found ROL!")
else()
  message(FATAL_ERROR "ROL library not found!")
endif()

add_library(rol::rol UNKNOWN IMPORTED)
set_target_properties(rol::rol PROPERTIES IMPORTED_LOCATION "${ROL_LIB}")
add_executable(ROLtest mytest.cpp)
target_link_libraries(ROLtest PRIVATE rol::rol) 

add_executable(ROLtest mytest.cpp)
```
The file in mytest then simply tests that it can find a 
basic ROL type:  
```
#include "ROL_StdVector.hpp"
typedef double RealT;

int main() {
  int dim = 100; 

  ROL::StdVector<RealT> x(dim);
  x.randomize(-1.0, 1.0);
  std::cout << "What's the norm of x? " << x.norm() << std::endl;
}
```
Compiling and executing should then be a matter of: 
```
$ cmake . -DROL_DIR=/Users/../path/to/rol/build/install
$ make
$ ./ROLtest
What's the norm of x? "some random number"
```


## Python Interface
The python pip package rol-python is the official python interface to ROL.
It is available at https://pypi.org/project/rol-python/ and can be installed via
```
pip install rol-python
```


## Copyright and License
See COPYRIGHT and LICENSE.


## Questions?
Contact team or developers:

* ROL Team     (GitHub handle: @sandialabs/rol)
* Drew Kouri   (GitHub handle: [dpkouri](https://github.com/dpkouri) or dpkouri@sandia.gov)
* Denis Ridzal (GitHub handle: [dridzal](https://github.com/dridzal) or dridzal@sandia.gov)

