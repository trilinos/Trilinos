#include <Kokkos_Core.hpp>

#include <iostream>

#ifdef KOKKOS_ENABLE_CUDA
  typedef Kokkos::Cuda      ExecSpace ;
  typedef Kokkos::CudaSpace  MemSpace ;
#elif defined(KOKKOS_ENABLE_OPENMP)
  typedef Kokkos::OpenMP    ExecSpace ;
  typedef Kokkos::OpenMP     MemSpace ;
#else
  typedef Kokkos::Serial    ExecSpace ;
  typedef Kokkos::HostSpace  MemSpace ;
#endif

typedef Kokkos::View<double**[2], Kokkos::LayoutRight, MemSpace> Mesh2D;

struct StructuredMeshKernel {
  Mesh2D inputmesh, outputmesh;

  Mesh2D::HostMirror host_inputmesh, host_outputmesh;

  int NumIntervalsInput, NumIntervalsOutput, RefineMultiplier, nrepeat;
  int NumNodesInput, NumNodesOutput;

  double input_meshsize, output_meshsize;

  double elapsed_time;

  KOKKOS_INLINE_FUNCTION void operator() (const int& i) const {
    for (int j = 0; j < NumIntervalsInput; ++j) {
      
      // refine a single grid cell
      for ( int i1 = 0; i1 <= RefineMultiplier ; ++i1 ) {
        for ( int j1 = 0; j1 <= RefineMultiplier ; ++j1 ) {
          outputmesh(RefineMultiplier * i + i1, RefineMultiplier * j + j1, 0) = 
            inputmesh(i, j, 0) + i1 * output_meshsize;
          outputmesh(RefineMultiplier * i + i1, RefineMultiplier * j + j1, 1) = 
            inputmesh(i, j, 1) + j1 * output_meshsize;
        }
      }
    }
  }

  void init (const int& n, const int &r) {
    NumIntervalsInput = n;
    RefineMultiplier  = r;
    nrepeat = 100;
    NumIntervalsOutput = RefineMultiplier * NumIntervalsInput;
    
    input_meshsize  = 1.0 / (double) NumIntervalsInput;
    output_meshsize = 1.0 / (double) NumIntervalsOutput;
    
    NumNodesInput  = 1 + NumIntervalsInput;
    NumNodesOutput = 1 + NumIntervalsOutput;
    
    inputmesh  = Mesh2D("inputmesh",  NumNodesInput,  NumNodesInput);
    outputmesh = Mesh2D("outputmesh", NumNodesOutput, NumNodesOutput);
    
    host_inputmesh  = Kokkos::create_mirror_view(inputmesh); 

    // Initialize input mesh on host
    for (int i = 0; i < NumNodesInput; ++i) {
      for (int j = 0; j < NumNodesInput; ++j) {
        host_inputmesh(i,j,0) = i * input_meshsize;
        host_inputmesh(i,j,1) = j * input_meshsize;
      }
    }

    host_outputmesh  = Kokkos::create_mirror_view(outputmesh); 

    Kokkos::deep_copy(inputmesh, host_inputmesh);
  }

  void run_lambda () {
    // Timer products
    struct timeval begin,end;
    
    gettimeofday(&begin,NULL);
    
    for ( int repeat = 0; repeat < nrepeat; repeat++) {
      Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,NumIntervalsInput), KOKKOS_LAMBDA (const int i) {
          for (int j = 0; j < NumIntervalsInput; ++j) {
            
            // refine a single grid cell
            for ( int i1 = 0; i1 <= RefineMultiplier ; ++i1 ) {
              for ( int j1 = 0; j1 <= RefineMultiplier ; ++j1 ) {
                outputmesh(RefineMultiplier * i + i1, RefineMultiplier * j + j1, 0) = 
                  inputmesh(i, j, 0) + i1 * output_meshsize;
                outputmesh(RefineMultiplier * i + i1, RefineMultiplier * j + j1, 1) = 
                  inputmesh(i, j, 1) + j1 * output_meshsize;
              }
            }
          }
        });
    }
    
    gettimeofday(&end,NULL);

    // Calculate time
    elapsed_time = 1.0*(end.tv_sec-begin.tv_sec) +
      1.0e-6*(end.tv_usec-begin.tv_usec);
  }

  void run_functor () {
    // Timer products
    struct timeval begin,end;
    
    gettimeofday(&begin,NULL);
    
    for ( int repeat = 0; repeat < nrepeat; repeat++) {
      Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,NumIntervalsInput), *this);
    }

    gettimeofday(&end,NULL);

    // Calculate time
    elapsed_time = 1.0*(end.tv_sec-begin.tv_sec) +
      1.0e-6*(end.tv_usec-begin.tv_usec);
  }

  void postprocess() {

    // Calculate bandwidth.
    double Gbytes = 1.0e-9 * double(sizeof(double) * ( NumNodesOutput * NumNodesOutput )) ;
    
    // Print results (problem size, time and bandwidth in GB/s)
    printf("  M( %d ) N( %d ) nrepeat ( %d ) problem( %g MB ) time( %g s ) bandwidth( %g GB/s )\n",
           NumNodesOutput , NumNodesInput, nrepeat, Gbytes * 1000, elapsed_time, Gbytes * nrepeat/ elapsed_time );

    Kokkos::deep_copy(host_outputmesh, outputmesh);

    for (int i = 0; i < NumIntervalsOutput; ++i) { // process rows
      for (int j = 0; j < NumIntervalsOutput; ++j) { // process cols
        double error = fabs( host_outputmesh( i, j, 0 ) - i * output_meshsize) +
                       fabs( host_outputmesh( i, j, 1 ) - j * output_meshsize);
        EXPECT_NEAR(error, 0.0, 1.e-12);
      }
    }
  }
  
};
