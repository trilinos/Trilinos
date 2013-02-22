// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include <string>
#include <iostream>
#include <cstdlib>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "TestStochastic.hpp"

#include "Stokhos_Host_CrsMatrix.hpp"
#include "Stokhos_Host_BlockCrsMatrix.hpp"
#include "Stokhos_Host_StochasticProductTensor.hpp"
#include "Stokhos_Host_CrsProductTensor.hpp"

// Algorithms
enum SG_Alg { ORIG_MAT_FREE, PROD_CRS };
const int num_sg_alg = 2;
const SG_Alg sg_alg_values[] = { ORIG_MAT_FREE, PROD_CRS };
const char *sg_alg_names[] = { "Original Matrix-Free", "Product CRS" };

std::vector<double>
run_test(const size_t num_cpu, const size_t num_core_per_cpu, 
	 const size_t num_threads_per_core, 
	 const size_t p, const size_t d, const size_t nGrid, const size_t nIter,
	 const bool symmetric, SG_Alg sg_alg, 
	 const std::vector<double>& perf1 = std::vector<double>())
{
  typedef double Scalar;
  typedef KokkosArray::Host Device;
  const size_t num_core = num_cpu * num_core_per_cpu;
  const size_t num_threads_per_cpu = num_core_per_cpu * num_threads_per_core;
  KokkosArray::Host::initialize( num_cpu , num_threads_per_cpu );

  std::vector<int> var_degree( d , p );

  std::vector<double> perf;
  if (sg_alg == PROD_CRS)
    perf = 
      unit_test::test_product_tensor_matrix<Scalar,Stokhos::CrsProductTensor<Scalar,Device>,Device>(var_degree , nGrid , nIter , symmetric );
  else if (sg_alg == ORIG_MAT_FREE)
    perf =
      unit_test::test_original_matrix_free_vec<Scalar,Device>( 
	var_degree , nGrid , nIter , true , symmetric );

  KokkosArray::Host::finalize();

  double speed_up;
  if (perf1.size() > 0)
    speed_up = perf1[1] / perf[1];
  else
    speed_up = perf[1] / perf[1];
  double efficiency = speed_up / num_core;

  std::cout << num_core << " , " 
	    << nGrid << " , "
	    << d << " , " 
	    << p << " , "
	    << perf[1] << " , "
	    << perf[2] << " , "
	    << speed_up << " , "
	    << 100.0 * efficiency << " , "
	    << std::endl;

  return perf;
}

int main(int argc, char *argv[])
{
  bool success = true;

  try {
    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    int p = 3;
    CLP.setOption("p", &p, "Polynomial order");
    int d = 4;
    CLP.setOption("d", &d, "Stochastic dimension");
    int nGrid = 64;
    CLP.setOption("n", &nGrid, "Number of spatial grid points in each dimension");
    int nIter = 1;
    CLP.setOption("niter", &nIter, "Number of iterations");
    int n_thread_per_core = 1;
    CLP.setOption("nthread", &n_thread_per_core, "Number of threads per core to use");
    int n_hyperthreads = 2;
    CLP.setOption("nht", &n_hyperthreads, "Number of hyperthreads per core available");
    SG_Alg sg_alg = PROD_CRS;
    CLP.setOption("alg", &sg_alg, num_sg_alg, sg_alg_values, sg_alg_names,
		  "SG Mat-Vec Algorithm");
    bool symmetric = true;
    CLP.setOption("symmetric", "asymmetric", &symmetric, "Use symmetric PDF");
    CLP.parse( argc, argv );
    
    // Detect number of CPUs and number of cores
    const size_t num_cpu  = 
      KokkosArray::Host::detect_gang_capacity();
    const size_t num_core_per_cpu = 
      KokkosArray::Host::detect_gang_worker_capacity() / 
      n_hyperthreads;
    
    // Print header
    std::cout << std::endl
	      << "\"#nCore\" , "
	      << "\"#nGrid\" , "
	      << "\"#Variable\" , "
	      << "\"PolyDegree\" , "
	      << "\"" << sg_alg_names[sg_alg] << " MXV Time\" , "
	      << "\"" << sg_alg_names[sg_alg] << " MXV GFLOPS\" , "
	      << "\"" << sg_alg_names[sg_alg] << " MXV Speedup\" , "
	      << "\"" << sg_alg_names[sg_alg] << " MXV Efficiency\" , "
	      << std::endl ;
    
    // Do a serial run to base speedup & efficiency from
    const std::vector<double> perf1 = 
      run_test(1, 1, 1, p, d, nGrid, nIter, symmetric, sg_alg);
    
    // First do 1 core per cpu
    for (size_t n=2; n<=num_cpu; ++n) {
      const std::vector<double> perf = 
	run_test(n, 1, 1, p, d, nGrid, nIter, symmetric, sg_alg, perf1);
    }
    
    // Now do all cpus, increasing number of cores
    for (size_t n=2; n<=num_core_per_cpu; ++n) {
      const std::vector<double> perf = 
	run_test(num_cpu, n, 1, p, d, nGrid, nIter, symmetric, sg_alg, perf1);
    }

    // Now do all cpus, all cores, with nthreads/core
    const std::vector<double> perf = 
      run_test(num_cpu, num_core_per_cpu, n_thread_per_core, p, d, nGrid, 
	       nIter, symmetric, sg_alg, perf1);
    

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
    
  if (!success)
    return -1;
  return 0 ;
}

