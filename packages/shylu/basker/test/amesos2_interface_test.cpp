#include "test_util.hpp"
#include "shylubasker_decl.hpp"
#include "shylubasker_def.hpp"
#include <Kokkos_Core.hpp>

int main(int argc, char* argv[])
{

  typedef long           Int;
  typedef double         Entry;
  typedef Kokkos::OpenMP Exe_Space;

  if(argc < 2)
    {
      std::cout <<"Inputs nthreads matrix.mtx (optional rhs.mtx)"
		<< std::endl;
      return -1;
    }
  
  Int nthreads = atoi(argv[1]);
  std::string mname = std::string(argv[2]);
  std::string vname;
  

  //Load inital information
  //Matrix
  Int m, n, nnz; 
  m = n = nnz = 0;
  Int *col_ptr, *row_idx;
  Entry *val;
 
  readMatrix(mname, m, n, nnz, &col_ptr, &row_idx, &val);
  
  //RHS
  Int vn, vm;
  vn = vm = 0;
  Entry *x, *xhat, *y;
  vn = n;
  x = new Entry[n]();
  xhat = new Entry[n]();
  if(argc == 4)
    {
      vname = std::string(argv[3]);
      readVector(vname, vm, y);
    }
  else
    {
      vm = m;
      y = new Entry[m]();
      for(Int i = 0; i < vm; i++)
	{
	  y[i] = (Entry) i;
	}
    }
  
  //Starting up Kokkos
  Exe_Space::initialize(nthreads);
  std::cout << "Kokkos Settings" << std::endl;
  std::cout << "hwloc aval: " 
	    << Kokkos::hwloc::available()<< std::endl;
  std::cout << "numa count: " 
	    << Kokkos::hwloc::get_available_numa_count() 
	    << std::endl;
  std::cout << "thrd numa:  " 
	    << Kokkos::hwloc::get_available_cores_per_numa() 
	    << std::endl;
 
  //Start Basker
  {
    int result = 0;
    BaskerNS::Basker<Int, Entry, Exe_Space> mybaker;
    //---Options
    mybasker.Options.same_pattern       = BASKER_FALSE;
    mybasker.Options.verbose            = BASKER_FALSE;
    mybasker.Options.verbose_matrix_out = BASKER_FALSE;
    mybasker.Options.realloc            = BASKER_TRUE;
    mybasker.Options.transpose          = BASKER_FALSE;
    mybasker.Options.symmetric          = BASKER_FALSE;
    mybasker.Options.AtA                = BASKER_TRUE;
    mybasker.Options.A_plu_At           = BASKER_TRUE;
    mybasker.Options.matching           = BASKER_TRUE;
    mybasker.Options.matching_type      = BASKER_MATCHING_BN;
    mybasker.Options.btf                = BASKER_TRUE;
    mybasker.Options.btf_max_percent    = BASKER_BTF_MAX_PERCENT;
    mybasker.Options.btf_large          = BASKER_BTF_LARGE;
    mybasker.Options.no_pivot           = BASKER_FALSE;
    //mybasker.Options.pivot_tol          = .001;
    //mybasker.Options.pivot_bias         = .001;
    //mybasker.Options.btf_prune_size      = 2;

    std::cout << "Setting Threads:" << nthreads << std::endl;
    mybasker.SetThreads(nthreads);
    mybasker.Symbolic(m,n,nnz,col_ptr,row_idx,val);
    mybasker.Factor(m,n,nnz,col_ptr,row_idx,val);
    mybasker.DEBUG_PRINT();
    mybasker.Solve(y,x);
    
    //mybasker.GetPerm()
    
    
    mybasker.Finalize();
  }
  
  Kokkos::finalize();

}//end main
