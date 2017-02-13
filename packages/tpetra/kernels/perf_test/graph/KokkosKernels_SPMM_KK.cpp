/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels: Linear Algebra and Graph Kernels
//                 Copyright 2016 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#include <iostream>
#include "KokkosKernels_SPGEMM.hpp"
//#include <Kokkos_Sparse_CrsMatrix.hpp>
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_GraphColor.hpp"
#include "KokkosKernels_IOUtils.hpp"

typedef int size_type;
//typedef size_t size_type;
typedef int idx;
typedef double wt;

int mkl_sort_option = 7;
int mkl_keep_output = 1;
int check_output = 0;
#define TRANPOSEFIRST false
#define TRANPOSESECOND false
int vector_size = -1;
bool calculate_read_write_cost = false;

enum MEMSPACE{HBM, DDR4}; //GPUS GPU vs DDR4
MEMSPACE amemspace = HBM; //DEFAULT
MEMSPACE bmemspace = HBM; //DEFAULT
MEMSPACE cmemspace = HBM; //DEFAULT
MEMSPACE workmemspace = HBM; //DEFAULT

char *coloring_input_file = NULL;
char *coloring_output_file = NULL;
char spgemm_step = '0';

namespace MyKokkosSparse{

template <typename OrdinalType, typename Device, typename MemoryTraits, typename SizeType>
class StaticCrsGraph {

public:
	typedef OrdinalType                                            data_type;
	typedef typename Device::execution_space                    execution_space;
	typedef Device                       device_type;
	typedef SizeType                                            size_type;

	typedef Kokkos::View<const size_type* , device_type >  row_map_type;
	typedef Kokkos::View<data_type*  , device_type >  entries_type;

	entries_type entries;
	row_map_type row_map;

	//! Construct an empty view.
	StaticCrsGraph () : entries(), row_map() {}

	//! Copy constructor (shallow copy).
	StaticCrsGraph (const StaticCrsGraph& rhs) : entries (rhs.entries), row_map (rhs.row_map)
	{}

	template<class EntriesType, class RowMapType>
	StaticCrsGraph (const EntriesType& entries_,const RowMapType& row_map_) : entries (entries_), row_map (row_map_)
	{}

	/** \brief  Assign to a view of the rhs array.
	 *          If the old view is the last view
	 *          then allocated memory is deallocated.
	 */
	StaticCrsGraph& operator= (const StaticCrsGraph& rhs) {
		entries = rhs.entries;
		row_map = rhs.row_map;
		return *this;
	}

	/**  \brief  Destroy this view of the array.
	 *           If the last view then allocated memory is deallocated.
	 */
	~StaticCrsGraph() {}
	KOKKOS_INLINE_FUNCTION
	size_type numRows() const {
		return (row_map.dimension_0 () != 0) ?
				row_map.dimension_0 () - static_cast<size_type> (1) :
				static_cast<size_type> (0);
	}
};




template <typename ScalarType, typename OrdinalType, typename Device, typename MemoryTraits, typename SizeType>
class CrsMatrix{
public:
	typedef typename Device::execution_space execution_space;
	typedef typename Device::memory_space memory_space;
	typedef Kokkos::Device<execution_space, memory_space> device_type;
	typedef ScalarType value_type;
	typedef OrdinalType ordinal_type;
	typedef MemoryTraits memory_traits;
	typedef SizeType size_type;

	typedef StaticCrsGraph<OrdinalType, Device, MemoryTraits, SizeType> StaticCrsGraphType;
	typedef typename StaticCrsGraphType::entries_type index_type;
	typedef typename index_type::non_const_value_type const_ordinal_type;
	typedef typename index_type::non_const_value_type non_const_ordinal_type;
	typedef typename StaticCrsGraphType::row_map_type row_map_type;
	typedef Kokkos::View<value_type*, Kokkos::LayoutRight, device_type, MemoryTraits> values_type;
	StaticCrsGraphType graph;
	values_type values;
	CrsMatrix () :
		numCols_ (0)
	{}
	CrsMatrix (const std::string& label,
			const OrdinalType& ncols,
			const values_type& vals,
			const StaticCrsGraphType& graph_) :
				graph (graph_),
				values (vals),
				numCols_ (ncols)
	{
	}

	//! The number of rows in the sparse matrix.
	KOKKOS_INLINE_FUNCTION ordinal_type numRows () const {
		return graph.numRows ();
	}

	//! The number of columns in the sparse matrix.
	KOKKOS_INLINE_FUNCTION ordinal_type numCols () const {
		return numCols_;
	}

	//! The number of stored entries in the sparse matrix.
	KOKKOS_INLINE_FUNCTION size_type nnz () const {
		return graph.entries.dimension_0 ();
	}
	ordinal_type numCols_;
};
}


template <typename ExecSpace, typename crsMat_t, typename crsMat_t2 = crsMat_t, typename crsMat_t3 = crsMat_t, typename TempMemSpace = ExecSpace, typename PersistentMemSpace = ExecSpace>
crsMat_t3 run_experiment(
    crsMat_t crsmat, crsMat_t2 crsmat2,
    int algorithm, int repeat , int chunksize, int multi_color_scale, int shmemsize, int teamsize, int use_dynamic_scheduling, int verbose);


template <typename myExecSpace, typename crsMat_t>
crsMat_t get_crsmat(idx *xadj, idx *adj, wt *ew, idx ne, idx nv, int algo){

    typedef typename crsMat_t::StaticCrsGraphType graph_t;
    typedef typename crsMat_t::row_map_type::non_const_type row_map_view_t;
    typedef typename crsMat_t::index_type::non_const_type   cols_view_t;
    typedef typename crsMat_t::values_type::non_const_type values_view_t;

    row_map_view_t rowmap_view("rowmap_view", nv+1);
    cols_view_t columns_view("colsmap_view", ne);
    values_view_t values_view("values_view", ne);

    KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(ne, ew, values_view);
    KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(ne, adj, columns_view);
    KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(nv+1, xadj, rowmap_view);

    idx ncols = 0;
    KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(ne, columns_view, ncols);
    ncols += 1;

    if (algo == 5) {
      //if algorithm is mkl_csrmultcsr convert to 1 base so that we dont dublicate the memory at the experiments/
      KokkosKernels::Experimental::Util::kk_a_times_x_plus_b< row_map_view_t, row_map_view_t,   int, int, myExecSpace>(nv + 1,  rowmap_view, rowmap_view,  1, 1);
      KokkosKernels::Experimental::Util::kk_a_times_x_plus_b<  cols_view_t, cols_view_t,  int, int, myExecSpace>(ne, columns_view, columns_view,  1, 1);
    }

    graph_t static_graph (columns_view, rowmap_view);
    crsMat_t crsmat("CrsMatrix", ncols, values_view, static_graph);
    return crsmat;
}

template <typename myExecSpace, typename in_crsMat_t, typename out_crsMat_t>
out_crsMat_t copy_crsmat(in_crsMat_t inputMat){

    typedef typename out_crsMat_t::StaticCrsGraphType graph_t;
    typedef typename out_crsMat_t::row_map_type::non_const_type row_map_view_t;
    typedef typename out_crsMat_t::index_type::non_const_type   cols_view_t;
    typedef typename out_crsMat_t::values_type::non_const_type values_view_t;


    typedef typename in_crsMat_t::StaticCrsGraphType in_graph_t;
    typedef typename in_crsMat_t::row_map_type::const_type in_row_map_view_t;
    typedef typename in_crsMat_t::index_type::const_type   in_cols_view_t;
    typedef typename in_crsMat_t::values_type::const_type in_values_view_t;


    const idx nv = inputMat.numRows();
    const idx ne = inputMat.graph.entries.dimension_0();

    row_map_view_t rowmap_view("rowmap_view", nv+1);
    cols_view_t columns_view("colsmap_view", ne);
    values_view_t values_view("values_view", ne);

    KokkosKernels::Experimental::Util::copy_vector<in_values_view_t , values_view_t, myExecSpace>(ne, inputMat.values, values_view);
    KokkosKernels::Experimental::Util::copy_vector<in_cols_view_t , cols_view_t, myExecSpace>(ne, inputMat.graph.entries, columns_view);
    KokkosKernels::Experimental::Util::copy_vector<in_row_map_view_t , row_map_view_t, myExecSpace>(nv+1, inputMat.graph.row_map, rowmap_view);

    idx ncols = 0;
    KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(ne, columns_view, ncols);
    ncols += 1;

    graph_t static_graph (columns_view, rowmap_view);
    out_crsMat_t crsmat("CrsMatrix", ncols, values_view, static_graph);
    return crsmat;
}

template <typename v1>
struct compare{
  v1 f,s;
  compare (v1 f_ , v1 s_): f(f_), s(s_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i, size_t &diff) const {

    if (f[i] - s[i] > 0.00001 || f[i] - s[i] < -0.00001) diff++;
  }

};




enum {
	CMD_USE_THREADS = 0
  , CMD_USE_NUMA
  , CMD_USE_CORE_PER_NUMA
  , CMD_USE_CUDA
  , CMD_USE_OPENMP
  , CMD_USE_CUDA_DEV
  , CMD_SPGEMM_ALGO
  , CMD_BIN_AMTX
  , CMD_BIN_RMTX
  , CMD_BIN_PMTX
  , CMD_MM_MODE
  , CMD_REPEAT
  , CMD_CHUNKSIZE
  , CMD_MULTICOLORSCALE
  , CMD_SHMEMSIZE
  , CMD_TEAMSIZE
  , CMD_DYNAMIC_SCHEDULE
  , CMD_MKL_SORT_OPTION
  , CMD_MKL_KEEP_OUTPUT
  , CMD_CHECK_OUTPUT
  , CMD_MEMSPACES
  , CMD_VERBOSE
  , CMD_ERROR
  , CMD_COUNT };

int main (int argc, char ** argv){



  int cmdline[ CMD_COUNT ] ;
  char *r_mtx_bin_file = NULL;
  char *a_mtx_bin_file = NULL;
  char *p_mtx_bin_file = NULL;

  for ( int i = 0 ; i < CMD_COUNT ; ++i ) cmdline[i] = 0 ;
  cmdline[ CMD_REPEAT ] = 1;
  cmdline[ CMD_CHUNKSIZE ] = -1;
  cmdline[ CMD_MULTICOLORSCALE ] = 1;
  cmdline[ CMD_SHMEMSIZE ] = 16128;
  cmdline[ CMD_TEAMSIZE ] = -1;
  cmdline[ CMD_VERBOSE ] = 0;
  cmdline[ CMD_SPGEMM_ALGO ] = 7;
  cmdline[ CMD_MKL_SORT_OPTION ] = 7;
  cmdline[ CMD_MKL_KEEP_OUTPUT ] = 1;

  for ( int i = 1 ; i < argc ; ++i ) {
    if ( 0 == strcasecmp( argv[i] , "threads" ) ) {
      cmdline[ CMD_USE_THREADS ] = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "openmp" ) ) {
      cmdline[ CMD_USE_OPENMP ] = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "repeat" ) ) {
      cmdline[ CMD_REPEAT ] = atoi( argv[++i] );
    }
    else if ( 0 == strcasecmp( argv[i] , "cores" ) ) {
      sscanf( argv[++i] , "%dx%d" ,
          cmdline + CMD_USE_NUMA ,
          cmdline + CMD_USE_CORE_PER_NUMA );
    }
    else if ( 0 == strcasecmp( argv[i] , "cuda" ) ) {
      cmdline[ CMD_USE_CUDA ] = 1 ;
    }
    else if ( 0 == strcasecmp( argv[i] , "chunksize" ) ) {
      cmdline[ CMD_CHUNKSIZE ] = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "teamsize" ) ) {
      cmdline[ CMD_TEAMSIZE ] = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "vectorsize" ) ) {
      vector_size = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "cuda-dev" ) ) {
      cmdline[ CMD_USE_CUDA ] = 1 ;
      cmdline[ CMD_USE_CUDA_DEV ] = atoi( argv[++i] ) ;
    }

    else if ( 0 == strcasecmp( argv[i] , "mmmode" ) ) {
      cmdline[ CMD_MM_MODE ] = atoi( argv[++i] ) ;
    }

    else if ( 0 == strcasecmp( argv[i] , "memspaces" ) ) {
      cmdline[ CMD_MEMSPACES ] = atoi( argv[++i] ) ;
      int memspaceinfo = cmdline[ CMD_MEMSPACES ];
      std::cout << "memspaceinfo:" << memspaceinfo << std::endl;

      if (memspaceinfo & 1){
    	  amemspace = HBM;
    	  std::cout << "Using HBM for A" << std::endl;
      }
      else {
    	  amemspace = DDR4;
    	  std::cout << "Using DDR4 for A" << std::endl;
      }
      memspaceinfo  = memspaceinfo >> 1;
      if (memspaceinfo & 1){
    	  bmemspace = HBM;
    	  std::cout << "Using HBM for B" << std::endl;
      }
      else {
    	  bmemspace = DDR4;
    	  std::cout << "Using DDR4 for B" << std::endl;
      }
      memspaceinfo  = memspaceinfo >> 1;
      if (memspaceinfo & 1){
    	  cmemspace = HBM;
    	  std::cout << "Using HBM for C" << std::endl;
      }
      else {
    	  cmemspace = DDR4;
    	  std::cout << "Using DDR4 for C" << std::endl;
      }
      memspaceinfo  = memspaceinfo >> 1;
      if (memspaceinfo & 1){
    	  workmemspace = HBM;
    	  std::cout << "Using HBM for work memory space" << std::endl;
      }
      else {
    	  workmemspace = DDR4;
    	  std::cout << "Using DDR4 for work memory space" << std::endl;
      }
      memspaceinfo  = memspaceinfo >> 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "CRWC" ) ) {
    	calculate_read_write_cost = true;
    }
    else if ( 0 == strcasecmp( argv[i] , "CIF" ) ) {
    	coloring_input_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "COF" ) ) {
    	coloring_output_file = argv[++i];
    }

    else if ( 0 == strcasecmp( argv[i] , "mcscale" ) ) {
      cmdline[ CMD_MULTICOLORSCALE ] = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "shmem" ) ) {
      cmdline[ CMD_SHMEMSIZE ] = atoi( argv[++i] ) ;
    }
    else if ( 0 == strcasecmp( argv[i] , "mklsort" ) ) {
      cmdline[ CMD_MKL_SORT_OPTION ] = atoi( argv[++i] ) ;
      mkl_sort_option = cmdline[ CMD_MKL_SORT_OPTION ];
    }
    else if ( 0 == strcasecmp( argv[i] , "mklsort" ) ) {
      cmdline[ CMD_MKL_KEEP_OUTPUT ] = atoi( argv[++i] ) ;
      mkl_keep_output = cmdline[ CMD_MKL_SORT_OPTION ];
    }
    else if ( 0 == strcasecmp( argv[i] , "checkoutput" ) ) {
      cmdline[ CMD_CHECK_OUTPUT ] = 1;
      check_output = cmdline[ CMD_MKL_SORT_OPTION ];
    }
    else if ( 0 == strcasecmp( argv[i] , "amtx" ) ) {
      a_mtx_bin_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "rmtx" ) ) {
      r_mtx_bin_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "pmtx" ) ) {
      p_mtx_bin_file = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "dynamic" ) ) {
      cmdline[ CMD_DYNAMIC_SCHEDULE ]  = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "verbose" ) ) {
      cmdline[ CMD_VERBOSE ]  = 1;
    }
    else if ( 0 == strcasecmp( argv[i] , "algorithm" ) ) {
      ++i;
      if ( 0 == strcasecmp( argv[i] , "MKL" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 1;
      }
      else if ( 0 == strcasecmp( argv[i] , "CUSPARSE" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 2;
      }
      else if ( 0 == strcasecmp( argv[i] , "CUSP" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 3;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKDEBUG" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 4;
      }
      else if ( 0 == strcasecmp( argv[i] , "MKL2" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 5;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKMEM2" ) ) {
    	  cmdline[ CMD_SPGEMM_ALGO ] = 6;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKMEM" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 7;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKSPEED" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 8;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKCOLOR" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 9;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKMULTICOLOR" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 10;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKMULTICOLOR2" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 11;
      }
      else if ( 0 == strcasecmp( argv[i] , "VIENNA" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 12;
      }
      else if ( 0 == strcasecmp( argv[i] , "KKMEMSPEED" ) ) {
        cmdline[ CMD_SPGEMM_ALGO ] = 13;
      }
      else if ( 0 == strcasecmp( argv[i] , "MULTIMEM" ) ) {
          cmdline[ CMD_SPGEMM_ALGO ] = 14;
      }

      else {
        cmdline[ CMD_ERROR ] = 1 ;
        std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
        std::cerr << "Options\n" << std::endl;
        std::cerr << "\t[Required] BACKEND: 'threads[numThreads]' | 'openmp [numThreads]' | 'cuda'" << std::endl;
        std::cerr << "\t[Required] INPUT MATRICES: 'amtx [path_to_a_matrix_market_file.mtx]'  'pmtx [path_to_p_matrix_market_file.mtx]' 'rmtx [path_to_r_matrix_market_file.mtx]'" << std::endl;
        std::cerr << "\t[Required] 'algorithm [KKMEM|KKSPEED|KKCOLOR|KKMULTICOLOR|KKMULTICOLOR2|MKL|CUSPARSE|CUSP|]'" << std::endl;
        std::cerr << "\t[Required] Multiplication mode: 'mmmode [0|1|2]' --> 0:AxA, 1:Rx(AP), 2:(RA)xP" << std::endl;
        std::cerr << "\tThe memory space used for each matrix: 'memspaces [0|1|....15]' --> Bits representing the use of HBM for Work, C, B, and A respectively. For example 12 = 1100, will store work arrays and C on HBM. A and B will be stored DDR. To use this enable multilevel memory in Kokkos, then compile SPGEMM executable with -DKOKKOSKERNELS_MULTILEVELMEM." << std::endl;
        std::cerr << "\t'CRWC': it will perform hypergraph analysis for memory accesses" << std::endl;
        std::cerr << "\t'CIF path_to_coloring_file': If coloring variants are used, colors will be read from this file." << std::endl;
        std::cerr << "\t'COF path_to_coloring_file': If coloring variants are used, first graph coloring will be performed, then writtent to this file." << std::endl;
        std::cerr << "\tLoop scheduling: 'dynamic': Use this for dynamic scheduling of the loops. (Better performance most of the time)" << std::endl;
        std::cerr << "\tVerbose Output: 'verbose'" << std::endl;


        return 0;
      }
    }
    else {
      cmdline[ CMD_ERROR ] = 1 ;
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
      std::cerr << "Options\n" << std::endl;
      std::cerr << "\t[Required] BACKEND: 'threads[numThreads]' | 'openmp [numThreads]' | 'cuda'" << std::endl;
      std::cerr << "\t[Required] INPUT MATRICES: 'amtx [path_to_a_matrix_market_file.mtx]'  'pmtx [path_to_p_matrix_market_file.mtx]' 'rmtx [path_to_r_matrix_market_file.mtx]'" << std::endl;
      std::cerr << "\t[Required] 'algorithm [KKMEM|KKSPEED|KKCOLOR|KKMULTICOLOR|KKMULTICOLOR2|MKL|CUSPARSE|CUSP|]'" << std::endl;
      std::cerr << "\t[Required] Multiplication mode: 'mmmode [0|1|2]' --> 0:AxA, 1:Rx(AP), 2:(RA)xP" << std::endl;
      std::cerr << "\tThe memory space used for each matrix: 'memspaces [0|1|....15]' --> Bits representing the use of HBM for Work, C, B, and A respectively. For example 12 = 1100, will store work arrays and C on HBM. A and B will be stored DDR. To use this enable multilevel memory in Kokkos, then compile SPGEMM executable with -DKOKKOSKERNELS_MULTILEVELMEM." << std::endl;
      std::cerr << "\t'CRWC': it will perform hypergraph analysis for memory accesses" << std::endl;
      std::cerr << "\t'CIF path_to_coloring_file': If coloring variants are used, colors will be read from this file." << std::endl;
      std::cerr << "\t'COF path_to_coloring_file': If coloring variants are used, first graph coloring will be performed, then writtent to this file." << std::endl;
      std::cerr << "\tLoop scheduling: 'dynamic': Use this for dynamic scheduling of the loops. (Better performance most of the time)" << std::endl;
      std::cerr << "\tVerbose Output: 'verbose'" << std::endl;

      return 0;
    }
  }

  if (a_mtx_bin_file == NULL){
    std::cerr << "Provide a mtx binary file" << std::endl ;
    std::cerr << "Options\n" << std::endl;
    std::cerr << "\t[Required] BACKEND: 'threads[numThreads]' | 'openmp [numThreads]' | 'cuda'" << std::endl;
    std::cerr << "\t[Required] INPUT MATRICES: 'amtx [path_to_a_matrix_market_file.mtx]'  'pmtx [path_to_p_matrix_market_file.mtx]' 'rmtx [path_to_r_matrix_market_file.mtx]'" << std::endl;
    std::cerr << "\t[Required] 'algorithm [KKMEM|KKSPEED|KKCOLOR|KKMULTICOLOR|KKMULTICOLOR2|MKL|CUSPARSE|CUSP|]'" << std::endl;
    std::cerr << "\t[Required] Multiplication mode: 'mmmode [0|1|2]' --> 0:AxA, 1:Rx(AP), 2:(RA)xP" << std::endl;
    std::cerr << "\tThe memory space used for each matrix: 'memspaces [0|1|....15]' --> Bits representing the use of HBM for Work, C, B, and A respectively. For example 12 = 1100, will store work arrays and C on HBM. A and B will be stored DDR. To use this enable multilevel memory in Kokkos, then compile SPGEMM executable with -DKOKKOSKERNELS_MULTILEVELMEM." << std::endl;
    std::cerr << "\t'CRWC': it will perform hypergraph analysis for memory accesses" << std::endl;
    std::cerr << "\t'CIF path_to_coloring_file': If coloring variants are used, colors will be read from this file." << std::endl;
    std::cerr << "\t'COF path_to_coloring_file': If coloring variants are used, first graph coloring will be performed, then writtent to this file." << std::endl;
    std::cerr << "\tLoop scheduling: 'dynamic': Use this for dynamic scheduling of the loops. (Better performance most of the time)" << std::endl;
    std::cerr << "\tVerbose Output: 'verbose'" << std::endl;
    return 0;
  }
  if (cmdline[ CMD_MM_MODE ] == 1 || cmdline[ CMD_MM_MODE ] == 2){
    if (r_mtx_bin_file == NULL){
      std::cerr << "Provide a r mtx binary file rmtx rmatrix_file" << std::endl ;
      std::cerr << "Options\n" << std::endl;
      std::cerr << "\t[Required] BACKEND: 'threads[numThreads]' | 'openmp [numThreads]' | 'cuda'" << std::endl;
      std::cerr << "\t[Required] INPUT MATRICES: 'amtx [path_to_a_matrix_market_file.mtx]'  'pmtx [path_to_p_matrix_market_file.mtx]' 'rmtx [path_to_r_matrix_market_file.mtx]'" << std::endl;
      std::cerr << "\t[Required] 'algorithm [KKMEM|KKSPEED|KKCOLOR|KKMULTICOLOR|KKMULTICOLOR2|MKL|CUSPARSE|CUSP|]'" << std::endl;
      std::cerr << "\t[Required] Multiplication mode: 'mmmode [0|1|2]' --> 0:AxA, 1:Rx(AP), 2:(RA)xP" << std::endl;
      std::cerr << "\tThe memory space used for each matrix: 'memspaces [0|1|....15]' --> Bits representing the use of HBM for Work, C, B, and A respectively. For example 12 = 1100, will store work arrays and C on HBM. A and B will be stored DDR. To use this enable multilevel memory in Kokkos, then compile SPGEMM executable with -DKOKKOSKERNELS_MULTILEVELMEM." << std::endl;
      std::cerr << "\t'CRWC': it will perform hypergraph analysis for memory accesses" << std::endl;
      std::cerr << "\t'CIF path_to_coloring_file': If coloring variants are used, colors will be read from this file." << std::endl;
      std::cerr << "\t'COF path_to_coloring_file': If coloring variants are used, first graph coloring will be performed, then writtent to this file." << std::endl;
      std::cerr << "\tLoop scheduling: 'dynamic': Use this for dynamic scheduling of the loops. (Better performance most of the time)" << std::endl;
      std::cerr << "\tVerbose Output: 'verbose'" << std::endl;
      return 0;
    }
    if (p_mtx_bin_file == NULL){
      std::cerr << "Provide a p mtx binary file pmtx rmatrix_file" << std::endl ;
      std::cerr << "Options\n" << std::endl;
      std::cerr << "\t[Required] BACKEND: 'threads[numThreads]' | 'openmp [numThreads]' | 'cuda'" << std::endl;
      std::cerr << "\t[Required] INPUT MATRICES: 'amtx [path_to_a_matrix_market_file.mtx]'  'pmtx [path_to_p_matrix_market_file.mtx]' 'rmtx [path_to_r_matrix_market_file.mtx]'" << std::endl;
      std::cerr << "\t[Required] 'algorithm [KKMEM|KKSPEED|KKCOLOR|KKMULTICOLOR|KKMULTICOLOR2|MKL|CUSPARSE|CUSP|]'" << std::endl;
      std::cerr << "\t[Required] Multiplication mode: 'mmmode [0|1|2]' --> 0:AxA, 1:Rx(AP), 2:(RA)xP" << std::endl;
      std::cerr << "\tThe memory space used for each matrix: 'memspaces [0|1|....15]' --> Bits representing the use of HBM for Work, C, B, and A respectively. For example 12 = 1100, will store work arrays and C on HBM. A and B will be stored DDR. To use this enable multilevel memory in Kokkos, then compile SPGEMM executable with -DKOKKOSKERNELS_MULTILEVELMEM." << std::endl;
      std::cerr << "\t'CRWC': it will perform hypergraph analysis for memory accesses" << std::endl;
      std::cerr << "\t'CIF path_to_coloring_file': If coloring variants are used, colors will be read from this file." << std::endl;
      std::cerr << "\t'COF path_to_coloring_file': If coloring variants are used, first graph coloring will be performed, then writtent to this file." << std::endl;
      std::cerr << "\tLoop scheduling: 'dynamic': Use this for dynamic scheduling of the loops. (Better performance most of the time)" << std::endl;
      std::cerr << "\tVerbose Output: 'verbose'" << std::endl;
      return 0;
    }
  }

  idx m = 0, nnzA = 0;
  idx *xadj, *adj;
  wt *ew;
#if defined( KOKKOS_HAVE_PTHREAD )

  if ( cmdline[ CMD_USE_THREADS ] ) {

    if ( cmdline[ CMD_USE_NUMA ] && cmdline[ CMD_USE_CORE_PER_NUMA ] ) {
      Kokkos::Threads::initialize( cmdline[ CMD_USE_THREADS ] ,
          cmdline[ CMD_USE_NUMA ] ,
          cmdline[ CMD_USE_CORE_PER_NUMA ] );
    }
    else {
      Kokkos::Threads::initialize( cmdline[ CMD_USE_THREADS ] );
    }

    if (cmdline[ CMD_SPGEMM_ALGO ] == 2 || cmdline[ CMD_SPGEMM_ALGO ] == 3){
      std::cerr << "CUSP and CUSPARSE cannot be run with PTHREADS" << std::endl ;
      return 0;
    }

    KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, a_mtx_bin_file);
    idx nv = m;
    idx ne = nnzA;
    Kokkos::Threads::print_configuration(std::cout);


    typedef Kokkos::Threads myExecSpace;
    typedef typename MyKokkosSparse::CrsMatrix<wt, idx, myExecSpace, void, size_type > crsMat_t;

    typedef typename crsMat_t::StaticCrsGraphType graph_t;
    typedef typename graph_t::row_map_type::non_const_type row_map_view_t;
    typedef typename graph_t::entries_type::non_const_type   cols_view_t;
    typedef typename crsMat_t::values_type::non_const_type values_view_t;

    row_map_view_t rowmap_view("rowmap_view", nv+1);
    cols_view_t columns_view("colsmap_view", ne);
    values_view_t values_view("values_view", ne);

    KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(ne, ew, values_view);
    KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(ne, adj, columns_view);
    KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(nv+1, xadj, rowmap_view);

    graph_t static_graph (columns_view, rowmap_view);
    crsMat_t crsmat("CrsMatrix", nv, values_view, static_graph);
    delete [] xadj;
    delete [] adj;
    delete [] ew;

    if (cmdline[ CMD_MM_MODE ] == 0){
      std::cout << "MULTIPLYING A*A" << std::endl;
      run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ], cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
    }
    else if (cmdline[ CMD_MM_MODE ] == 1){
      {
        std::cout << "MULTIPLYING A*P" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, p_mtx_bin_file);


        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;

        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ], cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
      }
      {
        std::cout << "MULTIPLYING R*(AP)" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, r_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;


        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ], cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
      }
    }
    else if (cmdline[ CMD_MM_MODE ] == 2){
      {
        std::cout << "MULTIPLYING R*A" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, r_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;


        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ], cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
      }
      {
        std::cout << "MULTIPLYING (RA)*P" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, p_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;


        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);

        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
      }
    }


    myExecSpace::finalize();
  }

#endif

#if defined( KOKKOS_HAVE_OPENMP )

  if ( cmdline[ CMD_USE_OPENMP ] ) {

	  if ( cmdline[ CMD_USE_NUMA ] && cmdline[ CMD_USE_CORE_PER_NUMA ] ) {
		  Kokkos::OpenMP::initialize( cmdline[ CMD_USE_OPENMP ] ,
				  cmdline[ CMD_USE_NUMA ] ,
				  cmdline[ CMD_USE_CORE_PER_NUMA ] );
	  }
	  else {
		  Kokkos::OpenMP::initialize( cmdline[ CMD_USE_OPENMP ] );
	  }
	  if (cmdline[ CMD_SPGEMM_ALGO ] == 2 || cmdline[ CMD_SPGEMM_ALGO ] == 3){
		  std::cerr << "CUSP and CUSPARSE cannot be run with OPENMP" << std::endl ;
		  return 0;
	  }

	  Kokkos::OpenMP::print_configuration(std::cout);

	  KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, a_mtx_bin_file);
	  idx nv = m;
	  idx ne = nnzA;


	  typedef Kokkos::OpenMP myExecSpace;
	  typedef typename MyKokkosSparse::CrsMatrix<wt, idx, myExecSpace, void, size_type > crsMat_t;
	  typedef typename crsMat_t::StaticCrsGraphType graph_t;
	  typedef typename crsMat_t::row_map_type::non_const_type row_map_view_t;
	  typedef typename crsMat_t::index_type::non_const_type   cols_view_t;
	  typedef typename crsMat_t::values_type::non_const_type values_view_t;
	  crsMat_t crsmat;

#ifdef KOKKOSKERNELS_MULTILEVELMEM
	  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> myHostExecSpace;
	  typedef typename MyKokkosSparse::CrsMatrix<wt, idx, myHostExecSpace, void, size_type > crsMat_host_t;
	  typedef typename crsMat_host_t::StaticCrsGraphType host_graph_t;
	  typedef typename crsMat_host_t::row_map_type::non_const_type host_row_map_view_t;
	  typedef typename crsMat_host_t::index_type::non_const_type   host_cols_view_t;
	  typedef typename crsMat_host_t::values_type::non_const_type host_values_view_t;
	  crsMat_host_t host_crsmat;
#endif







#ifdef KOKKOSKERNELS_MULTILEVELMEM
	  if (cmdline[ CMD_MM_MODE ] != 2){


		  if (amemspace == HBM){
			  crsmat = get_crsmat<myExecSpace, crsMat_t>(xadj, adj, ew, ne, nv, cmdline[ CMD_SPGEMM_ALGO ]);
		  }
		  else {
			  host_crsmat = get_crsmat<myExecSpace, crsMat_host_t>(xadj, adj, ew, ne, nv, cmdline[ CMD_SPGEMM_ALGO ]);
		  }

		  crsmat = get_crsmat<myExecSpace, crsMat_t>(xadj, adj, ew, ne, nv, cmdline[ CMD_SPGEMM_ALGO ]);

	  }
	  else {
		  if (bmemspace == HBM){
			  crsmat = get_crsmat<myExecSpace, crsMat_t>(xadj, adj, ew, ne, nv, cmdline[ CMD_SPGEMM_ALGO ]);
		  }
		  else {
			  host_crsmat = get_crsmat<myExecSpace, crsMat_host_t>(xadj, adj, ew, ne, nv, cmdline[ CMD_SPGEMM_ALGO ]);
		  }

	  }
#else
	  crsmat = get_crsmat<myExecSpace, crsMat_t>(xadj, adj, ew, ne, nv, cmdline[ CMD_SPGEMM_ALGO ]);
#endif
	  delete [] xadj;
	  delete [] adj;
	  delete [] ew;


	  //std::cout << "STARTUP MULTIPLYING A*A" << std::endl;
	  //run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat,  cmdline[ CMD_SPGEMM_ALGO ], cmdline[ CMD_REPEAT ], cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
	  //std::cout << "STARTUP DONE  A*A\n\n\n\n\n" << std::endl;

	  if (cmdline[ CMD_MM_MODE ] == 0){

		  run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>(
				  crsmat, crsmat,
				  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ], cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
	  }else if (cmdline[ CMD_MM_MODE ] == 1){
		  {
			  std::cout << "MULTIPLYING A*P" << std::endl;
			  idx m_ = 0, nnzA_ = 0;
			  idx *xadj_, *adj_;
			  wt *ew_;
			  KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m_, &nnzA_, &xadj_, &adj_, &ew_, p_mtx_bin_file);
			  crsMat_t crsmat2;

#ifdef KOKKOSKERNELS_MULTILEVELMEM
			  crsMat_host_t host_crsmat2;
			  if (bmemspace == HBM){
				  crsmat2 = get_crsmat<myExecSpace, crsMat_t>(xadj_, adj_, ew_, nnzA_, m_, cmdline[ CMD_SPGEMM_ALGO ]);
			  }
			  else {
				  host_crsmat2 = get_crsmat<myExecSpace, crsMat_host_t>(xadj_, adj_, ew_, nnzA_, m_, cmdline[ CMD_SPGEMM_ALGO ]);
			  }
#else
			  crsmat2 = get_crsmat<myExecSpace, crsMat_t>(xadj_, adj_, ew_, nnzA_, m_, cmdline[ CMD_SPGEMM_ALGO ]);

#endif
			  delete [] xadj_;
			  delete [] adj_;
			  delete [] ew_;
#ifdef KOKKOSKERNELS_MULTILEVELMEM

			  if (amemspace == HBM){
				  if (bmemspace == HBM){
					  if (cmemspace == HBM){
						  if (workmemspace == HBM){
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
							  (crsmat, crsmat2,
									  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
									  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
									  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
									  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  }
						  else {
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
							  (crsmat, crsmat2,
									  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
									  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
									  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
									  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  }
					  }
					  else{

						  if (workmemspace == HBM){
							  host_crsmat =
									  run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_host_t, myExecSpace, myExecSpace>
							  (crsmat, crsmat2,
									  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
									  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
									  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
									  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  }
						  else {
							  host_crsmat =
									  run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
							  (crsmat, crsmat2,
									  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
									  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
									  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
									  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  }

						  crsmat = copy_crsmat<myExecSpace, crsMat_host_t, crsMat_t>(host_crsmat);
						  host_crsmat = crsMat_host_t ();
					  }
				  }
				  else{
					  if (cmemspace == HBM){
						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_t, myExecSpace, myExecSpace>
						  (crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);

						  host_crsmat = copy_crsmat<myExecSpace, crsMat_t, crsMat_host_t>(crsmat);
						  crsmat = crsMat_t ();
					  }
					  else{
						  if (workmemspace == HBM)
							  host_crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
				  }
			  }
			  else {
				  if (bmemspace == HBM){
					  if (cmemspace == HBM){
						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
						  (host_crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
					  else{
						  if (workmemspace == HBM)
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (host_crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  crsmat = copy_crsmat<myExecSpace, crsMat_host_t, crsMat_t>(host_crsmat);
						  host_crsmat = crsMat_host_t ();
					  }
				  }
				  else{
					  if (cmemspace == HBM){
						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_t, myExecSpace, myExecSpace>
						  (host_crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  host_crsmat = copy_crsmat<myExecSpace, crsMat_t, crsMat_host_t>(crsmat);
						  crsmat = crsMat_t ();
					  }
					  else{
						  if (workmemspace == HBM)
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (host_crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
				  }
			  }
#else
			  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
			  (crsmat, crsmat2,
					  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
					  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
					  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
					  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
#endif



			  /*
        char *file = "AP.mtx";
        KokkosKernels::Experimental::Util::write_graph_bin(
            m, (idx) crsmat.graph.entries.dimension_0(),
            ( const idx *)crsmat.graph.row_map.ptr_on_device(),
            ( const idx *)crsmat.graph.entries.ptr_on_device(),
            ( const wt *)crsmat.values.ptr_on_device(),
            ( const char *)file);
			   */
		  }
		  {
			  std::cout << "MULTIPLYING R*(AP)" << std::endl;
			  idx m__ = 0, nnzA__ = 0;
			  idx *xadj__, *adj__;
			  wt *ew__;
			  KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m__, &nnzA__, &xadj__, &adj__, &ew__, r_mtx_bin_file);

			  crsMat_t crsmat2;
#ifdef KOKKOSKERNELS_MULTILEVELMEM
			  crsMat_host_t host_crsmat2;

			  if (amemspace == HBM){
				  crsmat2 = get_crsmat<myExecSpace, crsMat_t>(xadj__, adj__, ew__, nnzA__, m__, cmdline[ CMD_SPGEMM_ALGO ]);
			  }
			  else {
				  host_crsmat2 = get_crsmat<myExecSpace, crsMat_host_t>(xadj__, adj__, ew__, nnzA__, m__, cmdline[ CMD_SPGEMM_ALGO ]);
			  }
#else
			  crsmat2 = get_crsmat<myExecSpace, crsMat_t>(xadj__, adj__, ew__, nnzA__, m__, cmdline[ CMD_SPGEMM_ALGO ]);

#endif

			  delete [] xadj__;
			  delete [] adj__;
			  delete [] ew__;

			  //cmdline[ CMD_SPGEMM_ALGO ] = 1;
			  /* crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
        		(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ], cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
			   */
#ifdef KOKKOSKERNELS_MULTILEVELMEM

			  if (amemspace == HBM){
				  if (bmemspace == HBM){
					  if (cmemspace == HBM){

						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
						  (crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else

							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t,Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
					  else{
						  if (workmemspace == HBM)
							  host_crsmat =
									  run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat =
									  run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);


					  }
				  }
				  else{
					  if (cmemspace == HBM){
						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_t, myExecSpace, myExecSpace>
						  (crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);

					  }
					  else{

						  if (workmemspace == HBM)
							  host_crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
				  }
			  }

			  else {
				  if (bmemspace == HBM){
					  if (cmemspace == HBM){
						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
						  (host_crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
					  else{

						  if (workmemspace == HBM)
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (host_crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);

					  }
				  }
				  else{
					  if (cmemspace == HBM){

						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_t, myExecSpace, myExecSpace>
						  (host_crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);

					  }
					  else{
						  if (workmemspace == HBM)
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (host_crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
				  }
			  }
#else
			  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
		  (crsmat2, crsmat,
				  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
				  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
				  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
				  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
#endif

		  }

	  }
	  else if (cmdline[ CMD_MM_MODE ] == 2){
		  {
			  std::cout << "MULTIPLYING R*A" << std::endl;
			  idx m_ = 0, nnzA_ = 0;
			  idx *xadj_, *adj_;
			  wt *ew_;
			  KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m_, &nnzA_, &xadj_, &adj_, &ew_, r_mtx_bin_file);

			  crsMat_t crsmat2;

#ifdef KOKKOSKERNELS_MULTILEVELMEM
			  crsMat_host_t host_crsmat2;
			  if (amemspace == HBM){
				  crsmat2 = get_crsmat<myExecSpace, crsMat_t>(xadj_, adj_, ew_, nnzA_, m_, cmdline[ CMD_SPGEMM_ALGO ]);
			  }
			  else {
				  host_crsmat2 = get_crsmat<myExecSpace, crsMat_host_t>(xadj_, adj_, ew_, nnzA_, m_, cmdline[ CMD_SPGEMM_ALGO ]);
			  }
#else
			  crsmat2 = get_crsmat<myExecSpace, crsMat_t>(xadj_, adj_, ew_, nnzA_, m_, cmdline[ CMD_SPGEMM_ALGO ]);

#endif

			  delete [] xadj_;
			  delete [] adj_;
			  delete [] ew_;

			  //crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
			  //		(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
#ifdef KOKKOSKERNELS_MULTILEVELMEM
			  if (amemspace == HBM){
				  if (bmemspace == HBM){
					  if (cmemspace == HBM){

						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
						  (crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
					  else{
						  if (workmemspace == HBM)
							  host_crsmat =
									  run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat =
									  run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);

						  crsmat = copy_crsmat<myExecSpace, crsMat_host_t, crsMat_t>(host_crsmat);
						  host_crsmat = crsMat_host_t ();
					  }
				  }
				  else{
					  if (cmemspace == HBM){

						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_t, myExecSpace, myExecSpace>
						  (crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);

					  }
					  else{
						  if (workmemspace == HBM)
							  host_crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  crsmat = copy_crsmat<myExecSpace, crsMat_host_t, crsMat_t>(host_crsmat);
						  host_crsmat = crsMat_host_t ();
					  }
				  }
			  }


			  else {
				  if (bmemspace == HBM){
					  if (cmemspace == HBM){
						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
						  (host_crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  host_crsmat = copy_crsmat<myExecSpace, crsMat_t, crsMat_host_t>(crsmat);
						  crsmat = crsMat_t ();
					  }
					  else{
						  if (workmemspace == HBM)
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (host_crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat2, crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
				  }
				  else{
					  if (cmemspace == HBM){
						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_t, myExecSpace, myExecSpace>
						  (host_crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  host_crsmat = copy_crsmat<myExecSpace, crsMat_t, crsMat_host_t>(crsmat);
						  crsmat = crsMat_t ();
					  }
					  else{

						  if (workmemspace == HBM)
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (host_crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat2, host_crsmat,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
				  }
			  }
#else
			  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
		  (crsmat2, crsmat,
				  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
				  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
				  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
				  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
#endif


		  }
		  {
			  std::cout << "MULTIPLYING (RA)*P" << std::endl;
			  idx m_ = 0, nnzA_ = 0;
			  idx *xadj_, *adj_;
			  wt *ew_;
			  KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m_, &nnzA_, &xadj_, &adj_, &ew_, p_mtx_bin_file);

			  crsMat_t crsmat2;
#ifdef KOKKOSKERNELS_MULTILEVELMEM
			  crsMat_host_t host_crsmat2;
			  if (bmemspace == HBM){
				  crsmat2 = get_crsmat<myExecSpace, crsMat_t>(xadj_, adj_, ew_, nnzA_, m_, cmdline[ CMD_SPGEMM_ALGO ]);
			  }
			  else {
				  host_crsmat2 = get_crsmat<myExecSpace, crsMat_host_t>(xadj_, adj_, ew_, nnzA_, m_, cmdline[ CMD_SPGEMM_ALGO ]);
			  }
#else
			  crsmat2 = get_crsmat<myExecSpace, crsMat_t>(xadj_, adj_, ew_, nnzA_, m_, cmdline[ CMD_SPGEMM_ALGO ]);

#endif

			  delete [] xadj_;
			  delete [] adj_;
			  delete [] ew_;

			  /*
        crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
        		(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
			   */
#ifdef KOKKOSKERNELS_MULTILEVELMEM

			  if (amemspace == HBM){
				  if (bmemspace == HBM){
					  if (cmemspace == HBM){

						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
						  (crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
					  else{

						  if (workmemspace == HBM)
							  host_crsmat =
									  run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat =
									  run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);



					  }
				  }
				  else{
					  if (cmemspace == HBM){
						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_t, myExecSpace, myExecSpace>
						  (crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);

					  }
					  else{

						  if (workmemspace == HBM)
							  host_crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);

					  }
				  }
			  }

			  else {
				  if (bmemspace == HBM){
					  if (cmemspace == HBM){
						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
						  (host_crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
					  else{

						  if (workmemspace == HBM)
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (host_crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);

					  }
				  }
				  else{
					  if (cmemspace == HBM){
						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_t, myExecSpace, myExecSpace>
						  (host_crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);

					  }
					  else{
						  if (workmemspace == HBM)
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (host_crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
				  }
			  }
#else
			  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
		  (crsmat, crsmat2,
				  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
				  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
				  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
				  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
#endif
		  }
	  }
	  else if (cmdline[ CMD_MM_MODE ] == 3){
		  {
			  std::cout << "MULTIPLYING R*P" << std::endl;
			  idx m_ = 0, nnzA_ = 0;
			  idx *xadj_, *adj_;
			  wt *ew_;
			  KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m_, &nnzA_, &xadj_, &adj_, &ew_, r_mtx_bin_file);



			  crsMat_t crsmat2;
#ifdef KOKKOSKERNELS_MULTILEVELMEM
			  crsMat_host_t host_crsmat2;
			  if (amemspace == HBM){
				  crsmat2 = get_crsmat<myExecSpace, crsMat_t>(xadj_, adj_, ew_, nnzA_, m_, cmdline[ CMD_SPGEMM_ALGO ]);
			  }
			  else {
				  host_crsmat2 = get_crsmat<myExecSpace, crsMat_host_t>(xadj_, adj_, ew_, nnzA_, m_, cmdline[ CMD_SPGEMM_ALGO ]);
			  }
#else
			  crsmat2 = get_crsmat<myExecSpace, crsMat_t>(xadj_, adj_, ew_, nnzA_, m_, cmdline[ CMD_SPGEMM_ALGO ]);
#endif


			  delete [] xadj_;
			  delete [] adj_;
			  delete [] ew_;
			  crsmat = crsmat2;
#ifdef KOKKOSKERNELS_MULTILEVELMEM

			  host_crsmat= host_crsmat2;

			  KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m_, &nnzA_, &xadj_, &adj_, &ew_, p_mtx_bin_file);
			  if (bmemspace == HBM){
				  crsmat2 = get_crsmat<myExecSpace, crsMat_t>(xadj_, adj_, ew_, nnzA_, m_, cmdline[ CMD_SPGEMM_ALGO ]);
			  }
			  else {
				  host_crsmat2 = get_crsmat<myExecSpace, crsMat_host_t>(xadj_, adj_, ew_, nnzA_, m_, cmdline[ CMD_SPGEMM_ALGO ]);
			  }
#else
			  crsmat2 = get_crsmat<myExecSpace, crsMat_t>(xadj_, adj_, ew_, nnzA_, m_, cmdline[ CMD_SPGEMM_ALGO ]);

#endif
			  delete [] xadj_;
			  delete [] adj_;
			  delete [] ew_;
			  //crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
			  //		(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
#ifdef KOKKOSKERNELS_MULTILEVELMEM

			  if (amemspace == HBM){
				  if (bmemspace == HBM){
					  if (cmemspace == HBM){

						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
						  (crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
					  else{
						  if (workmemspace == HBM)
							  host_crsmat =
									  run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat =
									  run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);


					  }
				  }
				  else{
					  if (cmemspace == HBM){
						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_t, myExecSpace, myExecSpace>
						  (crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else

							  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
					  else{
						  if (workmemspace == HBM)

							  host_crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_host_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);

					  }
				  }
			  }

			  else {
				  if (bmemspace == HBM){
					  if (cmemspace == HBM){
						  if (workmemspace == HBM)

							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
						  (host_crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
					  else{
						  if (workmemspace == HBM)
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (host_crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat, crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);

					  }
				  }
				  else{
					  if (cmemspace == HBM){

						  if (workmemspace == HBM)
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_t, myExecSpace, myExecSpace>
						  (host_crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else
							  crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);


					  }
					  else{
						  if (workmemspace == HBM)
							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_host_t, myExecSpace, myExecSpace>
						  (host_crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
						  else

							  host_crsmat = run_experiment<myExecSpace, crsMat_host_t,crsMat_host_t,crsMat_host_t, Kokkos::HostSpace, Kokkos::HostSpace>
						  (host_crsmat, host_crsmat2,
								  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
								  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
								  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
								  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
					  }
				  }
			  }
#else
			  crsmat = run_experiment<myExecSpace, crsMat_t,crsMat_t,crsMat_t, myExecSpace, myExecSpace>
		  (crsmat, crsmat2,
				  cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],
				  cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ],
				  cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ],
				  cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
#endif

		  }

	  }


	  myExecSpace::finalize();
  }

#endif

#if defined( KOKKOS_HAVE_CUDA )
  if ( cmdline[ CMD_USE_CUDA ] ) {
    // Use the last device:

    if (cmdline[ CMD_SPGEMM_ALGO ] == 1){
      std::cerr << "MKL cannot be run with CUDA" << std::endl ;
      return 0;
    }
    Kokkos::HostSpace::execution_space::initialize();
    Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice( cmdline[ CMD_USE_CUDA_DEV ] ) );
    Kokkos::Cuda::print_configuration(std::cout);

    {
      //just warm up gpu.
    Kokkos::View<int *, Kokkos::Cuda> tmp_view("rowmap_view", 20000000);
    Kokkos::deep_copy(tmp_view, 2);
    }



    KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, a_mtx_bin_file);
    idx nv = m;
    idx ne = nnzA;

    typedef Kokkos::Cuda myExecSpace;
    typedef typename MyKokkosSparse::CrsMatrix<wt, idx, myExecSpace, void, size_type > crsMat_t;

    typedef typename crsMat_t::StaticCrsGraphType graph_t;
    typedef typename crsMat_t::row_map_type::non_const_type row_map_view_t;
    typedef typename crsMat_t::index_type::non_const_type   cols_view_t;
    typedef typename crsMat_t::values_type::non_const_type values_view_t;

    row_map_view_t rowmap_view("rowmap_view", nv+1);
    cols_view_t columns_view("colsmap_view", ne);
    values_view_t values_view("values_view", ne);


    {
      typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
      typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);
      typename values_view_t::HostMirror hv = Kokkos::create_mirror_view (values_view);

      for (idx i = 0; i <= nv; ++i){
        hr(i) = xadj[i];
      }

      for (idx i = 0; i < nnzA; ++i){
        hc(i) = adj[i];
        hv(i) = ew[i];
      }
      Kokkos::deep_copy (rowmap_view , hr);
      Kokkos::deep_copy (columns_view , hc);
      Kokkos::deep_copy (values_view , hv);


    }
    graph_t static_graph (columns_view, rowmap_view);
    crsMat_t crsmat("CrsMatrix", nv, values_view, static_graph);


    //n = m;


    if (cmdline[ CMD_MM_MODE ] == 0){
      std::cout << "MULTIPLYING A*A" << std::endl;
      run_experiment<Kokkos::Cuda, crsMat_t>(crsmat, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
    }
    else if (cmdline[ CMD_MM_MODE ] == 1){

      {
        std::cout << "MULTIPLYING A*P" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, p_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        {

          typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
          typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);
          typename values_view_t::HostMirror hv = Kokkos::create_mirror_view (values_view);

          for (idx i = 0; i <= m; ++i){
            hr(i) = xadj[i];
          }

          for (idx i = 0; i < nnzA; ++i){
            hc(i) = adj[i];
            hv(i) = ew[i];
          }
          Kokkos::deep_copy (rowmap_view , hr);
          Kokkos::deep_copy (columns_view , hc);
          Kokkos::deep_copy (values_view , hv);
        }

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;

        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);

      }
      {
        std::cout << "MULTIPLYING R*(AP)" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, r_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);

        {

          typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
          typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);
          typename values_view_t::HostMirror hv = Kokkos::create_mirror_view (values_view);

          for (idx i = 0; i <= m; ++i){
            hr(i) = xadj[i];
          }

          for (idx i = 0; i < nnzA; ++i){
            hc(i) = adj[i];
            hv(i) = ew[i];
          }
          Kokkos::deep_copy (rowmap_view , hr);
          Kokkos::deep_copy (columns_view , hc);
          Kokkos::deep_copy (values_view , hv);
        }


        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;

        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
      }
    }
    else if (cmdline[ CMD_MM_MODE ] == 2){
      {
        std::cout << "MULTIPLYING R*A" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, r_mtx_bin_file);

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);



        {

          typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
          typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);
          typename values_view_t::HostMirror hv = Kokkos::create_mirror_view (values_view);

          for (idx i = 0; i <= m; ++i){
            hr(i) = xadj[i];
          }

          for (idx i = 0; i < nnzA; ++i){
            hc(i) = adj[i];
            hv(i) = ew[i];
          }
          Kokkos::deep_copy (rowmap_view , hr);
          Kokkos::deep_copy (columns_view , hc);
          Kokkos::deep_copy (values_view , hv);
        }

        //KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        //KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        //KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;


        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat2, crsmat, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);

      }
      {
        std::cout << "MULTIPLYING (RA)*P" << std::endl;
        idx m = 0, nnzA = 0;
        idx *xadj, *adj;
        wt *ew;
        KokkosKernels::Experimental::Util::read_matrix<idx, wt> (&m, &nnzA, &xadj, &adj, &ew, p_mtx_bin_file);
        std::cout << 1 << std::endl;

        row_map_view_t rowmap_view("rowmap_view", m+1);
        cols_view_t columns_view("colsmap_view", nnzA);
        values_view_t values_view("values_view", nnzA);


        {

          typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view (rowmap_view);
          typename cols_view_t::HostMirror hc = Kokkos::create_mirror_view (columns_view);
          typename values_view_t::HostMirror hv = Kokkos::create_mirror_view (values_view);
          std::cout << 2 << std::endl;
          for (idx i = 0; i <= m; ++i){
            hr(i) = xadj[i];
          }
          std::cout << 3 << std::endl;
          for (idx i = 0; i < nnzA; ++i){
            hc(i) = adj[i];
            hv(i) = ew[i];
          }
          std::cout << 4 << std::endl;
          Kokkos::deep_copy (rowmap_view , hr);
          Kokkos::deep_copy (columns_view , hc);
          Kokkos::deep_copy (values_view , hv);
          std::cout << 5 << std::endl;
        }

        //KokkosKernels::Experimental::Util::copy_vector<wt * , values_view_t, myExecSpace>(nnzA, ew, values_view);
        //KokkosKernels::Experimental::Util::copy_vector<idx * , cols_view_t, myExecSpace>(nnzA, adj, columns_view);
        //KokkosKernels::Experimental::Util::copy_vector<idx * , row_map_view_t, myExecSpace>(m+1, xadj, rowmap_view);

        idx ncols = 0;
        KokkosKernels::Experimental::Util::view_reduce_max<cols_view_t, myExecSpace>(nnzA, columns_view, ncols);
        ncols += 1;


        graph_t static_graph (columns_view, rowmap_view);
        crsMat_t crsmat2("CrsMatrix2", ncols, values_view, static_graph);
        delete [] xadj;
        delete [] adj;
        delete [] ew;
        crsmat = run_experiment<myExecSpace, crsMat_t>(crsmat, crsmat2, cmdline[ CMD_SPGEMM_ALGO ],cmdline[ CMD_REPEAT ],cmdline[ CMD_CHUNKSIZE ], cmdline[ CMD_MULTICOLORSCALE ], cmdline[ CMD_SHMEMSIZE ], cmdline[ CMD_TEAMSIZE ], cmdline[ CMD_DYNAMIC_SCHEDULE ], cmdline[ CMD_VERBOSE]);
      }
    }




    Kokkos::Cuda::finalize();
    Kokkos::HostSpace::execution_space::finalize();
  }

#endif


  return 0;

}


template <typename crsMat_t, typename device>
bool is_same_matrix(crsMat_t output_mat1, crsMat_t output_mat2){

  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type   lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  size_t nrows1 = output_mat1.graph.row_map.dimension_0();
  size_t nentries1 = output_mat1.graph.entries.dimension_0() ;
  size_t nvals1 = output_mat1.values.dimension_0();

  size_t nrows2 = output_mat2.graph.row_map.dimension_0();
  size_t nentries2 = output_mat2.graph.entries.dimension_0() ;
  size_t nvals2 = output_mat2.values.dimension_0();


  lno_nnz_view_t h_ent1 (Kokkos::ViewAllocateWithoutInitializing("e1"), nentries1);
  scalar_view_t h_vals1 (Kokkos::ViewAllocateWithoutInitializing("v1"), nvals1);


  KokkosKernels::Experimental::Util::kk_sort_graph<typename graph_t::row_map_type,
    typename graph_t::entries_type,
    typename crsMat_t::values_type,
    lno_nnz_view_t,
    scalar_view_t,
    typename device::execution_space
    >(
    output_mat1.graph.row_map, output_mat1.graph.entries, output_mat1.values,
    h_ent1, h_vals1
  );

  lno_nnz_view_t h_ent2 (Kokkos::ViewAllocateWithoutInitializing("e1"), nentries2);
  scalar_view_t h_vals2 (Kokkos::ViewAllocateWithoutInitializing("v1"), nvals2);

  if (nrows1 != nrows2) return false;
  if (nentries1 != nentries2) return false;
  if (nvals1 != nvals2) return false;

  KokkosKernels::Experimental::Util::kk_sort_graph
      <typename graph_t::row_map_type,
      typename graph_t::entries_type,
      typename crsMat_t::values_type,
      lno_nnz_view_t,
      scalar_view_t,
      typename device::execution_space
      >(
      output_mat2.graph.row_map, output_mat2.graph.entries, output_mat2.values,
      h_ent2, h_vals2
    );

  bool is_identical = true;
  is_identical = KokkosKernels::Experimental::Util::kk_is_identical_view
      <typename graph_t::row_map_type, typename graph_t::row_map_type, typename lno_view_t::value_type,
      typename device::execution_space>(output_mat1.graph.row_map, output_mat2.graph.row_map, 0);
  if (!is_identical) return false;

  is_identical = KokkosKernels::Experimental::Util::kk_is_identical_view
      <lno_nnz_view_t, lno_nnz_view_t, typename lno_nnz_view_t::value_type,
      typename device::execution_space>(h_ent1, h_ent2, 0 );
  if (!is_identical) return false;

  is_identical = KokkosKernels::Experimental::Util::kk_is_identical_view
      <scalar_view_t, scalar_view_t, typename scalar_view_t::value_type,
      typename device::execution_space>(h_vals1, h_vals2, 0.000001);
  if (!is_identical) {
    std::cout << "Incorret values" << std::endl;
  }
  return true;
}


template <typename ExecSpace, typename crsMat_t, typename crsMat_t2 , typename crsMat_t3 , typename TempMemSpace , typename PersistentMemSpace >
crsMat_t3 run_experiment(
    crsMat_t crsMat, crsMat_t2 crsMat2,
    int algorithm, int repeat, int chunk_size ,int multi_color_scale, int shmemsize, int team_size, int use_dynamic_scheduling, int verbose){

  spgemm_step++;
  typedef typename crsMat_t3::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t3::StaticCrsGraphType::row_map_type::non_const_type lno_view_t;
  typedef typename crsMat_t3::StaticCrsGraphType::entries_type::non_const_type lno_nnz_view_t;

  lno_view_t row_mapC;
  lno_nnz_view_t entriesC;
  scalar_view_t valuesC;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle
      <lno_view_t,lno_nnz_view_t, scalar_view_t,
      ExecSpace, TempMemSpace,PersistentMemSpace > KernelHandle;

  KernelHandle kh;
  kh.set_team_work_size(chunk_size);
  kh.set_shmem_size(shmemsize);
  kh.set_suggested_team_size(team_size);
  kh.set_suggested_vector_size(vector_size);

  if (use_dynamic_scheduling){
    kh.set_dynamic_scheduling(true);
  }
  if (verbose){
    kh.set_verbose(true);
  }

  const idx m = crsMat.numRows();
  const idx n = crsMat2.numRows();
  const idx k = crsMat2.numCols();

  std::cout << "m:" << m << " n:" << n << " k:" << k << std::endl;
  if (n != crsMat.numCols()){
    std::cout << "crsMat.numCols():" << crsMat.numCols() << " crsMat2.numRows():" << crsMat2.numRows() << std::endl;
    exit(1);
  }

  lno_view_t row_mapC_ref;
  lno_nnz_view_t entriesC_ref;
  scalar_view_t valuesC_ref;
  crsMat_t3 Ccrsmat_ref;
  if (check_output)
  {
    std::cout << "Running a reference algorithm" << std::endl;
    row_mapC_ref = lno_view_t ("non_const_lnow_row", m + 1);
    entriesC_ref = lno_nnz_view_t ("");
    valuesC_ref = scalar_view_t ("");
    KernelHandle kh;
    kh.set_team_work_size(chunk_size);
    kh.set_shmem_size(shmemsize);
    kh.set_suggested_team_size(team_size);
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_MEMORY);

    if (use_dynamic_scheduling){
      kh.set_dynamic_scheduling(true);
    }

    Kokkos::Impl::Timer timer1;
    KokkosKernels::Experimental::Graph::spgemm_symbolic (
        &kh,
        m,
        n,
        k,
        crsMat.graph.row_map,
        crsMat.graph.entries,
        TRANPOSEFIRST,
        crsMat2.graph.row_map,
        crsMat2.graph.entries,
        TRANPOSESECOND,
        row_mapC_ref
    );

    ExecSpace::fence();
    double symbolic_time = timer1.seconds();

    Kokkos::Impl::Timer timer3;
    size_type c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
    if (c_nnz_size){
      entriesC_ref = lno_nnz_view_t (Kokkos::ViewAllocateWithoutInitializing("entriesC"), c_nnz_size);
      valuesC_ref = scalar_view_t (Kokkos::ViewAllocateWithoutInitializing("valuesC"), c_nnz_size);
    }

    KokkosKernels::Experimental::Graph::spgemm_numeric(
        &kh,
        m,
        n,
        k,
        crsMat.graph.row_map,
        crsMat.graph.entries,
        crsMat.values,
        TRANPOSEFIRST,

        crsMat2.graph.row_map,
        crsMat2.graph.entries,
        crsMat2.values,
        TRANPOSESECOND,
        row_mapC_ref,
        entriesC_ref,
        valuesC_ref
    );
    ExecSpace::fence();
    double numeric_time = timer3.seconds();

    typename crsMat_t3::StaticCrsGraphType static_graph (entriesC_ref, row_mapC_ref);
    crsMat_t3 Ccrsmat("CrsMatrixC", k, valuesC_ref, static_graph);
    Ccrsmat_ref = Ccrsmat;
  }


  switch (algorithm){
  case 1:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_MKL);
    break;
  case 2:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_CUSPARSE);
    break;
  case 3:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_CUSP);
    break;
  case 4:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_DEBUG);
    break;
  case 5:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_MKL2PHASE);
    kh.get_spgemm_handle()->set_mkl_sort_option(mkl_sort_option);
    break;
  case 6:
	  kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_MEMORY2);
	  break;
  case 7:
      kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_MEMORY);
      break;
  case 8:
      kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_SPEED);
      break;
  case 9:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_COLOR);
    break;

  case 10:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR);
    break;

  case 11:
      kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTICOLOR2);
      break;
  case 12:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_VIENNA);
    break;
  case 13:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_MEMSPEED);
    break;
  case 14:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_MULTIMEM);
    break;

  default:
    kh.create_spgemm_handle(KokkosKernels::Experimental::Graph::SPGEMM_KK_MEMORY);
    break;
  }

  kh.get_spgemm_handle()->set_multi_color_scale(multi_color_scale);
  kh.get_spgemm_handle()->mkl_keep_output = mkl_keep_output;
  kh.get_spgemm_handle()->mkl_convert_to_1base = false;
  kh.get_spgemm_handle()->set_read_write_cost_calc (calculate_read_write_cost);

  if (coloring_input_file){
	  kh.get_spgemm_handle()->coloring_input_file = std::string(&spgemm_step) + "_" + std::string(coloring_input_file);
  }
  if (coloring_output_file){
	  kh.get_spgemm_handle()->coloring_output_file = std::string(&spgemm_step) + "_" + std::string(coloring_output_file);
  }

  for (int i = 0; i < repeat; ++i){

    row_mapC = lno_view_t
              ("non_const_lnow_row",
                  m + 1);
    entriesC = lno_nnz_view_t ("");
    valuesC = scalar_view_t ("");

    Kokkos::Impl::Timer timer1;
    KokkosKernels::Experimental::Graph::spgemm_symbolic (
        &kh,
        m,
        n,
        k,
        crsMat.graph.row_map,
        crsMat.graph.entries,
        TRANPOSEFIRST,
        crsMat2.graph.row_map,
        crsMat2.graph.entries,
        TRANPOSESECOND,
        row_mapC
    );

    ExecSpace::fence();
    double symbolic_time = timer1.seconds();

    Kokkos::Impl::Timer timer3;
    size_type c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
    if (c_nnz_size){
      entriesC = lno_nnz_view_t (Kokkos::ViewAllocateWithoutInitializing("entriesC"), c_nnz_size);
      valuesC = scalar_view_t (Kokkos::ViewAllocateWithoutInitializing("valuesC"), c_nnz_size);
    }

    KokkosKernels::Experimental::Graph::spgemm_numeric(
        &kh,
        m,
        n,
        k,
        crsMat.graph.row_map,
        crsMat.graph.entries,
        crsMat.values,
        TRANPOSEFIRST,

        crsMat2.graph.row_map,
        crsMat2.graph.entries,
        crsMat2.values,
        TRANPOSESECOND,
        row_mapC,
        entriesC,
        valuesC
    );
    ExecSpace::fence();
    double numeric_time = timer3.seconds();

    std::cout
    << "mm_time:" << symbolic_time + numeric_time
    << " symbolic_time:" << symbolic_time
    << " numeric_time:" << numeric_time << std::endl;



  }

  std::cout << "row_mapC:" << row_mapC.dimension_0() << std::endl;
  std::cout << "entriesC:" << entriesC.dimension_0() << std::endl;
  std::cout << "valuesC:" << valuesC.dimension_0() << std::endl;
  KokkosKernels::Experimental::Util::print_1Dview(valuesC);
  KokkosKernels::Experimental::Util::print_1Dview(entriesC);


  typename crsMat_t3::StaticCrsGraphType static_graph (entriesC, row_mapC);
  crsMat_t3 Ccrsmat("CrsMatrixC", k, valuesC, static_graph);
  if (check_output){
    bool is_identical = is_same_matrix<crsMat_t3, typename crsMat_t3::device_type>(Ccrsmat_ref, Ccrsmat);
    if (!is_identical){
      std::cout << "Result is wrong." << std::endl;
      exit(1);
    }
  }


  return Ccrsmat;

}


