#ifndef __TSQR_Test_generateStack_hpp
#define __TSQR_Test_generateStack_hpp

#include <Tsqr_Blas.hpp>
#include <Tsqr_Lapack.hpp>
#include <Tsqr_Matrix.hpp>
#include <Tsqr_Util.hpp>
#include <Tsqr_Random_MatrixGenerator.hpp>
#include <Tsqr_RMessenger.hpp>

#include <algorithm>
#include <functional>
#include <sstream>
#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Test {

    /// Generate a (pseudo)random test problem consisting of nprocs
    /// different ncols by ncols upper triangular matrices, stacked
    /// vertically.
    ///
    /// \param generator [in/out] (Pseudo)random number generator,
    ///   that generates according to a normal(0,1) distribution.
    /// \param A_global [out] Matrix to fill.  Will set dimensions.
    /// \param singular_values [in] ncols singular values to use
    template< class Ordinal, class Scalar, class Generator >
    static void
    generateStack (Generator& generator,
		   Matrix< Ordinal, Scalar >& A_global,
		   const typename ScalarTraits< Scalar >::magnitude_type singular_values[],
		   const int nprocs,
		   const Ordinal ncols)
    {
      TSQR::Random::MatrixGenerator< Ordinal, Scalar, Generator > matGen (generator);
      const Ordinal nrows = nprocs * ncols;
      A_global.reshape (nrows, ncols);
      A_global.fill (Scalar(0));

      for (int p = 0; p < nprocs; ++p)
	{
	  Scalar* const curptr = A_global.get() + p*ncols;
	  MatView< Ordinal, Scalar > R_cur (ncols, ncols, curptr, nrows);
	  matGen.fill_random_R (ncols, R_cur.get(), nrows, singular_values);
	}
    }

    /// \brief Test problem for distributed-memory part of TSQR
    ///
    /// \param A_local [out] ncols by ncols upper triangular matrix
    ///   (on each MPI process)
    /// 
    /// \param A_global [out] Empty on all procs but Proc 0, where it
    ///   starts empty (not required, but this is more efficient) and
    ///   gets resized here
    ///
    /// \param ncols [in] Number of columns in the matrix to generate.
    ///   Number of rows in the matrix will be (# MPI processes)*ncols.
    ///
    /// \param generator [in/out] Normal(0,1) (pseudo)random number
    ///   generator
    ///
    /// \param messenger [in/out] MPI communicator object for Scalars
    ///
    template< class Ordinal, class Scalar, class Generator >
    void
    par_tsqr_test_problem (Generator& generator,
			   Matrix< Ordinal, Scalar >& A_local,
			   Matrix< Ordinal, Scalar >& A_global, 
			   const Ordinal ncols,
			   MessengerBase< Scalar >* const messenger)
    {
      const int nprocs = messenger->size();
      const int my_rank = messenger->rank();
      A_local.reshape (ncols, ncols);

      if (my_rank == 0)
	{
	  typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
	  typedef typename std::vector< magnitude_type >::size_type index_type;

	  std::vector< magnitude_type > singular_values (ncols);
	  singular_values[0] = magnitude_type(1);
	  for (index_type k = 1; k < ncols; ++k)
	    singular_values[k] = singular_values[k-1] / magnitude_type(2);

	  generateStack (generator, A_global, &singular_values[0], nprocs, ncols);
	  scatterStack (A_global, A_local, messenger);
	}
      else
	scatterStack (A_global, A_local, messenger);
    }


  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_generateStack_hpp
