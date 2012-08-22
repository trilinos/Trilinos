// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef __TSQR_Trilinos_TsqrTpetraTest_hpp
#define __TSQR_Trilinos_TsqrTpetraTest_hpp

#include "Teuchos_Tuple.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "Kokkos_SerialNode.hpp"
#ifdef HAVE_KOKKOSCLASSIC_TBB
#include "Kokkos_TBBNode.hpp"
#endif // HAVE_KOKKOSCLASSIC_TBB

#include "Teuchos_Time.hpp"
#include "AnasaziTpetraAdapter.hpp" // sic (not "-or")
#include "Tsqr_Random_NormalGenerator.hpp"

#include <sstream>
#include <stdexcept>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR { 
  namespace Trilinos { 
    namespace Test {
      using Teuchos::RCP;
      using Teuchos::Tuple;

      template< class S, class LO, class GO, class Node >
      class TpetraTsqrTest {
      public:
	typedef S scalar_type;
	typedef LO local_ordinal_type;
	typedef GO global_ordinal_type;
	typedef Node node_type;

	typedef typename TSQR::ScalarTraits< S >::magnitude_type magnitude_type;
	typedef TSQR::Trilinos::TsqrTpetraAdaptor< S, LO, GO, Node > adaptor_type;

	TpetraTsqrTest (const Tpetra::global_size_t nrowsGlobal,
			const size_t ncols,
			const Teuchos::RCP< const Teuchos::Comm<int> >& comm,
			const Teuchos::RCP< Node >& node,
			const Teuchos::ParameterList& params) :
	  results_ (magnitude_type(0), magnitude_type(0))
	{
	  using Teuchos::Tuple;
	  using Teuchos::Exceptions::InvalidParameter;
	  typedef typename adaptor_type::factor_output_type factor_output_type;
	  typedef Teuchos::SerialDenseMatrix< LO, S > matrix_type;

	  bool contiguousCacheBlocks = false;
	  try {
	    contiguousCacheBlocks = params.get<bool>("contiguousCacheBlocks");
	  } catch (InvalidParameter&) {
	    contiguousCacheBlocks = false;
	  }

	  triple_type testProblem = 
	    makeTestProblem (nrowsGlobal, ncols, comm, node, params);
	  // A is already filled in with the test problem.  A_copy and
	  // Q are just multivectors with the same layout as A.
	  // A_copy will be used for temporary storage, and Q will
	  // store the (explicitly represented) Q factor output.  R
	  // will store the R factor output.
	  RCP< MV > A = testProblem[0]; 
	  RCP< MV > A_copy = testProblem[1];
	  RCP< MV > Q = testProblem[2];
	  matrix_type R (ncols, ncols);

	  // Adaptor uses one of the multivectors only to reference
	  // the underlying communicator object.
	  adaptor_type adaptor (*A, params);
	  if (contiguousCacheBlocks)
	    adaptor.cacheBlock (*A, *A_copy);

	  factor_output_type factorOutput = 
	    adaptor.factor (*A_copy, R, contiguousCacheBlocks);
	  adaptor.explicitQ (*A_copy, factorOutput, *Q, contiguousCacheBlocks);
	  if (contiguousCacheBlocks)
	    {
	      // Use A_copy as temporary storage for un-cache-blocking
	      // Q.  Tpetra::MultiVector objects copy deeply.
	      *A_copy = *Q;
	      adaptor.unCacheBlock (*A_copy, *Q);
	    }
	  results_ = adaptor.verify (*A, *Q, R);
	}

	/// Return the residual error and the departure from
	/// orthogonality of the factorization output.
	std::pair< magnitude_type, magnitude_type > 
	getResults() const 
	{ 
	  return results_;
	}

      private:
	typedef Tpetra::MultiVector< S, LO, GO, Node >   MV;
	typedef Teuchos::Tuple< RCP< MV >, 3 >           triple_type;
	typedef Teuchos::RCP< const Teuchos::Comm<int> > comm_ptr;
	typedef Tpetra::Map< LO, GO, Node >              map_type;
	typedef Teuchos::RCP< const map_type >           map_ptr;
	typedef TSQR::Random::NormalGenerator< LO, S >   normalgen_type;
	typedef Teuchos::RCP< Node >                     node_ptr;

	/// Results of the factorization: residual error, and
	/// departure from orthgonality, each in the Frobenius norm,
	/// and each scaled by the Frobenius norm of the original
	/// matrix A.
	std::pair< magnitude_type, magnitude_type > results_;

	/// \brief Make a map for creating Tpetra test multivectors.
	///
	/// \param nrowsGlobal [in] Number of rows in the entire MultiVector
	/// \param comm [in] Communications handler object
	/// \param node [in] Kokkos node object
	///
	/// \return The desired map object
	///
	static map_ptr
	makeMap (const Tpetra::global_size_t nrowsGlobal,
		 const comm_ptr& comm,
		 const node_ptr& node)
	{
	  using Tpetra::createUniformContigMapWithNode;
	  return createUniformContigMapWithNode< LO, GO, Node > (nrowsGlobal, 
								 comm, node);
	}

	/// \brief Make a Tpetra test multivector for filling in.
	///
	static RCP< MV >
	makeMultiVector (const map_ptr& map,
			 const size_t ncols)
	{
	  // "false" means "don't fill with zeros"; we'll fill it in
	  // fillTpetraMultiVector().
	  return Teuchos::rcp (new MV (map, ncols, false));
	}

	/// \brief Fill in a given Tpetra test multivector with random values.
	///
	static void
	fillMultiVector (const RCP< MV >& mv,
			 const RCP< normalgen_type >& pGen)
	{
	  using TSQR::Trilinos::TpetraRandomizer;
	  typedef TpetraRandomizer< S, LO, GO, Node, normalgen_type > randomizer_type;

	  const LO ncols = mv->getNumVectors();
	  std::vector< S > singular_values (ncols);
	  if (ncols > 0)
	    {
	      singular_values[0] = S(1);
	      for (LO k = 1; k < ncols; ++k)
		singular_values[k] = singular_values[k-1] / S(2);
	    }
	  randomizer_type randomizer (*mv, pGen);
	  randomizer.randomMultiVector (*mv, &singular_values[0]);
	}

	static triple_type
	makeTestProblem (const Tpetra::global_size_t nrowsGlobal,
			 const size_t ncols,
			 const comm_ptr& comm,
			 const node_ptr& node,
			 const Teuchos::ParameterList& params)
	{
	  using TSQR::Trilinos::TpetraMessenger;
	  using TSQR::MessengerBase;

	  map_ptr map = makeMap (nrowsGlobal, comm, node);
	  RCP< MV > A = makeMultiVector (map, ncols);
	  RCP< MV > A_copy = makeMultiVector (map, ncols);
	  RCP< MV > Q = makeMultiVector (map, ncols);

	  // Fill A with the random test problem
	  RCP< normalgen_type > pGen (new normalgen_type);
	  fillMultiVector (A, pGen);

	  return Teuchos::tuple (A, A_copy, Q);
	}
      };

    } // namespace Test
  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrTpetraTest_hpp
