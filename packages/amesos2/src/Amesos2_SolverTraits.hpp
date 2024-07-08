// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/**
 * \file   Amesos2_SolverTraits.hpp
 * \author Eric Bavier <etbavie@muenster.srn.sandia.gov>
 * \date   Mon Jul 18 16:10:34 2011
 * 
 * \brief  Provides access to interesting solver traits
 * 
 * 
 */

#ifndef AMESOS2_SOLVERTRAITS_HPP
#define AMESOS2_SOLVERTRAITS_HPP

#include "Amesos2_Meta.hpp"
#include "Amesos2_MatrixAdapter.hpp"

namespace Amesos2 {

  /** \internal
   *  
   * \brief Provides traits about solvers.
   *
   * \internal Each concrete solver interface should specialize this
   * struct for that particular solver.
   */
  template <template <class,class> class ConcreteSolver>
  struct solver_traits {
    typedef Meta::nil_t supported_scalars;
  };

  
  ////////////////////////////
  // Related meta-functions //
  ////////////////////////////
  
  /** \internal
   *
   * \brief Check whether a solver supports a scalar type
   * 
   * A meta-function for checking a solver's support for a scalar.
   * The special-case is a Meta::nil_t list, which we are going to
   * interpret as the solver supporting any and all types.
   * If a solver's supported_scalars type-list is nil, this
   * meta-function returns true always.  Otherwise, it checks whether
   * the given type is in the solver's supported_scalars list.
   */

    /* SR: We will not use external initialization for the static const types.
     * Combined with template meta programming this fails in Intel compilers
     * 11-13. Moving all the initializations inside the declarations.
     */
  template <template <class,class> class ConcreteSolver,
	    typename Scalar>
  struct solver_supports_scalar {
    static const bool value =
              std::conditional_t<std::is_same_v<typename solver_traits<ConcreteSolver>::supported_scalars, Meta::nil_t>,
		       std::true_type,
		       Meta::type_list_contains<
			 typename solver_traits<ConcreteSolver>::supported_scalars,
			 Scalar> >::value;
  };

  template <template <class,class> class ConcreteSolver,
      typename Matrix>
  struct solver_supports_matrix {
    static const bool value = true;
  };

  // for kokkos adapter we only allow this for the specific solvers which
  // are using it. This is to avoid having ETI setup for all solvers. The
  // kokkos adapter is for testing UVM off and would become relic when Tpetra
  // switches to UVM off. To support this, solvers like Tacho override this
  // method and return true.
  template <template <class,class> class ConcreteSolver,
    typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  struct solver_supports_matrix<ConcreteSolver,
    KokkosSparse::CrsMatrix<Scalar, LocalOrdinal, ExecutionSpace>> {
    static const bool value = false;
  };

} // end namespace Amesos2

#endif	// AMESOS2_SOLVERTRAITS_HPP
