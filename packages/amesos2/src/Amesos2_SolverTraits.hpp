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
   */
  template <template <class,class> class ConcreteSolver,
	    typename Scalar>
  struct solver_supports_scalar {
    static const bool value;
  };

  /*
   * If a solver's supported_scalars type-list is nil, this
   * meta-function returns true always.  Otherwise, it checks whether
   * the given type is in the solver's supported_scalars list.
   */
  template <template <class,class> class ConcreteSolver, typename Scalar>
  const bool solver_supports_scalar<ConcreteSolver,Scalar>::value
  = Meta::if_then_else<Meta::is_same<typename solver_traits<ConcreteSolver>::supported_scalars, Meta::nil_t>::value,
		       Meta::true_type,
		       Meta::type_list_contains<
			 typename solver_traits<ConcreteSolver>::supported_scalars,
			 Scalar> >::type::value;


} // end namespace Amesos2

#endif	// AMESOS2_SOLVERTRAITS_HPP
