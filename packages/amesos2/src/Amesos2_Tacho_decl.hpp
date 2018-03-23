// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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
// Questions? Contact Sivasankaran Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef AMESOS2_TACHO_DECL_HPP
#define AMESOS2_TACHO_DECL_HPP

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"
#include "Amesos2_Tacho_FunctionMap.hpp"

#include "Tacho.hpp"
#include "Tacho_Solver.hpp"

namespace Amesos2 {

/** \brief Amesos2 interface to the Tacho package.
 *
 * See the \ref tacho_parameters "summary of Tacho parameters"
 * supported by this Amesos2 interface.
 *
 * \ingroup amesos2_solver_interfaces
 */
template <class Matrix,
          class Vector>
class TachoSolver : public SolverCore<Amesos2::TachoSolver, Matrix, Vector>
{
  friend class SolverCore<Amesos2::TachoSolver,Matrix,Vector>; // Give our base access
                                                          // to our private
                                                          // implementation funcs
public:

  /// Name of this solver interface.
  static const char* name;      // declaration. Initialization outside.

  typedef TachoSolver<Matrix,Vector>                                   type;
  typedef SolverCore<Amesos2::TachoSolver,Matrix,Vector>         super_type;

  // Since typedef's are not inheritted, go grab them
  typedef typename super_type::scalar_type                      scalar_type;
  typedef typename super_type::local_ordinal_type        local_ordinal_type;
  typedef typename super_type::global_size_type            global_size_type;

  typedef TypeMap<Amesos2::TachoSolver,scalar_type>                type_map;

  /*
   * The Tacho interface will need two other typedef's, which are:
   * - the tacho type that corresponds to scalar_type and
   * - the corresponding type to use for magnitude
   */
  typedef typename type_map::type                                tacho_type;
  typedef typename type_map::magnitude_type                  magnitude_type;

  typedef FunctionMap<Amesos2::TachoSolver,tacho_type>         function_map;

  // TODO - Not sure yet best place for organizing these typedefs
  typedef Tacho::ordinal_type                                  ordinal_type;
  typedef Tacho::size_type                                        size_type;
  typedef Kokkos::DefaultHostExecutionSpace                   HostSpaceType;
  typedef Kokkos::View<size_type*,HostSpaceType>            size_type_array;
  typedef Kokkos::View<ordinal_type*,HostSpaceType>      ordinal_type_array;
  typedef Kokkos::View<tacho_type*, HostSpaceType>         value_type_array;

  /// \name Constructor/Destructor methods
  //@{

  /**
   * \brief Initialize from Teuchos::RCP.
   *
   * \warning Should not be called directly!  Use instead
   * Amesos2::create() to initialize a Tacho interface.
   */
  TachoSolver(Teuchos::RCP<const Matrix> A,
          Teuchos::RCP<Vector>       X,
          Teuchos::RCP<const Vector> B);


  /// Destructor
  ~TachoSolver( );

  //@}

  /// Returns a short description of this Solver
  std::string description() const;

private:

  /**
   * \brief Performs pre-ordering on the matrix to increase efficiency.
   *
   * Tacho does not support pre-ordering, so this method does nothing.
   */
  int preOrdering_impl();


  /**
   * \brief Perform symbolic factorization of the matrix using Tacho.
   *
   * Called first in the sequence before numericFactorization.
   *
   * \throw std::runtime_error Tacho is not able to factor the matrix.
   */
  int symbolicFactorization_impl();


  /**
   * \brief Tacho specific numeric factorization
   *
   * \throw std::runtime_error Tacho is not able to factor the matrix
   */
  int numericFactorization_impl();

  /**
   * \brief Tacho specific solve.
   *
   * Uses the symbolic and numeric factorizations, along with the RHS
   * vector \c B to solve the sparse system of equations.  The
   * solution is placed in X.
   *
   * \throw std::runtime_error Tacho is not able to solve the system.
   *
   * \callgraph
   */
  int solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                 const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const;


  /**
   * \brief Determines whether the shape of the matrix is OK for this solver.
   */
  bool matrixShapeOK_impl() const;


  /**
   */
  void setParameters_impl(
    const Teuchos::RCP<Teuchos::ParameterList> & parameterList );


  /**
   * Hooked in by Amesos2::SolverCore parent class.
   *
   * \return a const Teuchos::ParameterList of all valid parameters for this
   * solver.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters_impl() const;


  /**
   * \brief Reads matrix data into internal structures
   *
   * \param [in] current_phase an indication of which solution phase this
   *                           load is being performed for.
   *
   * \return \c true if the matrix was loaded, \c false if not
   */
  bool loadA_impl(EPhase current_phase);


  /**
   * \brief can we optimize size_type and ordinal_type for straight pass through
   */
  bool do_optimization() const;

  // struct holds all data necessary to make a tacho factorization or solve call
  mutable struct TACHOData {
    typename Tacho::Solver<tacho_type,HostSpaceType> solver;

    // TODO: Implement the paramter options - confirm which we want and which have been implemented
    // int num_kokkos_threads;
    // int max_num_superblocks;
  } data_;

  // The following Arrays are persisting storage arrays for A, X, and B
  /// Stores the values of the nonzero entries for Tacho
  Teuchos::Array<tacho_type> nzvals_;
  /// Stores the location in \c Ai_ and Aval_ that starts row j
  Teuchos::Array<ordinal_type> colind_;
  /// Stores the row indices of the nonzero entries
  Teuchos::Array<size_type> rowptr_;

  // TODO: Decide handling for CUDA - how to fail
#ifdef KOKKOS_HAVE_OPENMP
  typedef Kokkos::OpenMP DeviceSpaceType;
#else
  typedef Kokkos::Serial DeviceSpaceType;
#endif
  typedef typename Tacho::Solver<tacho_type,DeviceSpaceType>::value_type_matrix
    solve_array_t;

  // used as an internal workspace - possibly we can store this better in TACHOData
  mutable solve_array_t workspace_;
};                              // End class Tacho


// Specialize solver_traits struct for Tacho
template <>
struct solver_traits<TachoSolver> {
#ifdef HAVE_TEUCHOS_COMPLEX
  typedef Meta::make_list4<float,
                           double,
                           std::complex<float>,
                           std::complex<double>
                           >supported_scalars;
#else
  typedef Meta::make_list2<float,
                           double
                           >supported_scalars;
#endif
};

} // end namespace Amesos2

#endif  // AMESOS2_TACHO_DECL_HPP
