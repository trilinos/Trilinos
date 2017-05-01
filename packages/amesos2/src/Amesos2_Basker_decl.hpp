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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/**
  \file   Amesos2_Basker_decl.hpp
  \author Joshua Dennis Booth <jdbooth@sandia.gov>
          Siva Rajamanickam <srajama@sandia.gov>

  \brief  Amesos2 Basker declarations.
*/


#ifndef AMESOS2_BASKER_DECL_HPP
#define AMESOS2_BASKER_DECL_HPP

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"
#include "Amesos2_Basker_FunctionMap.hpp"

//Note:  We got an error while being a class variable and mutable.  Need to comeback and fix!!


namespace Amesos2 {

/** \brief Amesos2 interface to the Baker package.
 *
 * See the \ref Basker_parameters "summary of Basker parameters"
 * supported by this Amesos2 interface.
 *
 * \ingroup amesos2_solver_interfaces
 */
template <class Matrix,class Vector>
class Basker : public SolverCore<Amesos2::Basker, Matrix, Vector>
{
  friend class SolverCore<Amesos2::Basker,Matrix,Vector>; // Give our base access
                                                          // to our private
                                                          // implementation funcs
public:

  /// Name of this solver interface.
  static const char* name;      // declaration. Initialization outside.


  typedef Basker<Matrix,Vector>                                        type;

  typedef SolverCore<Amesos2::Basker,Matrix,Vector>              super_type;

  // Since typedef's are not inheritted, go grab them
  typedef typename super_type::scalar_type                      scalar_type;
  typedef typename super_type::local_ordinal_type        local_ordinal_type;
  typedef typename super_type::global_ordinal_type      global_ordinal_type;
  typedef typename super_type::global_size_type            global_size_type;
  typedef typename super_type::node_type                          node_type;

  typedef TypeMap<Amesos2::Basker,scalar_type>                     type_map;

  typedef typename type_map::type                                  slu_type;
  typedef typename type_map::magnitude_type                  magnitude_type;

  typedef FunctionMap<Amesos2::Basker,slu_type>                function_map;

  typedef Matrix                                                matrix_type;


  Basker( Teuchos::RCP<const Matrix> A,
          Teuchos::RCP<Vector>       X,
          Teuchos::RCP<const Vector> B);
   ~Basker( );


private:

  /**
   * \brief Performs pre-ordering on the matrix to increase efficiency.
   *
   *   Come back to add support to Amesos for preordering
 */
  int preOrdering_impl();


  int symbolicFactorization_impl();


  /**
   * \brief Basker specific numeric factorization
   *
   * \throw std::runtime_error Basker is not able to factor the matrix
   */
  int numericFactorization_impl();


  /**
   * \brief Basker specific solve.
   *
   * Uses the symbolic and numeric factorizations, along with the RHS
   * vector \c B to solve the sparse system of equations.  The
   * solution is placed in X.
   *
   * \throw std::runtime_error Basker is not able to solve the system.
   *
   * \callgraph
   */
  int solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                 const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const;


  /**
   * \brief Determines whether the shape of the matrix is OK for this solver.
   */
  bool matrixShapeOK_impl() const;


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


  // Members
  int num_threads;

  // The following Arrays are persisting storage arrays for A, X, and B
  /// Stores the values of the nonzero entries for Basker
  Teuchos::Array<slu_type> nzvals_;
  /// Stores the location in \c Ai_ and Aval_ that starts row j
  Teuchos::Array<local_ordinal_type> rowind_;
  /// Stores the row indices of the nonzero entries
  Teuchos::Array<local_ordinal_type> colptr_;

  bool is_contiguous_;


  /// Persisting 1D store for X
  mutable Teuchos::Array<slu_type> xvals_;  local_ordinal_type ldx_;
  /// Persisting 1D store for B
  mutable Teuchos::Array<slu_type> bvals_;  local_ordinal_type ldb_;

    /*Handle for Basker object*/
 
#ifdef SHYLUBASKER
#ifdef HAVE_AMESOS2_KOKKOS
#ifdef KOKKOS_HAVE_OPENMP
  /*
  typedef typename node_type::device_type  kokkos_device;
  typedef typename kokkos_device::execution_space kokkos_exe;
  static_assert(std::is_same<kokkos_exe,Kokkos::OpenMP>::value,
  "Kokkos node type not support by experimental Basker Amesos2");
  */
  typedef Kokkos::OpenMP Exe_Space;
#elif defined(KOKKOS_HAVE_SERIAL)
  typedef Kokkos::Serial Exe_Space;
#else
#pragma message("Kokkos Node type not supported by Basker") 
#endif // OpenMP vs Serial
   ::BaskerNS::BaskerTrilinosInterface<local_ordinal_type, slu_type, Exe_Space>
       *basker;
#else
  #pragma message("HAVE SHYLUBASKER AND NOT KOKKOS! ERROR")
#endif
#else
   mutable ::Basker::Basker<local_ordinal_type,slu_type> basker;
 #endif


};                              // End class Basker


// Specialize solver_traits struct for Basker
// TODO
template <>
struct solver_traits<Basker> {
#ifdef HAVE_TEUCHOS_COMPLEX
  typedef Meta::make_list4<float,
                           double,
                           std::complex<float>,
                           std::complex<double> > supported_scalars;
#else
  typedef Meta::make_list2<float, double> supported_scalars;
#endif
};

} // end namespace Amesos2

#endif  // AMESOS2_BASKER_DECL_HPP
