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
  \file   Amesos2_MUMPS_decl.hpp
  \author Joshua Dennis Booth <jdbooth@sandia.gov>

  \brief  Amesos2 MUMPS declarations.
*/


#ifndef AMESOS2_MUMPS_DECL_HPP
#define AMESOS2_MUMPS_DECL_HPP

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"
#include "Amesos2_MUMPS_FunctionMap.hpp"


namespace Amesos2 {


/** \brief Amesos2 interface to the MUMPS package.
 *
 * See the \ref MUMPS_parameters "summary of MUMPS parameters"
 * supported by this MUMPS interface.
 *
 * Currently, special care is needed to build Trilinos with MUMPS
 * This is due to how Tribits deals with scalapack/blacs (outdated)
 * Therefore, the linking of blacs and scalapack needs to be done in the Trilinos_EXTRA_LINK_FILES CMake directive, e.g.,
 * -DTrilinos_EXTRA_LINK_FLAGS="-I/local/openmpi/mpif.h -lmpi -lmpiblacs -L/local/scalapack -lblas -llapack -lscalapack"
 * Additionally, ETI is best if ON, since MUMPS has limited supported types
 *
 *
 *
 * Serial MUMPS
 * Note: serial mumps is supported. 
 * However, one must link provided seqlib mpi that comes with mumps
 * -DTrilinos_EXTRA_LIN_FLAGS="-I/../libseq -L/.../libseq -lmpiseq
 *
 * \ingroup amesos2_solver_interfaces
 */
template <class Matrix,class Vector>
class MUMPS : public SolverCore<Amesos2::MUMPS, Matrix, Vector>
{
  friend class SolverCore<Amesos2::MUMPS,Matrix,Vector>; // Give our base access
                                                          // to our private
                                                          // implementation funcs
public:

  /// Name of this solver interface.
  static const char* name;      // declaration. Initialization outside.

  typedef MUMPS<Matrix,Vector>                                       type;
  typedef SolverCore<Amesos2::MUMPS,Matrix,Vector>             super_type;

  // Since typedef's are not inheritted, go grab them
  typedef typename super_type::scalar_type                      scalar_type;
  typedef typename super_type::local_ordinal_type        local_ordinal_type;
  typedef typename super_type::global_ordinal_type      global_ordinal_type;
  typedef typename super_type::global_size_type            global_size_type;

  typedef TypeMap<Amesos2::MUMPS,scalar_type>                    type_map;

  typedef typename type_map::type                                  slu_type;
  typedef typename type_map::magnitude_type                  magnitude_type;
  typedef typename type_map::MUMPS_STRUC_C                    MUMPS_STRUC_C;

  typedef FunctionMap<Amesos2::MUMPS,slu_type>               function_map;

  MUMPS(Teuchos::RCP<const Matrix> A,
          Teuchos::RCP<Vector>       X,
          Teuchos::RCP<const Vector> B);
   ~MUMPS( );


private:

  /**
   * \brief Performs pre-ordering on the matrix to increase efficiency.
   *
   *   Come back to add support to Amesos for preordering
 */
  int preOrdering_impl();


  int symbolicFactorization_impl();


  /**
   * \brief MUMPS specific numeric factorization
   *
   * \throw std::runtime_error MUMPS is not able to factor the matrix
   */
  int numericFactorization_impl();


  /**
   * \brief MUMPS specific solve.
   *
   * Uses the symbolic and numeric factorizations, along with the RHS
   * vector \c B to solve the sparse system of equations.  The
   * solution is placed in X.
   *
   * \throw std::runtime_error MUMPS is not able to solve the system.
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
   *Currently, the following MUMPS parameters/options are recognized and acted upon:
   *Please see MUMPS manual for details as parameters changed based on MUMPS version
   *<ul>
   *   <li> \c "ICNTL(1)" </li>
   *   <li> \c "ICNTL(2)" </li>
   *   <li> \c "ICNTL(3)" </li>
   *   <li> \c "ICNTL(4)" </li>
   *   <li> \c "ICNTL(6)" </li>
   *   <li> \c "ICNTL(9)" </li>
   *   <li> \c "ICNTL(11)" </li>
   *</ul>
   *
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


  int ConvertToTriplet();

  void MUMPS_ERROR() const;

#ifdef HAVE_MPI
  MPI_Comm MUMPSComm;
#endif


  bool MUMPS_MATRIX_LOAD;
  bool MUMPS_STRUCT;

  // The following Arrays are persisting storage arrays for A, X, and B
  /// Stores the values of the nonzero entries for MUMPS
  Teuchos::Array<slu_type> nzvals_;
  /// Stores the location in \c Ai_ and Aval_ that starts row j
  Teuchos::Array<local_ordinal_type> rowind_;
  /// Stores the row indices of the nonzero entries
  Teuchos::Array<local_ordinal_type> colptr_;

  /// Persisting 1D store for X
  mutable Teuchos::Array<slu_type> xvals_;  local_ordinal_type ldx_;
  /// Persisting 1D store for B
  mutable Teuchos::Array<slu_type> bvals_;  local_ordinal_type ldb_;

  mutable MUMPS_STRUC_C mumps_par;

  bool is_contiguous_;

};                              // End class MUMPS


template <class Matrix,class Vector>
class MUMPSNS : public SolverCore<Amesos2::MUMPS, Matrix, Vector>
{
  friend class SolverCore<Amesos2::MUMPS,Matrix,Vector>; // Give our base access
                                                          // to our private
                                                          // implementation funcs
public:


  MUMPSNS(Teuchos::RCP<const Matrix> A,
          Teuchos::RCP<Vector>       X,
                     Teuchos::RCP<const Vector> B);
  /*
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "This solver is not support for these types.");

  }

  */
   ~MUMPSNS( )
  {}


private:

  int preOrdering_impl()
  {
      
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "This solver is not support for these types.");

  }

  int symbolicFactorization_impl()
  {

    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "This solver is not support for these types.");

  }

  int numericFactorization_impl()
 {

    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "This solver is not support for these types.");

  }

  int solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                 const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
 {

    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "This solver is not support for these types.");

  }

  bool matrixShapeOK_impl() const
 {

    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "This solver is not support for these types.");

  }

  void setParameters_impl(
    const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
 {

    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "This solver is not support for these types.");

  }

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters_impl() const
 {

    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "This solver is not support for these types.");

  }


  bool loadA_impl(EPhase current_phase)
 {

    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "This solver is not support for these types.");

  }


  int ConvertToTriplet()
 {

    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "This solver is not support for these types.");

  }


  void MUMPS_ERROR() const
 {

    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "This solver is not support for these types.");

  }


#ifdef HAVE_MPI
  MPI_Comm MUMPSComm;
#endif

    bool is_contiguous_;

};



// Specialize solver_traits struct for MUMPS
// TODO
template <>
struct solver_traits<MUMPS> {
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

#endif  // AMESOS2_MUMPS_DECL_HPP
