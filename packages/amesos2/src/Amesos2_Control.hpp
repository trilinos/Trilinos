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
  \file   Amesos2_Control.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Sat Jan 16 08:52:02 2010

  \brief  Container class for control variables.
*/
#ifndef AMESOS2_CONTROL_HPP
#define AMESOS2_CONTROL_HPP

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

namespace Amesos2 {


struct Control {
  /// Default constructor.
  Control()
    : verbose_(0)
    , debug_(0)
    , useTranspose_(false)
    , addToDiag_("0.0")
    , addZeroToDiag_(false)
    , matrixProperty_(0)
    , rcond_threshold_(1e-12)
    , refactorize_(false)
    , scaleMethod_(0)
    , reindex_(false)
    { }

  /// Default destructor.
  ~Control() { };

  void setControlParameters(
    const Teuchos::RCP<Teuchos::ParameterList> & parameterList );

  /** \brief Sets the verbosity level.
   *
   * \internal Really should implement the Teuchos::VerboseObject
   * interface for this sort of thing.
   */ 
  int verbose_;

  /// Sets the level of debug output
  int debug_;


  /// When solving, use A^T instead of A
  bool useTranspose_;


  /**
   * \brief Add this value to the diagonal.
   *
   * This value is declared as a string, so that it may be later appropriately
   * interpreted as the scalar type of the matrix being manipulated.
   */
  std::string addToDiag_;


  /**
   * \brief Adds zero to diagonal of redistributed matrix.
   *
   * (some solvers choke on a matrix with a partly empty diag)
   */
  bool addZeroToDiag_;


  /**
   * \brief  Set the matrix property.
   *
   * Matrix property can be
   * - 0 : general asymmetric matrix;
   * - 1 : SPD;
   * - 2 : general symmetric matrix.
   * UNUSED - See bug #2331 and bug #2332
   */
  int matrixProperty_;


  /// If error is greater than \c this value, perform symbolic and
  ///  numeric factorization with full partial pivoting
  double rcond_threshold_; // if we refactorize, the factorization
  // may suffer in numeric quality.  We
  // compute rcond = min (abs (diag (U))) /
  // max (abs(diag (U))).  If this ratio is
  // <= rcond_threshold_, then the
  // "refactorization" is scrapped, and we
  // factor with full partial pivoting
  // instead.


  bool refactorize_; // if true, and if the Symbolic and Numeric
  // objects have already been created, then
  // attempt to "refactorize" (factor the matrix
  // with no changes to the pivot order since the
  // last call the klu_btf_factor).


  int scaleMethod_; // Most methods (klu, UMFPACK, Mumps, ...) can
  // scale the input matrix prior to
  // factorization.  This can improve pivoting,
  // reduces fill-in, and leads to a better
  // quality factorization.  The options are:
  // 0: no scaling
  // 1: use the default method for the specific
  // package
  // 2: use the method's 1st alternative (if it
  // has one)
  // 3: use the method's 2nd alternative, and so
  // on.
  //
  // Amesos2_Klu is, at present, the only code
  // which implements this


  /**
   * \brief Reindex A to standard indexing.
   *
   * If true, the Amesos2 class should reindex the matrix to standard
   * indexing (i.e. 0-(n-1)) At present, only Amesos2::Klu supports this
   * option.
   */
  bool reindex_;

};                              // end class Control


} // end namespace Amesos2

#endif	// AMESOS2_CONTROL_HPP
