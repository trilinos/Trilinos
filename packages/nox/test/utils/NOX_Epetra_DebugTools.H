//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_EPETRA_DEBUGTOOLS_H
#define NOX_EPETRA_DEBUGTOOLS_H

#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_TestCompare.H"

//----------------------------------------------------------------------------//
//----------------- DebugTools::ObjectOutput ----------------------------------//
//----------------------------------------------------------------------------//
//
// This class allows linear solver tolerance criteria to be used as a NOX
// convergence test.

namespace NOX {
  namespace Epetra {

class DebugTools
{

 public:

  enum NOX_EPETRA_OP_TYPE { MATRIX_FREE                 ,
                            FINITE_DIFFERENCE           ,
                            FINITE_DIFFERENCE_COLORNG   ,
                            CRS_MATRIX                  ,
                            NONE                           };

  enum FORMAT_TYPE        { ASCII         = 0           ,
                            MATRIX_MARKET = 1              };

  static NOX_EPETRA_OP_TYPE get_operator_type( const Epetra_Operator *                     );
  static NOX_EPETRA_OP_TYPE get_operator_type( const Teuchos::RCP<Epetra_Operator> );

  // Factory-like method to create an operator - needed RWH ?
  //static Epetra_Operator * create_nonlinear_operator( NOX_EPETRA_OP_TYPE );

  // Method to compute a matrix from an operator acting on unit Cartesian vect// spanning the space of problem dofs
  static Epetra_CrsMatrix * compute_matrix_using_operator( const Epetra_Operator * );

#ifdef HAVE_NOX_EPETRAEXT
  // I/O methods via EpetraExt
  static int readVector    ( std::string, const Epetra_Comm &, Epetra_Vector*&                    );
  static void writeVector  ( std::string, const Epetra_Vector &    , FORMAT_TYPE outFormat= ASCII, bool writeMap = true );
  static void writeMatrix  ( std::string, const Epetra_RowMatrix & , FORMAT_TYPE outFormat= ASCII );
  static void writeOperator( std::string, const Epetra_Operator &  , FORMAT_TYPE outFormat= ASCII );
  static void readMatrix   ( std::string,       Epetra_CrsMatrix*&                                );
#endif

};

class TestCompare : public NOX::TestCompare
{

 public:

  TestCompare(std::ostream& outputStream, const NOX::Utils& utilities);

  virtual ~TestCompare() {}

  // Method to compute a matrix from an operator acting on unit Cartesian vect// spanning the space of problem dofs

    //! Compare two Epetra_CrsMatrix's
    /*!
     * @param \c mat [in] Computed vector
     * @param \c mat_expected [in] Vector to be compared to
     * @param \c rtol [in] Relative tolerance \f$\epsilon_r\f$
     * @param \c atol [in] Absolute tolerance \f$\epsilon_a\f$
     * @param \c name [in] Name for this test when printing details
     * @return Returns 0 if
     * \f[
     *     \frac{|mat_ij - mat_expected_ij|}{\epsilon_a + \epsilon_r|mat_expected_ij|} < 1
     * \f]
     * for each $i$, $i$, and 1 otherwise.
     */
    virtual int testCrsMatrices(const Epetra_CrsMatrix & mat,
                                const Epetra_CrsMatrix & mat_expected,
                                double rtol, double atol, const std::string& name,
                                bool enforceStructure = true );

};

} // Epetra namespace
} // NOX namespace

#endif
