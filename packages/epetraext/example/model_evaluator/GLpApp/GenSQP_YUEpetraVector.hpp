/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
//@HEADER
*/

#ifndef GENSQP_YUEPETRAVECTOR_H
#define GENSQP_YUEPETRAVECTOR_H

#include "GenSQP_Vector.hpp"
#include "Epetra_MultiVector.h"


/** \class GenSQP::YUEpetraVector
    \brief The GenSQP::Vector / (y,u) Epetra_MultiVector adapter class.
      
    Holds a pointer to two Epetra_MultiVectors y_epetra_vec and u_epetra_vec and implements the member functions 
    of the GenSQP::Vector class.
    Common use: optimal control. The vector y_epetra_vec represents the state variables, the vector
    u_epetra_vec represents the control variables.
*/


namespace GenSQP {

class YUEpetraVector : public Vector {

private:
  
  Teuchos::RefCountPtr<Epetra_MultiVector>  y_epetra_vec_;
  Teuchos::RefCountPtr<Epetra_MultiVector>  u_epetra_vec_;

public:

  YUEpetraVector( const Teuchos::RefCountPtr<Epetra_MultiVector> &y_epetra_vec,
                  const Teuchos::RefCountPtr<Epetra_MultiVector> &u_epetra_vec );

  /** \name Overridden from Vector */
  //@{

  double innerProd( const Vector &x ) const;

  void linComb( const double &alpha, const Vector &x, const double &beta );

  void Scale( const double &alpha );

  void Set( const double &alpha );

  void Set( const double &alpha, const Vector &x );

  Teuchos::RefCountPtr<Vector> createVector() const;

  //@}
  
  /** Returns a reference counted pointer to the private y_epetra_vec data container ("state variables").
  */
  Teuchos::RefCountPtr<const Epetra_MultiVector> getYVector() const;

  /** Returns a reference counted pointer to the private u_epetra_vec data container ("control variables").
  */
  Teuchos::RefCountPtr<const Epetra_MultiVector> getUVector() const;


}; // class YUEpetraVector

} // namespace GenSQP

#endif
