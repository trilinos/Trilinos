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

#ifndef NOX_EPETRA_SCHURINTERFACE_H
#define NOX_EPETRA_SCHURINTERFACE_H

#include "NOX_Epetra.H"

//-----------------------------------------------------------------------------

// This is a new class that may eventually get moved into NOX::Epetra.  For now,
// this is simply used as a testbed for driving new coupling algorithms - RWH 10/23/2006

namespace NOX {

namespace Epetra {

class SchurInterface
{
public:

  SchurInterface() { };

  virtual ~SchurInterface() { };

  //! Apply the action of a matrix block on a vector
  virtual bool applyBlockAction( int rowBlock, int colBlock, const Epetra_Vector & x, Epetra_Vector & result ) = 0;

  //! Apply the action of a matrix block inverse on a vector
  virtual bool applyBlockInverseAction( int rowBlock, int colBlock, const Epetra_Vector & x, Epetra_Vector & result ) = 0;

  //! Get the communicator
  virtual const Epetra_Comm & getComm() = 0;

  //! Query availability of explicit operator
  virtual bool hasExplicitOperator( int rowBlock, int colBlock )
    { return false; }

  //! Access explicit operator
  virtual Teuchos::RCP<Epetra_Operator> getExplicitOperator( int rowBlock, int colBlock )
    { return Teuchos::null; }

};
} // Epetra
} // NOX

#endif
