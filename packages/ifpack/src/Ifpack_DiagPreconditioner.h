/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK_DIAG_PRECONDITIONER_H
#define IFPACK_DIAG_PRECONDITIONER_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_ConfigDefs.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"
class Epetra_BlockMap;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Comm;

//! Ifpack_DiagPreconditioner: a class for diagonal preconditioning.
/*
Ifpack_DiagPreconditioner: a class to wrap a vector as diagonal preconditioner. The preconditioner is simply defined by
\f[
z_i = D_i r_i,
\f]
where \f$r,z\f$ are the vector to be preconditioned and the preconditioned vector, and \f$D_i\f$ is the i-th element of the scaling vector.

\author Marzio Sala, ETHZ/D-INFK

\date Last updated on 17-Apr-06

 */
class Ifpack_DiagPreconditioner : public Epetra_Operator
{
  public:

    //! ctor
    Ifpack_DiagPreconditioner(const Epetra_Map& DomainMap,
                              const Epetra_Map& RangeMap,
                              const Epetra_Vector& diag);

    //! dtor
    ~Ifpack_DiagPreconditioner();

    int SetUseTranspose(bool UseTranspose_in)
    {
      UseTranspose_ = UseTranspose_in;
      return(0);
    }

    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    double NormInf() const
    {
      return(-1.0);
    }

    const char* Label() const
    {
      return("Ifpack_DiagPreconditioner");
    }

    bool UseTranspose() const
    {
      return(UseTranspose_);
    }

    bool HasNormInf() const
    {
      return(false);
    }

    const Epetra_Comm& Comm() const
    {
      return(diag_.Comm());
    }

    const Epetra_Map& OperatorDomainMap() const
    {
      return(RangeMap_);
    }

    const Epetra_Map& OperatorRangeMap() const
    {
      return(DomainMap_);
    }

    const Epetra_BlockMap& Map() const
    {
      return(diag_.Map());
    }

  private:
    bool UseTranspose_;
    const Epetra_Map& DomainMap_;
    const Epetra_Map& RangeMap_;
    const Epetra_Vector& diag_;
};

#endif
