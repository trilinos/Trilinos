/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef TIFPACK_DIAG_PRECONDITIONER_HPP
#define TIFPACK_DIAG_PRECONDITIONER_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_Vector.hpp"
class Tpetra_BlockMap;
class Tpetra_Map;
class Tpetra_MultiVector;
class Tpetra_Comm;

using namespace std;

//! Tifpack_DiagPreconditioner: a class for diagonal preconditioning.
/*
Tifpack_DiagPreconditioner: a class to wrap a vector as diagonal preconditioner. The preconditioner is simply defined by
\f[
z_i = D_i r_i,
\f]
where \f$r,z\f$ are the vector to be preconditioned and the preconditioned vector, and \f$D_i\f$ is the i-th element of the scaling vector.

\author Michael Heroux, ETHZ/D-INFK

\date Last updated on 17-Apr-06

 */
class Tifpack_DiagPreconditioner : public Tpetra_Operator
{
  public:

    //! ctor
    Tifpack_DiagPreconditioner(const Tpetra_Map& DomainMap,
                              const Tpetra_Map& RangeMap,
                              const Tpetra_Vector& diag);

    //! dtor
    ~Tifpack_DiagPreconditioner();

    int SetUseTranspose(bool UseTranspose_in)
    {
      UseTranspose_ = UseTranspose_in;
      return(0);
    }

    int Apply(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const;

    int ApplyInverse(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const;

    double NormInf() const
    {
      return(-1.0);
    }

    const char* Label() const
    {
      return("Tifpack_DiagPreconditioner");
    }

    bool UseTranspose() const
    {
      return(UseTranspose_);
    }

    bool HasNormInf() const
    {
      return(false);
    }

    const Tpetra_Comm& Comm() const
    {
      return(diag_.Comm());
    }

    const Tpetra_Map& OperatorDomainMap() const
    {
      return(RangeMap_);
    }

    const Tpetra_Map& OperatorRangeMap() const
    {
      return(DomainMap_);
    }

    const Tpetra_BlockMap& Map() const
    {
      return(diag_.Map());
    }

  private:
    bool UseTranspose_;
    const Tpetra_Map& DomainMap_;
    const Tpetra_Map& RangeMap_;
    const Tpetra_Vector& diag_;
};

#endif
