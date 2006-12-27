// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include <iostream>
#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

class PyOperator : public Epetra_Operator
{
public:
  PyOperator(const Epetra_Comm& Comm) :
    Comm_(Comm)
  {}

  ~PyOperator() {}

  int SetUseTranspose(bool UseTranspose) = 0;

  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

  virtual double NormInf() const = 0;

  virtual const char * Label() const
  {
    return("PySerialOperator");
  }

  virtual bool UseTranspose() const = 0;

  virtual bool HasNormInf() const = 0;

  virtual const Epetra_Comm & Comm() const
  {
    return(Comm_);
  }

  virtual const Epetra_Map & OperatorDomainMap() const = 0;

  virtual const Epetra_Map & OperatorRangeMap() const = 0;

  virtual const Epetra_Map & Map() const = 0;

private:

  const Epetra_Comm& Comm_;
};
