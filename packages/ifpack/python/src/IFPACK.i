// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//            PyTrilinos.Epetra: Python Interface to Epetra
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

%module(package="PyTrilinos") IFPACK

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

#include "Epetra_VbrMatrix.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Operator.h"
#include "Teuchos_ParameterList.hpp"

// Amesos includes
#include "Ifpack.h"
#include "Ifpack_Utils.h"
#include "Ifpack_Preconditioner.h"
%}

// Auto-documentation feature
%feature("autodoc", "1");

%rename (Factory) Ifpack;
%rename (ILU) Ifpack_ILU;

using namespace std;

// Epetra interface includes
%import "Epetra.i"
%import "AztecOO.i"

// Amesos interface includes
%include "Ifpack.h"
%include "Ifpack_Utils.h"
%include "Ifpack_Preconditioner.h"

// Extensions for Ifpack
%extend Ifpack_Preconditioner
{
  void SetInt(char* Name, int value)
  {
    Teuchos::ParameterList List;
    List.set(Name, (int)value);
    self->SetParameters(List);
  }

  void SetBool(char* Name, bool value)
  {
    Teuchos::ParameterList List;
    List.set(Name, (bool)value);
    self->SetParameters(List);
  }

  void SetDouble(char* Name, double value)
  {
    Teuchos::ParameterList List;
    List.set(Name, (double)value);
    self->SetParameters(List);
  }

  void SetString(char* Name, char* value)
  {
    Teuchos::ParameterList List;
    List.set(Name, value);
    self->SetParameters(List);
  }

  string __str__() {
    stringstream os;
    os << *self;
    return os.str();
  }

  void __del__()
  {
    delete self;
  }
}

