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

%module(package="PyTrilinos") ML

#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#ifndef HAVE_ML_EPETRA
#define HAVE_ML_EPETRA
#endif

#ifndef HAVE_ML_TEUCHOS
#define HAVE_ML_TEUCHOS
#endif

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ParameterList.hpp"
%}

// Auto-documentation feature
%feature("autodoc", "1");

// Epetra interface includes
using namespace std;
%import "Epetra.i"

// Amesos interface includes
%include "ml_MultiLevelPreconditioner.h"

// Extensions for ML
%extend ML_Epetra::MultiLevelPreconditioner
{
  void SetInt(char* Name, int value)
  {
    Teuchos::ParameterList List;
    List.set(Name, (int)value);
    self->SetParameterList(List);
  }

  void SetBool(char* Name, bool value)
  {
    Teuchos::ParameterList List;
    List.set(Name, (bool)value);
    self->SetParameterList(List);
  }

  void SetDouble(char* Name, double value)
  {
    Teuchos::ParameterList List;
    List.set(Name, (double)value);
    self->SetParameterList(List);
  }

  void SetString(char* Name, char* value)
  {
    Teuchos::ParameterList List;
    List.set(Name, value);
    self->SetParameterList(List);
  }
}
