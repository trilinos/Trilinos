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

%module(package="PyTrilinos") Amesos

#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

// Epetra includes
#include "Epetra_Object.h"
#include "Epetra_DistObject.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_DataAccess.h"

// Amesos includes
#include "Amesos_ConfigDefs.h"
#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#ifdef HAVE_AMESOS_LAPACK
#include "Amesos_Lapack.h"
#endif
#ifdef HAVE_AMESOS_KLU
#include "Amesos_Klu.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_SUPERLU
#include "Amesos_Superlu.h"
#endif
#ifdef HAVE_AMESOS_MUMPS
#include "Amesos_Mumps.h"
#endif
%}

// Auto-documentation feature
%feature("autodoc", "1");

// Rename directives for Epetra
%rename(Object              ) Epetra_Object;
%rename(Comm                ) Epetra_Comm;
%rename(SerialComm          ) Epetra_SerialComm;
%rename(BlockMap            ) Epetra_BlockMap;
%rename(Map                 ) Epetra_Map;
%rename(LinearProblem       ) Epetra_LinearProblem;
%rename(MultiVector         ) Epetra_MultiVector;
%rename(Operator            ) Epetra_Operator ;
%rename(RowMatrix           ) Epetra_RowMatrix;
%rename(CrsMatrix           ) Epetra_CrsMatrix;

// Rename directives for Amesos
%rename(BaseSolver          ) Amesos_BaseSolver;
%rename(Factory             ) Amesos;
%rename(Klu                 ) Amesos_Klu;
%rename(Lapack              ) Amesos_Lapack;
%rename(Umfpack             ) Amesos_Umfpack;
%rename(Superlu             ) Amesos_Superlu;
%rename(Superludist         ) Amesos_Superludist;
%rename(Mumps               ) Amesos_Mumps;
%rename(Dscpack             ) Amesos_Dscpack;

// SWIG library includes
%include "std_string.i"

// Epetra interface includes
using namespace std;
%include "Epetra_Object.h"
%include "Epetra_DistObject.h"
%include "Epetra_Comm.h"
%include "Epetra_SerialComm.h"
%include "Epetra_BlockMap.h"
%include "Epetra_Map.h"
%include "Epetra_MultiVector.h"
%include "Epetra_Operator.h"
%include "Epetra_RowMatrix.h"
%include "Epetra_CrsMatrix.h"
%include "Epetra_LinearProblem.h"
%include "Epetra_DataAccess.h"

// Amesos interface includes
%include "Amesos_config.h"
%include "Amesos.h"
%include "Amesos_BaseSolver.h"
#ifdef HAVE_AMESOS_LAPACK
%include "Amesos_Lapack.h"
#endif
#ifdef HAVE_AMESOS_KLU
%include "Amesos_Klu.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
%include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_SUPERLU
%include "Amesos_Superlu.h"
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
%include "Amesos_Superludist.h"
#endif
#ifdef HAVE_AMESOS_MUMPS
%include "Amesos_Mumps.h"
#endif

// Extensions for Epetra
%extend Epetra_Object {
  string __str__() {
    stringstream os;
    self->Print(os);                  // Put the output in os
    string s = os.str();              // Extract the string from os
    return s.substr(0,s.length()-1);  // Return the string minus trailing \n
  }
}

%extend Epetra_CrsMatrix 
{
  int InsertGlobalValue(int i, int j, double val) {
    double val2 = val;
    int j2 = j;
    return self->InsertGlobalValues(i, 1, &val2, &j2);
  }

  int ReplaceGlobalValue(int i, int j, double val) {
    double val2 = val;
    int j2 = j;
    return self->ReplaceGlobalValues(i, 1, &val2, &j2);
  }

  int ReplaceMyValue(int i, int j, double val) {
    double val2 = val;
    int j2 = j;
    return self->ReplaceMyValues(i, 1, &val2, &j2);
  }
}

%extend Epetra_MultiVector
{
  void Set(const int vector, const int element, const double value)
  {
    (*self)[vector][element] = value;
  }

  double Get(const int vector, const int element)
  {
    return((*self)[vector][element]);
  }
}

// Extensions for Amesos
%extend Amesos_BaseSolver 
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
    os << "*** Amesos_BaseSolver ***";
    return os.str();
  }

  void __del__()
  {
    delete self;
  }
}
