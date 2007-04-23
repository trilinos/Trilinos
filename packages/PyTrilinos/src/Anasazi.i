// -*- c++ -*-

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

// This documentation string will be the python help facility help
// string
%define %anasazi_docstring
"
PyTrilinos.Anasazi is the python interface to Trilinos package
Anasazi:

    http://software.sandia.gov/trilinos/packages/anasazi

Anasazi is a collection of eigensolver technologies.  Currently, this
python module is a stub, and none of the Anasazi classes are yet
supported.
"
%enddef

// Define the module name, its package and documentation string
%module(package   = "PyTrilinos",
	autodoc   = "1",
	docstring = %anasazi_docstring) Anasazi

%{
// System includes
#include <sstream>

// Configuration includes
#include "PyTrilinos_config.h"

// Epetra includes
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialSymDenseMatrix.h"
#include "Epetra_SerialDenseSVD.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_SerialDistributor.h"
#include "Epetra_MapColoring.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_JadMatrix.h"

// Epetra NumPy includes
#include "Epetra_NumPyIntSerialDenseMatrix.h"
#include "Epetra_NumPyIntSerialDenseVector.h"
#include "Epetra_NumPySerialDenseMatrix.h"
#include "Epetra_NumPySerialDenseVector.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyFEVector.h"

// Anasazi includes
#include "AnasaziVersion.cpp"
#include "AnasaziTypes.hpp"
#include "AnasaziOutputManager.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziEpetraAdapter.hpp"
%}

// General ignore directives
%ignore *::operator=;
%ignore *::print;

// Auto-documentation feature
%feature("autodoc", "1");

// C++ STL support
using namespace std;
%include "stl.i"

// Support for other Trilinos packages
%import "Teuchos.i"
%import "Epetra.i"

/////////////////////////////
// Anasazi Version support //
/////////////////////////////
%include "AnasaziVersion.cpp"
%pythoncode %{
__version__ = Anasazi_Version().split()[2]
%}

///////////////////////////
// Anasazi Types support //
///////////////////////////
%include "AnasaziTypes.hpp"

///////////////////////////////////
// Anasazi OutputManager support //
///////////////////////////////////
%include "AnasaziOutputManager.hpp"
%template (OutputManagerDouble)
          Anasazi::OutputManager<double>;

////////////////////////////////////////
// Anasazi BasicOutputManager support //
////////////////////////////////////////
%include "AnasaziBasicOutputManager.hpp"
%template (BasicOutputManagerDouble)
          Anasazi::BasicOutputManager<double>;

/////////////////////////////////
// Anasazi SortManager support //
/////////////////////////////////
%include "AnasaziSortManager.hpp"

///////////////////////////////
// Anasazi BasicSort support //
///////////////////////////////
%include "AnasaziBasicSort.hpp"

////////////////////////////////////
// Anasazi MultiVecTraits support //
////////////////////////////////////
%include "AnasaziMultiVecTraits.hpp"

//////////////////////////////
// Anasazi MultiVec support //
//////////////////////////////
%include "AnasaziMultiVec.hpp"
%template (MultiVecDouble)
          Anasazi::MultiVec<double>;

////////////////////////////////////
// Anasazi OperatorTraits support //
////////////////////////////////////
%include "AnasaziOperatorTraits.hpp"

//////////////////////////////
// Anasazi Operator support //
//////////////////////////////
%include "AnasaziOperator.hpp"
%template (OperatorDouble)
          Anasazi::Operator<double>;

///////////////////////////////////
// Anasazi EpetraAdapter support //
///////////////////////////////////
%include "AnasaziEpetraAdapter.hpp"
%template (SortManagerEpetra)
          Anasazi::SortManager<double, Epetra_MultiVector, Epetra_Operator>;
%template (BasicSortEpetra)
          Anasazi::BasicSort<double, Epetra_MultiVector, Epetra_Operator>;
%template (MultiVecTraitsEpetra)
          Anasazi::MultiVecTraits<double, Epetra_MultiVector>;
%template (OperatorTraitsEpetra)
          Anasazi::OperatorTraits<double, Epetra_MultiVector, Epetra_Operator>;

//////////////////////////////
// Generic python interface //
//////////////////////////////
%pythoncode %{

%}
