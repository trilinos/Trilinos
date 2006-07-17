// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//            PyTrilinos.Galeri: Python Interface to Galeri
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

// This documentation string will be the python help facility help
// string
%define GALERI_DOCSTRING
"Galeri: Matrix Generation Package.
The Galeri module allows an easy creation of Epetra_Map's and
Epetra_CrsMatrix's.  Use the python help() facility for local documentation on
classes and methods, or see the on-line documentation for more in-depth
information. Also give a look to the examples in galeri/python/example.

The most important classes of the IFPACK module are:
- Galeri.CreateMap()
- Galeri.CreateCrsMatrix()
- Galeri.ReadHB()

Galeri requires the Epetra and Teuchos modules of PyTrilinos.
"

%enddef

// Define the module name, its package and documentation string
%module(package="PyTrilinos", docstring=GALERI_DOCSTRING) Galeri

%{
// System includes
#include <sstream>

// Epetra includes
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_LinearProblem.h"

// Epetra Python utility code
#include "Epetra_PyRowMatrix.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"

// Teuchos includes
#include "Teuchos_ParameterList.hpp"

// Teuchos Python utility code
#include "PyTeuchos_Utils.h"

// Galeri includes
#include "Galeri_ConfigDefs.h"
#include "Galeri_Version.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_VbrMatrices.h"
#include "Galeri_Utils.h"
#include "Galeri_ReadHB.h"
// Only needed to let SWIG define SWIGTYPE_p_Epetra_Vector, otherwise
// undefined. 
void Galeri_Dummy(Epetra_Vector* xxx)
{
  cout << "Don't use this..." << endl;
}

%}

%import "Epetra.i"

%feature("autodoc", "1");

%include "std_string.i"
using namespace std;

// typemaps
%typemap(in) (Teuchos::ParameterList& List)
{
  $1 = CreateList($input);
  if ($1 == 0)
  {
    PyErr_SetString(PyExc_ValueError, "Expecting a dictionary");
    return NULL;
  }
}

%typemap(freearg) (Teuchos::ParameterList& List)
{
  if ($1) delete($1);
}

%include "Galeri_Version.h"
%include "Galeri_Maps.h"
%include "Galeri_CrsMatrices.h"
%include "Galeri_VbrMatrices.h"
%include "Galeri_Utils.h"
void Galeri_Dummy(Epetra_Vector* xxx);

namespace Galeri {
%typemap(argout) (Epetra_Map       *& OutMap,
                  Epetra_CrsMatrix *& OutMatrix,
                  Epetra_Vector    *& OutX,
                  Epetra_Vector    *& OutB,
                  Epetra_Vector    *& OutXexact) 
{
  PyObject *oMap, *oMatrix, *oX, *oB, *oXexact;
    oMap    = SWIG_NewPointerObj((void*)(*$1), SWIGTYPE_p_Epetra_Map, 1);
    oMatrix = SWIG_NewPointerObj((void*)(*$2),
                                 SWIGTYPE_p_Epetra_CrsMatrix, 1);
    oX      = SWIG_NewPointerObj((void*)(*$3), SWIGTYPE_p_Epetra_Vector, 1);
    oB      = SWIG_NewPointerObj((void*)(*$4),
                                 SWIGTYPE_p_Epetra_Vector,    1);
    oXexact = SWIG_NewPointerObj((void*)(*$5),
                                 SWIGTYPE_p_Epetra_Vector,    1);
    $result = Py_BuildValue("(OOOOO)",oMap,oMatrix,oX,oB,oXexact);
}

%typemap(in,numinputs=0) Epetra_Map *&OutMap(Epetra_Map* _map) {
  $1 = &_map;
}

%typemap(in,numinputs=0) Epetra_CrsMatrix *&OutMatrix(Epetra_CrsMatrix* _matrix) {
  $1 = &_matrix;
}

%typemap(in,numinputs=0) Epetra_Vector *&OutX(Epetra_Vector* _x) {
  $1 = &_x;
}

%typemap(in,numinputs=0) Epetra_Vector *&OutB(Epetra_Vector* _B) {
  $1 = &_B;
}

%typemap(in,numinputs=0) Epetra_Vector *&OutXexact(Epetra_Vector* _xexact) {
  $1 = &_xexact;
}
void ReadHB(char* data_file, const Epetra_Comm& comm, 
            Epetra_Map*& OutMap,  Epetra_CrsMatrix*& OutMatrix, 
            Epetra_Vector*& OutX, Epetra_Vector*& OutB,
            Epetra_Vector*& OutXexact) throw(int);
}

%pythoncode %{
__version__ = Galeri_Version().split()[2]
%}
