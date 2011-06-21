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

%define %galeri_docstring
"
PyTrilinos.Galeri is the python interface to the Trilinos example
matrix package Galeri:

    http://trilinos.sandia.gov/packages/galeri

The purpose of Galeri is to provide a set of example matrices
distributed across a set of common processor decompositions for
testing purposes.

Galeri provides the following user-level functions:

    * CreateCartesianCoordinates  - Create cartesian coordinates
    * GetNeighboursCartesian2d    - Get neighbor nodes from a 2D grid
    * GetNeighboursCartesian3d    - Get neighbor nodes from a 3D grid
    * PrintStencil2D              - Output a 2D stencil
    * CreateMap                   - Create a Map object
    * CreateCrsMatrix             - Create a specified CrsMatrix
    * CreateVbrMatrix             - Create a specified VbrMatrix
    * ReadHB                      - Read a problem definition from an HB file

For examples of usage, please consult the following scripts in the
example subdirectory of the PyTrilinos package:

    * exGaleri.py
    * exGaleri_ReadHB.py
    * exAztecOO.py
    * exIFPACK.py
    * exMLAPI_Simple.py
"
%enddef

%module(package   = "PyTrilinos",
	autodoc   = "1",
	docstring = %galeri_docstring) Galeri

%{
// System includes
#include <sstream>

// Configuration includes
#include "PyTrilinos_config.h"
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif

// Epetra includes
#ifdef HAVE_EPETRA
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_MapColoring.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_InvOperator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_SerialDenseSVD.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDistributor.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_OffsetIndex.h"
#include "Epetra_Time.h"

// Epetra Python utility code
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyFEVector.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_NumPySerialDenseMatrix.h"
#include "Epetra_NumPySerialSymDenseMatrix.h"
#include "Epetra_NumPySerialDenseVector.h"
#include "Epetra_NumPyIntSerialDenseMatrix.h"
#include "Epetra_NumPyIntSerialDenseVector.h"
#endif

// Teuchos includes
#ifdef HAVE_TEUCHOS
#include "Teuchos_PythonParameter.h"
#endif

// Galeri includes
#include "Galeri_Version.h"
#include "Galeri_Utils.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_VbrMatrices.h"
#include "Galeri_ReadHB.h"
%}

// Include PyTrilinos configuration
%include "PyTrilinos_config.h"

// Standard exception handling
%include "exception.i"

// Turn on autodocumentation
%feature("autodoc", "1");

// Include Galeri documentation
%include "Galeri_dox.i"

// Typemap support for STL
%include "stl.i"

// Trilinos package imports
#ifdef HAVE_TEUCHOS
%import "Teuchos.i"
#endif
#ifdef HAVE_EPETRA
%import "Epetra.i"
#endif

// General exception handling
%feature("director:except")
{
  if ($error != NULL)
  {
    throw Swig::DirectorMethodException();
  }
}

%exception
{
  try
  {
    $action
  }
  catch(Swig::DirectorException &e)
  {
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

////////////////////////////
// Galeri_Version support //
////////////////////////////
%include "Galeri_Version.h"
%pythoncode
%{
__version__ = Galeri_Version().split()[2]
%}

//////////////////////////
// Galeri_Utils support //
//////////////////////////
%newobject Galeri::CreateCartesianCoordinates;
%include "Galeri_Utils.h"

/////////////////////////
// Galeri_Maps support //
/////////////////////////
%newobject Galeri::CreateMap;
%include "Galeri_Maps.h"

////////////////////////////////
// Galeri_CrsMatrices support //
////////////////////////////////
%newobject Galeri::CreatCrsMatrix;
%include "Galeri_CrsMatrices.h"

////////////////////////////////
// Galeri_VbrMatrices support //
////////////////////////////////
%newobject Galeri::CreateVbrMatrix;
%include "Galeri_VbrMatrices.h"

///////////////////////////
// Galeri_ReadHB support //
///////////////////////////
#ifdef HAVE_EPETRA
%feature("autodoc",
"ReadHB(str filename, Epetra.Comm comm) -> (Epetra.Map map, Epetra.CrsMatrix A,
                                           Epetra.Vector x, Epetra.Vector b,
                                           Epetra.Vector exact)

Given an HB filename and an Epetra communicator, return a tuple of
Epetra objects generated from the problem described in the file.  This
tuple contains an Epetra.Map, Epetra.CrsMatrix, and Epetra.Vectors for
the solution, right-hand side and the exact solution.")
Galeri::ReadHB;
%include "Galeri_ReadHB.h"
#endif

// Turn off the exception handling
%exception;
