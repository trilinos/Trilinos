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

%define %epetra_docstring
"
PyTrilinos.Epetra is the python interface to the Trilinos linear
algebra services package Epetra:

    http://trilinos.sandia.gov/packages/epetra

The purpose of Epetra is to provide fundamental linear algebra
services to the rest of Trilinos.  These services include parallel
decomposition and communication, vectors and multivectors, graphs,
operators, and dense and sparse matrices.  Note that the C++ version
of Epetra uses the prefix 'Epetra_' which has been stripped from the
python version.

Epetra provides the following user-level classes (by category):

    * Communicators: PyComm, SerialComm, MpiComm (if built with mpi
      support)
    * Data distribution maps: Map, BlockMap, LocalMap
    * Vectors: Vector, MultiVector, IntVector
    * Graphs: CrsGraph, FECrsGraph
    * Operators and matrices: Operator, RowMatrix, CrsMatrix,
      FECrsMatrix, VbrMatrix
    * Serial dense objects: SerialDenseVector, SerialDenseMatrix,
      SerialDenseOperator, SerialDenseSolver, IntSerialDenseVector,
      IntSerialDenseMatrix
    * Aggregates: LinearProblem
    * Utilities: Import, Export, Time, MapColoring, Util

For examples of usage, please consult the following scripts in the
example subdirectory of the PyTrilinos package:

    * exEpetra.py
    * exEpetra_Comm.py
    * exEpetra_ImportExport.py
    * exEpetra_CrsMatrix_Easy.py
    * exEpetra_CrsMatrix_Efficient.py
    * exEpetra_FECrsMatrix_Easy.py

The Epetra module has been designed to use and interoperate with the
numpy module, which provides multidimensional array support.  Epetra
class constructors or methods that expect C arrays in C++ can
typically accept numpy arrays in python (or any python sequence that
numpy can convert to an array).  Similarly, methods that return C
arrays in C++ will return numpy arrays in python.  Also, certain
Epetra classes that represent contiguous blocks of homogeneous data
have been given the attributes of numpy arrays using multiple
inheritance: Vector, MultiVector, IntVector, SerialDenseVector,
SerialDenseMatrix, IntSerialDenseVector and IntSerialDenseMatrix.
"
%enddef

%module(package   = "PyTrilinos",
	directors = "1",
	docstring = %epetra_docstring) Epetra

%{
// System includes
#include <sstream>

// Configuration includes
#include "PyTrilinos_config.h"
#include "Epetra_ConfigDefs.h"

// Import the numpy interface
#include "NumPyImporter.h"
%}

// Global swig features
%feature("autodoc", "1");
%feature("compactdefaultargs");

// Configuration support
%include "PyTrilinos_config.h"
%rename(FormatStdout) Epetra_FormatStdout;
%include "Epetra_ConfigDefs.h"

// Include Epetra documentation
%include "Epetra_dox.i"

// SWIG library includes
using std::string;
%include "stl.i"

// SWIG NumPy interface file
%include "numpy.i"

// Epetra interface includes
%include "Epetra_Base.i"              // Base classes and utility classes
%include "Epetra_SerialDense.i"       // SerialDense classes
%include "Epetra_Comm.i"              // Communicator classes
%include "Epetra_Maps.i"              // Map classes
%include "Epetra_Vectors.i"           // Vectors and MultiVectors
%include "Epetra_Graphs.i"            // Graph classes
%include "Epetra_Operators.i"         // Operator and matrix classes

// Turn off the exception handling
%exception;
