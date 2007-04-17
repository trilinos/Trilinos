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
"""
Epetra is a module under the PyTrilinos package.  PyTrilinos is a
python interface to select components of the Trilinos set of solver
packages developed at Sandia National Laboratories.  Documentation for
Trilinos can be found in the source directory or online at
http://software.sandia.gov/Trilinos.

The Epetra module offers a suite of linear algebra services,
supporting both serial and parallel data distribution, and dense and
sparse storage.  The special Epetra.PyComm function returns either a
SerialComm or MpiComm communicator, depending on the environment
Epetra is imported under.

Note that in the python interface, the 'Epetra_' prefix has been
stripped from all Epetra classes, but (unlike the C++ implementation)
Epetra resides in its own namespace.  Use the python help() facility
for local documentation on classes and methods, or see the online
documentation for more in-depth information.

Epetra provides the following user-level classes:

    Communicators: PyComm, SerialComm, MpiComm (if built with mpi
        support)
    Data distribution maps: Map, BlockMap, LocalMap
    Vectors: Vector, MultiVector, IntVector
    Graphs: CrsGraph, FECrsGraph
    Operators and matrices: Operator, RowMatrix, CrsMatrix,
        FECrsMatrix, VbrMatrix
    Serial dense objects: SerialDenseVector, SerialDenseMatrix,
        SerialDenseOperator, SerialDenseSolver, IntSerialDenseVector,
        IntSerialDenseMatrix
    Aggregates: LinearProblem
    Utilities: Import, Export, Time, MapColoring

The Epetra module has been designed to use and interoperate with the
numpy module, which provides multidimensional array support.  Epetra
class constructors or methods that expect C arrays in C++ can
typically accept numpy arrays in python (or any python sequence that
numpy can convert to an array).  Similarly, methods that return C
arrays in C++ will return numpy arrays in python.  Also, certain
Epetra classes represent contiguous blocks of homogeneous data.  These
classes have been given the attributes of numpy arrays using multiple
inheritance, and include the Vector, MultiVector, IntVector,
SerialDenseVector, SerialDenseMatrix, IntSerialDenseVector and
IntSerialDenseMatrix classes.
"""
%enddef

%module(package   = "PyTrilinos",
	directors = "1",
	docstring = %epetra_docstring) Epetra

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

// Configuration includes
#include "PyTrilinos_config.h"

// Import the numpy interface
#include "NumPyImporter.h"
%}

// Global swig features
%feature("autodoc", "1");
%feature("compactdefaultargs");

// SWIG library includes
using namespace std;
%include "stl.i"
%include "exception.i"

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
