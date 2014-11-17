// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
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
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif
#include "Epetra_ConfigDefs.h"
#include "Epetra_DLLExportMacro.h"

// Import the numpy interface
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
%}

// Global swig features
%feature("autodoc", "1");
%feature("compactdefaultargs");

// SWIG library includes
using std::string;
%include "stl.i"

// SWIG NumPy interface file
%include "numpy.i"

// PyTrilinos configuration support
%include "PyTrilinos_config.h"

// Epetra configuration support

// The Epetra configuration header file %include'd below contains many
// 'using' statements, which result in 'Warning(315): Nothing known
// about ...' messages in newer versions of SWIG.  To get SWIG to
// learn about these std functions, etc., would require providing an
// include filepath, which would be difficult to do portably.
// Instead, I provide forward declarations of the entities causing
// warnings, after instructing SWIG to ignore (not wrap) them.

%ignore std::FILE;
%ignore std::istream;
%ignore std::ostream;
%ignore std::abort;
%ignore std::abs;
%ignore std::asin;
%ignore std::atof;
%ignore std::atoi;
%ignore std::ceil;
%ignore std::cerr;
%ignore std::cos;
%ignore std::cout;
%ignore std::endl;
%ignore std::cout;
%ignore std::endl;
%ignore std::exit;
%ignore std::fabs;
%ignore std::fclose;
%ignore std::gets;
%ignore std::fgets;
%ignore std::floor;
%ignore std::flush;
%ignore std::fopen;
%ignore std::fprintf;
%ignore std::free;
%ignore std::malloc;
%ignore std::memcpy;
%ignore std::pow;
%ignore std::rand;
%ignore std::realloc;
%ignore std::sin;
%ignore std::sprintf;
%ignore std::sqrt;
%ignore std::sscanf;
%ignore std::strchr;
%ignore std::strcmp;
%ignore std::strcpy;
%ignore std::strlen;
%ignore std::strtok;

namespace std
{
struct FILE;
class istream;
class ostream;
void abort;
int abs(int);
long abs(long);
float abs(float);
double abs(double);
long double abs(long double);
float asin(float);
double asin(double);
long double asin(long double);
double atof;
int atoi;
float ceil(float);
double ceil(double);
long double ceil(long double);
extern ostream cerr;
float cos(float);
double cos(double);
long double cos(long double);
extern ostream cout;
ostream & endl(ostream&);
void exit(int);
float fabs(float);
double fabs(double);
long double fabs(long double);
int fclose;
char* gets;
char* fgets;
float floor(float);
double floor(double);
long double floor(long double);
int flush;
FILE* fopen;
int fprintf;
void free(void*);
void* malloc;
void* memcpy;
float pow(float);
double pow(double);
long double pow(long double);
int rand;
void* realloc;
float sin(float);
double sin(double);
long double sin(long double);
int sprintf;
float sqrt(float);
double sqrt(double);
long double sqrt(long double);
int sscanf;
char* strchr;
int strcmp;
char* strcpy;
size_t strlen;
char* strtok;
}

%rename(FormatStdout) Epetra_FormatStdout;
%include "Epetra_config.h"
%include "Epetra_ConfigDefs.h"
%include "Epetra_DLLExportMacro.h"

// Teuchos::RCP<> support
%include "Teuchos_RCP.i"

// Include Epetra documentation
%include "Epetra_dox.i"

// Epetra interface includes
%include "Epetra_Base.i"              // Base classes and utility classes
%include "Epetra_SerialDense.i"       // SerialDense classes
%include "Epetra_Comm.i"              // Communicator classes
%include "Epetra_Maps.i"              // Map classes
%include "Epetra_Dist.i"              // Distributed object support
%include "Epetra_Vectors.i"           // Vectors and MultiVectors
%include "Epetra_Graphs.i"            // Graph classes
%include "Epetra_Operators.i"         // Operator and matrix classes

// Turn off the exception handling
%exception;
