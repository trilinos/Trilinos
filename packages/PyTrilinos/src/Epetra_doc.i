// -*- C++ -*-
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

// The purpose of this file is to provide manually-written
// documentation strings for PyTrilinos.Epetra functions, classes and
// methods, when the doxygen- or swig-generated documentation is
// either insufficient or inaccurate.

// Epetra.Object class

%feature("docstring")
Epetra_Object
"The base Epetra class.
    
The Epetra_Object class provides capabilities common to all Epetra
objects, such as a label that identifies an object instance, constant
definitions, enum types.  In C++, it supports a ``Print()`` method
that takes an output stream as an argument.  In the python
implementation for this and all derived classes, this method takes an
optional file object argument whose default value is standard out."

%feature("docstring")
Epetra_Object::__str__
"Returns the results of ``Print()`` in a string, so that
the ``print`` command will work on ``Epetra`` objects.  The
``Print()`` methods are designed to run correctly in parallel, so do
not execute ``print`` on an Epetra object conditionally on the
processor number.  For example, do not do

  ``if comm.MyPID() == 0: print epetra_obj``

or it will hang your code."

// Epetra.Util class

%feature("docstring")
Epetra_Util
"Epetra Util Wrapper Class.

The Epetra.Util class is a collection of useful functions that cut
across a broad set of other classes.  A random number generator is
provided, along with methods to set and retrieve the random-number
seed.

The random number generator is a multiplicative linear congruential
generator, with multiplier 16807 and modulus 2^31 - 1. It is based on
the algorithm described in 'Random Number Generators: Good Ones Are
Hard To Find', S. K. Park and K. W. Miller, Communications of the ACM,
vol. 31, no. 10, pp. 1192-1201.

The C++ Sort() method is not supported in python.

A static function is provided for creating a new Epetra.Map object
with 1-to-1 ownership of entries from an existing map which may have
entries that appear on multiple processors.

Epetra.Util is a serial interface only.  This is appropriate since the
standard utilities are only specified for serial execution (or shared
memory parallel)."

%feature("autodoc",
"Broadcast(self, numpy.ndarray myObj, int root)

Argument myObj must be a numpy array, so that the Broadcast can be
performed in-place.  Its scalar data type must be int, long or double.
In C++, this routine has an integer error return code.  In python, a
non-zero return code is converted to an exception.")
Epetra_Comm::Broadcast;

%feature("docstring")
Epetra_Comm::GatherAll
"Argument myObj can be a numpy array or any sequence that can be
converted to a numpy array.  Its scalar data type must be int, long or
double.  The return argument is a numpy array of the same type.  In
C++, this routine has an integer error return code.  In python, a
non-zero return code is converted to an exception."

%feature("docstring")
Epetra_Comm::SumAll
"Argument myObj can be a numpy array or any sequence that can be
converted to a numpy array.  Its scalar data type must be int, long or
double.  The return argument is a numpy array of the same type.  In
C++, this routine has an integer error return code.  In python, a
non-zero return code is converted to an exception."

%feature("docstring")
Epetra_Comm::MaxAll
"Argument myObj can be a numpy array or any sequence that can be
converted to a numpy array.  Its scalar data type must be int, long or
double.  The return argument is a numpy array of the same type.  In
C++, this routine has an integer error return code.  In python, a
non-zero return code is converted to an exception."

%feature("docstring")
Epetra_Comm::MinAll
"Argument myObj can be a numpy array or any sequence that can be
converted to a numpy array.  Its scalar data type must be int, long or
double.  The return argument is a numpy array of the same type.  In
C++, this routine has an integer error return code.  In python, a
non-zero return code is converted to an exception."

%feature("docstring")
Epetra_Comm::ScanSum
"Argument myObj can be a numpy array or any sequence that can be
converted to a numpy array.  Its scalar data type must be int, long or
double.  The return argument is a numpy array of the same type.  In
C++, this routine has an integer error return code.  In python, a
non-zero return code is converted to an exception."

// Epetra.BlockMap class

%feature("autodoc",
"
__init__(self, int numGlobalElements, int elementSize, int indexBase,
     Comm comm) -> BlockMap

BlockMap constructor with implicit local elements and constant element
size.  Arguments are:

     numGlobalElements  - Total number of elements over all processors.
                          Specify -1 to have the constructor compute
                          the number of global elements
     elementSize        - The number of degrees of freedom associated
                          with every element.
     indexBase          - The base integer value for indexed array
                          references.  Typically this is 0 for C/C++ and 1
                          for Fortran, but it can be set to any integer
                          value.
     comm               - The Epetra.Comm communicator. This communicator
                          can in turn be queried for processor rank and
                          size information.")
Epetra_BlockMap::Epetra_BlockMap(int, int, int, const Epetra_Comm &);

%feature("autodoc",
"
__init__(self, int numGlobalElements, int numMyElements, int elementSize,
     int indexBase, Comm comm) -> BlockMap

BlockMap constructor with specified number of local elements and
constant element size.  Arguments are:

     numGlobalElements  - Total number of elements over all processors.
                          Specify -1 to have the constructor compute
                          the number of global elements
     numMyElements      - Number of local elements on this processor.
     elementSize        - The number of degrees of freedom associated
                          with every element.
     indexBase          - The base integer value for indexed array
                          references.  Typically this is 0 for C/C++ and 1
                          for Fortran, but it can be set to any integer
                          value.
     comm               - The Epetra.Comm communicator. This communicator
                          can in turn be queried for processor rank and
                          size information.")
Epetra_BlockMap::Epetra_BlockMap(int, int, int, int, const Epetra_Comm &);

%feature("autodoc",
"
__init__(self, int numGlobalElements, PySequence myGlobalElements,
     int elementSize, int indexBase, Comm comm) -> BlockMap

BlockMap constructor with specified list of local elements and
constant element size.  Arguments are:

     numGlobalElements  - Total number of elements over all processors.
                          Specify -1 to have the constructor compute
                          the number of global elements
     myGlobalElements   - A sequence of integers specifying the global
                          element indexes on this processor.
     elementSize        - The number of degrees of freedom associated
                          with every element.
     indexBase          - The base integer value for indexed array
                          references.  Typically this is 0 for C/C++ and 1
                          for Fortran, but it can be set to any integer
                          value.
     comm               - The Epetra.Comm communicator. This communicator
                          can in turn be queried for processor rank and
                          size information.")
Epetra_BlockMap::Epetra_BlockMap(int, int, const int*, int, int,
				 const Epetra_Comm &);

%feature("autodoc",
"
__init__(self, int numGlobalElements, PySequence myGlobalElements,
     PySequence elementsSizes, int indexBase, Comm comm) -> BlockMap

BlockMap constructor with specified list of local elements and
specified list of element sizes.  Arguments are:

     numGlobalElements  - Total number of elements over all processors.
                          Specify -1 to have the constructor compute
                          the number of global elements
     myGlobalElements   - A sequence of integers specifying the global
                          element indexes on this processor.
     elementSizes       - A sequence of integers specifying the number of
                          degrees of freedom associated with each element
                          on this processor.
     indexBase          - The base integer value for indexed array
                          references.  Typically this is 0 for C/C++ and 1
                          for Fortran, but it can be set to any integer
                          value.
     comm               - The Epetra.Comm communicator. This communicator
                          can in turn be queried for processor rank and
                          size information.")
Epetra_BlockMap::Epetra_BlockMap(int, int, const int*, int, const int*, int,
				 const Epetra_Comm &);

%feature("autodoc",
"
__init__(self, BlockMap map) -> BlockMap

BlockMap copy constructor.")
Epetra_BlockMap::Epetra_BlockMap(const Epetra_BlockMap &);

%feature("docstring")
Epetra_BlockMap::RemoteIDList
"``GIDList`` is a sequence of integer global IDs, and the return
argument is the three-tuple ``(PIDList, LIDList, sizeList)``, which
are ``numpy.ndarray`` objects of integers representing the processor
IDs, local IDs and element sizes, respectively."

%feature("docstring")
Epetra_BlockMap::FindLocalElementID
"Returns a tuple containing the local ID of the element that contains
the given local pointID, and the offset of the point in that element."

%feature("docstring")
Epetra_BlockMap::MyGlobalElements
"Returns a numpy array of integers specifying the list of global IDs on
the processor."

%feature("docstring")
Epetra_BlockMap::FirstPointInElementList
"Returns a numpy array of integer first local point numbers for all of
the local elements."

%feature("docstring")
Epetra_BlockMap::ElementSizeList
"Returns a numpy array of integer sizes for each local element."

%feature("docstring")
Epetra_BlockMap::PointToElementList
"Returns a numpy array of integers such that for each local point, it
indicates the local element ID that the point belongs to."

// Epetra.Map class

%feature("autodoc",
"
__init__(self, int numGlobalElements, int indexBase, Comm comm) -> Map

Map constructor with implicit number of elements per processor.
Arguments are:

     numGlobalElements  - Total number of elements over all processors.
                          Specify -1 to have the constructor compute
                          the number of global elements
     indexBase          - The base integer value for indexed array
                          references.  Typically this is 0 for C/C++ and 1
                          for Fortran, but it can be set to any integer
                          value.
     comm               - The Epetra.Comm communicator. This communicator
                          can in turn be queried for processor rank and
                          size information.")
Epetra_Map::Epetra_Map(int, int, const Epetra_Comm &);

%feature("autodoc",
"
__init__(self, int numGlobalElements, int numMyElements, int indexBase,
     Comm comm) -> Map

Map constructor with specified number of elements per processor.
Arguments are:

     numGlobalElements  - Total number of elements over all processors.
                          Specify -1 to have the constructor compute
                          the number of global elements
     numMyElements      - Number of local elements on this processor.
     indexBase          - The base integer value for indexed array
                          references.  Typically this is 0 for C/C++ and 1
                          for Fortran, but it can be set to any integer
                          value.
     comm               - The Epetra.Comm communicator. This communicator
                          can in turn be queried for processor rank and
                          size information.")
Epetra_Map::Epetra_Map(int, int, int, const Epetra_Comm &);

%feature("autodoc",
"
__init__(self, int numGlobalElements, PySequence myGlobalElements,
     int indexBase, Comm comm) -> Map

Map constructor with specified list of global element IDs for each
processor.  Arguments are:

     numGlobalElements  - Total number of elements over all processors.
                          Specify -1 to have the constructor compute
                          the number of global elements
     myGlobalElements   - A sequence of integers specifying the global
                          element indexes on this processor.
     indexBase          - The base integer value for indexed array
                          references.  Typically this is 0 for C/C++ and 1
                          for Fortran, but it can be set to any integer
                          value.
     comm               - The Epetra.Comm communicator. This communicator
                          can in turn be queried for processor rank and
                          size information.")
Epetra_Map::Epetra_Map(int, int, const int*, int, const Epetra_Comm &);

%feature("autodoc",
"
__init__(self, Map map) -> Map

Map copy constructor.")
Epetra_Map::Epetra_Map(const Epetra_Map &);

// Epetra.NumPyMultiVector class

%feature("docstring")
Epetra_NumPyMultiVector::ExtractCopy
"Return a numpy.ndarray that is a copy of the MultiVector."

%feature("docstring")
Epetra_NumPyMultiVector::ExtractView
"Return a numpy.ndarray that is a view of the MultiVector."

%feature("docstring")
Epetra_NumPyMultiVector::Dot
"Return a numpy.ndarray of the dot products of the MultiVector and a."

%feature("docstring")
Epetra_NumPyMultiVector::Norm1
"Return a numpy.ndarray of the L-1 norms of MultiVector."

%feature("docstring")
Epetra_NumPyMultiVector::Norm2
"Return a numpy.ndarray of the the L-2 norms of MultiVector."

%feature("docstring")
Epetra_NumPyMultiVector::NormInf
"Return a numpy.ndarray of the L-infinity norms of MultiVector."

%feature("docstring")
Epetra_NumPyMultiVector::NormWeighted
"Return a numpy.ndarray of the weighted norms of MultiVector."

%feature("docstring")
Epetra_NumPyMultiVector::MinValue
"Return a numpy.ndarray of the minimum values in MultiVector."

%feature("docstring")
Epetra_NumPyMultiVector::MaxValue
"Return a numpy.ndarray of the maximum values in MultiVector."

%feature("docstring")
Epetra_NumPyMultiVector::MeanValue
"Return a numpy.ndarray of the mean values of the MultiVector."

// Epetra.NumPyVector class

%feature("docstring")
Epetra_NumPyVector::ExtractCopy
"Return a numpy.ndarray that is a copy of the Vector."

%feature("docstring")
Epetra_NumPyVector::ExtractView
"Return a numpy.ndarray that is a view of the Vector."

%feature("docstring")
Epetra_NumPyVector::Dot
"Return the dot product of the Vector and a."

%feature("docstring")
Epetra_NumPyVector::Norm1
"Return the L-1 norm of Vector."

%feature("docstring")
Epetra_NumPyVector::Norm2
"Return the the L-2 norm of Vector."

%feature("docstring")
Epetra_NumPyVector::NormInf
"Return the L-infinity norm of Vector."

%feature("docstring")
Epetra_NumPyVector::NormWeighted
"Return the weighted norm of Vector."

%feature("docstring")
Epetra_NumPyVector::MinValue
"Return the minimum values in Vector."

%feature("docstring")
Epetra_NumPyVector::MaxValue
"Return the maximum values in Vector."

%feature("docstring")
Epetra_NumPyVector::MeanValue
"Return the mean value of the Vector."

%feature("docstring")
Epetra_NumPyVector::ReplaceGlobalValues
"Replace global values at specified index (and offset)"

%feature("docstring")
Epetra_NumPyVector::ReplaceMyValues
"Replace local values at specified index (and offset)"

%feature("docstring")
Epetra_NumPyVector::SumIntoGlobalValues
"Sum into global values at specified indices (and offset)"

%feature("docstring")
Epetra_NumPyVector::SumIntoMyValues
"Sum into local values at specified indices (and offset)"

// Epetra.NumPyFEVector class

%feature("docstring")
Epetra_NumPyFEVector::ExtractCopy
"Return a numpy.ndarray that is a copy of the FEVector."

%feature("docstring")
Epetra_NumPyFEVector::ExtractView
"Return a numpy.ndarray that is a view of the FEVector."

%feature("docstring")
Epetra_NumPyFEVector::Dot
"Return the dot product of the FEVector and a."

%feature("docstring")
Epetra_NumPyFEVector::Norm1
"Return the L-1 norm of FEVector."

%feature("docstring")
Epetra_NumPyFEVector::Norm2
"Return the the L-2 norm of FEVector."

%feature("docstring")
Epetra_NumPyFEVector::NormInf
"Return the L-infinity norm of FEVector."

%feature("docstring")
Epetra_NumPyFEVector::NormWeighted
"Return the weighted norm of FEVector."

%feature("docstring")
Epetra_NumPyFEVector::MinValue
"Return the minimum values in FEVector."

%feature("docstring")
Epetra_NumPyFEVector::MaxValue
"Return the maximum values in FEVector."

%feature("docstring")
Epetra_NumPyFEVector::MeanValue
"Return the mean value of the FEVector."

%feature("docstring")
Epetra_NumPyFEVector::ReplaceGlobalValues
"Replace global values at specified index (and offset)"

%feature("docstring")
Epetra_NumPyFEVector::SumIntoGlobalValues
"Sum into global values at specified indices (and offset)"

// Epetra.IntVector class

%feature("docstring")
Epetra_NumPyIntVector::ExtractCopy
"Return a numpy.ndarray that is a copy of the IntVector."

%feature("docstring")
Epetra_NumPyIntVector::ExtractView
"Return a numpy.ndarray that is a view of the IntVector."

%feature("docstring")
Epetra_NumPyIntVector::Values
"Return a numpy.ndarray that is a view of the IntVector."

// Epetra.CrsGraph class

%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, int numIndicesPerRow,
    bool staticProfile=False) -> CrsGraph

  Constructor with implicit column map and constant indices per row.
  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - Map describing distribution of rows across processors
    numIndicesPerRow  - Integer number of indices per row
    staticProfile     - Static profile flag")
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess, const Epetra_BlockMap&,
				 int, bool);

%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, PySequence
    numIndicesPerRow, bool staticProfile=False) -> CrsGraph

  Constructor with implicit column map and variable indices per row.
  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - Map describing distribution of rows across processors
    numIndicesPerRow  - Sequence of integers representing the number of indices
                        per row
    staticProfile     - Static profile flag")
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess, const Epetra_BlockMap&,
				 int, const int*, bool);

%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, BlockMap colMap,
    int numIndicesPerRow, bool staticProfile=False) -> CrsGraph

  Constructor with specified column map and constant indices per row.
  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - Map describing distribution of rows across processors
    colMap            - Map describing distribution of columns across processors
    numIndicesPerRow  - Integer number of indices per row
    staticProfile     - Static profile flag")
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess, const Epetra_BlockMap&,
				 const Epetra_BlockMap&, int, bool);

%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, BlockMap colMap,
    PySequence numIndicesPerRow, bool staticProfile=False) -> CrsGraph

  Constructor with specified column map and variable indices per row.
  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - Map describing distribution of rows across processors
    colMap            - Map describing distribution of columns across processors
    numIndicesPerRow  - Sequence of integers representing the number of indices
                        per row
    staticProfile     - Static profile flag")
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess, const Epetra_BlockMap&,
				 const Epetra_BlockMap&, int, const int*, bool);

%feature("autodoc",
"
__init__(self, CrsGraph graph) -> CrsGraph

  Copy constructor.  Arguments:

    graph - Source graph for copy constructor")
Epetra_CrsGraph::Epetra_CrsGraph(const Epetra_CrsGraph&);

%feature("autodoc",
"InsertGlobalIndices(self, int globalRow, PySequence indices) -> int

Insert a sequence of global indices into the set of nonzero columns
for the specified global row.  Argument indices can be a numpy array
of integers or any python sequence that can be converted to a numpy
array of integers.  The integers represent global IDs that are to be
inserted into the graph.  An integer error/warning code is returned.
")
Epetra_CrsGraph::InsertGlobalIndices;

%feature("autodoc",
"RemoveGlobalIndices(self, int globalRow, PySequence indices) -> int

Remove a sequence of global indices from the set of nonzero columns
for the specified global row.  Argument indices can be a numpy array
of integers or any python sequence that can be converted to a numpy
array of integers.  The integers represent global IDs that are to be
removed from the graph.  An integer error/warning code is returned.
")
Epetra_CrsGraph::RemoveGlobalIndices(int, int, int*);

%feature("autodoc",
"InsertMyIndices(self, int localRow, PySequence indices) -> int

Insert a sequence of local indices into the set of nonzero columns for
the specified local row.  Argument indices can be a numpy array of
integers or any python sequence that can be converted to a numpy array
of integers.  The integers represent local IDs that are to be inserted
into the graph.  An integer error/warning code is returned.
")
Epetra_CrsGraph::InsertMyIndices;

%feature("autodoc",
"RemoveMyIndices(self, int localRow, PySequence indices) -> int

Remove a sequence of local indices from the set of nonzero columns for
the specified local row.  Argument indices can be a numpy array of
integers or any python sequence that can be converted to a numpy array
of integers.  The integers represent local IDs that are to be removed
from the graph.  An integer error/warning code is returned.
")
Epetra_CrsGraph::RemoveMyIndices(int, int, int*);

// Epetra.CrsMatrix class

%feature("autodoc",
"__init__(self, Epetra_DataAccess CV, Map rowMap, int numEntriesPerRow, 
    bool staticProfile=False) -> CrsMatrix

  CrsMatrix constructor with implicit column map and constant number
  of entries per row.  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - describes distribution of rows across processors
    numEntriesPerRow  - constant number of entries per row
    staticProfile     - static profile flag
"
)
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess, const Epetra_Map&, int,
				   bool);

%feature("autodoc",
"__init__(self, Epetra_DataAccess CV, Map rowMap, PySequence numEntriesPerRow, 
    bool staticProfile=False) -> CrsMatrix

  CrsMatrix constructor with implicit column map and variable number
  of entries per row.  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - describes distribution of rows across processors
    numEntriesPerRow  - variable number of entries per row
    staticProfile     - static profile flag
"
)
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess, const Epetra_Map&,
				   const int*, int, bool);


%feature("autodoc",
"__init__(self, Epetra_DataAccess CV, Map rowMap, Map colMap, int numEntriesPerRow, 
    bool staticProfile=False) -> CrsMatrix

  CrsMatrix constructor with specified column map and constant number
  of entries per row.  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - describes distribution of rows across processors
    colMap            - describes distribution of columns across processors
    numEntriesPerRow  - constant number of entries per row
    staticProfile     - static profile flag
"
)
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess, const Epetra_Map&,
				   const Epetra_Map&, int, bool);


%feature("autodoc",
"__init__(self, Epetra_DataAccess CV, Map rowMap, Map colMap, PySequence
    numEntriesPerRow, bool staticProfile=False) -> CrsMatrix

  CrsMatrix constructor with specified column map and variable number
  of entries per row.  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - describes distribution of rows across processors
    colMap            - describes distribution of columns across processors
    numEntriesPerRow  - variable number of entries per row
    staticProfile     - static profile flag
"
)
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess, const Epetra_Map&,
				   const Epetra_Map&, const int*, int, bool);


%feature("autodoc",
"__init__(self, Epetra_DataAccess CV, CrsGraph graph) -> CrsMatrix

  CrsMatrix constructor with CrsGraph.  Arguments:

    CV     - Epetra.Copy or Epetra.View
    graph  - CrsGraph describing structure of matrix
"
)
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess, const Epetra_CrsGraph&);


%feature("autodoc",
"__init__(self, CrsMatrix matrix) -> CrsMatrix

  CrsMatrix copy constructor.  Argument:

    matrix  - source CrsMatrix
"
)
Epetra_CrsMatrix::Epetra_CrsMatrix(const Epetra_CrsMatrix&);

%feature("autodoc",
"ExtractGlobalRowCopy(self, int globalRow) -> (numpy.ndarray,numpy.ndarray)

  Returns a two-tuple of numpy arrays of the same size; the first is
  an array of integers that represent the nonzero columns on the
  matrix; the second is an array of doubles that represent the values
  of the matrix entries.  The input argument is a global row index."
)
Epetra_CrsMatrix::ExtractGlobalRowCopy;

%feature("autodoc",
"ExtractMyRowCopy(self, int myRow) -> (numpy.ndarray,numpy.ndarray)

  Returns a two-tuple of numpy arrays of the same size; the first is
  an array of integers that represent the nonzero columns on the
  matrix; the second is an array of doubles that represent the values
  of the matrix entries.  The input argument is a local row index."
)
Epetra_CrsMatrix::ExtractMyRowCopy;

%feature("autodoc",
"InsertGlobalValues(self, int globalRow, PySequence values, PySequence
    indices) -> int

  Arguments:

    globalRow - global row index
    values    - a sequence of doubles that represent the values to insert
    indices   - a sequence of integers that represent the indices to insert
")
Epetra_CrsMatrix::InsertGlobalValues(int,double*,int,int*,int);

%feature("autodoc",
"InsertGlobalValues(self, PySequence rows, PySequence cols, PySequence
    values) -> int

  Arguments:

    rows    - a sequence of integers that represent the row indices to insert
    cols    - a sequence of integers that represent the column indices to
              insert
    values  - a sequence of doubles that represent the values to insert
")
Epetra_CrsMatrix::InsertGlobalValues(PyObject*,PyObject*,PyObject*);

%feature("autodoc",
"InsertMyValues(self, int myRow, PySequence values, PySequence indices) -> int

  Arguments:

    myRow     - local row index
    values    - a sequence of doubles that represent the values to insert
    indices   - a sequence of integers that represent the indices to insert
")
Epetra_CrsMatrix::InsertMyValues(int,double*,int,int*,int);

%feature("autodoc",
"InsertMyValues(self, PySequence rows, PySequence cols, PySequence
    values) -> int

  Arguments:

    rows    - a sequence of integers that represent the row indices to insert
    cols    - a sequence of integers that represent the column indices to
              insert
    values  - a sequence of doubles that represent the values to insert
")
Epetra_CrsMatrix::InsertMyValues(PyObject*,PyObject*,PyObject*);

%feature("autodoc",
"ReplaceGlobalValues(self, int globalRow, PySequence values, PySequence
    indices) -> int

  Arguments:

    globalRow - global row index
    values    - a sequence of doubles that represent the values to replace
    indices   - a sequence of integers that represent the indices to replace
")
Epetra_CrsMatrix::ReplaceGlobalValues(int,double*,int,int*,int);

%feature("autodoc",
"ReplaceGlobalValues(self, PySequence rows, PySequence cols, PySequence
    values) -> int

  Arguments:

    rows    - a sequence of integers that represent the row indices to replace
    cols    - a sequence of integers that represent the column indices to
              replace
    values  - a sequence of doubles that represent the values to replace
")
Epetra_CrsMatrix::ReplaceGlobalValues(PyObject*,PyObject*,PyObject*);

%feature("autodoc",
"ReplaceMyValues(self, int myRow, PySequence values, PySequence indices) -> int

  Arguments:

    myRow     - local row index
    values    - a sequence of doubles that represent the values to replace
    indices   - a sequence of integers that represent the indices to replace
")
Epetra_CrsMatrix::ReplaceMyValues(int,double*,int,int*,int);

%feature("autodoc",
"ReplaceMyValues(self, PySequence rows, PySequence cols, PySequence
    values) -> int

  Arguments:

    rows    - a sequence of integers that represent the row indices to replace
    cols    - a sequence of integers that represent the column indices to
              replace
    values  - a sequence of doubles that represent the values to replace
")
Epetra_CrsMatrix::ReplaceMyValues(PyObject*,PyObject*,PyObject*);

%feature("autodoc",
"SumIntoGlobalValues(self, int globalRow, PySequence values, PySequence
    indices) -> int

  Arguments:

    globalRow - global row index
    values    - a sequence of doubles that represent the values to sum into
    indices   - a sequence of integers that represent the indices to sum into
")
Epetra_CrsMatrix::SumIntoGlobalValues(int,double*,int,int*,int);

%feature("autodoc",
"SumIntoGlobalValues(self, PySequence rows, PySequence cols, PySequence
    values) -> int

  Arguments:

    rows    - a sequence of integers that represent the row indices to sum into
    cols    - a sequence of integers that represent the column indices to
              sum into
    values  - a sequence of doubles that represent the values to sum into
")
Epetra_CrsMatrix::SumIntoGlobalValues(PyObject*,PyObject*,PyObject*);

%feature("autodoc",
"SumIntoMyValues(self, int myRow, PySequence values, PySequence indices) -> int

  Arguments:

    myRow     - local row index
    values    - a sequence of doubles that represent the values to sum into
    indices   - a sequence of integers that represent the indices to sum into
")
Epetra_CrsMatrix::SumIntoMyValues(int,double*,int,int*,int);

%feature("autodoc",
"SumIntoMyValues(self, PySequence rows, PySequence cols, PySequence
    values) -> int

  Arguments:

    rows    - a sequence of integers that represent the row indices to sum into
    cols    - a sequence of integers that represent the column indices to
              sum into
    values  - a sequence of doubles that represent the values to sum into
")
Epetra_CrsMatrix::SumIntoValues(PyObject*,PyObject*,PyObject*);

%feature("autodoc",
"__setitem__(self, PyTuple index, double val)

The __setitem__() method is called when square-bracket indexing is
used to set a value of the matrix.  For example, the last line of::

    comm = Epetra.SerialComm()
    m = Epetra.CrsMatrix(9,0,comm)
    m[0,0] = 3.14

calls::

    m.__setitem__((0,0), 3.14)

Thus, argument 'index' is a tuple filled with whatever indices you
give the square-bracket operator when setting.  For __setitem__(),
this raises an IndexError unless 'index' is a two-tuple of integers.
Argument 'val' must be convertible to a double.  Under the covers,
__setitem__() calls ReplaceGlobalValues() or InsertGlobalValues() as
necessary, so the indices are expected to be global IDs.  Note that if
you use __setitem__() to insert a new matrix element, you will need to
call FillComplete() again, whether or not you have called it before.
")
Epetra_CrsMatrix::__setitem__;

%feature("autodoc",
"
__getitem__(self, PyTuple index) -> double
__getitem__(self, int row) -> numpy.ndarray

The __getitem__() method is called when square-bracket indexing is
used to get a value from the matrix.  For example, the last two lines
of::

    comm = Epetra.SerialComm()
    m = Epetra.CrsMatrix(9,0,comm)
    m.InsertGlobalValues(0, [0.0, 1.0, 2.0], [0,1,2])
    diag = m[0,0]
    row  = m[0]

call::

    m.__getitem__((0,0))
    m.__getitem__(0)

The __getitem__() method behaves according to the following table:

                    FillComplete()    #    
    Index               called      procs  Return value
    --------------  --------------  -----  ---------------------------
    single integer       true        any   numpy array of doubles
    single integer       false        1    numpy array of doubles
    single integer       false       >1    raise IndexError
    two integers         either      any   double

You should provide global IDs as the integer indices if FillComplete()
has been called.  If not, you should provide local IDs.  If you
reference a matrix element that is off-processor, __getitem__() will
raise an IndexError.

Under the covers, __getitem__() will call ExtractGlobalRowView() if
FillComplete() has been called, or ExtractMyRowView() if it has not.
If either of these return a non-zero return code, this is converted to
a python RuntimeError.  The resulting data is copied to the output
array.
")
Epetra_CrsMatrix::__getitem__;

// Epetra.VbrMatrix class

%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, int
    numBlockEntriesPerRow) -> VbrMatrix

  VbrMatrix constructor with implicit column map and constant number
  of block entries per row.
")
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess,
			           const Epetra_BlockMap&,
				   int);

%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, PySequence
    numBlockEntriesPerRow) -> VbrMatrix

  VbrMatrix constructor with implicit column map and variable number
  of block entries per row.

")
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess,
				   const Epetra_BlockMap&,
				   int*, int);

%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, BlockMap colMap,
    int numBlockEntriesPerRow) -> VbrMatrix

  VbrMatrix constructor with specified column map and constant number
  of block entries per row.

")
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess,
			           const Epetra_BlockMap&,
			           const Epetra_BlockMap&,
				   int);

%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, BlockMap colMap,
    PySequence numBlockEntriesPerRow) -> VbrMatrix

  VbrMatrix constructor with specified column map and variable number
  of block entries per row.

")
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess,
				   const Epetra_BlockMap&,
				   const Epetra_BlockMap&,
				   int*, int);

%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, CrsGraph graph) -> VbrMatrix

  CrsGraph constructor.
")
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess,
				   const Epetra_CrsGraph&);

%feature("autodoc",
"
__init__(self, VbrMatrix matrix) -> VbrMatrix

  Copy constructor.
")
Epetra_VbrMatrix::Epetra_VbrMatrix(const Epetra_VbrMatrix&);

// Epetra.Operator class

%feature("docstring")
Epetra_Operator
"
For cross-language polymorphism to work in python, you must call this
constructor::

    from PyTrilinos import Epetra
    class MyOperator(Epetra.Operator):
        def __init__(self):
            Epetra.Operator.__init__(self)

Other than that, the Epetra.Operator class is much more forgiving than
its C++ counterpart.  Often, you can override just the Label() and
Apply() methods.
"

%feature("autodoc",
"Apply(self, MultiVector x, MultiVector y) -> int

In C++, the Apply() method is pure virtual, thus intended to be
overridden by derived classes.  In python, cross-language polymorphism
is supported, and you are expected to derive classes from this base
class and redefine the Apply() method.  C++ code (e.g., AztecOO
solvers) can call back to your Apply() method as needed.  You must
support two arguments, labeled here MultiVector x and MultiVector y.
These will be converted from Epetra_MultiVector C++ objects to
numpy-hybrid Epetra.MultiVector objects before they are passed to you.
Thus, it is legal to use slice indexing and other numpy features to
compute y from x.

If application of your operator is successful, return 0; else return
some non-zero error code.

It is strongly suggested that you prevent Apply() from raising any
exceptions.  Accidental errors can be prevented by wrapping your code
in a try block:

    try:
        # Your code goes here...
    except Exception, e:
        print 'A python exception was raised by method Apply:'
        print e
        return -1

By returning a -1, you inform the calling routine that Apply() was
unsuccessful.
")
Epetra_Operator::Apply;

%feature("autodoc",
"ApplyInverse(self, MultiVector x, MultiVector y) -> int

In C++, the ApplyInverse() method is pure virtual, thus intended to be
overridden by derived classes.  In python, cross-language polymorphism
is supported, and you are expected to derive classes from this base
class and redefine the ApplyInverse() method.  C++ code (e.g., AztecOO
solvers) can call back to your ApplyInverse() method as needed.  You
must support two arguments, labeled here MultiVector x and MultiVector
y.  These will be converted from Epetra_MultiVector C++ objects to
numpy-hybrid Epetra.MultiVector objects before they are passed to you.
Thus, it is legal to use slice indexing and other numpy features to
compute y from x.

If application of your operator is successful, return 0; else return
some non-zero error code.

It is strongly suggested that you prevent ApplyInverse() from raising
any exceptions.  Accidental errors can be prevented by wrapping your
code in a try block:

    try:
        # Your code goes here...
    except Exception, e:
        print 'A python exception was raised by method ApplyInverse:'
        print e
        return -1

By returning a -1, you inform the calling routine that ApplyInverse()
was unsuccessful.
")
Epetra_Operator::ApplyInverse;

// Epetra.RowMatrix class

%feature("autodoc",
"NumMyRowEntries(int myRow, numpy.ndarray numEntries) -> int

In C++, numEntries in an int&.  In python, it is provided to you as a
numpy array of length one so that you can set its value in-place using
numEntries[0] = ....")
Epetra_RowMatrix::NumMyRowEntries;

%feature("autodoc",
"ExtractMyRowCopy(int myRow, int length, numpy.ndarray numEntries,
    numpy.ndarray values, numpy.ndarray indices) -> int

In C++, numEntries in an int&.  In python, it is provided to you as a
numpy array of length one so that you can set its value in-place using
numEntries[0] = ....

Arguments values and indices are double* and int*, respectively, in
C++.  In python, these are provided to you as numpy arrays of the
given length, so that you may alter their entries in-place.")
Epetra_RowMatrix::ExtractMyRowCopy;

%feature("autodoc",
"ExtractDiagonalCopy(Vector diagonal) -> int

Argument diagonal is provided to you as a numpy-hybrid Epetra.Vector,
giving you access to the numpy interface in addition to the
Epetra_Vector C++ interface.")
Epetra_RowMatrix::ExtractDiagonalCopy;

%feature("autodoc",
"Multiply(bool useTranspose, MultiVector x, MultiVector y) -> int

In C++, arguments x and y are Epetra_MultiVectors.  In python, they
are provided to you as numpy-hybrid Epetra.MultiVectors, giving you
access to the numpy interface in addition to the Epetra_MultiVector
C++ interface.")
Epetra_RowMatrix::Multiply;

%feature("autodoc",
"Solve((bool upper, bool trans, bool unitDiagonal, MultiVector x,
    MultiVector y) -> int

In C++, arguments x and y are Epetra_MultiVectors.  In python, they
are provided to you as numpy-hybrid Epetra.MultiVectors, giving you
access to the numpy interface in addition to the Epetra_MultiVector
C++ interface.")
Epetra_RowMatrix::Solve;

%feature("autodoc",
"InvRowSum(Vector x) -> int

Argument x is provided to you as a numpy-hybrid Epetra.Vector, giving
you access to the numpy interface in addition to the Epetra_Vector C++
interface.")
Epetra_RowMatrix::InvRowSum;

%feature("autodoc",
"LeftScale(Vector x) -> int

Argument x is provided to you as a numpy-hybrid Epetra.Vector, giving
you access to the numpy interface in addition to the Epetra_Vector C++
interface.")
Epetra_RowMatrix::LeftScale;

%feature("autodoc",
"InvColSums(Vector x) -> int

Argument x is provided to you as a numpy-hybrid Epetra.Vector, giving
you access to the numpy interface in addition to the Epetra_Vector C++
interface.")
Epetra_RowMatrix::InvColSums;

%feature("autodoc",
"InvRowSums(Vector x) -> int

Argument x is provided to you as a numpy-hybrid Epetra.Vector, giving
you access to the numpy interface in addition to the Epetra_Vector C++
interface.")
Epetra_RowMatrix::InvRowSums;

%feature("autodoc",
"RightScale(Vector x) -> int

Argument x is provided to you as a numpy-hybrid Epetra.Vector, giving
you access to the numpy interface in addition to the Epetra_Vector C++
interface.")
Epetra_RowMatrix::RightScale;
