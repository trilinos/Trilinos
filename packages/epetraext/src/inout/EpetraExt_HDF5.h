//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef EPETRAEXT_HDF5_H
#define EPETRAEXT_HDF5_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

#include "EpetraExt_ConfigDefs.h"
#ifdef HAVE_EPETRAEXT_HDF5

#include <string>
#include "hdf5.h"
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Comm;
class Epetra_IntVector;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_RowMatrix;
class Epetra_CrsMatrix;
class Epetra_VbrMatrix;
namespace Teuchos {
  class ParameterList;
}
namespace EpetraExt {
  class Handle;
  template<class T>
  class DistArray;
}

namespace EpetraExt 
{
/*! 
\brief class HDF5: A class for storing Epetra objects in parallel binary files

<P><B>Introduction</B>

<P>The HDF5 class reads and writes data using the HDF5 parallel binary data
format. HDF5 has the following advantages:
- the file format is binary and portable. Being binary, its size is
  very compact;
- the file format can be read and written in parallel;
- MATLAB contains a built-in HDF5 data reader, making it easy to
  interface Trilinos applications with MATLAB and vice-versa;
- the HDF5 library may already be installed on your system; otherwise,
  it is easy to build from the source tarball on the official NCSA
  HDF5 web page;
- HDF5 objects are like a file system, and they can contain
  directories and files. This means that more data can be stored in
  the same file.  One can also append new data to an existing file, or
  read only one data set;
- It is very easy to add metadata to your data; for example a string
  describing the machine where the code was executed, the parameter
  list used for the preconditioner or the nonlinear solver, the
  convergence history, and so on;
- The number of processors reading a given object does not have to be
  the same as the number used to write the object.

In order to use HDF5, make sure to set the CMake Boolean option
EpetraExt_USING_HDF5 when building Trilinos.  You may have to tell
CMake where to find the HDF5 library and header files, if they are not
installed in their default locations.

The class supports input (I), output (O), or both (I/O) for the
following distributed Epetra objects:
- Epetra_Map (I/O);
- Epetra_BlockMap (I/O);
- Epetra_CrsGraph (I/O);
- Epetra_RowMatrix (O);
- Epetra_CrsMatrix (I/O);
- Epetra_IntVector (I/O);
- Epetra_MultiVector (I/O).

The class also supports some non-distributed types:
- Teuchos::ParameterList (I/O);
- basic data types like scalar int and double variables (I/O);
- and arrays of arbitrary data type are supported (I/O).

The HDF5 class assumes that non-distributed data types have the same
value on all processors.

This class also provides utility methods:
- CreateGroup() creates the specified group in the file;
- IsContained() returns \c true is the specified group is already
  contained in the file;
- WriteComment() allows to write any string as comment for a group;
- a method to write and read a distributed array of arbitrary type.

By using these methods, as well as the other methods to write
non-distributed types, one can read and write any serial or
distributed object.

<P><B>Data Model</B>

The HDF5 library itself can be used to define very general data
formats.  Our HDF5 class, instead, is structured around the concept of
\e groups. A \e group is an entity, for example a scalar value, an
Epetra_Map, or a Teuchos::ParameterList. Within each group, different
datasets describe the content of the group. For example, an
Epetra_MultiVector is specified by datasets \c NumVectors and \c
Values, which contain the number of vectors and the numerical values,
respectively. The \c comment of each group is a character string that
must match the class name.

Our HDF5 class has the following limitations:
- Objects stored in a file cannot be deleted.  If you want to do that,
  you should read the content of the file, then create a new file and
  store in it the information to keep;
- it is not possible to overwrite distributed objects.


<P><B>Errors</B>

When an error occurs, a EpetraExt::Exception is thrown.  The Print()
method of the Exception class returns a description of what went
wrong.


<P><B>Example of usage</B>

First, one must create an HDF5 class, then either Open() or Create()
the file:
\code
// Comm is an Epetra_Comm communicator wrapper object.
EpetraExt::HDF5 HDF5(Comm);
HDF5.Create("myfile.h5");
\endcode

Writing commands might be as follows:
\code
Epetra_Map* Map = <create map here>
Epetra_Map* BlockMap = <create block map here>
Epetra_RowMatrix* Matrix = <create matrix here>
Epetra_MultiVector* LHS = <...>
Epetra_MultiVector* RHS = <...>

// write a map, whose group name contains the number of processors
HDF5.Write("map-" + toString(Comm.NumProc()), *Map);
HDF5.Write("matrix", *Matrix);
HDF5.Write("LHS", LHS);
HDF5.Write("RHS", RHS);
\endcode

To write a Teuchos::ParameterList, simply do
\code
HDF5.Write("EpetraExtList", EpetraExtList);

\endcode
The file can contain metadata as well. Each metadatum is defined by a
group and a dataset name. A group may contain more than one dataset,
and may be a new group or an already existing group. For example, to
specify the numerical quadrature formula used to assemble the matrix,
do as follows:
\code
HDF5.Write("matrix", "quadrature order", 3);
\endcode
Alternatively, datasets may be assigned to a new group, let's say \c
"my parameters":
\code
HDF5.Write("my parameters", "latitude", 12);
HDF5.Write("my parameters", "longitude", 67);
HDF5.Write("my parameters", "angle", 12.3);

vector<int> iarray(3); 
iarray[0] = 0, iarray[1] = 1; iarray[2] = 2;
HDF5.Write("my parameters", "int array", H5T_NATIVE_INT, 3, &iarray[0]);

vector<double> darray(3); 
darray[0] = 0.1, darray[1] = 1.1; darray[2] = 2.1;
HDF5.Write("my parameters", "double array", H5T_NATIVE_DOUBLE, 3, &darray[0]);
\endcode

Note that all non-distributed datasets must have the same value on all
processors.

Reading data is as easy as writing. Let us consider how to read an
Epetra_CrsMatrix.  Other Epetra objects having a similar behavior. The
ReadCrsMatrixProperties() method can be used to query for some matrix
properties without reading the whole matrix:
\code
int NumGlobalRows, NumGlobalCols, NumGlobalNonzeros;
int NumGlobalDiagonals, MaxNumEntries;
double NormOne, NormInf;

ReadCrsMatrixProperties(GroupName, NumGlobalRows, NumGlobalCols,
                        NumGlobalNonzeros, NumGlobalDiagonals, MaxNumEntries,
                        NormOne, NormInf);
\endcode

The above call is not required, and can be skipped.  In order to read
the Epetra_CrsMatrix, do as follows:
\code
Epetra_CrsMatrix* NewMatrix = NULL;
HDF5.Read("matrix", NewMatrix);
\endcode

In this case, \c NewMatrix is based on a linear map. If the matrix's
DomainMap() and RangeMap() are known and non-trivial, one can use
\code
HDF5.Read("matrix", DomainMap, RangeMap, NewMatrix);
\endcode

Reading metadata looks like:
\code
HDF5.Read("my parameters", "latitude", new_latitude);
HDF5.Read("my parameters", "longitude", new_longitude);
HDF5.Read("my parameters", "int array", H5T_NATIVE_INT, 3, &new_iarray[0]);
HDF5.Read("my parameters", "double array", H5T_NATIVE_DOUBLE, 3, &new_darray[0]);
\endcode

To analyze the content of the file, one can use 
\c "h5dump filename.h5" or "h5dump filename.h5 -H".


<P><B>MATLAB Interface</B>

MATLAB provides built-in functions to read, write, and query HDF5
files: \c hdf5read, \c hdf5write, and \c hdf5info, respectively.  For
example, to read the above \c Epetra_CrsMatrix into MATLAB as a sparse
matrix, do as follows:
\code
NumGlobalRows = double(hdf5read('myfile.h5', '/matrix/NumGlobalRows/'));
NumGlobalCols = double(hdf5read('myfile.h5', '/matrix/NumGlobalCols/'));
ROW = double(hdf5read('myfile.h5', '/matrix/ROW/'));
COL = double(hdf5read('myfile.h5', '/matrix/COL/'));
VAL = hdf5read('myfile.h5', '/matrix/VAL/');
A = sparse(ROW + 1, COL + 1, VAL, NumGlobalRows, NumGlobalCols);
\endcode
The use of \c double() is required by Matlab's \c sparse function,
since it does not accept \c int32 data.

To dump a MATLAB matrix (in this case, \c A) to the file \c matlab.h5,
do as follows:
\code
n = 10;
A = speye(n, n);
[ROW,COL,VAL] = find(A);
hdf5write('matlab.h5', '/speye/__type__',           'Epetra_RowMatrix');
hdf5write('matlab.h5', '/speye/NumGlobalRows',      int32(n), 'WriteMode', 'append');
hdf5write('matlab.h5', '/speye/NumGlobalCols',      int32(n), 'WriteMode', 'append');
hdf5write('matlab.h5', '/speye/NumGlobalNonzeros',  int32(n), 'WriteMode', 'append');
hdf5write('matlab.h5', '/speye/NumGlobalDiagonals', int32(n), 'WriteMode', 'append');
hdf5write('matlab.h5', '/speye/MaxNumEntries',      int32(1), 'WriteMode', 'append');
hdf5write('matlab.h5', '/speye/NormOne',            1.0,      'WriteMode', 'append');
hdf5write('matlab.h5', '/speye/NormInf',            1.0,      'WriteMode', 'append');
hdf5write('matlab.h5', '/speye/ROW', int32(ROW - 1), 'WriteMode', 'append');
hdf5write('matlab.h5', '/speye/COL', int32(COL - 1), 'WriteMode', 'append');
hdf5write('matlab.h5', '/speye/VAL', VAL,            'WriteMode', 'append');
\endcode
Note that the \c __type__ specification must reflect the Epetra class name.

To dump a MATLAB dense array (in this case, \c x) to the file \c matlab.h5, do as follows:
\code
n = 10;
x = [zeros(n,1), rand(n, 1)]';
hdf5write('matlab.h5', '/x/__type__',    'Epetra_MultiVector');
hdf5write('matlab.h5', '/x/GlobalLength',int32(n), 'WriteMode', 'append');
hdf5write('matlab.h5', '/x/NumVectors',  int32(2), 'WriteMode', 'append');
hdf5write('matlab.h5', '/x/Values',      x,        'WriteMode', 'append');
\endcode
Note that MATLAB vectors must be stored as \e row vectors.  

You can also write a Map from MATLAB for use by Epetra.  The following
example shows how to define an Epetra_Map that distributes data over
two processors:
\code
IndexBase = 0;
NumMyElements = [5 5];
n = 10;
MyGlobalElements = [5 6 7 8 9 0 1 2 3 4];
hdf5write('matlab.h5', '/map-2/__type__',          'Epetra_Map');
hdf5write('matlab.h5', '/map-2/NumGlobalElements', int32(n),                'WriteMode', 'append');
hdf5write('matlab.h5', '/map-2/IndexBase',         int32(IndexBase),        'WriteMode', 'append');
hdf5write('matlab.h5', '/map-2/NumProc',           int32(2),                'WriteMode', 'append');
hdf5write('matlab.h5', '/map-2/NumMyElements',     int32(NumMyElements),    'WriteMode', 'append');
hdf5write('matlab.h5', '/map-2/MyGlobalElements',  int32(MyGlobalElements), 'WriteMode', 'append');
\endcode


\author Marzio Sala, D-INFK/ETHZ

\date Last updated on 16-Mar-06.  Documentation revised by Mark Hoemmen (05 Oct 2011) for spelling, grammar, and clarity.

\todo 
- all distributed objects are assumed in local state, this is not necessary (just easier)
- Epetra_VbrMatrix input has to be done, now it is considered as an Epetra_RowMatrix

*/
class HDF5 
{
  public: 
    // @{ \name Constructor and destructor.

    //! Constructor
    HDF5(const Epetra_Comm& Comm); 

    //! Destructor
    ~HDF5() 
    {
      if (IsOpen())
        Close();
    }

    // @}
    // @{ \name Basic operations
    
    //! Create a new file.
    void Create(const std::string FileName);

    //! Open specified file with given access type.
    void Open(const std::string FileName, int AccessType = H5F_ACC_RDWR);

    //! Close the file.
    void Close()
    {
      H5Fclose(file_id_);
      IsOpen_ = false;
    }
    
    //! Flush the content to the file
    void Flush()
    {
      H5Fflush(file_id_, H5F_SCOPE_GLOBAL);
    }

    //! Return \c true if a file has already been opened using Open()/Create()
    bool IsOpen() const
    {
      return(IsOpen_);
    }

    //! Create group \c GroupName.
    void CreateGroup(const std::string& GroupName)
    {
      hid_t group_id = H5Gcreate(file_id_, GroupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group_id);
    }

    //! Return \c true if \c Name is contained in the database.
    bool IsContained(std::string Name, std::string GroupName = "");

    // @}
    // @{ \name basic non-distributed data types
    
    //! Write an integer in group \c GroupName using the given \c DataSetName.
    void Write(const std::string& GroupName, const std::string& DataSetName, int data);

    //! Read an integer from group \c /GroupName/DataSetName
    void Read(const std::string& GroupName, const std::string& DataSetName, int& data);

    //! Write a double in group \c GroupName using the given \c DataSetName.
    void Write(const std::string& GroupName, const std::string& DataSetName, double data);

    //! Read a double from group \c /GroupName/DataSetName
    void Read(const std::string& GroupName, const std::string& DataSetName, double& data);

    //! Write a string in group \c GroupName using the given \c DataSetName.
    void Write(const std::string& GroupName, const std::string& DataSetName, const std::string& data);

    //! Read a string from group \c /GroupName/DataSetName
    void Read(const std::string& GroupName, const std::string& DataSetName, std::string& data);

    //! Read the serial array \c data, of type \c type, from group \c GroupName, using the dataset name \c DataSetName.
    void Read(const std::string& GroupName, const std::string& DataSetName,
              const hid_t type, const int Length, void* data);

    //! Write the serial array \c data, of type \c type, to group \c GroupName, using the dataset name \c DataSetName
    void Write(const std::string& GroupName, const std::string& DataSetName,
                         const hid_t type, const int Length, 
                         const void* data);

    //! Associate string \c Comment with group \c GroupName.
    void WriteComment(const std::string& GroupName, std::string Comment)
    {
      H5Gset_comment(file_id_, GroupName.c_str(), Comment.c_str());
    }

    //! Read the string associated with group \c GroupName.
    void ReadComment(const std::string& GroupName, std::string& Comment)
    {
      char comment[128];
      H5Gget_comment(file_id_, GroupName.c_str(), 128, comment);
      Comment = comment;
    }

    // @}
    // @{ \name Distributed arrays
    
    //! Write the distributed array \c data, of type \c type, to group \c GroupName, using dataset name \c DataSetName
    void Write(const std::string& GroupName, const std::string& DataSetName, int MySize, int GlobalSize, hid_t type, const void* data);

    //! Read the distributed array \c data, of type \c type, from group \c GroupName, using dataset name \c DataSetName
    void Read(const std::string& GroupName, const std::string& DataSetName,
              int MySize, int GlobalSize,
              const hid_t type, void* data);

    // @}
    // @{ \name Epetra_Map/Epetra_BlockMap

    //! Write a Map to group \c GroupName.
    void Write(const std::string& GroupName, const Epetra_Map& Map);

    //! Read a map from \c GroupName.
    void Read(const std::string& GroupName, Epetra_Map*& Map);

    //! Read basic properties of specified Epetra_Map.
    void ReadMapProperties(const std::string& GroupName, 
                           int& NumGlobalElements,
                           int& IndexBase,
                           int& NumProc);

    //! Read a block map from \c GroupName.
    void Read(const std::string& GroupName, Epetra_BlockMap*& Map);

    //! Write a block map to group \c GroupName.
    void Write(const std::string& GroupName, const Epetra_BlockMap& Map);

    //! Read basic properties of specified Epetra_BlockMap.
    void ReadBlockMapProperties(const std::string& GroupName, 
                                int& NumGlobalElements,
                                int& NumGlobalPoints,
                                int& IndexBase,
                                int& NumProc);

    // @}
    // @{ \name Epetra_CrsGraph

    //! Read a vector from group \c GroupName, assuming linear distribution.
    void Read(const std::string& GroupName, Epetra_CrsGraph*& Graph);

    //! Read a vector from group \c GroupName using the given map.
    void Read(const std::string& GroupName, const Epetra_Map& DomainMap, 
              const Epetra_Map& RangeMap, Epetra_CrsGraph*& Graph);

    //! Write a distributed vector to group \c GroupName.
    void Write(const std::string& GroupName, const Epetra_CrsGraph& Graph);

    //! Read basic properties of specified Epetra_CrsGraph.
    void ReadCrsGraphProperties(const std::string& GroupName, 
                                int& NumGlobalRows,
                                int& NumGlobalCols,
                                int& NumGlobalNonzeros,
                                int& NumGlobalDiagonals,
                                int& MaxNumIndices);

    // @}
    // @{ \name Epetra_IntVector

    //! Write a distributed vector to group \c GroupName.
    void Write(const std::string& GroupName, const Epetra_IntVector& x);

    //! Read a vector from group \c GroupName, assuming linear distribution.
    void Read(const std::string& GroupName, Epetra_IntVector*& X);

    //! Read a vector from group \c GroupName using the given map.
    void Read(const std::string& GroupName, const Epetra_Map& Map, Epetra_IntVector*& X);

    //! Read basic properties of specified Epetra_IntVector.
    void ReadIntVectorProperties(const std::string& GroupName, int& GlobalLength);

    // @}
    // @{ \name Epetra_MultiVector

    /// \brief Write a distributed vector (or its transpose) to group \c GroupName.  
    ///
    /// Write the transpose if writeTranspose is true.
    void Write(const std::string& GroupName, const Epetra_MultiVector& x, bool writeTranspose = false);

    /// \brief Read a vector (or its transpose) from group \c GroupName.
    /// 
    /// This method assumes a linear distribution.  Read the transpose
    /// if writeTranspose is true.
    void Read(const std::string& GroupName, Epetra_MultiVector*& X,
              bool writeTranspose = false, const int& indexBase = 0);

    /// \brief Read a vector from group \c GroupName using the given map. 
    ///
    /// Read the transpose if writeTranspose is true.
    void Read(const std::string& GroupName, const Epetra_Map& Map, Epetra_MultiVector*& X,
              bool writeTranspose = false);

    //! Read basic properties of specified Epetra_MultiVector.
    void ReadMultiVectorProperties(const std::string& GroupName, 
                                   int& GlobalLength,
                                   int& NumVectors);

    // @}
    // @{ \name Epetra_RowMatrix/Epetra_CrsMatrix

    //! Write a distributed RowMatrix to group \c GroupName.
    void Write(const std::string& GroupName, const Epetra_RowMatrix& Matrix);

    //! Read a square matrix from group \c GroupName, assuming linear distribution.
    void Read(const std::string& GroupName, Epetra_CrsMatrix*& A);

    //! Read a matrix from group \c GroupName with given range and domain maps.
    void Read(const std::string& GroupName, 
              const Epetra_Map& DomainMap, 
              const Epetra_Map& RangeMap, 
              Epetra_CrsMatrix*& A);

    //! Read basic properties of specified Epetra_CrsMatrix.
    void ReadCrsMatrixProperties(const std::string& GroupName, 
                                 int& NumGlobalRows,
                                 int& NumGlobalCols,
                                 int& NumNonzeros,
                                 int& NumGlobalDiagonals,
                                 int& MaxNumEntries,
                                 double& NormOne,
                                 double& NormInf);

    // @}
    // @{ \name Teuchos::ParameterList

    //! Write a parameter list to group \c GroupName.
    void Write(const std::string& GroupName, const Teuchos::ParameterList& List);

    //! Read a parameter list from group \c GroupName.
    void Read(const std::string& GroupName, Teuchos::ParameterList& List);

    // @}
    // @{ \name EpetraExt::DistArray<int>

    //! Write an EpetraExt::DistArray<int> to group \c GroupName.
    void Write(const std::string& GroupName, const DistArray<int>& array);

    //! Read an EpetraExt::DistArray<int> from group \c GroupName.
    void Read(const std::string& GroupName, DistArray<int>*& array);

    //! Read an EpetraExt::DistArray<int> from group \c GroupName.
    void Read(const std::string& GroupName, const Epetra_Map& Map, DistArray<int>*& array);

    //! Read the global number of elements and type for a generic handle object
    void ReadIntDistArrayProperties(const std::string& GroupName, 
                                    int& GlobalLength,
                                    int& RowSize);

    // @}
    // @{ \name EpetraExt::DistArray<double>

    //! Write an EpetraExt::DistArray<int> to group \c GroupName.
    void Write(const std::string& GroupName, const DistArray<double>& array);

    //! Read an EpetraExt::DistArray<int> from group \c GroupName.
    void Read(const std::string& GroupName, DistArray<double>*& array);

    //! Read an EpetraExt::DistArray<int> from group \c GroupName.
    void Read(const std::string& GroupName, const Epetra_Map& Map, DistArray<double>*& array);

    //! Read the global number of elements and type for a generic handle object
    void ReadDoubleDistArrayProperties(const std::string& GroupName, 
                                       int& GlobalLength,
                                       int& RowSize);
    // @}
    // @}
    // @{ \name Generic distributed object

    //! Write an Epetra_DistObject to group \c GroupName.
    void Write(const std::string& GroupName, const Handle& List);

    //! Read an Epetra_DistObject from group \c GroupName.
    void Read(const std::string& GroupName, Handle& List);

    //! Read the global number of elements and type for a generic handle object
    void ReadHandleProperties(const std::string& GroupName, 
                              std::string& Type,
                              int& NumGlobalElements);

    // @}
  private:
    // @{ \name Private Data

    //! Returns a reference to this object's communicator.
    const Epetra_Comm& Comm() const
    {
      return(Comm_);
    }

    //! This object's communicator.
    const Epetra_Comm& Comm_; 
    //! FileName currently open.
    std::string FileName_;
    //! If \c true, a file is currently open.
    bool IsOpen_;

    //! file ID for HDF5.
    hid_t       file_id_;
    hid_t	plist_id_;
    herr_t	status;

    // @}
};
}
#endif
#endif /* EPETRAEXT_HDF5_H */
