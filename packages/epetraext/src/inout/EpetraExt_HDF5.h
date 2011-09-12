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

#include "EpetraExt_ConfigDefs.h"
#ifdef HAVE_EPETRAEXT_HDF5

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

<P>Class HDF5 allows to read and write using the HDF5 parallel binary data
format. HDF5 has the following advantages:
- the file format is binary and portable. Being binary, its size is very compact;
- the file format can be read and written in parallel;
- MATLAB contains a built-in HDF5 data reader, making it easy to interface Trilinos
  application with MATLAB, and vice-versa;
- the library is easy to build: just download the tarball from the official NCSA HDF5 
  web page, configure/make/make install, and everything should work. While building
  Trilinos, then you just have to specify the location of the header files using
  --with-incdirs, the location of the library using --with-ldflags, and the HDF5 library 
  using --with-libs. Note that you might need to add "-lz" to the list of libraries; 
  this library is typically already installed on your system;
- HDF5 are like a file system, and they can contain directories and files. This means
  that more data can be stored in the same file. Besides, one can append new data to
  an existing file, or read only one data set;
- It is very easy to add "meta-data" to your data; for example a string containing
  information about the machine where the code was executed, the parameter list used
  for the preconditioner or the nonlinear solver, the convergence history, and so on;
- The number of processors reading a given data does not have to coincide with that
  used to write the data.

The class supports I/O for the following distributed Epetra objects:
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
It is meant that the non-distributed data types have the same value on all processors.

Some utility methods are provided:
- CreateGroup() creates the specified group in the file;
- IsContained() returns \c true is the specified group is already contained in the file;
- WriteComment() allows to write any string as comment for a group;
- a method to write and read a distributed array of arbitrary type is provided. 

By using these methods, as well as the other methods to write non-distributed types, one 
can read and write any serial or distributed object.

<P><B>Data Model</B>

The HDF5 library itself can be used to define very general data formats; this class, 
instead, is only structured around the concept of \e groups. A \e group is an entity,
like for example a scalar value, an Epetra_Map, or a Teuchos::ParameterList. Within each
group, different datasets describe the content of the group. For example, an 
Epetra_MultiVector is specified by datasets \c NumVectors and \c Values, which contain 
the number of vectors, and the numerical values, respectively. The \c comment of each group 
is a character string that must match the class name.

The HDF5 class has the following limitations:
- Objects stored in a file cannot be deleted; if you want to do that, you should read the 
  content of the file, then create a new file and store on it the information to keep;
- it is not possible to overwrite distributed objects.


<P><B>Errors</B>

When an error occurs, a EpetraExt::Exception is thrown. Method Print() of the Exception class
gives a description of what went wrong.


<P><B>Example of usage</B>

First, one has to create an HDF5 class, then either Open() or Create() the file:
\code
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
The file can contain meta-data as well. Each meta-data is defined by a group,
and a dataset name. A group can contain more than one dataset, and can be a new group
or an already existing group. For example, to specify the numerical quadrature formula
used to assemble the matrix, one can do as follows:
\code
HDF5.Write("matrix", "quadrature order", 3);
\endcode
Alternatively, datasets can be assigned to a new group, let's say \c "my parameters":
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

Note that all non-distributed datasets are supposed to have the same value
on all processors.

Reading data is a easy as writing. Let us consider how to read an Epetra_CrsMatrix, 
other Epetra objects having a similar bejavior. Method ReadCrsMatrixProperties()
can be used to query for some matrix properties, without reading the whole matrix:
\code
int NumGlobalRows, NumGlobalCols, NumGlobalNonzeros;
int NumGlobalDiagonals, MaxNumEntries;
double NormOne, NormInf;

ReadCrsMatrixProperties(GroupName, NumGlobalRows, NumGlobalCols,
                        NumGlobalNonzeros, NumGlobalDiagonals, MaxNumEntries,
                        NormOne, NormInf);
\endcode

The above call is not required, and can be skipped. The actual reading is:
\code
Epetra_CrsMatrix* NewMatrix = 0;
HDF5.Read("matrix", NewMatrix);
\endcode

In this case, \c NewMatrix is based on a linear map. If the DomainMap() and RangeMap() are
known and non-trivial, one can use
\code
HDF5.Read("matrix", DomainMap, RangeMap, NewMatrix);
\endcode

Reading meta-data looks like:
\code
HDF5.Read("my parameters", "latitude", new_latitude);
HDF5.Read("my parameters", "longitude", new_longitude);
HDF5.Read("my parameters", "int array", H5T_NATIVE_INT, 3, &new_iarray[0]);
HDF5.Read("my parameters", "double array", H5T_NATIVE_DOUBLE, 3, &new_darray[0]);
\endcode

To analyze the content of the file, one can use 
\c "h5dump filename.h5" or "h5dump filename.h5 -H".


<P><B>MATLAB Interface</B>

Reading HDF5 files from MATLAB is very, since the built-in functions 
\c hdf5read, \c hdf5write, \c hdf5info. For example, to read the above \c Matrix from
MATLAB, one can do:
\code
NumGlobalRows = double(hdf5read('myfile.h5', '/matrix/NumGlobalRows/'));
NumGlobalCols = double(hdf5read('myfile.h5', '/matrix/NumGlobalCols/'));
ROW = double(hdf5read('myfile.h5', '/matrix/ROW/'));
COL = double(hdf5read('myfile.h5', '/matrix/COL/'));
VAL = hdf5read('myfile.h5', '/matrix/VAL/');
A = sparse(ROW + 1, COL + 1, VAL, NumGlobalRows, NumGlobalCols);
\endcode
The use of \c double() is required by \c sparse, which does not accept \c int32 data.

To dump on file \c matlab.h5 a MATLAB matrix (in this case, \c A), one can proceed as 
follows:
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
Note that \c __type__ specification, that must reflect the Epetra class name.

To dump on file \c matlab.h5 a MATLAB dense array (in this case, \c x), one can proceed as follows:
\code
n = 10;
x = [zeros(n,1), rand(n, 1)]';
hdf5write('matlab.h5', '/x/__type__',    'Epetra_MultiVector');
hdf5write('matlab.h5', '/x/GlobalLength',int32(n), 'WriteMode', 'append');
hdf5write('matlab.h5', '/x/NumVectors',  int32(2), 'WriteMode', 'append');
hdf5write('matlab.h5', '/x/Values',      x,        'WriteMode', 'append');
\endcode
Note that MATLAB vectors must be stored as \e row vectors.  

To write a Map from MATLAB follows a very similar pattern. The following example shows how to define a map that can be used with two processors:
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

\date Last updated on 16-Mar-06.

\todo 
- all distributed objects are assumed in local state, this is not necessary (just easier)
- Epetra_VbrMatrix input has to be done, now it is considered as an Epetra_RowMatrix

*/
class HDF5 
{
  public: 
    // @{ \name Constructor and destructor.
    //! ctor
    HDF5(const Epetra_Comm& Comm); 

    //! dtor
    ~HDF5() 
    {
      if (IsOpen())
        Close();
    }

    // @}
    // @{ \name Basic operations
    
    //! Creates a new file.
    void Create(const string FileName);

    //! Opens specified file with given access type.
    void Open(const string FileName, int AccessType = H5F_ACC_RDWR);

    //! Closes the file.
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

    //! Returns \c true is a file has already been open using Open()/Create()
    bool IsOpen() const
    {
      return(IsOpen_);
    }

    //! Creates group \c GroupName.
    void CreateGroup(const string& GroupName)
    {
      hid_t group_id = H5Gcreate(file_id_, GroupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group_id);
    }

    //! Returns \c true is \c Name is contained in the database.
    bool IsContained(const string Name);

    // @}
    // @{ \name basic non-distributed data types
    
    //! Writes an integer in group \c GroupName using intentified \c DataSetName.
    void Write(const string& GroupName, const string& DataSetName, int data);

    //! Reads an integer from group \c /GroupName/DataSetName
    void Read(const string& GroupName, const string& DataSetName, int& data);

    //! Writes a double in group \c GroupName using intentified \c DataSetName.
    void Write(const string& GroupName, const string& DataSetName, double data);

    //! Reads a double from group \c /GroupName/DataSetName
    void Read(const string& GroupName, const string& DataSetName, double& data);

    //! Writes a string in group \c GroupName using intentified \c DataSetName.
    void Write(const string& GroupName, const string& DataSetName, const string& data);

    //! Reads a string from group \c /GroupName/DataSetName
    void Read(const string& GroupName, const string& DataSetName, string& data);

    //! Reads serial array \c data, of type \c type, from group \c GroupNameusing dataset name \c DataSetName
    void Read(const string& GroupName, const string& DataSetName,
              const int type, const int Length, void* data);

    //! Writes serial array \c data, of type \c type, to group \c GroupNameusing dataset name \c DataSetName
    void Write(const string& GroupName, const string& DataSetName,
                         const int type, const int Length, 
                         void* data);

    //! Associates string \c Comment with group \c GroupName.
    void WriteComment(const string& GroupName, string Comment)
    {
      H5Gset_comment(file_id_, GroupName.c_str(), Comment.c_str());
    }

    //! Reads the string associated with group \c GroupName.
    void ReadComment(const string& GroupName, string& Comment)
    {
      char comment[128];
      H5Gget_comment(file_id_, GroupName.c_str(), 128, comment);
      Comment = comment;
    }

    // @}
    // @{ \name Distributed arrays
    
    //! Writes distributed array \c data, of type \c type, to group \c GroupNameusing dataset name \c DataSetName
    void Write(const string& GroupName, const string& DataSetName, int MySize, int GlobalSize, int type, const void* data);

    //! Reads distributed array \c data, of type \c type, from group \c GroupNameusing dataset name \c DataSetName
    void Read(const string& GroupName, const string& DataSetName,
              int MySize, int GlobalSize,
              const int type, void* data);

    // @}
    // @{ \name Epetra_Map/Epetra_BlockMap

    //! Writes a Map to group \c GroupName.
    void Write(const string& GroupName, const Epetra_Map& Map);

    //! Reads a map from \c GroupName.
    void Read(const string& GroupName, Epetra_Map*& Map);

    //! Reads basic properties of specified Epetra_Map.
    void ReadMapProperties(const string& GroupName, 
                           int& NumGlobalElements,
                           int& IndexBase,
                           int& NumProc);

    //! Reads a block map from \c GroupName.
    void Read(const string& GroupName, Epetra_BlockMap*& Map);

    //! Writes a block map to group \c GroupName.
    void Write(const string& GroupName, const Epetra_BlockMap& Map);

    //! Reads basic properties of specified Epetra_BlockMap.
    void ReadBlockMapProperties(const string& GroupName, 
                                int& NumGlobalElements,
                                int& NumGlobalPoints,
                                int& IndexBase,
                                int& NumProc);

    // @}
    // @{ \name Epetra_CrsGraph

    //! Reads a vector from group \c GroupName, assumes linear distribution.
    void Read(const string& GroupName, Epetra_CrsGraph*& Graph);

    //! Reads a vector from group \c GroupName using given map.
    void Read(const string& GroupName, const Epetra_Map& DomainMap, 
              const Epetra_Map& RangeMap, Epetra_CrsGraph*& Graph);

    //! Writes a distributed vector to group \c GroupName.
    void Write(const string& GroupName, const Epetra_CrsGraph& Graph);

    //! Reads basic properties of specified Epetra_CrsGraph.
    void ReadCrsGraphProperties(const string& GroupName, 
                                int& NumGlobalRows,
                                int& NumGlobalCols,
                                int& NumGlobalNonzeros,
                                int& NumGlobalDiagonals,
                                int& MaxNumIndices);

    // @}
    // @{ \name Epetra_IntVector

    //! Writes a distributed vector to group \c GroupName.
    void Write(const string& GroupName, const Epetra_IntVector& x);

    //! Reads a vector from group \c GroupName, assumes linear distribution.
    void Read(const string& GroupName, Epetra_IntVector*& X);

    //! Reads a vector from group \c GroupName using given map.
    void Read(const string& GroupName, const Epetra_Map& Map, Epetra_IntVector*& X);

    //! Reads basic properties of specified Epetra_IntVector.
    void ReadIntVectorProperties(const string& GroupName, int& GlobalLength);

    // @}
    // @{ \name Epetra_MultiVector

    //! Writes a distributed vector to group \c GroupName, write transpose if writeTranspose set to true.
    void Write(const string& GroupName, const Epetra_MultiVector& x, bool writeTranspose = false);

    //! Reads a vector from group \c GroupName, assumes linear distribution. Read transpose if writeTranspose set to true.
    void Read(const string& GroupName, Epetra_MultiVector*& X,
              bool writeTranspose = false, const int& indexBase = 0);

    //! Reads a vector from group \c GroupName using given map. Read transpose if writeTranspose set to true.
    void Read(const string& GroupName, const Epetra_Map& Map, Epetra_MultiVector*& X,
              bool writeTranspose = false);

    //! Reads basic properties of specified Epetra_MultiVector.
    void ReadMultiVectorProperties(const string& GroupName, 
                                   int& GlobalLength,
                                   int& NumVectors);

    // @}
    // @{ \name Epetra_RowMatrix/Epetra_CrsMatrix

    //! Writes a distributed RowMatrix to group \c GroupName.
    void Write(const string& GroupName, const Epetra_RowMatrix& Matrix);

    //! Reads a square matrix from group \c GroupName, assumes linear distribution.
    void Read(const string& GroupName, Epetra_CrsMatrix*& A);

    //! Reads a matrix from group \c GroupName with given range and domain maps.
    void Read(const string& GroupName, 
              const Epetra_Map& DomainMap, 
              const Epetra_Map& RangeMap, 
              Epetra_CrsMatrix*& A);

    //! Reads basic properties of specified Epetra_CrsMatrix.
    void ReadCrsMatrixProperties(const string& GroupName, 
                                 int& NumGlobalRows,
                                 int& NumGlobalCols,
                                 int& NumNonzeros,
                                 int& NumGlobalDiagonals,
                                 int& MaxNumEntries,
                                 double& NormOne,
                                 double& NormInf);

    // @}
    // @{ \name Teuchos::ParameterList

    //! Writes a parameter list to group \c GroupName.
    void Write(const string& GroupName, const Teuchos::ParameterList& List);

    //! Reads a parameter list from group \c GroupName.
    void Read(const string& GroupName, Teuchos::ParameterList& List);

    // @}
    // @{ \name EpetraExt::DistArray<int>

    //! Writes an EpetraExt::DistArray<int> to group \c GroupName.
    void Write(const string& GroupName, const DistArray<int>& array);

    //! Reads an EpetraExt::DistArray<int> from group \c GroupName.
    void Read(const string& GroupName, DistArray<int>*& array);

    //! Reads an EpetraExt::DistArray<int> from group \c GroupName.
    void Read(const string& GroupName, const Epetra_Map& Map, DistArray<int>*& array);

    //! Reads the global number of elements and type for a generic handle object
    void ReadIntDistArrayProperties(const string& GroupName, 
                                    int& GlobalLength,
                                    int& RowSize);

    // @}
    // @{ \name EpetraExt::DistArray<double>

    //! Writes an EpetraExt::DistArray<int> to group \c GroupName.
    void Write(const string& GroupName, const DistArray<double>& array);

    //! Reads an EpetraExt::DistArray<int> from group \c GroupName.
    void Read(const string& GroupName, DistArray<double>*& array);

    //! Reads an EpetraExt::DistArray<int> from group \c GroupName.
    void Read(const string& GroupName, const Epetra_Map& Map, DistArray<double>*& array);

    //! Reads the global number of elements and type for a generic handle object
    void ReadDoubleDistArrayProperties(const string& GroupName, 
                                       int& GlobalLength,
                                       int& RowSize);
    // @}
    // @}
    // @{ \name Generic distributed object

    //! Writes an Epetra_DistObject to group \c GroupName.
    void Write(const string& GroupName, const Handle& List);

    //! Reads an Epetra_DistObject from group \c GroupName.
    void Read(const string& GroupName, Handle& List);

    //! Reads the global number of elements and type for a generic handle object
    void ReadHandleProperties(const string& GroupName, 
                              string& Type,
                              int& NumGlobalElements);

    // @}
  private:
    // @{ \name Private Data

    //! Returns a reference to the communicator of \c this object.
    const Epetra_Comm& Comm() const
    {
      return(Comm_);
    }

    //! Communicator of this object.
    const Epetra_Comm& Comm_; 
    //! FileName currently open.
    string FileName_;
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
