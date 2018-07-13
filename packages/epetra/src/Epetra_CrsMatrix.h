/*
//@HEADER
// ************************************************************************
//
//               Epetra: Linear Algebra Services Package
//                 Copyright 2011 Sandia Corporation
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
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_CRSMATRIX_H
#define EPETRA_CRSMATRIX_H

#include "Epetra_ConfigDefs.h"
#include "Epetra_DistObject.h"
#include "Epetra_CompObject.h"
#include "Epetra_BLAS.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Map.h"

#ifdef Epetra_ENABLE_CASK
#include "cask.h"
#endif

class Epetra_Map;
class Epetra_Import;
class Epetra_Export;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_IntSerialDenseVector;


// Define this to see a complete dump a an Epetra_CrsMatrix::Multiply(...) call
//#define EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY

#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
extern bool Epetra_CrsMatrixTraceDumpMultiply;
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY

//! Epetra_CrsMatrix: A class for constructing and using real-valued double-precision sparse compressed row matrices.

/*! The Epetra_CrsMatrix class is a sparse compressed row matrix object. This matrix can be
  used in a parallel setting, with data distribution described by Epetra_Map attributes.
  The structure or graph of the matrix is defined by an Epetra_CrsGraph attribute.

 In addition to coefficient access, the primary operations provided by Epetra_CrsMatrix are matrix
 times vector and matrix times multi-vector multiplication.

 Epetra_CrsMatrix matrices can be square or rectangular.

<b>Creating and filling Epetra_CrsMatrix objects</b>

Constructing Epetra_CrsMatrix objects is a multi-step process.  The basic steps are as follows:
<ol>
  <li> Create Epetra_CrsMatrix instance, including storage, via one of the constructors:
    <ul>
     <li>Constructor that accepts one Epetra_Map object, a row-map defining the distribution of matrix rows.
     <li>Constructor that accepts two Epetra_Map objects. (The second map is a column-map, and describes the set
       of column-indices that appear in each processor's portion of the matrix. Generally these are
       overlapping sets -- column-indices may appear on more than one processor.)
     <li>Constructor that accepts an Epetra_CrsGraph object, defining the non-zero structure of the matrix.
    </ul>
  Note that the constructors which accept Epetra_Map arguments also accept an argument that gives an
  estimate of the number of nonzeros per row. This allows storage to be pre-allocated and can improve
  the performance of the data input methods. The estimate need not be accurate, as additional storage is
  allocated automatically when needed. However, a more accurate estimate helps performance by reducing
  the amount of extra memory allocation.
  <li> Enter values via one or more Insert/Replace/SumInto functions.
  <li> Complete construction by calling FillComplete.
</ol>

Note that, even after a matrix is constructed (FillComplete has been called), it is possible to update existing
matrix entries.  It is \e not possible to create new entries.

<b>Epetra_Map attributes</b>

Epetra_CrsMatrix objects have four Epetra_Map attributes, which are held by the Epetra_CrsGraph attribute.

The Epetra_Map attributes can be obtained via these accessor methods:
<ul>
 <li>RowMap() Describes the numbering and distribution of the rows of the matrix. The row-map exists and is valid
 for the entire life of the matrix. The set of matrix rows is defined by the row-map and may not be changed. Rows
 may not be inserted or deleted by the user. The only change that may be made is that the user can replace the
 row-map with a compatible row-map (which is the same except for re-numbering) by calling the ReplaceRowMap() method.
 <li>ColMap() Describes the set of column-indices that appear in the rows in each processor's portion of the matrix.
 Unless provided by the user at construction time, a valid column-map doesn't exist until FillComplete() is called.
 <li>RangeMap() Describes the range of the matrix operator. e.g., for a matrix-vector product operation, the result
   vector's map must be compatible with the range-map of this matrix. The range-map is usually the same as the row-map.
   The range-map is set equal to the row-map at matrix creation time, but may be specified by the user when
   FillComplete() is called.
 <li>DomainMap() Describes the domain of the matrix operator. The domain-map can be specified by the user when
 FillComplete() is called. Until then, it is set equal to the row-map.
</ul>

It is important to note that while the row-map and the range-map are often the same, the column-map and the domain-map
are almost never the same. The set of entries in a distributed column-map almost always form overlapping sets, with
entries being associated with more than one processor. A domain-map, on the other hand, must be a 1-to-1 map, with
entries being associated with only a single processor.

<b>Local versus Global Indices</b>

Epetra_CrsMatrix has query functions IndicesAreLocal() and IndicesAreGlobal(), which are used to determine whether the
underlying Epetra_CrsGraph attribute's column-indices have been transformed into a local index space or not. (This
transformation occurs when the method Epetra_CrsGraph::FillComplete() is called, which happens when the
method Epetra_CrsMatrix::FillComplete() is called.) The state of the indices in the
graph determines the behavior of many Epetra_CrsMatrix methods. If an Epetra_CrsMatrix instance is constructed using
one of the constructors that does not accept a pre-existing Epetra_CrsGraph object, then an Epetra_CrsGraph attribute
is created internally and its indices remain untransformed (IndicesAreGlobal()==true) until Epetra_CrsMatrix::FillComplete()
is called. The query function Epetra_CrsMatrix::Filled() returns true if Epetra_CrsMatrix::FillComplete() has been
called.

Note the following method characteristics:

<ul>
 <li>InsertGlobalValues() may only be used to insert new nonzeros in the matrix if indices are global.
 <li>SumIntoGlobalValues() may be used regardless of whether indices are global or local, but can only be used
  to update matrix locations that already exist; it can never be used to establish new nonzero locations.
 <li>ReplaceGlobalValues() may also be used only to update matrix locations that already exist, and works
  regardless of whether indices are local or global.
 <li>SumIntoMyValues() and ReplaceMyValues() may only be used if indices are local.
 <li>Multiply() may only be used after FillComplete() has been called.
</ul>

Most methods have preconditions documented, check documentation for specific methods not mentioned here.

<b> Counting Floating Point Operations </b>

Each Epetra_CrsMatrix object keeps track of the number
of \e serial floating point operations performed using the specified object as the \e this argument
to the function.  The Flops() function returns this number as a double precision number.  Using this
information, in conjunction with the Epetra_Time class, one can get accurate parallel performance
numbers.  The ResetFlops() function resets the floating point counter.

\warning A Epetra_Map is required for the Epetra_CrsMatrix constructor.

*/

class EPETRA_LIB_DLL_EXPORT Epetra_CrsMatrix: public Epetra_DistObject, public Epetra_CompObject, public Epetra_BLAS, public virtual Epetra_RowMatrix {
 public:

   //! @name Constructors/Destructor
  //@{
  //! Epetra_CrsMatrix constructor with variable number of indices per row.
  /*! Creates a Epetra_CrsMatrix object and allocates storage.

  \param CV - (In) An Epetra_DataAccess enumerated type set to Copy or View.
  \param RowMap - (In) An Epetra_Map defining the numbering and distribution of matrix rows.
  \param NumEntriesPerRow - (In) An integer array of length NumRows
  such that NumEntriesPerRow[i] indicates the (approximate if StaticProfile=false) number of entries in the ith row.
  \param StaticProfile - (In) Optional argument that indicates whether or not NumIndicesPerRow should be interpreted as an exact
  count of nonzeros, or should be used as an approximation.  By default this value is false, allowing the profile to be determined
  dynamically.  If the user sets it to true, then the memory allocation for the Epetra_CrsGraph object will be done in one large
  block, saving on memory fragmentation and generally improving the performance of matrix multiplication and solve kernels.
  */
  Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const int* NumEntriesPerRow, bool StaticProfile = false);

  //! Epetra_CrsMatrix constructor with fixed number of indices per row.
  /*! Creates a Epetra_CrsMatrix object and allocates storage.

  \param CV - (In) An Epetra_DataAccess enumerated type set to Copy or View.
  \param RowMap - (In) An Epetra_Map defining the numbering and distribution of matrix rows.
  \param NumEntriesPerRow - (In) An integer that indicates the (approximate) number of entries in the each row.
  Note that it is possible to use 0 for this value and let fill occur during the insertion phase.
  \param StaticProfile - (In) Optional argument that indicates whether or not NumIndicesPerRow should be interpreted as an exact
  count of nonzeros, or should be used as an approximation.  By default this value is false, allowing the profile to be determined
  dynamically.  If the user sets it to true, then the memory allocation for the Epetra_CrsGraph object will be done in one large
  block, saving on memory fragmentation and generally improving the performance of matrix multiplication and solve kernels.

  */
  Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, int NumEntriesPerRow, bool StaticProfile = false);

  //! Epetra_CrsMatrix constructor with variable number of indices per row.
  /*! Creates a Epetra_CrsMatrix object and allocates storage.

  \param CV - (In) An Epetra_DataAccess enumerated type set to Copy or View.
  \param RowMap - (In) An Epetra_Map defining the numbering and distribution of matrix rows.
  \param ColMap - (In) An Epetra_Map defining the set of column-indices that appear in each processor's
  locally owned matrix rows.
  \param NumEntriesPerRow - (In) An integer array of length NumRows
  such that NumEntriesPerRow[i] indicates the (approximate if StaticProfile=false) number of entries in the ith row.
  \param StaticProfile - (In) Optional argument that indicates whether or not NumIndicesPerRow should be interpreted as an exact
  count of nonzeros, or should be used as an approximation.  By default this value is false, allowing the profile to be determined
  dynamically.  If the user sets it to true, then the memory allocation for the Epetra_CrsGraph object will be done in one large
  block, saving on memory fragmentation and generally improving the performance of matrix multiplication and solve kernels.
  */
  Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, const int* NumEntriesPerRow, bool StaticProfile = false);

  //! Epetra_CrsMatrix constuctor with fixed number of indices per row.
  /*! Creates a Epetra_CrsMatrix object and allocates storage.

  \param CV - (In) An Epetra_DataAccess enumerated type set to Copy or View.
  \param RowMap - (In) An Epetra_Map defining the numbering and distribution of matrix rows.
  \param ColMap - (In) An Epetra_Map defining the set of column-indices that appear in each processor's
  locally owned matrix rows.
  \param NumEntriesPerRow - (In) An integer that indicates the (approximate if StaticProfile=false) number of entries in the each row.
  Note that it is possible to use 0 for this value and let fill occur during the insertion phase.
  \param StaticProfile - (In) Optional argument that indicates whether or not NumIndicesPerRow should be interpreted as an exact
  count of nonzeros, or should be used as an approximation.  By default this value is false, allowing the profile to be determined
  dynamically.  If the user sets it to true, then the memory allocation for the Epetra_CrsGraph object will be done in one large
  block, saving on memory fragmentation and generally improving the performance of matrix multiplication and solve kernels.

  */
  Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, int NumEntriesPerRow, bool StaticProfile = false);

  //! Construct a matrix using an existing Epetra_CrsGraph object.
  /*! Allows the nonzero structure from another matrix, or a structure that was
    constructed independently, to be used for this matrix.
    \param CV - (In) An Epetra_DataAccess enumerated type set to Copy or View.
    \param Graph - (In) A Epetra_CrsGraph object, constructed directly or extracted from another Epetra matrix object.
  */

  Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_CrsGraph& Graph);

//! Epetra CrsMatrix constructor that also fuses Import and FillComplete().
  /*!
    A common use case is to create an empty destination Epetra_CrsMatrix,
    redistribute from a source CrsMatrix (by an Import or Export
    operation), then call FillComplete() on the destination
    CrsMatrix.  This constructor fuses these three cases, for an
    Import redistribution.

    Fusing redistribution and FillComplete() exposes potential
    optimizations.  For example, it may make constructing the column
    map faster, and it may avoid intermediate unoptimized storage in
    the destination Epetra_CrsMatrix.  These optimizations may improve
    performance for specialized kernels like sparse matrix-matrix
    multiply, as well as for redistributing data after doing load
    balancing.

    The resulting matrix is fill complete (in the sense of
    Filled()) and has optimized storage (in the sense of
    StorageOptimized()).  It the DomainMap is taken from the SourceMatrix,
    the RangeMap is presumed to be RowImporter.TargetMap() if not specified

    \param SourceMatrix [in] The source matrix from which to
    import.  The source of an Import must have a nonoverlapping
    distribution.

    \param RowImporter [in] The Import instance containing a
    precomputed redistribution plan.  The source Map of the
    Import must be the same as the row Map of sourceMatrix.

    \param DomainMap [in] The new domainMap for the new matrix. If not specified,
    then the DomainMap of the SourceMatrix is used.

    \param RangeMap [in] The new rangeMap for the new matrix. If not specified,
    then RowImporter.TargetMap() is used.

    \param RestrictCommunicator [in] Restricts the resulting communicator to active
    processes only.
  */
  Epetra_CrsMatrix(const Epetra_CrsMatrix & SourceMatrix, const Epetra_Import & RowImporter, const Epetra_Map * DomainMap=0, const Epetra_Map * RangeMap=0, bool RestrictCommunicator = false);

  //! Epetra CrsMatrix constructor that also fuses Import and FillComplete().
  /*!
    A common use case is to create an empty destination Epetra_CrsMatrix,
    redistribute from a source CrsMatrix (by an Import or Export
    operation), then call FillComplete() on the destination
    CrsMatrix.  This constructor fuses these three cases, for an
    Import redistribution.

    Fusing redistribution and FillComplete() exposes potential
    optimizations.  For example, it may make constructing the column
    map faster, and it may avoid intermediate unoptimized storage in
    the destination Epetra_CrsMatrix.  These optimizations may improve
    performance for specialized kernels like sparse matrix-matrix
    multiply, as well as for redistributing data after doing load
    balancing.

    The resulting matrix is fill complete (in the sense of
    Filled()) and has optimized storage (in the sense of
    StorageOptimized()).  It the DomainMap is taken from the SourceMatrix,
    the RangeMap is presumed to be RowImporter.TargetMap() if not specified

    \param SourceMatrix [in] The source matrix from which to
    import.  The source of an Import must have a nonoverlapping
    distribution.

    \param RowImporter [in] The Import instance containing a
    precomputed redistribution plan.  The source Map of the
    Import must be the same as the row Map of sourceMatrix.

    \param DomainImporter [in] The Import instance containing a
    precomputed redistribution plan (for the domain maps).
    The source Map of the Import must be the same as the domain
    Map of sourceMatrix.

    \param DomainMap [in] The new domainMap for the new matrix.

    \param RangeMap [in] The new rangeMap for the new matrix.

    \param RestrictCommunicator [in] Restricts the resulting communicator to active
    processes only.
  */
  Epetra_CrsMatrix(const Epetra_CrsMatrix & SourceMatrix, const Epetra_Import & RowImporter, const Epetra_Import * DomainImporter, const Epetra_Map * DomainMap, const Epetra_Map * RangeMap, bool RestrictCommunicator);

  //! Epetra CrsMatrix constructor that also fuses Ex[prt and FillComplete().
  /*!
    A common use case is to create an empty destination Epetra_CrsMatrix,
    redistribute from a source CrsMatrix (by an Import or Export
    operation), then call FillComplete() on the destination
    CrsMatrix.  This constructor fuses these three cases, for an
    Import redistribution.

    Fusing redistribution and FillComplete() exposes potential
    optimizations.  For example, it may make constructing the column
    map faster, and it may avoid intermediate unoptimized storage in
    the destination Epetra_CrsMatrix.  These optimizations may improve
    performance for specialized kernels like sparse matrix-matrix
    multiply, as well as for redistributing data after doing load
    balancing.

    The resulting matrix is fill complete (in the sense of
    Filled()) and has optimized storage (in the sense of
    StorageOptimized()).  It the DomainMap is taken from the SourceMatrix,
    the RangeMap is presumed to be RowImporter.TargetMap() if not specified

    \param SourceMatrix [in] The source matrix from which to
    import.  The source of an Import must have a nonoverlapping
    distribution.

    \param RowExporter [in] The Export instance containing a
    precomputed redistribution plan.  The source Map of the
    Import must be the same as the row Map of sourceMatrix.

    \param DomainMap [in] The new domainMap for the new matrix. If not specified,
    then the DomainMap of the SourceMatrix is used.

    \param RangeMap [in] The new rangeMap for the new matrix. If not specified,
    then RowExporter.TargetMap() is used.

    \param RestrictCommunicator [in] Restricts the resulting communicator to active
    processes only.
  */
  Epetra_CrsMatrix(const Epetra_CrsMatrix & SourceMatrix, const Epetra_Export & RowExporter, const Epetra_Map * DomainMap=0, const Epetra_Map * RangeMap=0, bool RestrictCommunicator = false);

  //! Epetra CrsMatrix constructor that also fuses Ex[prt and FillComplete().
  /*!
    A common use case is to create an empty destination Epetra_CrsMatrix,
    redistribute from a source CrsMatrix (by an Import or Export
    operation), then call FillComplete() on the destination
    CrsMatrix.  This constructor fuses these three cases, for an
    Import redistribution.

    Fusing redistribution and FillComplete() exposes potential
    optimizations.  For example, it may make constructing the column
    map faster, and it may avoid intermediate unoptimized storage in
    the destination Epetra_CrsMatrix.  These optimizations may improve
    performance for specialized kernels like sparse matrix-matrix
    multiply, as well as for redistributing data after doing load
    balancing.

    The resulting matrix is fill complete (in the sense of
    Filled()) and has optimized storage (in the sense of
    StorageOptimized()).  It the DomainMap is taken from the SourceMatrix,
    the RangeMap is presumed to be RowImporter.TargetMap() if not specified

    \param SourceMatrix [in] The source matrix from which to
    import.  The source of an Import must have a nonoverlapping
    distribution.

    \param RowExporter [in] The Export instance containing a
    precomputed redistribution plan.  The source Map of the
    Import must be the same as the row Map of sourceMatrix.

    \param DomainExporter [in] The Export instance containing a
    precomputed redistribution plan (for the domain map.
    The source Map of the Import must be the same as the domain Map
    of sourceMatrix.

    \param DomainMap [in] The new domainMap for the new matrix.

    \param RangeMap [in] The new rangeMap for the new matrix.

    \param RestrictCommunicator [in] Restricts the resulting communicator to active
    processes only.
  */
  Epetra_CrsMatrix(const Epetra_CrsMatrix & SourceMatrix, const Epetra_Export & RowExporter, const Epetra_Export * DomainExporter, const Epetra_Map * DomainMap, const Epetra_Map * RangeMap, bool RestrictCommunicator);


  //! Copy constructor.
  Epetra_CrsMatrix(const Epetra_CrsMatrix& Matrix);

  //! Epetra_CrsMatrix Destructor
  virtual ~Epetra_CrsMatrix();
  //@}

  //! @name Insertion/Replace/SumInto methods
  //@{

  //! Assignment operator
  Epetra_CrsMatrix& operator=(const Epetra_CrsMatrix& src);

  //! Initialize all values in the matrix with constant value.
  /*!
    \param ScalarConstant - (In) Value to use.

    \return Integer error code, set to 0 if successful.
    \pre None.
    \post All values in \e this set to ScalarConstant.
  */
  int PutScalar(double ScalarConstant);

  //! Multiply all values in the matrix by a constant value (in place: A <- ScalarConstant * A).
  /*!
    \param ScalarConstant - (In) Value to use.

    \return Integer error code, set to 0 if successful.
    \pre None.
    \post All values of \e this have been multiplied by ScalarConstant.
  */
  int Scale(double ScalarConstant);

  //! Insert a list of elements in a given global row of the matrix.
  /*!
    This method is used to construct a matrix for the first time.  It cannot
    be used if the matrix structure has already been fixed (via a call to FillComplete()).
    If multiple values are inserted for the same matrix entry, the values are initially
    stored separately, so memory use will grow as a result.  However, when FillComplete is called
    the values will be summed together and the additional memory will be released.

    For example, if the values 2.0, 3.0 and 4.0 are all inserted in Row 1, Column 2, extra storage
    is used to store each of the three values separately.  In this way, the insert process does not
    require any searching and can be faster.  However, when FillComplete() is called, the values
    will be summed together to equal 9.0 and only a single entry will remain in the matrix for
    Row 1, Column 2.

    \param GlobalRow - (In) Row number (in global coordinates) to put elements.
    \param NumEntries - (In) Number of entries.
    \param Values - (In) Values to enter.
    \param Indices - (In) Global column indices corresponding to values.

    \return Integer error code, set to 0 if successful. Note that if the
    allocated length of the row has to be expanded, a positive warning code
    will be returned.

    \warning This method may not be called once FillComplete() has been called.

    \pre IndicesAreLocal()==false && IndicesAreContiguous()==false
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  virtual int InsertGlobalValues(int GlobalRow, int NumEntries, const double* Values, const int* Indices);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  virtual int InsertGlobalValues(long long GlobalRow, int NumEntries, const double* Values, const long long* Indices);
#endif

  //! Insert a list of elements in a given global row of the matrix.
  /*!
    This method is used to construct a matrix for the first time.  It cannot
    be used if the matrix structure has already been fixed (via a call to FillComplete()).
    If multiple values are inserted for the same matrix entry, the values are initially
    stored separately, so memory use will grow as a result.  However, when FillComplete is called
    the values will be summed together and the additional memory will be released.

    For example, if the values 2.0, 3.0 and 4.0 are all inserted in Row 1, Column 2, extra storage
    is used to store each of the three values separately.  In this way, the insert process does not
    require any searching and can be faster.  However, when FillComplete() is called, the values
    will be summed together to equal 9.0 and only a single entry will remain in the matrix for
    Row 1, Column 2.

    \param GlobalRow - (In) Row number (in global coordinates) to put elements.
    \param NumEntries - (In) Number of entries.
    \param Values - (In) Values to enter.
    \param Indices - (In) Global column indices corresponding to values.

    \return Integer error code, set to 0 if successful. Note that if the
    allocated length of the row has to be expanded, a positive warning code
    will be returned.

    \warning This method may not be called once FillComplete() has been called.

    \pre IndicesAreLocal()==false && IndicesAreContiguous()==false
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  virtual int InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  virtual int InsertGlobalValues(long long GlobalRow, int NumEntries, double* Values, long long* Indices);
#endif

  //! Replace specified existing values with this list of entries for a given global row of the matrix.
  /*!
    \param GlobalRow - (In) Row number (in global coordinates) to put elements.
    \param NumEntries - (In) Number of entries.
    \param Values - (In) Values to enter.
    \param Indices - (In) Global column indices corresponding to values.

    \return Integer error code, set to 0 if successful. Note that if a value
    is not already present for the specified location in the matrix, the
    input value will be ignored and a positive warning code will be returned.

    \pre IndicesAreLocal()==false && IndicesAreContiguous()==false
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  virtual int ReplaceGlobalValues(int GlobalRow, int NumEntries, const double* Values, const int* Indices);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  virtual int ReplaceGlobalValues(long long GlobalRow, int NumEntries, const double* Values, const long long* Indices);
#endif

  //! Add this list of entries to existing values for a given global row of the matrix.
  /*!
    \param GlobalRow - (In) Row number (in global coordinates) to put elements.
    \param NumEntries - (In) Number of entries.
    \param Values - (In) Values to enter.
    \param Indices - (In) Global column indices corresponding to values.

    \return Integer error code, set to 0 if successful. Note that if a value
    is not already present for the specified location in the matrix, the
    input value will be ignored and a positive warning code will be returned.

    \pre IndicesAreLocal()==false && IndicesAreContiguous()==false
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  virtual int SumIntoGlobalValues(int GlobalRow, int NumEntries, const double* Values, const int* Indices);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  virtual int SumIntoGlobalValues(long long GlobalRow, int NumEntries, const double* Values, const long long* Indices);
#endif

  //! Insert a list of elements in a given local row of the matrix.
  /*!
    \param MyRow - (In) Row number (in local coordinates) to put elements.
    \param NumEntries - (In) Number of entries.
    \param Values - (In) Values to enter.
    \param Indices - (In) Local column indices corresponding to values.

    \return Integer error code, set to 0 if successful. Note that if the
    allocated length of the row has to be expanded, a positive warning code
    will be returned.
    \pre IndicesAreGlobal()==false && (IndicesAreContiguous()==false || CV_==View)
    \post The given local row of the matrix has been updated as described above.

  */
  int InsertMyValues(int MyRow, int NumEntries, const double* Values, const int* Indices);

  //! Insert a list of elements in a given local row of the matrix.
  /*!
    \param MyRow - (In) Row number (in local coordinates) to put elements.
    \param NumEntries - (In) Number of entries.
    \param Values - (In) Values to enter.
    \param Indices - (In) Local column indices corresponding to values.

    \return Integer error code, set to 0 if successful. Note that if the
    allocated length of the row has to be expanded, a positive warning code
    will be returned.
    \pre IndicesAreGlobal()==false && (IndicesAreContiguous()==false || CV_==View)
    \post The given local row of the matrix has been updated as described above.

  */
  int InsertMyValues(int MyRow, int NumEntries, double* Values, int* Indices);

  //! Replace current values with this list of entries for a given local row of the matrix.
  /*!
    \param MyRow - (In) Row number (in local coordinates) to put elements.
    \param NumEntries - (In) Number of entries.
    \param Values - (In) Values to enter.
    \param Indices - (In) Local column indices corresponding to values.

    \return Integer error code, set to 0 if successful. Note that if a value
    is not already present for the specified location in the matrix, the
    input value will be ignored and a positive warning code will be returned.
    \pre IndicesAreLocal()==true
    \post MyRow contains the given list of Values at the given Indices.
  */
  int ReplaceMyValues(int MyRow, int NumEntries, const double* Values, const int* Indices);

  //! Add this list of entries to existing values for a given local row of the matrix.
  /*!
    \param MyRow - (In) Row number (in local coordinates) to put elements.
    \param NumEntries - (In) Number of entries.
    \param Values - (In) Values to enter.
    \param Indices - (In) Local column indices corresponding to values.

    \return Integer error code, set to 0 if successful. Note that if the
    allocated length of the row has to be expanded, a positive warning code
    will be returned.
    \pre IndicesAreLocal()==true
    \post The given Values at the given Indices have been summed into the
    entries of MyRow.
  */
  int SumIntoMyValues(int MyRow, int NumEntries, const double* Values, const int* Indices);

  //! Replaces diagonal values of the matrix with those in the user-provided vector.
  /*! This routine is meant to allow replacement of {\bf existing} diagonal values.
    If a diagonal value does not exist for a given row, the corresponding value in
    the input Epetra_Vector will be ignored and the return code will be set to 1.

    The Epetra_Map associated with the input Epetra_Vector must be compatible with
    the RowMap of the matrix.

    \param Diagonal - (In) New values to be placed in the main diagonal.

    \return Integer error code, set to 0 if successful, set to 1 on the calling processor if one or more diagonal entries not present in matrix.
    \pre Filled()==true
    \post Diagonal values have been replaced with the values of Diagonal.
  */
  int ReplaceDiagonalValues(const Epetra_Vector& Diagonal);

  //@}

  //! @name Transformation methods
  //@{


  //! Signal that data entry is complete.  Perform transformations to local index space.
  /* This version of FillComplete assumes that the domain and range
     distributions are identical to the matrix row distributions.
    \param OptimizeDataStorage - (In) If true, storage will be packed for optimal performance.  Depending
           on how the matrix was constructed, optimizing the storage may have no impact on performance
     or one-time memory use, or may have a large impact.  If the user was careful in allocating memory
     for the matrix by setting StaticProfile to true in the matrix constructor, then no extra storage
     will be allocated in attempting to optimize storage.  If the user did not set StaticProfile to true,
     then optimizing the storage will temporarily use additional memory, will have a noticeable impact
     on performance and ultimately reduce the storage associated with the matrix.

     By default storage will be optimized.  If you cannot tolerate the increased temporary memory use,
     should set this value to false.

     \return error code, 0 if successful. Returns a positive warning code of 3
        if the matrix is rectangular (meaning that the other overloading of
        FillComplete should have been called, with differen domain-map and
        range-map specified).
  */
  int FillComplete(bool OptimizeDataStorage = true);

  //! Signal that data entry is complete.  Perform transformations to local index space.
  /* This version of FillComplete requires the explicit specification of the domain
     and range distribution maps.  These maps are used for importing and exporting vector
     and multi-vector elements that are needed for distributed matrix computations.  For
     example, to compute y = Ax in parallel, we would specify the DomainMap as the distribution
     of the vector x and the RangeMap as the distribution of the vector y.
    \param DomainMap - (In) Map that describes the distribution of vector and multi-vectors in the
    matrix domain.
    \param RangeMap - (In) Map that describes the distribution of vector and multi-vectors in the
    matrix range.

    \param OptimizeDataStorage - (In) If true, storage will be packed for optimal performance.  Depending
           on how the matrix was constructed, optimizing the storage may have no impact on performance
     or one-time memory use, or may have a large impact.  If the user was careful in allocating memory
     for the matrix by setting StaticProfile to true in the matrix constructor, then no extra storage
     will be allocated in attempting to optimize storage.  If the user did not set StaticProfile to true,
     then optimizing the storage will temporarily use additional memory, will have a noticeable impact
     on performance and ultimately reduce the storage associated with the matrix.

     By default storage will be optimized.  If you cannot tolerate the increased temporary memory use,
     should set this value to false.

    \return error code, 0 if successful. positive warning code of 2 if it is detected that the
    matrix-graph got out of sync since this matrix was constructed (for instance if
    graph.FillComplete() was called by another matrix that shares the graph)

    \post IndicesAreLocal()==true
    */
  int FillComplete(const Epetra_Map& DomainMap, const Epetra_Map& RangeMap, bool OptimizeDataStorage = true);

  //! Make consecutive row index sections contiguous, minimize internal storage used for constructing graph.
  /*! After construction and during initialization (when values are being added), the matrix coefficients
    for each row are managed as separate segments of memory. This method moves the values for all rows
    into one large contiguous array and eliminates internal storage that is not needed after matrix construction. Calling this
    method can have a significant impact on memory costs and machine performance.

    If this object was constructed in View mode then this method can't make non-contiguous values contiguous and will
    return a warning code of 1 if the viewed data isn't already contiguous.

    \note A call to this method will also call the OptimizeStorage method for the associated Epetra_CrsGraph object.  If
    the storage for this graph has already been optimized this additional call will have no effect.

    \return Integer error code, set to 0 if successful.

    \pre Filled()==true.
    \pre If CV=View when the graph was constructed, then this method will be effective \only if the indices of the graph were already contiguous.  In this case, the indices are left untouched and internal storage for the graph is minimized.

    \post StorageOptimized()==true, if successful.
    \post Graph().StorageOptimized()==true, if successful.

  */
  int OptimizeStorage();


  //! Eliminates memory that is used for construction.  Make consecutive row index sections contiguous.
  int MakeDataContiguous() {EPETRA_CHK_ERR(OptimizeStorage()); return(0);}
  //@}

  //! @name Extraction methods
  //@{

  //! Returns a copy of the specified global row in user-provided arrays.
  /*!
    \param GlobalRow - (In) Global row to extract.
    \param ILength - (In) Length of Values and Indices.
    \param NumEntries - (Out) Number of nonzero entries extracted.
    \param Values - (Out) Extracted values for this row.
    \param Indices - (Out) Extracted global column indices for the corresponding values.

    \return Integer error code, set to 0 if successful, non-zero if global row is not owned by calling process
or if the number of entries in this row exceed the Length parameter.
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int ExtractGlobalRowCopy(int GlobalRow, int Length, int& NumEntries, double* Values, int* Indices) const;
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  int ExtractGlobalRowCopy(long long GlobalRow, int Length, int& NumEntries, double* Values, long long* Indices) const;
#endif

  //! Returns a copy of the specified local row in user-provided arrays.
  /*!
    \param MyRow - (In) Local row to extract.
    \param Length - (In) Length of Values and Indices.
    \param NumEntries - (Out) Number of nonzero entries extracted.
    \param Values - (Out) Extracted values for this row.
    \param Indices - (Out) Extracted local column indices for the corresponding values.

    \return Integer error code, set to 0 if successful.

    \pre IndicesAreLocal()==true
  */
  int ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, double* Values, int* Indices) const;

  //! Returns a copy of the specified global row values in user-provided array.
  /*!
    \param GlobalRow - (In) Global row to extract.
    \param Length - (In) Length of Values.
    \param NumEntries - (Out) Number of nonzero entries extracted.
    \param Values - (Out) Extracted values for this row.

    \return Integer error code, set to 0 if successful.
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int ExtractGlobalRowCopy(int GlobalRow, int Length, int& NumEntries, double* Values) const;
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  int ExtractGlobalRowCopy(long long GlobalRow, int Length, int& NumEntries, double* Values) const;
#endif

  //! Returns a copy of the specified local row values in user-provided array.
  /*!
    \param MyRow - (In) Local row to extract.
    \param Length - (In) Length of Values.
    \param NumEntries - (Out) Number of nonzero entries extracted.
    \param Values - (Out) Extracted values for this row.

    \return Integer error code, set to 0 if successful.
  */
  int ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, double* Values) const;

  //! Returns a copy of the main diagonal in a user-provided vector.
  /*!
    \param Diagonal - (Out) Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
    \pre Filled()==true
    \post Unchanged.
  */
  int ExtractDiagonalCopy(Epetra_Vector& Diagonal) const;

  //! Returns a view of the specified global row values via pointers to internal data.
  /*!
    \param GlobalRow - (In) Global row to view.
    \param NumEntries - (Out) Number of nonzero entries extracted.
    \param Values - (Out) Extracted values for this row.
    \param Indices - (Out) Extracted global column indices for the corresponding values.

    \return Integer error code, set to 0 if successful. Returns -1 of row not on this processor.
    Returns -2 if matrix is not in global form (if FillComplete() has already been called).

    \pre IndicesAreGlobal()==true
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int ExtractGlobalRowView(int GlobalRow, int& NumEntries, double*& Values, int*& Indices) const;
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  int ExtractGlobalRowView(long long GlobalRow, int& NumEntries, double*& Values, long long*& Indices) const;
#endif

  //! Returns a view of the specified local row values via pointers to internal data.
  /*!
    \param MyRow - (In) Local row to view.
    \param NumEntries - (Out) Number of nonzero entries extracted.
    \param Values - (Out) Extracted values for this row.
    \param Indices - (Out) Extracted local column indices for the corresponding values.

    \return Integer error code, set to 0 if successful. Returns -1 of row not on this processor.
    Returns -2 if matrix is not in local form (if FillComplete() has \e not been called).

    \pre IndicesAreLocal()==true
  */
  int ExtractMyRowView(int MyRow, int& NumEntries, double*& Values, int*& Indices) const;

  //! Returns a view of the specified global row values via pointers to internal data.
  /*!
    \param GlobalRow - (In) Global row to extract.
    \param NumEntries - (Out) Number of nonzero entries extracted.
    \param Values - (Out) Extracted values for this row.

    \return Integer error code, set to 0 if successful.
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int ExtractGlobalRowView(int GlobalRow, int& NumEntries, double*& Values) const;
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  int ExtractGlobalRowView(long long GlobalRow, int& NumEntries, double*& Values) const;
#endif

  //! Returns a view of the specified local row values via pointers to internal data.
  /*!
    \param MyRow - (In) Local row to extract.
    \param NumEntries - (Out) Number of nonzero entries extracted.
    \param Values - (Out) Extracted values for this row.

    \return Integer error code, set to 0 if successful.
  */
  int ExtractMyRowView(int MyRow, int& NumEntries, double*& Values) const;
  //@}

  //! @name Computational methods
  //@{

  //! Returns the result of a Epetra_CrsMatrix multiplied by a Epetra_Vector x in y.
  /*!
    \param TransA - (In) If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param x - (In) An Epetra_Vector to multiply by.
    \param y - (Out) An Epetra_Vector containing result.

    \return Integer error code, set to 0 if successful.
    \pre Filled()==true
    \post Unchanged.
  */
  int Multiply(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const;
  int Multiply1(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const;


  //! Returns the result of a Epetra_CrsMatrix multiplied by a Epetra_MultiVector X in Y.
  /*!
    \param TransA - (In) If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param X - (In) An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Y - (Out) An Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
    \pre Filled()==true
    \post Unchanged.
  */
  int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  int Multiply1(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of a local solve using the Epetra_CrsMatrix on a Epetra_Vector x in y.
  /*! This method solves a triangular system of equations asynchronously on each processor.
    \param Upper - (In) If true, solve Uy = x, otherwise solve Ly = x.
    \param Trans - (In) If true, solve transpose problem.
    \param UnitDiagonal - (In) If true, assume diagonal is unit (whether it's stored or not).
    \param x - (In) An Epetra_Vector to solve for.
    \param y - (Out) An Epetra_Vector containing result.

    \return Integer error code, set to 0 if successful.
    \pre Filled()==true
    \post Unchanged.
  */
  int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_Vector& x, Epetra_Vector& y) const;

  //! Returns the result of a local solve using the Epetra_CrsMatrix a Epetra_MultiVector X in Y.
  /*! This method solves a triangular system of equations asynchronously on each processor.
    \param Upper - (In) If true, solve Uy = x, otherwise solve Ly = x.
    \param Trans - (In) If true, solve transpose problem.
    \param UnitDiagonal - (In) If true, assume diagonal is unit (whether it's stored or not).
    \param X - (In) An Epetra_MultiVector of dimension NumVectors to solve for.
    \param Y - (Out) An Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
    \pre Filled()==true
    \post Unchanged.
  */
  int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Computes the inverse of the sum of absolute values of the rows of the Epetra_CrsMatrix, results returned in x.
  /*! The vector x will return such that x[i] will contain the inverse of the sum of the absolute values of the entries in the
          ith row of the \e this matrix.  Using the resulting vector from this function as input to LeftScale()
                will make the infinity norm of the resulting matrix exactly 1.
                \warning The NormInf() method will not properly calculate the infinity norm for a matrix that has entries that are
                replicated on multiple processors.  In this case, if the rows are fully replicated, NormInf() will return a
                value equal to the maximum number of processors that any individual row of the matrix is replicated on.
    \param x - (Out) An Epetra_Vector containing the inverse of the row sums of the \e this matrix.
                \warning When rows are fully replicated on multiple processors, it is assumed that the distribution of x is
                the same as the rows (RowMap())of \e this.  When multiple processors contain partial sums for individual entries, the
                distribution of x is assumed to be the same as the RangeMap() of \e this.  When each row of \e this is
                uniquely owned, the distribution of x can be that of the RowMap() or the RangeMap().

    \return Integer error code, set to 0 if successful.
    \pre Filled()==true
    \post Unchanged.
  */
  int InvRowSums(Epetra_Vector& x) const;

  //! Computes the inverse of the max of absolute values of the rows of the Epetra_CrsMatrix, results returned in x.
  /*! The vector x will return such that x[i] will contain the inverse of max of the absolute values of the entries in the ith
          row of the \e this matrix.
                \warning This method will not work when multiple processors contain partial sums for individual entries.
    \param x - (Out) An Epetra_Vector containing the inverse of the row maxs of the \e this matrix.
    \warning When rows are fully replicated on multiple processors, it is assumed that the distribution of x is
                the same as the rows (RowMap())of \e this.  When each row of \e this is uniquely owned, the distribution of
    x can be that of the RowMap() or the RangeMap().

    \return Integer error code, set to 0 if successful.
    \pre Filled()==true
    \post Unchanged.
  */
  int InvRowMaxs(Epetra_Vector& x) const;

  //! Scales the Epetra_CrsMatrix on the left with a Epetra_Vector x.
  /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
    and j denotes the column number of A.
    \param x - (In) An Epetra_Vector to scale with.

    \return Integer error code, set to 0 if successful.
    \pre Filled()==true
    \post The matrix will be scaled as described above.
  */
  int LeftScale(const Epetra_Vector& x);

  //! Computes the inverse of the sum of absolute values of the columns of the Epetra_CrsMatrix, results returned in x.
  /*! The vector x will return such that x[j] will contain the inverse of the sum of the absolute values of the
    entries in the jth column of the \e this matrix.   Using the resulting vector from this function as input to
    RightScale() will make the one norm of the resulting matrix exactly 1.
                \warning The NormOne() method will not properly calculate the one norm for a matrix that has entries that are
                replicated on multiple processors.  In this case, if the columns are fully replicated, NormOne() will return a
                value equal to the maximum number of processors that any individual column of the matrix is repliated on.

    \param x - (Out) An Epetra_Vector containing the column sums of the \e this matrix.
                \warning When columns are fully replicated on multiple processors, it is assumed that the distribution of x is
                the same as the columns (ColMap()) of \e this.  When multiple processors contain partial sums for entries, the
                distribution of x is assumed to be the same as the DomainMap() of \e this.  When each column of \e this is
                uniquely owned, the distribution of x can be that of the ColMap() or the DomainMap().

    \return Integer error code, set to 0 if successful.
    \pre Filled()==true
    \post Unchanged.
  */
  int InvColSums(Epetra_Vector& x) const;

  //! Computes the max of absolute values of the columns of the Epetra_CrsMatrix, results returned in x.
  /*! The vector x will return such that x[j] will contain the inverse of max of the absolute values of the entries
    in the jth row of the \e this matrix.
                \warning This method will not work when multiple processors contain partial sums for individual entries.
    \param x - (Out) An Epetra_Vector containing the column maxs of the \e this matrix.
                \warning When columns are fully replicated on multiple processors, it is assumed that the distribution of x is
                the same as the columns (ColMap()) of \e this.  When each column of \e this is
                uniquely owned, the distribution of x can be that of the ColMap() or the DomainMap().

    \return Integer error code, set to 0 if successful.
    \pre Filled()==true
    \post Unchanged.
  */
  int InvColMaxs(Epetra_Vector& x) const;

  //! Scales the Epetra_CrsMatrix on the right with a Epetra_Vector x.
  /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
    and j denotes the global column number of A.
    \param x - (In) The Epetra_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
    \pre Filled()==true
    \post The matrix will be scaled as described above.
  */
  int RightScale(const Epetra_Vector& x);
  //@}

  //! @name Matrix Properties Query Methods
  //@{


  //! If FillComplete() has been called, this query returns true, otherwise it returns false.
  bool Filled() const {return(Graph_.Filled());}

  //! If OptimizeStorage() has been called, this query returns true, otherwise it returns false.
  bool StorageOptimized() const {return(StorageOptimized_);}

  //! If matrix indices has not been transformed to local, this query returns true, otherwise it returns false.
  bool IndicesAreGlobal() const {return(Graph_.IndicesAreGlobal());}

  //! If matrix indices has been transformed to local, this query returns true, otherwise it returns false.
  bool IndicesAreLocal() const {return(Graph_.IndicesAreLocal());}

  //! If matrix indices are packed into single array (done in OptimizeStorage()) return true, otherwise false.
  bool IndicesAreContiguous() const {return(Graph_.IndicesAreContiguous());}

  //! If matrix is lower triangular in local index space, this query returns true, otherwise it returns false.
  bool LowerTriangular() const {return(Graph_.LowerTriangular());}

  //! If matrix is upper triangular in local index space, this query returns true, otherwise it returns false.
  bool UpperTriangular() const {return(Graph_.UpperTriangular());}

  //! If matrix has no diagonal entries in global index space, this query returns true, otherwise it returns false.
  bool NoDiagonal() const {return(Graph_.NoDiagonal());}

  //@}

  //! @name Attribute access functions
  //@{

  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f]
     \warning The NormInf() method will not properly calculate the infinity norm for a matrix that has entries that are
     replicated on multiple processors.  */
  double NormInf() const;

  //! Returns the one norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_1\f$ such that
     \f[\| A \|_1= \max_{1\lej\len} \sum_{i=1}^m |a_{ij}| \f].
     \warning The NormOne() method will not properly calculate the one norm for a matrix that has entries that are
     replicated on multiple processors.
  */
  double NormOne() const;

  //! Returns the frobenius norm of the global matrix.
  /* Returns the quantity \f[ \| A \|_{Frobenius} = \sqrt{\sum_{i=1}^m \sum_{j=1}^n\|a_{ij}\|^2}\f]
     \warning the NormFrobenius() method will not properly calculate the frobenius norm for a matrix that
     has entries which are replicated on multiple processors. In that case, the returned
     norm will be larger than the true norm.
   */
  double NormFrobenius() const;

  //! Returns the number of nonzero entries in the global matrix.
  /*
    Note that if maps are defined such that some nonzeros appear on
    multiple processors, then those nonzeros will be counted multiple times.
    If the user wishes to assemble a matrix from overlapping submatrices,
    they can use Epetra_FECrsMatrix.
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int  NumGlobalNonzeros() const {
      if(RowMap().GlobalIndicesInt())
        return (int) NumGlobalNonzeros64();
      throw "Epetra_CrsMatrix::NumGlobalNonzeros: GlobalIndices not int.";
    }
#endif
  long long NumGlobalNonzeros64() const {return(Graph_.NumGlobalNonzeros64());}

  //! Returns the number of global matrix rows.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int  NumGlobalRows() const {
      if(RowMap().GlobalIndicesInt())
        return (int) NumGlobalRows64();
      throw "Epetra_CrsMatrix::NumGlobalRows: GlobalIndices not int.";
    }
#endif
  long long NumGlobalRows64() const {return(Graph_.NumGlobalRows64());}

  //! Returns the number of global matrix columns.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int  NumGlobalCols() const {
      if(RowMap().GlobalIndicesInt())
        return (int) NumGlobalCols64();
      throw "Epetra_CrsMatrix::NumGlobalCols: GlobalIndices not int.";
    }
#endif
  long long NumGlobalCols64() const {return(Graph_.NumGlobalCols64());}

  //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int  NumGlobalDiagonals() const {
      if(RowMap().GlobalIndicesInt())
        return (int) NumGlobalDiagonals64();
      throw "Epetra_CrsMatrix::NumGlobalDiagonals: GlobalIndices not int.";
    }
#endif
  long long NumGlobalDiagonals64() const {return(Graph_.NumGlobalDiagonals64());}

  //! Returns the number of nonzero entries in the calling processor's portion of the matrix.
  int NumMyNonzeros() const {return(Graph_.NumMyNonzeros());}

  //! Returns the number of matrix rows owned by the calling processor.
  int NumMyRows() const {return(Graph_.NumMyRows());}

  //! Returns the number of entries in the set of column-indices that appear on this processor.
  /*! The set of column-indices that appear on this processor is the union of column-indices that
    appear in all local rows. The size of this set isn't available until FillComplete() has been called.
    \pre Filled()==true
  */
  int NumMyCols() const {return(Graph_.NumMyCols());}

  //! Returns the number of local nonzero diagonal entries, based on global row/column index comparisons.
  /*!
    \pre Filled()==true
  */
  int NumMyDiagonals() const {return(Graph_.NumMyDiagonals());}

  //! Returns the current number of nonzero entries in specified global row on this processor.
  int NumGlobalEntries(long long Row) const {return(Graph_.NumGlobalIndices(Row));}

  //! Returns the allocated number of nonzero entries in specified global row on this processor.
  int NumAllocatedGlobalEntries(int Row) const{return(Graph_.NumAllocatedGlobalIndices(Row));}

  //! Returns the maximum number of nonzero entries across all rows on this processor.
  /*!
    \pre Filled()==true
  */
  int MaxNumEntries() const {return(Graph_.MaxNumIndices());}

  //! Returns the maximum number of nonzero entries across all rows on all processors.
  /*!
    \pre Filled()==true
  */
  int GlobalMaxNumEntries() const {return(Graph_.GlobalMaxNumIndices());}

  //! Returns the current number of nonzero entries in specified local row on this processor.
  int NumMyEntries(int Row) const {return(Graph_.NumMyIndices(Row));}

  //! Returns the allocated number of nonzero entries in specified local row on this processor.
  int NumAllocatedMyEntries(int Row) const {return(Graph_.NumAllocatedMyIndices(Row));}

  //! Returns the index base for row and column indices for this graph.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  //! Index base for this map.
  int  IndexBase() const {
    if(RowMap().GlobalIndicesInt())
      return (int) IndexBase64();
    throw "Epetra_CrsMatrix::IndexBase: GlobalIndices not int.";
  }
#endif
  long long  IndexBase64() const {return(Graph_.IndexBase64());};


  //! Returns true if the graph associated with this matrix was pre-constructed and therefore not changeable.
  bool StaticGraph() {return(StaticGraph_);}

  //! Returns a reference to the Epetra_CrsGraph object associated with this matrix.
  const Epetra_CrsGraph& Graph() const {return(Graph_);}

  //! Returns the Epetra_Map object associated with the rows of this matrix.
  const Epetra_Map& RowMap() const {return((Epetra_Map &)Graph_.RowMap());}

  //! Replaces the current RowMap with the user-specified map object.
  /** Replaces the current RowMap with the user-specified map object, but only
      if currentmap->PointSameAs(newmap) is true. This is a collective function.
      Returns 0 if map is replaced, -1 if not.

      \pre RowMap().PointSameAs(newmap)==true
  */
  int ReplaceRowMap(const Epetra_BlockMap& newmap);

  //! Returns true if we have a well-defined ColMap, and returns false otherwise.
  /*! \pre We have a well-defined ColMap if a) a ColMap was passed in at construction,
    or b) the MakeColMap function has been called. (Calling either of the FillComplete functions
    will result in MakeColMap being called.)
  */
  bool HaveColMap() const {return(Graph_.HaveColMap());}

  //! Replaces the current ColMap with the user-specified map object.
  /** Replaces the current ColMap with the user-specified map object, but only
      if no entries have been inserted into the matrix (both IndicesAreLocal()
      and IndicesAreGlobal() are false) or currentmap->PointSameAs(newmap) is true.
      This is a collective function.
      Returns 0 if map is replaced, -1 if not.

      \pre (IndicesAreLocal()==false && IndicesAreGlobal()==false) || ColMap().PointSameAs(newmap)==true
  */
  int ReplaceColMap(const Epetra_BlockMap& newmap);

  //! Replaces the current DomainMap & Importer with the user-specified map object.
  /** Replaces the current DomainMap and Importer with the user-specified map object, but only
      if the matrix has been FillCompleted, Importer's TargetMap matches the ColMap
      and Importer's SourceMap matches the DomainMap (assuming the importer isn't null).  If an Importer
      is passed in, Epetra_CrsMatrix will copy it.

      Returns 0 if map/importer is replaced, -1 if not.

      \pre (!NewImporter && ColMap().PointSameAs(NewDomainMap)) || (NewImporter && ColMap().PointSameAs(NewImporter->TargetMap()) && NewDomainMap.PointSameAs(NewImporter->SourceMap()))

  */
  int ReplaceDomainMapAndImporter(const Epetra_Map & NewDomainMap, const Epetra_Import * NewImporter);

  //! Remove processes owning zero rows from the Maps and their communicator.
  /** Remove processes owning zero rows from the Maps and their communicator.
     \warning This method is ONLY for use by experts.

     \warning We make NO promises of backwards compatibility.
     This method may change or disappear at any time.

     \param newMap [in] This <i>must</i> be the result of calling
     the removeEmptyProcesses() method on the row Map.  If it
     is not, this method's behavior is undefined.  This pointer
     will be null on excluded processes.
  */
  int RemoveEmptyProcessesInPlace(const Epetra_BlockMap * NewMap);

  //! Returns the Epetra_Map object that describes the set of column-indices that appear in each processor's locally owned matrix rows.
  /*!Note that if the matrix was constructed with only a row-map, then until FillComplete() is called, this method returns
    a column-map that is a copy of the row-map. That 'initial' column-map is replaced with a computed column-map (that
    contains the set of column-indices appearing in each processor's local portion of the matrix) when FillComplete() is
    called.

    \pre HaveColMap()==true
  */
  const Epetra_Map& ColMap() const {return((Epetra_Map &) Graph_.ColMap());}

  //! Returns the Epetra_Map object associated with the domain of this matrix operator.
  /*!
    \pre Filled()==true
  */
  const Epetra_Map& DomainMap() const {return((Epetra_Map &)Graph_.DomainMap());}

  //! Returns the Epetra_Map object associated with the range of this matrix operator.
  /*!
    \pre Filled()==true
  */
  const Epetra_Map& RangeMap() const  {return((Epetra_Map &)Graph_.RangeMap());}

  //! Returns the Epetra_Import object that contains the import operations for distributed operations.
  const Epetra_Import* Importer() const {return(Graph_.Importer());}

  //! Returns the Epetra_Export object that contains the export operations for distributed operations.
  const Epetra_Export* Exporter() const {return(Graph_.Exporter());}

  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm& Comm() const {return(Epetra_DistObject::Comm());}
  //@}

  //! @name Local/Global ID methods
  //@{
  //! Returns the local row index for given global row index, returns -1 if no local row for this global row.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int LRID( int GRID_in) const {return(Graph_.LRID(GRID_in));}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  int LRID( long long GRID_in) const {return(Graph_.LRID(GRID_in));}
#endif

#if defined(EPETRA_NO_32BIT_GLOBAL_INDICES) && defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  // default implementation so that no compiler/linker error in case neither 32 nor 64
  // bit indices present.
  int LRID(long long GRID_in) const {return(Graph_.LRID(GRID_in));}
#endif

  //! Returns the global row index for give local row index, returns IndexBase-1 if we don't have this local row.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    int  GRID(int LRID_in) const {
      if(RowMap().GlobalIndicesInt())
        return (int) GRID64(LRID_in);
      throw "Epetra_CrsMatrix::GRID: GlobalIndices not int.";
    }
#endif
  long long GRID64( int LRID_in) const {return(Graph_.GRID64(LRID_in));}

  //! Returns the local column index for given global column index, returns -1 if no local column for this global column.
  /*!
    \pre HaveColMap()==true (If HaveColMap()==false, returns -1)
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int LCID( int GCID_in) const {return(Graph_.LCID(GCID_in));}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  int LCID( long long GCID_in) const {return(Graph_.LCID(GCID_in));}
#endif

  //! Returns the global column index for give local column index, returns IndexBase-1 if we don't have this local column.
  /*!
    \pre HaveColMap()==true (If HaveColMap()==false, returns -1)
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    int  GCID(int LCID_in) const {
      if(RowMap().GlobalIndicesInt())
        return (int) GCID64(LCID_in);
      throw "Epetra_CrsMatrix::GCID: GlobalIndices not int.";
    }
#endif
  long long GCID64( int LCID_in) const {return(Graph_.GCID64(LCID_in));}

  //! Returns true if the GRID passed in belongs to the calling processor in this map, otherwise returns false.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  bool MyGRID(int GRID_in) const {return(Graph_.MyGRID(GRID_in));}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  bool MyGRID(long long GRID_in) const {return(Graph_.MyGRID(GRID_in));}
#endif

  //! Returns true if the LRID passed in belongs to the calling processor in this map, otherwise returns false.
  bool MyLRID(int LRID_in) const {return(Graph_.MyLRID(LRID_in));}

  //! Returns true if the GCID passed in belongs to the calling processor in this map, otherwise returns false.
  /*!
    \pre HaveColMap()==true (If HaveColMap()==false, returns -1)
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  bool MyGCID(int GCID_in) const {return(Graph_.MyGCID(GCID_in));}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  bool MyGCID(long long GCID_in) const {return(Graph_.MyGCID(GCID_in));}
#endif

  //! Returns true if the LRID passed in belongs to the calling processor in this map, otherwise returns false.
  /*!
    \pre HaveColMap()==true (If HaveColMap()==false, returns -1)
   */
  bool MyLCID(int LCID_in) const {return(Graph_.MyLCID(LCID_in));}

  //! Returns true of GID is owned by the calling processor, otherwise it returns false.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  bool MyGlobalRow(int GID) const {return(Graph_.MyGlobalRow(GID));}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  bool MyGlobalRow(long long GID) const {return(Graph_.MyGlobalRow(GID));}
#endif
  //@}


  //! @name I/O Methods
  //@{

  //! Print method
  virtual void Print(std::ostream& os) const;
  //@}

  //! @name Additional methods required to support the Epetra_Operator interface
  //@{

  //! Returns a character string describing the operator
  const char* Label() const {return(Epetra_Object::Label());}

  //! If set true, transpose of this operator will be applied.
  /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
    affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface
    does not support transpose use, this method should return a value of -1.

    \param UseTranspose - (In) If true, multiply by the transpose of operator, otherwise just use operator.

    \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose_in) {UseTranspose_ = UseTranspose_in; return(0);}

  //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
  /*!
    \param X - (In) An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Y -(Out) An Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
    \pre Filled()==true
    \post Unchanged.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(Epetra_CrsMatrix::Multiply(Epetra_CrsMatrix::UseTranspose(), X, Y));}

  //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
  /*! In this implementation, we use several existing attributes to determine how virtual
    method ApplyInverse() should call the concrete method Solve().  We pass in the UpperTriangular(),
    the Epetra_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
    if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
    be accounted for when doing a triangular solve.

    \param X - (In) An Epetra_MultiVector of dimension NumVectors to solve for.
    \param Y - (Out) An Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
    \pre Filled()==true
    \post Unchanged.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(Solve(UpperTriangular(), Epetra_CrsMatrix::UseTranspose(), NoDiagonal(), X, Y));}

  //! Returns true because this class can compute an Inf-norm.
  bool HasNormInf() const {return(true);}

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);}

  //! Returns the Epetra_Map object associated with the domain of this matrix operator.
  const Epetra_Map& OperatorDomainMap() const
        {
          if (UseTranspose()) return(RangeMap());
          else return(DomainMap());
        }

  //! Returns the Epetra_Map object associated with the range of this matrix operator.
  const Epetra_Map& OperatorRangeMap() const
        {
          if (UseTranspose()) return(DomainMap());
          else return(RangeMap());
        }

  //@}
  //! @name Additional methods required to implement Epetra_RowMatrix interface
  //@{

  //! Return the current number of values stored for the specified local row.
  /*! Similar to NumMyEntries() except NumEntries is returned as an argument
    and error checking is done on the input value MyRow.
    \param MyRow - (In) Local row.
    \param NumEntries - (Out) Number of nonzero values.

    \return Integer error code, set to 0 if successful.
    \pre None.
    \post Unchanged.
  */
  int NumMyRowEntries(int MyRow, int& NumEntries) const;

  //! Map() method inherited from Epetra_DistObject
  const Epetra_BlockMap& Map() const { return Epetra_DistObject::Map(); }

  //! Returns the Epetra_Map object associated with the rows of this matrix.
  const Epetra_Map& RowMatrixRowMap() const {return(RowMap());}

  //! Returns the Epetra_Map object associated with columns of this matrix.
  const Epetra_Map& RowMatrixColMap() const {return(ColMap());}

  //! Returns the Epetra_Import object that contains the import operations for distributed operations.
  const Epetra_Import* RowMatrixImporter() const {return(Importer());}

  //@}

  //! @name Inlined Operator Methods
  //@{

  //! Inlined bracket operator for fast access to data. (Const and Non-const versions)
  /*! No error checking and dangerous for optimization purposes.
    \param Loc - (In) Local row.

    \return reference to pointer to locally indexed Loc row in matrix.
  */
  inline double* operator[] (int Loc) {
    if (StorageOptimized()){ int * ind = Graph().IndexOffset(); return(All_Values_+ind[Loc]);}
    else return Values_[Loc];}
  inline double* operator[] (int Loc) const {
    if (StorageOptimized()){ int * ind = Graph().IndexOffset(); return(All_Values_+ind[Loc]);}
    else return Values_[Loc];}
  //@}

  //! @name Expert-only methods:  These methods are intended for experts only and have some risk of changing in the future, since they rely on underlying data structure assumptions
  //@{
  //! Returns internal data pointers associated with Crs matrix format.
  /*! Returns data pointers to facilitate optimized code within external packages.

    \param IndexOffset - (Out) Extracted array of indices into Values[] and Indices[]. Local
                          row k is stored in Values[IndexOffset[k]:IndexOffset[k+1]-1] and
                          Indices[IndexOffset[k]:IndexOffset[k+1]-1].
    \param Values - (Out) Extracted values for all local rows.
    \param Indices - (Out) Extracted local column indices for the corresponding values.

    \return Integer error code, set to 0 if successful. Returns -1 if FillComplete has not been
                performed or Storage has not been Optimized.

    \warning This method is intended for expert only, its use may require user code modifications in future versions of Epetra.
  */
  int ExtractCrsDataPointers(int *& IndexOffset, int *& Indices, double *& Values_in) const {
    if (StorageOptimized()) {
      IndexOffset = Graph().IndexOffset();
      Indices = Graph().All_Indices();
      Values_in  = All_Values();
      return (0);
    }
    else { IndexOffset = 0; Indices = 0; Values_in  = 0; return (-1);} }


    //! Returns a reference to the Epetra_IntSerialDenseVector used to hold the local IndexOffsets (CRS rowptr)
    /*!
       \warning This method is intended for experts only, its use may require user code modifications in future versions of Epetra.
    */
    Epetra_IntSerialDenseVector& ExpertExtractIndexOffset();

    //! Returns a reference to the Epetra_IntSerialDenseVector used to hold the local All_Indices (CRS colind)
    /*!
       \warning This method is intended for experts only, its use may require user code modifications in future versions of Epetra.
    */
    Epetra_IntSerialDenseVector& ExpertExtractIndices();

    //! Returns a reference to the double* used to hold the values array
    /*!
       \warning This method is intended for experts only, its use may require user code modifications in future versions of Epetra.
    */
    double *& ExpertExtractValues() {return All_Values_;}


    //! Performs a FillComplete on an object that aready has filled CRS data
    /*! Performs a lightweight FillComplete on an object that already has filled IndexOffsets, All_Indices and All_Values.
      This routine is needed to support the EpetraExt::MatrixMatrix::Multiply and should not be called by users.
       \warning Epetra_CrsMatrix will assume ownership of the Importer/Exporter you pass in.  You should not deallocate it afterwards.
       \warning This method is intended for expert developer use only, and should never be called by user code.
    */
    int ExpertStaticFillComplete(const Epetra_Map & DomainMap,const Epetra_Map & RangeMap, const Epetra_Import * Importer=0, const Epetra_Export * Exporter=0, int NumMyDiagonals=-1);


    //! Makes sure this matrix has a unique CrsGraphData object
    /*! This routine is needed to support the EpetraExt::MatrixMatrix::Multiply and should not be called by users.
      \warning This method is intended for expert developer use only, and should never be called by user code.
    */
    int ExpertMakeUniqueCrsGraphData();

  //! Forces FillComplete() to locally order ghostnodes associated with each remote processor in ascending order.
  /*! To be compliant with AztecOO, FillComplete() already locally orders ghostnodes such that
      information received from processor k has a lower local numbering than information received
      from processor j if k is less than j.  SortGhostsAssociatedWithEachProcessor(True) further
      forces FillComplete() to locally number all ghostnodes received from processor k in ascending
      order. That is, the local numbering of b is less than c if the global numbering of b is less
      than c and if both b and c are owned by the same processor. This is done to be compliant with
      some limited block features within ML. In particular, some ML features require that a block
      structure of the matrix be maintained even within the ghost variables. Always returns 0.
  */

  int SortGhostsAssociatedWithEachProcessor(bool Flag)  {Graph_.SortGhostsAssociatedWithEachProcessor(Flag); return(0);}
  //@}

  //! @name Deprecated methods:  These methods still work, but will be removed in a future version
  //@{

  //! Use ColMap() instead.
  const Epetra_Map& ImportMap() const {return((Epetra_Map&) Graph_.ImportMap());}

  //! Use FillComplete() instead.
  int TransformToLocal();

  //! Use FillComplete(const Epetra_Map& DomainMap, const Epetra_Map& RangeMap) instead.
  int TransformToLocal(const Epetra_Map* DomainMap, const Epetra_Map* RangeMap);

  //@}


 protected:
  bool Allocated() const {return(Allocated_);}
  int SetAllocated(bool Flag) {Allocated_ = Flag; return(0);}
  double** Values() const {
    if (StorageOptimized()) throw ReportError("This method: double** Values() cannot be called when StorageOptimized()==true", -1);
    else return(Values_);}
  double* All_Values() const {
    if (!StorageOptimized()) throw ReportError("This method: double* All_Values()cannot be called when StorageOptimized()==false", -1);
    else return(All_Values_);}
  double* Values(int LocalRow) const {
    if (StorageOptimized())
      if (Graph().StorageOptimized())
  return(All_Values_+Graph().IndexOffset()[LocalRow]);
      else throw ReportError("This method: double* Values()cannot be called when StorageOptimized()==true and Graph().StorageOptimized()==false", -1);
    else return(Values_[LocalRow]);}

  void InitializeDefaults();
  int Allocate();

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int InsertValues(int LocalRow, int NumEntries, double* Values, int* Indices);
  int InsertValues(int LocalRow, int NumEntries, const double* Values, const int* Indices);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  int InsertValues(int LocalRow, int NumEntries, double* Values, long long* Indices);
  int InsertValues(int LocalRow, int NumEntries, const double* Values, const long long* Indices);
#endif

  int InsertOffsetValues(long long GlobalRow, int NumEntries, double *Values, int *Indices);
  int InsertOffsetValues(long long GlobalRow, int NumEntries, const double *Values, const int *Indices);
  int ReplaceOffsetValues(long long GlobalRow, int NumEntries, const double *Values, const int *Indices);
  int SumIntoOffsetValues(long long GlobalRow, int NumEntries, const double *Values, const int *Indices);
  void UpdateImportVector(int NumVectors) const;
  void UpdateExportVector(int NumVectors) const;
  void GeneralMV(double * x, double * y) const;
  void GeneralMTV(double * x, double * y) const;
  void GeneralMM(double ** X, int LDX, double ** Y, int LDY, int NumVectors) const;
  void GeneralMTM(double ** X, int LDX, double ** Y, int LDY, int NumVectors) const;
  void GeneralSV(bool Upper, bool Trans, bool UnitDiagonal, double * x, double * y) const;
  void GeneralSM(bool Upper, bool Trans, bool UnitDiagonal, double ** X, int LDX, double ** Y, int LDY, int NumVectors) const;

  void SetStaticGraph(bool Flag) {StaticGraph_ = Flag;}

  int CheckSizes(const Epetra_SrcDistObject& A);

  int CopyAndPermute(const Epetra_SrcDistObject& Source,
                     int NumSameIDs,
                     int NumPermuteIDs,
                     int* PermuteToLIDs,
                     int* PermuteFromLIDs,
                     const Epetra_OffsetIndex * Indexor,
                     Epetra_CombineMode CombineMode = Zero);
  int CopyAndPermuteCrsMatrix(const Epetra_CrsMatrix& A,
                              int NumSameIDs,
                              int NumPermuteIDs,
                              int* PermuteToLIDs,
                              int* PermuteFromLIDs,
                              const Epetra_OffsetIndex * Indexor,
                              Epetra_CombineMode CombineMode);
  int CopyAndPermuteRowMatrix(const Epetra_RowMatrix& A,
                              int NumSameIDs,
                              int NumPermuteIDs,
                              int* PermuteToLIDs,
                              int* PermuteFromLIDs,
                              const Epetra_OffsetIndex * Indexor,
                              Epetra_CombineMode CombineMode);

  int PackAndPrepare(const Epetra_SrcDistObject& Source,
                     int NumExportIDs,
                     int* ExportLIDs,
                     int& LenExports,
                     char*& Exports,
                     int& SizeOfPacket,
                     int* Sizes,
                     bool& VarSizes,
                     Epetra_Distributor& Distor);

  int UnpackAndCombine(const Epetra_SrcDistObject& Source,
                       int NumImportIDs,
                       int* ImportLIDs,
                       int LenImports,
                       char* Imports,
                       int& SizeOfPacket,
                       Epetra_Distributor& Distor,
                       Epetra_CombineMode CombineMode,
                       const Epetra_OffsetIndex * Indexor);

  //! Sort column entries, row-by-row, in ascending order.
  int SortEntries();

  //! If SortEntries() has been called, this query returns true, otherwise it returns false.
  bool Sorted() const {return(Graph_.Sorted());}

  //! Add entries that have the same column index. Remove redundant entries from list.
  int MergeRedundantEntries();

  //! If MergeRedundantEntries() has been called, this query returns true, otherwise it returns false.
  bool NoRedundancies() const {return(Graph_.NoRedundancies());}

  void DeleteMemory();

  Epetra_CrsGraph Graph_;
  bool Allocated_;
  bool StaticGraph_;
  bool UseTranspose_;
  bool constructedWithFilledGraph_;
  bool matrixFillCompleteCalled_;
  bool StorageOptimized_;

  double** Values_;
  int* Values_alloc_lengths_;
  double* All_Values_;
  mutable double NormInf_;
  mutable double NormOne_;
  mutable double NormFrob_;

  int NumMyRows_;
  mutable Epetra_MultiVector* ImportVector_;
  mutable Epetra_MultiVector* ExportVector_;

  Epetra_DataAccess CV_;

  bool squareFillCompleteCalled_;
#ifdef Epetra_ENABLE_CASK
  caskHandle_t cask;
#endif
 private:

  // These are the pre-5.0 versions of solve.  They are still faster that generic 5.0 solves, so we keep them around
  int Solve1(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_Vector& x, Epetra_Vector& y) const;
  int Solve1(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

private:

  template<typename int_type>
  int TInsertGlobalValues(int_type Row, int NumEntries, const double* values, const int_type* Indices);

  template<typename int_type>
  int TInsertGlobalValues(int_type Row, int NumEntries, double* values, int_type* Indices);

  template<typename int_type>
  int InsertValues(int Row, int NumEntries, const double* values, const int_type* Indices);

  template<typename int_type>
  int InsertValues(int Row, int NumEntries, double* values, int_type* Indices);

  template<typename int_type>
  int TReplaceGlobalValues(int_type Row, int NumEntries, const double * srcValues, const int_type *Indices);

  template<typename int_type>
  int TSumIntoGlobalValues(int_type Row, int NumEntries, const double * srcValues, const int_type *Indices);

  template<typename int_type>
  int ExtractGlobalRowCopy(int_type Row, int Length, int & NumEntries, double * values, int_type * Indices) const;

  template<typename int_type>
  int ExtractGlobalRowCopy(int_type Row, int Length, int & NumEntries, double * values) const;

  template<typename int_type>
  int ExtractGlobalRowView(int_type Row, int & NumEntries, double *& values, int_type *& Indices) const;

  template<typename int_type>
  int ExtractGlobalRowView(int_type Row, int & NumEntries, double *& values) const;

  template<typename int_type>
    int TCopyAndPermuteCrsMatrix(const Epetra_CrsMatrix& A,
                              int NumSameIDs,
                              int NumPermuteIDs,
                              int* PermuteToLIDs,
                              int* PermuteFromLIDs,
                              const Epetra_OffsetIndex * Indexor,
                              Epetra_CombineMode CombineMode);

  template<typename int_type>
    int TCopyAndPermuteRowMatrix(const Epetra_RowMatrix& A,
                              int NumSameIDs,
                              int NumPermuteIDs,
                              int* PermuteToLIDs,
                              int* PermuteFromLIDs,
                              const Epetra_OffsetIndex * Indexor,
                              Epetra_CombineMode CombineMode);

  template<typename int_type>
    int TUnpackAndCombine(const Epetra_SrcDistObject& Source,
                       int NumImportIDs,
                       int* ImportLIDs,
                       int LenImports,
                       char* Imports,
                       int& SizeOfPacket,
                       Epetra_Distributor& Distor,
                       Epetra_CombineMode CombineMode,
                       const Epetra_OffsetIndex * Indexor);

  // Used for fused[import|export] constructors
  template<class TransferType>
  void FusedTransfer(const Epetra_CrsMatrix & SourceMatrix,
          const TransferType & RowTransfer,
          const TransferType* DomainTransfer,
          const Epetra_Map * DomainMap,
          const Epetra_Map * RangeMap,
          bool RestrictCommunicator);

 public:

  void FusedImport(const Epetra_CrsMatrix & SourceMatrix,
        const Epetra_Import & RowImporter,
        const Epetra_Map * DomainMap,
        const Epetra_Map * RangeMap,
        bool RestrictCommunicator);

  void FusedExport(const Epetra_CrsMatrix & SourceMatrix,
        const Epetra_Export & RowExporter,
        const Epetra_Map * DomainMap,
        const Epetra_Map * RangeMap,
        bool RestrictCommunicator);

  void FusedImport(const Epetra_CrsMatrix & SourceMatrix,
        const Epetra_Import & RowImporter,
        const Epetra_Import * DomainImporter,
        const Epetra_Map * DomainMap,
        const Epetra_Map * RangeMap,
        bool RestrictCommunicator);

  void FusedExport(const Epetra_CrsMatrix & SourceMatrix,
        const Epetra_Export & RowExporter,
        const Epetra_Export * DomainExporter,
        const Epetra_Map * DomainMap,
        const Epetra_Map * RangeMap,
        bool RestrictCommunicator);


};
#endif /* EPETRA_CRSMATRIX_H */
