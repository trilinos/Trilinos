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

#ifndef EPETRA_FEVECTOR_H
#define EPETRA_FEVECTOR_H

#include <Epetra_CombineMode.h>
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>

#include <vector>

class Epetra_IntSerialDenseVector;
class Epetra_LongLongSerialDenseVector;
class Epetra_SerialDenseVector;

  // TODO this file needs to be changed for long long


/** Epetra Finite-Element Vector. This class inherits Epetra_MultiVector
  and thus provides all Epetra_MultiVector functionality.

  The added functionality provided by Epetra_FEVector is the ability to
  perform finite-element style vector assembly. It accepts sub-vector
  contributions, such as those that would come from element-load vectors, etc.
  These sub-vectors need not be owned by the local processor. In other
  words, the user can assemble overlapping data (e.g., corresponding to shared
  finite-element nodes). When the user is finished assembling their vector
  data, they then call the method Epetra_FEVector::GlobalAssemble() which
  gathers the overlapping data (all non-local data that was input on each
  processor) into the data-distribution specified by the map with which the
  Epetra_FEVector was constructed.
*/

class EPETRA_LIB_DLL_EXPORT Epetra_FEVector : public Epetra_MultiVector {
 public:
   /** Constructor that requires a map specifying a non-overlapping
      data layout.
      @param Map Map describing a non-overlapping distribution for the
        underlying Epetra_MultiVector into which this Epetra_FEVector will
        funnel data.

      @param numVectors Optional argument, default value is 1. (See the
        documentation for Epetra_MultiVector for the meaning of this
        argument.

      @param ignoreNonLocalEntries Optional argument, default value is
        false.  Under certain special circumstances it is desirable to
        have non-local contributions ignored rather than saving them
        for the GlobalAssemble step.
   */
   Epetra_FEVector(const Epetra_BlockMap& Map,
                   int numVectors=1,
		   bool ignoreNonLocalEntries=false);

  //! Set multi-vector values from two-dimensional array.
  /*!
    \param In 
    Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In 
    Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param In
    A - Pointer to an array of double precision numbers.  The first vector starts at A.
    The second vector starts at A+MyLDA, the third at A+2*MyLDA, and so on.
    \param In
    MyLDA - The "Leading Dimension", or stride between vectors in memory.
    \warning This value refers to the stride on the calling processor.  Thus it is a
    local quantity, not a global quantity.
    \param In 
    NumVectors - Number of vectors in multi-vector.

    \return Integer error code, set to 0 if successful.

    See Detailed Description section for further discussion.
  */
  Epetra_FEVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, 
		     double *A, int MyLDA, int NumVectors,
        bool ignoreNonLocalEntries=false);
  
  //! Set multi-vector values from array of pointers.
  /*!
    \param In 
    Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In 
    Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param In
    ArrayOfPointers - An array of pointers such that ArrayOfPointers[i] points to the memory
    location containing ith vector to be copied.
    \param In 
    NumVectors - Number of vectors in multi-vector.

    \return Integer error code, set to 0 if successful.

    See Detailed Description section for further discussion.
  */
  Epetra_FEVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, 
		     double **ArrayOfPointers, int NumVectors,
            bool ignoreNonLocalEntries=false);

  /** Copy constructor. */
  Epetra_FEVector(const Epetra_FEVector& source);

   /** Destructor */
   virtual ~Epetra_FEVector();

   /** Accumulate values into the vector, adding them to any values that
       already exist for the specified indices.
   */
   int SumIntoGlobalValues(int numIDs,
                           const int* GIDs, const double* values,
                           int vectorIndex=0);
   int SumIntoGlobalValues(int numIDs,
                           const long long* GIDs, const double* values,
                           int vectorIndex=0);

   /** Accumulate values into the vector, adding them to any values that
       already exist for the specified GIDs.

       @param GIDs List of global ids. Must be the same length as the
       accompanying list of values.

       @param values List of coefficient values. Must be the same length as
       the accompanying list of GIDs.
   */
   int SumIntoGlobalValues(const Epetra_IntSerialDenseVector& GIDs,
			   const Epetra_SerialDenseVector& values,
                           int vectorIndex=0);
   int SumIntoGlobalValues(const Epetra_LongLongSerialDenseVector& GIDs,
			   const Epetra_SerialDenseVector& values,
                           int vectorIndex=0);

   /** Copy values into the vector overwriting any values that already exist
        for the specified indices.
    */
   int ReplaceGlobalValues(int numIDs, const int* GIDs, const double* values,
                           int vectorIndex=0);
   int ReplaceGlobalValues(int numIDs, const long long* GIDs, const double* values,
                           int vectorIndex=0);

   /** Copy values into the vector, replacing any values that
       already exist for the specified GIDs.

       @param GIDs List of global ids. Must be the same length as the
       accompanying list of values.

       @param values List of coefficient values. Must be the same length as
       the accompanying list of GIDs.
   */
   int ReplaceGlobalValues(const Epetra_IntSerialDenseVector& GIDs,
			   const Epetra_SerialDenseVector& values,
                           int vectorIndex=0);
   int ReplaceGlobalValues(const Epetra_LongLongSerialDenseVector& GIDs,
			   const Epetra_SerialDenseVector& values,
                           int vectorIndex=0);

   int SumIntoGlobalValues(int numIDs,
                           const int* GIDs,
			   const int* numValuesPerID,
			   const double* values,
                           int vectorIndex=0);
   int SumIntoGlobalValues(int numIDs,
                           const long long* GIDs,
			   const int* numValuesPerID,
			   const double* values,
                           int vectorIndex=0);

   int ReplaceGlobalValues(int numIDs, const int* GIDs,
			   const int* numValuesPerID,
			   const double* values,
                           int vectorIndex=0);
   int ReplaceGlobalValues(int numIDs, const long long* GIDs,
			   const int* numValuesPerID,
			   const double* values,
                           int vectorIndex=0);

   /** Gather any overlapping/shared data into the non-overlapping partitioning
      defined by the Map that was passed to this vector at construction time.
      Data imported from other processors is stored on the owning processor
      with a "sumInto" or accumulate operation.
      This is a collective method -- every processor must enter it before any
      will complete it.

      Optimization for power-users:
      The optional parameter 'reuse_map_and_exporter' defaults to false.
      By default, a map that describes the non-local data is re-created at
      each call to GlobalAssemble, along with an exporter used to do the
      communication.
      This is expensive. If you know that the layout of your nonlocal data has
      not changed since your previous call to GlobalAssemble, you can set this
      flag to true and it will reuse the previously created map and exporter
      rather than creating new ones.
   */
   int GlobalAssemble(Epetra_CombineMode mode = Add,
                      bool reuse_map_and_exporter = false);

   /** Set whether or not non-local data values should be ignored.
    */
   void setIgnoreNonLocalEntries(bool flag) {
     ignoreNonLocalEntries_ = flag;
   }

   Epetra_FEVector& operator=(const Epetra_FEVector& source);

 protected:
  int inputValues(int numIDs,
                  const int* GIDs, const double* values,
                  bool suminto,
                  int vectorIndex);

  int inputValues(int numIDs,
                  const int* GIDs, const int* numValuesPerID,
		  const double* values,
                  bool suminto,
                  int vectorIndex);

  int inputNonlocalValue(int GID, double value, bool suminto,
                         int vectorIndex);

  int inputNonlocalValues(int GID, int numValues, const double* values,
			  bool suminto, int vectorIndex);

  //! Allocate the Map, Export object, and MultiVector for nonlocal data.
  void createNonlocalMapAndExporter();

  //! Deallocate the Map, Export object, and MultiVector for nonlocal data.
  void destroyNonlocalMapAndExporter();

  //! Make all the nonlocal multivector entries zero.
  void zeroNonlocalData();

  /// \brief Deallocate storage for nonlocal data.
  ///
  /// This does not include whatever \c destroyNonlocalMapAndExporter() 
  /// deallocates.
  void destroyNonlocalData();

  long long myFirstID_;
  int myNumIDs_;

  std::vector<int> nonlocalIDs_;
  std::vector<int> nonlocalElementSize_;

  /// \brief Array of arrays (one per column) of nonlocal coefficients.
  ///
  /// Entry j in this array is an array of nonlocal coefficients for
  /// column j of the multivector.  \c nonlocalCoefs_[j] in turn
  /// contains the data for all the nonlocal block entries.  Entries
  /// within a block are stored in column-major order (Fortran style),
  /// in contiguous chunks of \c Map().MaxElementSize() entries.
  ///
  /// In summary: for k = 0 to numNonlocalIDs_ - 1,
  /// nonlocalCoefs_[j][j*M + i] contains the entries of nonlocal
  /// block entry k, for i = 0 to nonlocalElementSize_[k] - 1.
  std::vector<std::vector<double> > nonlocalCoefs_;

  //! Map describing distribution of nonlocal data.
  Epetra_BlockMap* nonlocalMap_;

  //! Export object that sums nonlocal data into a nonoverlapping distribution.
  Epetra_Export* exporter_;

  //! Multivector that holds nonlocal data; source for the Export operation.
  Epetra_MultiVector* nonlocalVector_;

  bool ignoreNonLocalEntries_;
};

#endif

