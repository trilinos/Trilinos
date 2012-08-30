/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#ifndef _fei_Vector_hpp_
#define _fei_Vector_hpp_

#include <fei_iosfwd.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_Reducer.hpp>

namespace fei {
  /** Abstract representation of an algebraic multi-vector. This representation may
      be used with an overlapping data decomposition. The data distribution is
      defined by the fei::VectorSpace object returned by the method
      getVectorSpace().
      This representation does not require that data be accessed only on the
      'owning' processor. Data for any indices that are either shared or owned by
      the local processor may be passed to, or accessed from, the vector on the
      local processor. In most cases the underlying library-specific vector will
      have a non-overlapping data decomposition (each equation uniquely owned by
      a single processor). Overlapping data (shared by local processor but the
      equation is owned by another processor) may be assembled into this abstract
      vector locally, and will be moved into the underlying non-overlapping
      vector on the correct processor when the gatherFromOverlap() method is
      called. Conversely, if the user wants to retrieve overlapping data from
      the vector locally for an equation that resides on another processor, that
      data is not guaranteed to be available until the scatterToOverlap() method
      is called. The scatterToOverlap() method does communication necessary to
      populate shared-but-not-owned data in the fei::Vector from data in the
      underlying algebraic vector.

      From the point of view of fei::Vector, there are two types of data: owned
      and shared-but-not-owned.

      Data Input (passing user data into the vector):<br>
      When locally-owned data is input, fei::Vector relays it immediately to the
      underlying algebraic vector.
      When shared-but-not-owned data is input, fei::Vector holds it in temporary
      storage. When gatherToOverlap() is called, fei::Vector moves it to the
      owning processor and then relays it to the underlying algebraic vector. At
      that point the temporary storage is deleted.

      Data Access (retrieving data from the vector):<br>
      When locally-owned data is accessed, fei::Vector retrieves it from the
      underlying algebraic vector directly.
      In order to access shared-but-not-owned data (overlapped data), it is 
      necessary first to call the method scatterToOverlap(). This method does the
      communication necessary to re-create and populate temporary storage with
      the shared data by retrieving that data from the underlying algebraic
      vector on the owning processor and sending it to the sharing processors.
  */
  class Vector {
  public:
    /** Vector Factory interface */
    class Factory {
    public:
      /** Usual virtual destructor */
      virtual ~Factory(){}

      /** Produce an instance of a Vector using a VectorSpace. */
      virtual fei::SharedPtr<fei::Vector>
	createVector(fei::SharedPtr<fei::VectorSpace> vecSpace,
		     int numVectors=1) = 0;

      /** Produce an instance of a Vector using a VectorSpace. */
      virtual fei::SharedPtr<fei::Vector>
	createVector(fei::SharedPtr<fei::VectorSpace> vecSpace,
		     bool isSolutionVector,
		     int numVectors=1) = 0;

      /** Produce an instance of a Vector using a MatrixGraph. */
      virtual fei::SharedPtr<fei::Vector>
	createVector(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
		     int numVectors=1) = 0;

      /** Produce an instance of a Vector using a MatrixGraph. */
      virtual fei::SharedPtr<fei::Vector>
	createVector(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
		     bool isSolutionVector,
		     int numVectors=1) = 0;
    };

    /** Virtual destructor. */
    virtual ~Vector(){}

    /** Return an implementation-dependent name describing the run-time type
	of this object.
    */
    virtual const char* typeName() const = 0;

    /** Set a specified scalar throughout the vector. */
    virtual int putScalar(double scalar) = 0;

    /** Accumulate values into the vector, adding them to any values that
	already exist for the specified indices.
    */
    virtual int sumIn(int numValues, const int* indices, const double* values,
		      int vectorIndex=0) = 0;

    /** Copy values into the vector overwriting any values that already exist
	for the specified indices.
    */
    virtual int copyIn(int numValues, const int* indices, const double* values,
		       int vectorIndex=0) = 0;

    /** Retrieve a copy of values from the vector for the specified indices.
	Note that if the specified indices are not local in the underlying
	non-overlapping data decomposition, these values are not guaranteed to
	be correct until after the scatterToOverlap() method has been called.
    */
    virtual int copyOut(int numValues, const int* indices, double* values,
			int vectorIndex=0) const = 0;

    /** Update 'this' = b*'this' + a*x
     */
    virtual int update(double a,
		       const fei::Vector* x,
		       double b) = 0;

    /** Scatter data from the underlying non-overlapping data decomposition to
	the overlapping data decomposition. In other words, update values for
	shared indices from underlying uniquely owned data.
    */
    virtual int scatterToOverlap() = 0;

    /** perform initial communication to establish message sizes that will
      be needed for exchanging shared-node data.
      Called from within gatherFromOverlap usually, doesn't usually need to
      be explicitly called by client code. (Power users only...)
    */
    virtual void setCommSizes() = 0;

    /** Gather data from the overlapping data decomposition into the underlying
	non-overlapping data decomposition.
    */
    virtual int gatherFromOverlap(bool accumulate = true) = 0;

    /** Query for the VectorSpace object associated with this vector. */
    virtual fei::SharedPtr<fei::VectorSpace> getVectorSpace() const = 0;

    /** Set the VectorSpace object associated with this vector. */
    virtual void setVectorSpace(fei::SharedPtr<fei::VectorSpace> vecSpace) = 0;

    /** Sum field data into the vector, adding it to any data that may already
        be present at the specified locations.
        If the specified fieldID doesn't exist at one or more of the specified
        IDs, then the corresponding positions in the data array will simply
        not be used.
    */
    virtual int sumInFieldData(int fieldID,
			       int idType,
			       int numIDs,
			       const int* IDs,
			       const double* data,
			       int vectorIndex=0) = 0;

    /** Copy field data into the vector, overwriting any data that may already
        be present at the specified locations. 
        If the specified fieldID doesn't exist at one or more of the specified
        IDs, then the corresponding positions in the data array will simply
        not be used.
    */
    virtual int copyInFieldData(int fieldID,
				int idType,
				int numIDs,
				const int* IDs,
				const double* data,
				int vectorIndex=0) = 0;

    virtual int copyInFieldDataLocalIDs(int fieldID,
        int idType,
        int numIDs,
        const int* localIDs,
        const double* data,
        int vectorIndex=0) = 0;

    /** Copy field data out of the vector into the user-allocated data array.
      If the specified fieldID doesn't exist at one or more of the specified
      IDs, then the corresponding positions in the data array will simply not
      be referenced.
    */
    virtual int copyOutFieldData(int fieldID,
				 int idType,
				 int numIDs,
				 const int* IDs,
				 double* data,
				 int vectorIndex=0) = 0;

    /** Write the vector's contents into the specified file.
	@param filename Text name of the file to be created or overwritten.
	If in a parallel environment, each processor will take turns writing
	into the file.
	@param matrixMarketFormat Optional argument, defaults to true. If true
	the contents of the file will be MatrixMarket real array format. If not
	true, the contents of the file will contain the vector's global
	dimension on the first line, and all following lines will contain a
	space-separated pair with global index first and coefficient value
	second.
	@return error-code 0 if successful, -1 if some error occurs such as
	failure to open file.
     */
    virtual int writeToFile(const char* filename,
			    bool matrixMarketFormat=true) = 0;

    /** Write the vector's contents to the specified ostream.
	@param ostrm ostream to be written to.
	@param matrixMarketFormat Optional argument, defaults to true. If true
	the contents of the vector will be written in MatrixMarket real array
	format. If not true, the stream will be given the vector's global
	dimension on the first line, and all following lines will contain a
	space-separated pair with global index first and coefficient value
	second.
	@return error-code 0 if successful, -1 if some error occurs such as
	failure to open file.
     */
    virtual int writeToStream(FEI_OSTREAM& ostrm,
			      bool matrixMarketFormat=true) = 0;

  };//class Vector
}//namespace fei

#ifndef _fei_ostream_ops_hpp_
#include <fei_ostream_ops.hpp>
#endif

#endif // _fei_Vector_hpp_
