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


#ifndef _fei_Vector_Impl_hpp_
#define _fei_Vector_Impl_hpp_

#include <fei_macros.hpp>
#include <fei_VectorTraits.hpp>

#include <fei_VectorTraits_CSVec.hpp>
#include <fei_VectorTraits_FillableVec.hpp>
#include <fei_VectorTraits_LinSysCore.hpp>
#include <fei_VectorTraits_LinProbMgr.hpp>
#include <fei_VectorTraits_FEData.hpp>
#include <snl_fei_FEVectorTraits.hpp>
#include <snl_fei_FEVectorTraits_FED.hpp>
#include <fei_SharedIDs.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_Reducer.hpp>
#include <fei_Logger.hpp>
#include <fei_Vector.hpp>
#include <fei_Vector_core.hpp>
#include <fei_iosfwd.hpp>

#undef fei_file
#define fei_file "fei_Vector_Impl.hpp"

#include <fei_ErrMacros.hpp>

namespace fei {

  /** To be used for vector data assembly, including local assembly of shared
      data. Provides operations for gathering the overlapped data (locally-
      stored shared data) to a non-overlapped data distribution (e.g., send
      shared data to owning processor) and vice-versa for scattering non-
      overlapped data to the overlapped distribution.

      When shared data that is not locally-owned is assembled into this
      object, it will be held locally until the gatherFromOverlap()
      operation is performed. When data that is locally-owned is assembled
      into this object, it will be passed directly to the underlying
      algebraic (non-overlapping) vector.

      When non-locally-owned shared data is requested from this vector object,
      the operation is only guaranteed to succeed if scatterToOverlap() has
      been called. If non-local data has recently been input to this vector,
      but the vector has not been 'synchronized' using gatherFromOverlap() and
      scatterToOverlap(), then that data may be used to answer the request but
      the values may be erroneous due to not including contributions from
      other processors.
  */
  template<typename T>
  class Vector_Impl : public fei::Vector, public fei::Vector_core {
  public:

    /** Constructor that takes a VectorSpace */
    Vector_Impl(fei::SharedPtr<fei::VectorSpace> vecSpace,
	   T* vector, int numLocalEqns,
	   bool isSolutionVector=false,
           bool deleteVector=false);

    /** Destructor */
    virtual ~Vector_Impl();

    /** Return a name describing the run-time type
	of this object.
    */
    const char* typeName() const { return(fei::VectorTraits<T>::typeName()); }

    /** Update 'this' = b*'this' + a*x
     */
    int update(double a,
	       const fei::Vector* x,
	       double b);

    /** Use data in the underlying non-overlapping decomposition to update
	any shared data in the overlapping decomposition.

	If any data is already held for the shared positions, that data will
	be replaced by the data from the 'owning' processor.
    */
    int scatterToOverlap();

    void setCommSizes();

    /** Move any shared data from the overlapping decomposition to the
	underlying non-overlapping decomposition.
    */
    int gatherFromOverlap(bool accumulate = true);

    /** Set a specified scalar throughout the vector. */
    int putScalar(double scalar);

    /** Sum values into the vector, adding to any
	that may already exist at the specified indices.
    */
    int sumIn(int numValues, const int* indices, const double* values,
	      int vectorIndex=0);

    /** Copy values into the vector, overwriting any that may already exist
	at the specified indices.
    */
    int copyIn(int numValues, const int* indices, const double* values,
	       int vectorIndex=0);

    /** Obtain the VectorSpace associated with this vector.
     */
    fei::SharedPtr<fei::VectorSpace> getVectorSpace() const
      { return(get_vector_space()); }

    /** Set the VectorSpace associated with this vector.
     */
    void setVectorSpace(fei::SharedPtr<fei::VectorSpace> vecSpace)
    {
      set_vector_space( vecSpace );
    }

    /** Sum field data into the vector, adding to any coefficients that may
	already exist at the specified locations.
        If the specified fieldID doesn't exist at one or more of the specified
        IDs, then the corresponding positions in the data array will simply
        not be used.
    */
    int sumInFieldData(int fieldID,
		       int idType,
		       int numIDs,
		       const int* IDs,
		       const double* data,
		       int vectorIndex=0);

    /** Copy field data into the vector, overwriting any coefficients that may
	already exist at the specified locations.
        If the specified fieldID doesn't exist at one or more of the specified
        IDs, then the corresponding positions in the data array will simply
        not be used.
    */
    int copyInFieldData(int fieldID,
			int idType,
			int numIDs,
			const int* IDs,
			const double* data,
			int vectorIndex=0);

    int copyInFieldDataLocalIDs(int fieldID,
			int idType,
			int numIDs,
			const int* localIDs,
			const double* data,
			int vectorIndex=0);

    /** Copy field data out of the vector, into the caller-allocated data
	array.
        If the specified fieldID doesn't exist at one or more of the specified
        IDs, then the corresponding positions in the data array will simply
        not be referenced.
    */
    int copyOutFieldData(int fieldID,
			 int idType,
			 int numIDs,
			 const int* IDs,
			 double* data,
			 int vectorIndex=0);

    int writeToFile(const char* filename,
		    bool matrixMarketFormat=true);

    int writeToStream(FEI_OSTREAM& ostrm,
		      bool matrixMarketFormat=true);

    /** Set the library-specific underlying vector object that this
	snl_fei::Vector is filtering data in and out of.
    */
    void setUnderlyingVector(T* vec)
      {
	vector_ = vec;
      }

    /** Get the library-specific underlying vector object that this
	snl_fei::Vector is filtering data in and out of.
    */
    T* getUnderlyingVector() { return( vector_ ); }
    const T* getUnderlyingVector() const { return( vector_ ); }

    int copyOut(int numValues,
		const int* indices,
		double* values,
		int vectorIndex=0) const;

    /** please ignore
     */
    int copyOut_FE(int nodeNumber, int dofOffset, double& value);

  private:
    int giveToUnderlyingVector(int numValues,
			       const int* indices,
			       const double* values,
			       bool sumInto=true,
			       int vectorIndex=0);

    int copyOutOfUnderlyingVector(int numValues,
				  const int* indices,
				  double* values,
				  int vectorIndex=0) const;

    int sumIntoFEVector(int blockID,
			int connOffset,
			int numNodes,
			const int* nodeNumbers,
			const int* numIndicesPerNode,
      const int* dof_ids,
			const double* values);

    T* vector_;
    bool isSolution_;
    bool deleteVector_;

    int localProc_;
    int numProcs_;
    std::string dbgprefix_;
  };//class Vector_Impl

} //namespace fei

//----------------------------------------------------------------------------
template<typename T>
fei::Vector_Impl<T>::Vector_Impl(fei::SharedPtr<fei::VectorSpace> vecSpace,
			   T* vector, int numLocalEqns,
			   bool isSolutionVector,
                           bool deleteVector)
  : Vector_core(vecSpace, numLocalEqns),
    vector_(vector),
    isSolution_(isSolutionVector),
    deleteVector_(deleteVector),
    localProc_(0),
    numProcs_(1),
    dbgprefix_("VecImpl: ")
{
  if (strcmp(snl_fei::FEVectorTraits<T>::typeName(), "unsupported")) {
    setFEVector(true);
  }
  else {
    setFEVector(false);
  }

  localProc_ = fei::localProc(vecSpace->getCommunicator());
  numProcs_ = fei::numProcs(vecSpace->getCommunicator());

  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<" ctor, numLocalEqns="<<numLocalEqns
       <<", typeName: "<<typeName()<<FEI_ENDL;
  }

  std::vector<int> idTypes;
  vecSpace->getIDTypes(idTypes);
  std::vector<int> eqns;
  std::vector<double> zeros;
  for(size_t i=0; i<idTypes.size(); ++i) {
    int idType = idTypes[i];
    fei::SharedIDs<int>& sharedIDs = vecSpace->getSharedIDs(idType);
    const fei::SharedIDs<int>::map_type& idMap = sharedIDs.getSharedIDs();
    fei::SharedIDs<int>::map_type::const_iterator
      iter = idMap.begin(), iterEnd = idMap.end();
    for(; iter!=iterEnd; ++iter) {
      int ID = iter->first;
      int eqn;
      vecSpace->getGlobalIndex(idType, ID, eqn);
      int ndof = vecSpace->getNumDegreesOfFreedom(idType, ID);
      eqns.resize(ndof);
      zeros.resize(ndof, 0.0);
      for(int j=0; j<ndof; ++j) eqns[j] = eqn+j;
      if (!isSolutionVector) {
        sumIn(ndof, &eqns[0], &zeros[0]);
      }
      else {
        copyIn(ndof, &eqns[0], &zeros[0]);
      }
    }
  }

  setCommSizes();
  std::vector<CSVec*>& remoteVecs = remotelyOwned();
  for(size_t i=0; i<remoteVecs.size(); ++i) {
    remoteVecs[i]->clear();
  }
}

//----------------------------------------------------------------------------
template<typename T>
fei::Vector_Impl<T>::~Vector_Impl()
{
  if (deleteVector_) delete vector_;
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::putScalar(double scalar)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"putScalar("<<scalar<<")"<<FEI_ENDL;
  }

  if (haveFEVector()) {
    if (scalar != 0.0) return(-1);
    CHK_ERR( snl_fei::FEVectorTraits<T>::reset(vector_) );
  }
  else {
    CHK_ERR( fei::VectorTraits<T>::setValues(vector_, firstLocalOffset(), scalar) );
  }
  for(unsigned p=0; p<remotelyOwned().size(); ++p) {
    fei::set_values(*(remotelyOwned()[p]), scalar);
  }
  return(0);
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::update(double a,
			       const fei::Vector* x,
			       double b)
{
  const fei::Vector_Impl<T>* sx = dynamic_cast<const fei::Vector_Impl<T>* >(x);
  if (sx != 0) {
    const T* tx = sx->getUnderlyingVector();
    return( fei::VectorTraits<T>::update(vector_, a, tx, b) );
  }
  else {
    return( -1 );
  }
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::scatterToOverlap()
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"scatterToOverlap"<<FEI_ENDL;
  }

  return( Vector_core::scatterToOverlap() );
}

//----------------------------------------------------------------------------
template<typename T>
void fei::Vector_Impl<T>::setCommSizes()
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"setCommSizes"<<FEI_ENDL;
  }

  Vector_core::setCommSizes();
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::gatherFromOverlap(bool accumulate)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"gatherFromOverlap"<<FEI_ENDL;
  }

  return( Vector_core::gatherFromOverlap(accumulate) );
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::sumIn(int numValues,
			      const int* indices, const double* values,
			      int vectorIndex)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"sumIn(n="<<numValues<<")"<<FEI_ENDL;
  }

  return( giveToVector(numValues, indices, values, true, vectorIndex) );
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::copyIn(int numValues,
			       const int* indices, const double* values,
			       int vectorIndex)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"copyIn(n="<<numValues<<")"<<FEI_ENDL;
  }

  return( giveToVector(numValues, indices, values, false, vectorIndex) );
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::giveToUnderlyingVector(int numValues,
					       const int* indices,
					       const double* values,
					       bool sumInto,
					       int vectorIndex)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"giveToUnderlying(";
    for(int i=0; i<numValues; ++i) {
      os << "{"<<indices[i]<<","<<values[i]<<"} ";
    }
    os<<")"<<FEI_ENDL;
  }

  int err = fei::VectorTraits<T>::putValuesIn(vector_, firstLocalOffset(),
					     numValues, indices, values,
					     sumInto, isSolution_, vectorIndex);
  if (err < 0) {
    return(err);
  }
  return(0);
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::copyOutOfUnderlyingVector(int numValues,
						  const int* indices,
						  double* values,
						  int vectorIndex) const
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"copyOutOfUnderlying(n="<<numValues<<")"<<FEI_ENDL;
  }

  return( fei::VectorTraits<T>::copyOut(vector_, firstLocalOffset(),
					     numValues, indices, values,
					     isSolution_, vectorIndex) );
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::sumInFieldData(int fieldID,
				       int idType,
				       int numIDs,
				       const int* IDs,
				       const double* data,
				       int vectorIndex)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"sumInFieldData(n="<<numIDs<<")"<<FEI_ENDL;
  }

  return( assembleFieldData(fieldID, idType, numIDs, IDs, data, true, vectorIndex));
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::copyInFieldData(int fieldID,
					int idType,
					int numIDs,
					const int* IDs,
					const double* data,
					int vectorIndex)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"copyInFieldData(n="<<numIDs<<")"<<FEI_ENDL;
  }

  return(assembleFieldData(fieldID, idType, numIDs, IDs, data, false, vectorIndex));
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::copyInFieldDataLocalIDs(int fieldID,
					int idType,
					int numIDs,
					const int* localIDs,
					const double* data,
					int vectorIndex)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"copyInFieldDataLocalIDs(n="<<numIDs<<")"<<FEI_ENDL;
  }

  return(assembleFieldDataLocalIDs(fieldID, idType, numIDs, localIDs, data, false, vectorIndex));
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::copyOut_FE(int nodeNumber, int dofOffset, double& value)
{
  return( snl_fei::FEVectorTraits<T>::copyOut(vector_, nodeNumber, dofOffset, value) );
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::copyOutFieldData(int fieldID,
					 int idType,
					 int numIDs,
					 const int* IDs,
					 double* data,
					 int vectorIndex)
{
  return( Vector_core::copyOutFieldData(fieldID, idType, numIDs, IDs, data,
					 vectorIndex));
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::writeToFile(const char* filename,
				    bool matrixMarketFormat)
{
  return( Vector_core::writeToFile(filename, matrixMarketFormat) );
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::writeToStream(FEI_OSTREAM& ostrm,
				      bool matrixMarketFormat)
{
  return( Vector_core::writeToStream(ostrm, matrixMarketFormat) );
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::copyOut(int numValues,
				const int* indices,
				double* values,
				int vectorIndex) const
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"copyOut(n="<<numValues<<")"<<FEI_ENDL;
  }

  return( Vector_core::copyOut(numValues, indices, values, vectorIndex) );
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Vector_Impl<T>::sumIntoFEVector(int blockID,
					int connOffset,
					int numNodes,
					const int* nodeNumbers,
					const int* numIndicesPerNode,
          const int* dof_ids,
					const double* values)
{
  return( snl_fei::FEVectorTraits<T>::sumInElemVector(vector_, blockID, connOffset,
					      numNodes, nodeNumbers,
					      numIndicesPerNode, dof_ids, values) );
}

#undef fei_file
#define fei_file "unknown_file"

#endif // _fei_Vector_Impl_hpp_

