
#include "Petra_TSF_MultiVector.h"

namespace Petra_TSF {
//=============================================================================

// Petra_BlockMap Constructor
template<class scalarType>
MultiVector<scalarType>::MultiVector(const Petra_BlockMap& Map, int NumVectors)
  : Petra_RDP_MultiVector(Map, NumVectors){}
//==========================================================================

// Copy Constructor

template<class scalarType>
MultiVector<scalarType>::MultiVector(const Petra_TSF::MultiVector<scalarType>& Source)
  :Petra_RDP_MultiVector(Source) {}
//==========================================================================

// This constructor copies in or makes view of a standard Fortran array

template<class scalarType>
MultiVector<scalarType>::MultiVector(Petra_DataAccess CV, const Petra_BlockMap& Map, 
						     double *A, int MyLDA, int NumVectors)
  : Petra_RDP_MultiVector(CV, Map, A, MyLDA, NumVectors) {}

//==========================================================================

// This constructor copies in or makes view of a C/C++ array of pointer

template<class scalarType>
MultiVector::MultiVector(Petra_DataAccess CV, const Petra_BlockMap& Map, 
						     double **ArrayOfPointers, int NumVectors)
  : Petra_RDP_MultiVector(CV, Map, ArrayOfPointers, NumVectors) {}

//==========================================================================

// This constructor copies or makes view of selected vectors, specified in Indices, 
// from an existing MultiVector

template<class scalarType>
MultiVector::MultiVector(Petra_DataAccess CV, 
			 const Petra_RDP_MultiVector& Source, 
			 int *Indices, int NumVectors)
  : Petra_RDP_MultiVector(CV, Source, Indices, NumVectors) {}

//==========================================================================

// This interface copies or makes view of a range of vectors from an existing MultiVector

template<class scalarType>
MultiVector::MultiVector(Petra_DataAccess CV, 
						     const Petra_RDP_MultiVector& Source, 
						     int StartIndex, int NumVectors)
  : Petra_RDP_MultiVector(CV, Source, StartIndex, NumVectors){}

//=========================================================================
int MultiVector::CloneZero( TSF_MultiVector * & A, const int NumVectors) {

  Petra_TSF::MultiVector * A1 = new Petra_TSF::MultiVector(Map(), NumVectors); assert(A1);
  A = dynamic_cast<TSF_MultiVector *>(A1);
  return(0);

}
//=========================================================================
int MultiVector::CloneCopy( TSF_MultiVector * & A) {

  Petra_TSF::MultiVector * A1 = new Petra_TSF::MultiVector(*this); assert(A1);
  A = dynamic_cast<TSF_MultiVector *>(A1);
  return(0);

}
//=========================================================================
int MultiVector::CloneView( TSF_MultiVector * & A, int * Indices, int NumVectors) {

  Petra_TSF::MultiVector * A1 = new Petra_TSF::MultiVector(View, *this, Indices, NumVectors); assert(A1);
  A = dynamic_cast<TSF_MultiVector *>(A1);
  return(0);

}

//=========================================================================
int MultiVector::Dot(const TSF_MultiVector & A, 
				double *Result) {
  const Petra_TSF::MultiVector *A1 = dynamic_cast<const Petra_TSF::MultiVector *>(&A); assert(A1);
  return(Petra_RDP_MultiVector::Dot(*A1, Result));

}
//=========================================================================
int MultiVector::Abs(const TSF_MultiVector& A) {
  
  const Petra_TSF::MultiVector *A1 = dynamic_cast<const Petra_TSF::MultiVector *>(&A); assert(A1);
  return(Petra_RDP_MultiVector::Abs(*A1));

}
//=========================================================================
int MultiVector::Reciprocal(const TSF_MultiVector& A) {

  const Petra_TSF::MultiVector *A1 = dynamic_cast<const Petra_TSF::MultiVector *>(&A); assert(A1);
  return(Petra_RDP_MultiVector::Reciprocal(*A1));

}

  //=========================================================================
  int MultiVector::Scale (double ScalarA, const TSF_MultiVector& A) {

  const Petra_TSF::MultiVector *A1 = dynamic_cast<const Petra_TSF::MultiVector *>(&A); assert(A1);
  return(Petra_RDP_MultiVector::Scale(ScalarA, *A1));

  }

  //=========================================================================
  int MultiVector::Update(double ScalarA, const TSF_MultiVector& A, double Scalar) {

  const Petra_TSF::MultiVector *A1 = dynamic_cast<const Petra_TSF::MultiVector *>(&A); assert(A1);
  return(Petra_RDP_MultiVector::Update(ScalarA, *A1, Scalar));
 }

//=========================================================================
int MultiVector::Update(double ScalarA, const TSF_MultiVector& A, 
				  double ScalarB, const TSF_MultiVector& B, double Scalar) {

  const Petra_TSF::MultiVector *A1 = dynamic_cast<const Petra_TSF::MultiVector *>(&A); assert(A1);
  const Petra_TSF::MultiVector *B1 = dynamic_cast<const Petra_TSF::MultiVector *>(&B); assert(B1);
  return(Petra_RDP_MultiVector::Update(ScalarA, *A1, ScalarB, *B1, Scalar));
 }


//=========================================================================
int  MultiVector::NormWeighted (const TSF_MultiVector& Weights, double* Result)  {

  const Petra_TSF::MultiVector *A1 = dynamic_cast<const Petra_TSF::MultiVector *>(&Weights); assert(A1);
  return(Petra_RDP_MultiVector::NormWeighted(*A1, Result));
}

  //=========================================================================
  int  MultiVector::Multiply (char TransA, char TransB, double ScalarAB, 
					const TSF_MultiVector& A, 
					const TSF_MultiVector& B,
					double Scalar ) {

  const Petra_TSF::MultiVector *A1 = dynamic_cast<const Petra_TSF::MultiVector *>(&A); assert(A1);
  const Petra_TSF::MultiVector *B1 = dynamic_cast<const Petra_TSF::MultiVector *>(&B); assert(B1);
  return(Petra_RDP_MultiVector::Multiply(TransA, TransB, ScalarAB, *A1, *B1, Scalar));
 }


//=========================================================================
int MultiVector::Multiply(double ScalarAB, const TSF_MultiVector& A, 
					const TSF_MultiVector& B,
					double Scalar) {

  const Petra_TSF::MultiVector *A1 = dynamic_cast<const Petra_TSF::MultiVector *>(&A); assert(A1);
  const Petra_TSF::MultiVector *B1 = dynamic_cast<const Petra_TSF::MultiVector *>(&B); assert(B1);
  return(Petra_RDP_MultiVector::Multiply(ScalarAB, *A1, *B1, Scalar));
  
}
//=========================================================================
int MultiVector::ReciprocalMultiply(double ScalarAB, const TSF_MultiVector& A, 
						  const TSF_MultiVector& B,
						  double Scalar) {
  
  const Petra_TSF::MultiVector *A1 = dynamic_cast<const Petra_TSF::MultiVector *>(&A); assert(A1);
  const Petra_TSF::MultiVector *B1 = dynamic_cast<const Petra_TSF::MultiVector *>(&B); assert(B1);
  return(Petra_RDP_MultiVector::ReciprocalMultiply(ScalarAB, *A1, *B1, Scalar));
}
} // namespace Petra_TSF
