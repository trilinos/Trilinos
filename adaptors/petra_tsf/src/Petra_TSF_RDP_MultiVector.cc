
#include "Petra_TSF_RDP_MultiVector.h"


//=============================================================================

// Petra_BlockMap Constructor

Petra_TSF_RDP_MultiVector::Petra_TSF_RDP_MultiVector(const Petra_BlockMap& Map, int NumVectors)
  : Petra_RDP_MultiVector(Map, NumVectors){}
//==========================================================================

// Copy Constructor

Petra_TSF_RDP_MultiVector::Petra_TSF_RDP_MultiVector(const Petra_RDP_MultiVector& Source)
  :Petra_RDP_MultiVector(Source) {}
//==========================================================================

// This constructor copies in or makes view of a standard Fortran array

Petra_TSF_RDP_MultiVector::Petra_TSF_RDP_MultiVector(Petra_DataAccess CV, const Petra_BlockMap& Map, 
						     double *A, int MyLDA, int NumVectors)
  : Petra_RDP_MultiVector(CV, Map, A, MyLDA, NumVectors) {}

//==========================================================================

// This constructor copies in or makes view of a C/C++ array of pointer

Petra_TSF_RDP_MultiVector::Petra_TSF_RDP_MultiVector(Petra_DataAccess CV, const Petra_BlockMap& Map, 
						     double **ArrayOfPointers, int NumVectors)
  : Petra_RDP_MultiVector(CV, Map, ArrayOfPointers, NumVectors) {}

//==========================================================================

// This constructor copies or makes view of selected vectors, specified in Indices, 
// from an existing MultiVector

Petra_TSF_RDP_MultiVector::Petra_TSF_RDP_MultiVector(Petra_DataAccess CV, 
						     const Petra_RDP_MultiVector& Source, 
						     int *Indices, int NumVectors)
  : Petra_RDP_MultiVector(CV, Source, Indices, NumVectors) {}

//==========================================================================

// This interface copies or makes view of a range of vectors from an existing MultiVector

Petra_TSF_RDP_MultiVector::Petra_TSF_RDP_MultiVector(Petra_DataAccess CV, 
						     const Petra_RDP_MultiVector& Source, 
						     int StartIndex, int NumVectors)
  : Petra_RDP_MultiVector(CV, Source, StartIndex, NumVectors){}

//=========================================================================
Petra_TSF_RDP_MultiVector::~Petra_TSF_RDP_MultiVector(){}
//=========================================================================
int Petra_TSF_RDP_MultiVector::CloneZero( TSF_RDP_MultiVector * & A, const int NumVectors) {

  Petra_TSF_RDP_MultiVector * A1 = new Petra_TSF_RDP_MultiVector(Map(), NumVectors); assert(A1);
  A = dynamic_cast<TSF_RDP_MultiVector *>(A1);
  return(0);

}
//=========================================================================
int Petra_TSF_RDP_MultiVector::CloneCopy( TSF_RDP_MultiVector * & A) {

  Petra_TSF_RDP_MultiVector * A1 = new Petra_TSF_RDP_MultiVector(*this); assert(A1);
  A = dynamic_cast<TSF_RDP_MultiVector *>(A1);
  return(0);

}
//=========================================================================
int Petra_TSF_RDP_MultiVector::CloneView( TSF_RDP_MultiVector * & A, int * Indices, int NumVectors) {

  Petra_TSF_RDP_MultiVector * A1 = new Petra_TSF_RDP_MultiVector(View, *this, Indices, NumVectors); assert(A1);
  A = dynamic_cast<TSF_RDP_MultiVector *>(A1);
  return(0);

}

//=========================================================================
int Petra_TSF_RDP_MultiVector::Dot(const TSF_RDP_MultiVector & A, 
				double *Result) {
  const Petra_TSF_RDP_MultiVector *A1 = dynamic_cast<const Petra_TSF_RDP_MultiVector *>(&A); assert(A1);
  return(Petra_RDP_MultiVector::Dot(*A1, Result));

}
//=========================================================================
int Petra_TSF_RDP_MultiVector::Abs(const TSF_RDP_MultiVector& A) {
  
  const Petra_TSF_RDP_MultiVector *A1 = dynamic_cast<const Petra_TSF_RDP_MultiVector *>(&A); assert(A1);
  return(Petra_RDP_MultiVector::Abs(*A1));

}
//=========================================================================
int Petra_TSF_RDP_MultiVector::Reciprocal(const TSF_RDP_MultiVector& A) {

  const Petra_TSF_RDP_MultiVector *A1 = dynamic_cast<const Petra_TSF_RDP_MultiVector *>(&A); assert(A1);
  return(Petra_RDP_MultiVector::Reciprocal(*A1));

}

  //=========================================================================
  int Petra_TSF_RDP_MultiVector::Scale (double ScalarA, const TSF_RDP_MultiVector& A) {

  const Petra_TSF_RDP_MultiVector *A1 = dynamic_cast<const Petra_TSF_RDP_MultiVector *>(&A); assert(A1);
  return(Petra_RDP_MultiVector::Scale(ScalarA, *A1));

  }

  //=========================================================================
  int Petra_TSF_RDP_MultiVector::Update(double ScalarA, const TSF_RDP_MultiVector& A, double Scalar) {

  const Petra_TSF_RDP_MultiVector *A1 = dynamic_cast<const Petra_TSF_RDP_MultiVector *>(&A); assert(A1);
  return(Petra_RDP_MultiVector::Update(ScalarA, *A1, Scalar));
 }

//=========================================================================
int Petra_TSF_RDP_MultiVector::Update(double ScalarA, const TSF_RDP_MultiVector& A, 
				  double ScalarB, const TSF_RDP_MultiVector& B, double Scalar) {

  const Petra_TSF_RDP_MultiVector *A1 = dynamic_cast<const Petra_TSF_RDP_MultiVector *>(&A); assert(A1);
  const Petra_TSF_RDP_MultiVector *B1 = dynamic_cast<const Petra_TSF_RDP_MultiVector *>(&B); assert(B1);
  return(Petra_RDP_MultiVector::Update(ScalarA, *A1, ScalarB, *B1, Scalar));
 }


//=========================================================================
int  Petra_TSF_RDP_MultiVector::NormWeighted (const TSF_RDP_MultiVector& Weights, double* Result)  {

  const Petra_TSF_RDP_MultiVector *A1 = dynamic_cast<const Petra_TSF_RDP_MultiVector *>(&Weights); assert(A1);
  return(Petra_RDP_MultiVector::NormWeighted(*A1, Result));
}

  //=========================================================================
  int  Petra_TSF_RDP_MultiVector::Multiply (char TransA, char TransB, double ScalarAB, 
					const TSF_RDP_MultiVector& A, 
					const TSF_RDP_MultiVector& B,
					double Scalar ) {

  const Petra_TSF_RDP_MultiVector *A1 = dynamic_cast<const Petra_TSF_RDP_MultiVector *>(&A); assert(A1);
  const Petra_TSF_RDP_MultiVector *B1 = dynamic_cast<const Petra_TSF_RDP_MultiVector *>(&B); assert(B1);
  return(Petra_RDP_MultiVector::Multiply(TransA, TransB, ScalarAB, *A1, *B1, Scalar));
 }


//=========================================================================
int Petra_TSF_RDP_MultiVector::Multiply(double ScalarAB, const TSF_RDP_MultiVector& A, 
					const TSF_RDP_MultiVector& B,
					double Scalar) {

  const Petra_TSF_RDP_MultiVector *A1 = dynamic_cast<const Petra_TSF_RDP_MultiVector *>(&A); assert(A1);
  const Petra_TSF_RDP_MultiVector *B1 = dynamic_cast<const Petra_TSF_RDP_MultiVector *>(&B); assert(B1);
  return(Petra_RDP_MultiVector::Multiply(ScalarAB, *A1, *B1, Scalar));
  
}
//=========================================================================
int Petra_TSF_RDP_MultiVector::ReciprocalMultiply(double ScalarAB, const TSF_RDP_MultiVector& A, 
						  const TSF_RDP_MultiVector& B,
						  double Scalar) {
  
  const Petra_TSF_RDP_MultiVector *A1 = dynamic_cast<const Petra_TSF_RDP_MultiVector *>(&A); assert(A1);
  const Petra_TSF_RDP_MultiVector *B1 = dynamic_cast<const Petra_TSF_RDP_MultiVector *>(&B); assert(B1);
  return(Petra_RDP_MultiVector::ReciprocalMultiply(ScalarAB, *A1, *B1, Scalar));
}
