
//@HEADER
// ***********************************************************************
// 
//                Komplex: Complex Linear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#include "Epetra_Map.h"
#include "Epetra_Util.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Komplex_LinearProblem.h"
#include "Teuchos_Assert.hpp"
//==============================================================================
Komplex_LinearProblem::Komplex_LinearProblem(double c0r, double c0i, const Epetra_RowMatrix & A0,
					     double c1r, double c1i, const Epetra_RowMatrix & A1,
					     const Epetra_MultiVector & Xr, const Epetra_MultiVector & Xi,
					     const Epetra_MultiVector & Br, const Epetra_MultiVector & Bi) {
  bool firstTime = true;
  TEUCHOS_TEST_FOR_EXCEPT(ProcessValues(c0r, c0i, A0, c1r, c1i, A1, Xr, Xi, Br, Bi, firstTime)!=0);

  TEUCHOS_TEST_FOR_EXCEPT(KomplexMatrix_->FillComplete(*KomplexMatrixDomainMap_,*KomplexMatrixRangeMap_)!=0);

  KomplexProblem_ = Teuchos::rcp(new Epetra_LinearProblem(KomplexMatrix_.get(), KomplexLHS_.get(), KomplexRHS_.get()));

}
//==============================================================================
int Komplex_LinearProblem::UpdateValues(double c0r, double c0i, const Epetra_RowMatrix & A0,
					     double c1r, double c1i, const Epetra_RowMatrix & A1,
					     const Epetra_MultiVector & Xr, const Epetra_MultiVector & Xi,
					     const Epetra_MultiVector & Br, const Epetra_MultiVector & Bi) {
  bool firstTime = false;
  KomplexMatrix_->PutScalar(0.0);
  TEUCHOS_TEST_FOR_EXCEPT(ProcessValues(c0r, c0i, A0, c1r, c1i, A1, Xr, Xi, Br, Bi, firstTime)!=0);
  return(0);
}
//==============================================================================
int Komplex_LinearProblem::ProcessValues(double c0r, double c0i, const Epetra_RowMatrix & A0,
					 double c1r, double c1i, const Epetra_RowMatrix & A1,
					 const Epetra_MultiVector & Xr, const Epetra_MultiVector & Xi,
					 const Epetra_MultiVector & Br, const Epetra_MultiVector & Bi,
					 bool firstTime) {

  TEUCHOS_TEST_FOR_EXCEPT(TestMaps(A0, A1, Xr, Xi, Br, Bi)!=0);
  TEUCHOS_TEST_FOR_EXCEPT(InitMatrixAccess(A0, A1)!=0);
  TEUCHOS_TEST_FOR_EXCEPT(ConstructKomplexMaps(A0.OperatorDomainMap(), A0.OperatorRangeMap(), A0.RowMatrixRowMap())!=0);

  if (firstTime) {
    KomplexMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*KomplexMatrixRowMap_,0));
    KomplexLHS_ = Teuchos::rcp(new Epetra_MultiVector(*KomplexMatrixDomainMap_, Xr.NumVectors()));
    KomplexRHS_ = Teuchos::rcp(new Epetra_MultiVector(*KomplexMatrixRangeMap_, Xr.NumVectors()));
  }

  int NumMyRows = A0.NumMyRows();
  // Process X and B values
  for (int j=0; j<Xr.NumVectors(); j++) {
    double *localKX = &((*KomplexLHS_)[j][0]);
    double *localKB = &((*KomplexRHS_)[j][0]);
    for (int i=0; i< NumMyRows; i++) {
      localKX[2*i] = Xr[j][i];
      localKX[2*i+1] = Xi[j][i];
      localKB[2*i] = Br[j][i];
      localKB[2*i+1] = Bi[j][i];
    }
  }

  int NumIndices0;
  int NumIndices1;
  int * Indices0;
  int * Indices1;
  double * Values0;
  double * Values1;

  std::vector<int> KIndices0(MaxNumMyEntries0_);
  std::vector<double> KValues0(MaxNumMyEntries0_);
  std::vector<int> KIndices1(MaxNumMyEntries1_);
  std::vector<double> KValues1(MaxNumMyEntries1_);
  
  int * A0_ColGids = A0.RowMatrixColMap().MyGlobalElements();
  int * A1_ColGids = A1.RowMatrixColMap().MyGlobalElements();

  // For each row of A0 and A1 find the real and imag parts and place entries:
  // A_r -A_i
  // A_i  A_r

  for ( int i=0; i<NumMyRows; i++) {
    // Get ith row
    TEUCHOS_TEST_FOR_EXCEPT(GetRow(i, A0, A1, NumIndices0, Values0, Indices0, NumIndices1, Values1, Indices1)<0);
    int globalRow = A0.RowMatrixRowMap().GID(i);

    // Process real part due to c0r*A0
    if (c0r!=0.0) {
      // Upper Left
      for (int j=0; j<NumIndices0; j++) {
	KIndices0[j] = 2*A0_ColGids[Indices0[j]];
	KValues0[j] = c0r*Values0[j];
      }
      TEUCHOS_TEST_FOR_EXCEPT(PutRow(2*globalRow, NumIndices0, &KValues0[0], &KIndices0[0], firstTime)<0);

      // Lower Right
      for (int j=0; j<NumIndices0; j++) {
	KIndices0[j]++;
	// KValues0[j] = c0r*Values0[j];
      }
      TEUCHOS_TEST_FOR_EXCEPT(PutRow(2*globalRow+1, NumIndices0, &KValues0[0], &KIndices0[0], firstTime)<0);
    }
    // Process imag part due to c0i*A0
    if (c0i!=0.0) {
      // Upper Right (negate values)
      for (int j=0; j<NumIndices0; j++) {
	KIndices0[j] = 2*A0_ColGids[Indices0[j]]+1;
	KValues0[j] = - c0i*Values0[j];
      }
      TEUCHOS_TEST_FOR_EXCEPT(PutRow(2*globalRow, NumIndices0, &KValues0[0],  &KIndices0[0], firstTime)<0);
      // Lower Left
      for (int j=0; j<NumIndices0; j++) {
	KIndices0[j]--;
	KValues0[j] *= -1.0;
      }
      TEUCHOS_TEST_FOR_EXCEPT(PutRow(2*globalRow+1, NumIndices0, &KValues0[0], &KIndices0[0], firstTime)<0);
    }
    // Process real part due to c1r*A1
    if (c1r!=0.0) {
      // Upper Left
      for (int j=0; j<NumIndices1; j++) {
	KIndices1[j] = 2*A1_ColGids[Indices1[j]];
	KValues1[j] = c1r*Values1[j];
      }
      TEUCHOS_TEST_FOR_EXCEPT(PutRow(2*globalRow, NumIndices1, &KValues1[0], &KIndices1[0], firstTime)<0);

      // Lower Right
      for (int j=0; j<NumIndices1; j++) {
	KIndices1[j]++;
	//KValues1[j] = c1r*Values1[j];
      }
      TEUCHOS_TEST_FOR_EXCEPT(PutRow(2*globalRow+1, NumIndices1, &KValues1[0], &KIndices1[0], firstTime)<0);
    }
    // Process imag part due to c1i*A1
    if (c1i!=0.0) {
      // Upper Right (negate values)
      for (int j=0; j<NumIndices1; j++) {
	KIndices1[j] = 2*A1_ColGids[Indices1[j]]+1;
	KValues1[j] = - c1i*Values1[j];
      }
      TEUCHOS_TEST_FOR_EXCEPT(PutRow(2*globalRow, NumIndices1, &KValues1[0], &KIndices1[0], firstTime)<0);
      // Lower Left
      for (int j=0; j<NumIndices1; j++) {
	KIndices1[j]--;
	KValues1[j] *= -1.0;
      }
      TEUCHOS_TEST_FOR_EXCEPT(PutRow(2*globalRow+1, NumIndices1, &KValues1[0], &KIndices1[0], firstTime)<0);
    }
  }
  return(0);
}
//==============================================================================
Komplex_LinearProblem::~Komplex_LinearProblem(){

}
//==============================================================================
int Komplex_LinearProblem::TestMaps (const Epetra_RowMatrix & A0, const Epetra_RowMatrix & A1,
				      const Epetra_MultiVector & Xr, const Epetra_MultiVector & Xi,
				      const Epetra_MultiVector & Br, const Epetra_MultiVector & Bi){

  TEUCHOS_TEST_FOR_EXCEPT(!A0.RowMatrixRowMap().SameAs(A1.RowMatrixRowMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!A0.OperatorDomainMap().SameAs(A1.OperatorDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!A0.OperatorRangeMap().SameAs(A1.OperatorRangeMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!A0.OperatorDomainMap().SameAs(Xr.Map()));
  TEUCHOS_TEST_FOR_EXCEPT(!A0.OperatorRangeMap().SameAs(Br.Map()));
  TEUCHOS_TEST_FOR_EXCEPT(!A0.OperatorDomainMap().SameAs(Xi.Map()));
  TEUCHOS_TEST_FOR_EXCEPT(!A0.OperatorRangeMap().SameAs(Bi.Map()));

  // Test number of vectors also
  TEUCHOS_TEST_FOR_EXCEPT(Xr.NumVectors()!=Xi.NumVectors());
  TEUCHOS_TEST_FOR_EXCEPT(Xr.NumVectors()!=Br.NumVectors());
  TEUCHOS_TEST_FOR_EXCEPT(Xr.NumVectors()!=Bi.NumVectors());

  return(0);

}
//==============================================================================
int Komplex_LinearProblem::ConstructKomplexMaps(const Epetra_Map & A0DomainMap, const Epetra_Map & A0RangeMap,
						const Epetra_Map & A0RowMap) {
    
  // Always make a RowMap
  TEUCHOS_TEST_FOR_EXCEPT(MakeKomplexMap(A0RowMap, KomplexMatrixRowMap_)<0);

  // Remaining maps could be the same as the row map, so we check
  if (!A0RowMap.SameAs(A0RangeMap)) {
    TEUCHOS_TEST_FOR_EXCEPT(MakeKomplexMap(A0RangeMap, KomplexMatrixRangeMap_)<0);
  }
  else KomplexMatrixRangeMap_ = KomplexMatrixRowMap_;

  if (!A0RowMap.SameAs(A0DomainMap)) {
    TEUCHOS_TEST_FOR_EXCEPT(MakeKomplexMap(A0DomainMap, KomplexMatrixDomainMap_)<0);
  }
  else KomplexMatrixDomainMap_ = KomplexMatrixRowMap_;

  /*cout << "A0RowMap = " << A0RowMap << endl;
    cout << "KomplexMatrixRowMap_ = " << *KomplexMatrixRowMap_ << endl;
    cout << "A0RangeMap = " << A0RangeMap << endl;
    cout << "KomplexMatrixRangeMap_ = " << *KomplexMatrixRangeMap_ << endl;
    cout << "A0DomainMap = " << A0DomainMap << endl;
    cout << "KomplexMatrixDomainMap_ = " << *KomplexMatrixDomainMap_ << endl;
  */
  return(0);
}
//==============================================================================
int Komplex_LinearProblem::MakeKomplexMap(const Epetra_Map & Map, Teuchos::RCP<Epetra_Map> & KMap) {

  if (Map.LinearMap()) {
    KMap = Teuchos::rcp(new Epetra_Map(-1, Map.NumMyElements()*2, Map.IndexBase(), Map.Comm()));
  }
  else {
    Epetra_IntSerialDenseVector KMapGIDs(Map.NumMyElements()*2);
    int * myGlobalElements = Map.MyGlobalElements();
    for (int i=0; i<Map.NumMyElements(); i++) {
      KMapGIDs[2*i] = 2*myGlobalElements[i];
      KMapGIDs[2*i+1] = 2*myGlobalElements[i]+1;
    }
    KMap = Teuchos::rcp(new Epetra_Map(-1, KMapGIDs.Length(), KMapGIDs.Values(), Map.IndexBase(), Map.Comm()));
  }
  return(0);
}
//==============================================================================
int Komplex_LinearProblem::InitMatrixAccess(const Epetra_RowMatrix & A0, const Epetra_RowMatrix & A1) {

  MaxNumMyEntries0_ = A0.MaxNumEntries();
  MaxNumMyEntries1_ = A1.MaxNumEntries(); 

  // Cast to CrsMatrix, if possible.  Can save some work.
  CrsA0_ = dynamic_cast<const Epetra_CrsMatrix *>(&A0);
  CrsA1_ = dynamic_cast<const Epetra_CrsMatrix *>(&A1);
  A0A1AreCrs_ = (CrsA0_!=0 && CrsA1_!=0); // Pointers are non-zero if cast worked
  if (!A0A1AreCrs_) {
    Indices0_.resize(MaxNumMyEntries0_);
    Values0_.resize(MaxNumMyEntries0_);
    Indices1_.resize(MaxNumMyEntries1_);
    Values1_.resize(MaxNumMyEntries1_);
  }

  return(0);
}
//==============================================================================
int Komplex_LinearProblem::GetRow(int Row, const Epetra_RowMatrix & A0, const Epetra_RowMatrix & A1,
				  int & NumIndices0, double * & Values0, int * & Indices0,
				  int & NumIndices1, double * & Values1, int * & Indices1) {

  if (A0A1AreCrs_) { // View of current row
    TEUCHOS_TEST_FOR_EXCEPT(CrsA0_->ExtractMyRowView(Row, NumIndices0, Values0, Indices0)<0); 
    TEUCHOS_TEST_FOR_EXCEPT(CrsA1_->ExtractMyRowView(Row, NumIndices1, Values1, Indices1)<0); 
  }
  else { // Copy of current row
    TEUCHOS_TEST_FOR_EXCEPT(A0.ExtractMyRowCopy(Row, MaxNumMyEntries0_, NumIndices0, &Values0_[0], &Indices0_[0])<0); 
    TEUCHOS_TEST_FOR_EXCEPT(A1.ExtractMyRowCopy(Row, MaxNumMyEntries1_, NumIndices1, &Values1_[0], &Indices1_[0])<0); 
    Values0 = &Values0_[0];
    Indices0 = &Indices0_[0];
    Values1 = &Values1_[0];
    Indices1 = &Indices1_[0];
  } 
  return(0);
}
//==============================================================================
int Komplex_LinearProblem::PutRow(int Row, int & NumIndices, double * Values, int * Indices, bool firstTime) {

  if (firstTime) {
    TEUCHOS_TEST_FOR_EXCEPT(KomplexMatrix_->InsertGlobalValues(Row, NumIndices, Values, Indices)<0);
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT(KomplexMatrix_->SumIntoGlobalValues(Row, NumIndices, Values, Indices)<0);
  }

  return(0);
}
//==============================================================================
int Komplex_LinearProblem::ExtractSolution(Epetra_MultiVector & Xr, Epetra_MultiVector & Xi) {

  int NumMyRows = Xr.MyLength();
  // Process X and B values
  for (int j=0; j<Xr.NumVectors(); j++) {
    double *localKX = &((*KomplexLHS_)[j][0]);
    double *localXr = &(Xr[j][0]);
    double *localXi = &(Xi[j][0]);
    for (int i=0; i< NumMyRows; i++) {
      localXr[i] = localKX[2*i];
      localXi[i] = localKX[2*i+1];
    }
  }

  return(0);
}
