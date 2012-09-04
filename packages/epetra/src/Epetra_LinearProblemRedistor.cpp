
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

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrixTransposer.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_LinearProblemRedistor.h"
#include "Epetra_Util.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
//=============================================================================
Epetra_LinearProblemRedistor::Epetra_LinearProblemRedistor(Epetra_LinearProblem * OrigProblem, 
																													 const Epetra_Map & RedistMap)
  : OrigProblem_(OrigProblem),
		NumProc_(0),
    RedistProblem_(0),
    RedistMap_((Epetra_Map *) &RedistMap),
		Transposer_(0),
    Replicate_(false),
    ConstructTranspose_(false),
    MakeDataContiguous_(false),
    MapGenerated_(false),
    RedistProblemCreated_(false),
		ptr_(0)
{
}
//=============================================================================
Epetra_LinearProblemRedistor::Epetra_LinearProblemRedistor(Epetra_LinearProblem * OrigProblem, 
																													 int NumProc,
																													 bool Replicate)
  : OrigProblem_(OrigProblem),
		NumProc_(NumProc),
    RedistProblem_(0),
    RedistMap_(0),
		Transposer_(0),
    Replicate_(Replicate),
    ConstructTranspose_(false),
    MakeDataContiguous_(false),
    MapGenerated_(false),
    RedistProblemCreated_(false),
		ptr_(0)
{
}
//=============================================================================
Epetra_LinearProblemRedistor::Epetra_LinearProblemRedistor(const Epetra_LinearProblemRedistor& Source)
  : OrigProblem_(Source.OrigProblem_),
		NumProc_(Source.NumProc_),
    RedistProblem_(Source.RedistProblem_),
    RedistMap_(Source.RedistMap_),
		Transposer_(Source.Transposer_),
    Replicate_(Source.Replicate_),
    ConstructTranspose_(Source.ConstructTranspose_),
    MakeDataContiguous_(Source.MakeDataContiguous_),
    RedistProblemCreated_(Source.RedistProblemCreated_),
		ptr_(0)
{

	if (RedistProblem_!=0) RedistProblem_ = new Epetra_LinearProblem(*Source.RedistProblem_);
	if (Transposer_!=0) Transposer_ = new Epetra_RowMatrixTransposer(*Source.Transposer_);
}
//=========================================================================
Epetra_LinearProblemRedistor::~Epetra_LinearProblemRedistor(){

	if (ptr_!=0) {delete [] ptr_; ptr_=0;}
	if (MapGenerated_ && RedistMap_!=0) {delete RedistMap_; RedistMap_=0;}
	if (RedistProblem_!=0) {
		 // If no tranpose, then we must delete matrix (otherwise the transposer must).
		if (!ConstructTranspose_ && RedistProblem_->GetMatrix()!=0)
			delete RedistProblem_->GetMatrix();
		if (RedistProblem_->GetLHS()!=0) delete RedistProblem_->GetLHS();
		if (RedistProblem_->GetRHS()!=0) delete RedistProblem_->GetRHS();
		delete RedistProblem_; RedistProblem_=0;
	}
	if (RedistExporter_!=0) {delete RedistExporter_; RedistExporter_=0;}
	if (Transposer_!=0) {delete Transposer_; Transposer_ = 0;}

}

//=========================================================================
int Epetra_LinearProblemRedistor::GenerateRedistMap() {


  if (MapGenerated_) return(0);

	const Epetra_Map & SourceMap = OrigProblem_->GetMatrix()->RowMatrixRowMap();
	const Epetra_Comm & Comm = SourceMap.Comm();
	int IndexBase = SourceMap.IndexBase();

	int NumProc = Comm.NumProc();

	if (NumProc_!=NumProc) return(-1); // Right now we are only supporting redistribution to all processors.

	// Build a list of contiguous GIDs that will be used in either case below.
	int NumMyRedistElements = 0;
	if ((Comm.MyPID()==0) || Replicate_) NumMyRedistElements = SourceMap.NumGlobalElements64();
	
	// Now build the GID list to broadcast SourceMapGIDs to all processors that need it (or just to PE 0).
	int * ContigIDs = 0;
	if (NumMyRedistElements>0) ContigIDs = new int[NumMyRedistElements];
	for (int i=0; i<NumMyRedistElements; i++) ContigIDs[i] = IndexBase + i;
	
	// Case 1: If the map for the input matrix is not a linear contiguous map, then we have to collect
	// the map indices in order to construct a compatible map containing all GIDs.
	if (!SourceMap.LinearMap()) {

		// First generate a linear map of the same distribution as RowMatrixRowMap
		Epetra_Map SourceLinearMap(-1, SourceMap.NumMyElements(), IndexBase, Comm);
		// Generate Int vector containing GIDs of SourceMap
		Epetra_IntVector SourceMapGIDs(View, SourceLinearMap, SourceMap.MyGlobalElements());
		// Now Build target map for SourceMapGIDs and Importer to get IDs
		Epetra_Map GIDsTargetMap(-1, NumMyRedistElements, ContigIDs, IndexBase, Comm);
		if (NumMyRedistElements>0) delete [] ContigIDs;

		Epetra_Import GIDsImporter(GIDsTargetMap, SourceMap);
		Epetra_IntVector TargetMapGIDs(GIDsTargetMap);

		// Now Send SourceMap GIDs to PE 0, and all processors if Replicate is true
		EPETRA_CHK_ERR(TargetMapGIDs.Import(SourceMapGIDs, GIDsImporter, Insert));

		// Finally, create RedistMap containing all GIDs of SourceMap on PE 0, or all PEs of Replicate is true
		RedistMap_ = new Epetra_Map(-1, NumMyRedistElements, TargetMapGIDs.Values(), IndexBase, Comm);// FIXME long long

	}
	// Case 2: If the map has contiguous IDs then we can simply build the map right away using the list
	// of contiguous GIDs directly
	else
		RedistMap_ = new Epetra_Map(-1, NumMyRedistElements, ContigIDs, IndexBase, Comm);// FIXME long long

	MapGenerated_ = true;

  return(0);
}

//=========================================================================
int Epetra_LinearProblemRedistor::CreateRedistProblem(const bool ConstructTranspose, 
																											const bool MakeDataContiguous, 
																											Epetra_LinearProblem *& RedistProblem) {

	if (RedistProblemCreated_) EPETRA_CHK_ERR(-1);  // This method can only be called once

	Epetra_RowMatrix * OrigMatrix = OrigProblem_->GetMatrix();
	Epetra_MultiVector * OrigLHS = OrigProblem_->GetLHS();
	Epetra_MultiVector * OrigRHS = OrigProblem_->GetRHS();
		
	if (OrigMatrix==0) EPETRA_CHK_ERR(-2); // There is no matrix associated with this Problem


	if (RedistMap_==0) {
		EPETRA_CHK_ERR(GenerateRedistMap());
	}

	RedistExporter_ = new Epetra_Export(OrigProblem_->GetMatrix()->RowMatrixRowMap(), *RedistMap_);

	RedistProblem_ = new Epetra_LinearProblem();
	Epetra_CrsMatrix * RedistMatrix;

	// Check if the tranpose should be create or not
	if (ConstructTranspose) {
		Transposer_ = new Epetra_RowMatrixTransposer(OrigMatrix);
		EPETRA_CHK_ERR(Transposer_->CreateTranspose(MakeDataContiguous, RedistMatrix, RedistMap_));
	}
	else {
		// If not, then just do the redistribution based on the the RedistMap
		RedistMatrix = new Epetra_CrsMatrix(Copy, *RedistMap_, 0);
		// need to do this next step until we generalize the Import/Export ops for CrsMatrix
		Epetra_CrsMatrix * OrigCrsMatrix = dynamic_cast<Epetra_CrsMatrix *>(OrigMatrix); 
		EPETRA_CHK_ERR(RedistMatrix->Export(*OrigCrsMatrix, *RedistExporter_, Add));
		EPETRA_CHK_ERR(RedistMatrix->FillComplete());
	}

	RedistProblem_->SetOperator(RedistMatrix);

	// Now redistribute the RHS and LHS if non-zero

	Epetra_MultiVector * RedistLHS = 0;
	Epetra_MultiVector * RedistRHS = 0;

	int ierr = 0;

	if (OrigLHS!=0) {
		RedistLHS = new Epetra_MultiVector(*RedistMap_, OrigLHS->NumVectors());
		EPETRA_CHK_ERR(RedistLHS->Export(*OrigLHS, *RedistExporter_, Add));
	}
	else ierr = 1;

	if (OrigRHS!=0) {
		RedistRHS = new Epetra_MultiVector(*RedistMap_, OrigLHS->NumVectors());
		EPETRA_CHK_ERR(RedistRHS->Export(*OrigRHS, *RedistExporter_, Add));
	}
	else ierr ++;

	RedistProblem_->SetLHS(RedistLHS);
	RedistProblem_->SetRHS(RedistRHS);

	RedistProblemCreated_ = true;

  return(ierr);
}

//=========================================================================
int Epetra_LinearProblemRedistor::UpdateRedistProblemValues(Epetra_LinearProblem * ProblemWithNewValues) {
	
	if (!RedistProblemCreated_) EPETRA_CHK_ERR(-1);  // This method can only be called after CreateRedistProblem()

	Epetra_RowMatrix * OrigMatrix = ProblemWithNewValues->GetMatrix();
	Epetra_MultiVector * OrigLHS = ProblemWithNewValues->GetLHS();
	Epetra_MultiVector * OrigRHS = ProblemWithNewValues->GetRHS();
		
	if (OrigMatrix==0) EPETRA_CHK_ERR(-2); // There is no matrix associated with this Problem


	Epetra_CrsMatrix * RedistMatrix = dynamic_cast<Epetra_CrsMatrix *>(RedistProblem_->GetMatrix());

	// Check if the tranpose should be create or not
	if (ConstructTranspose_) {
		EPETRA_CHK_ERR(Transposer_->UpdateTransposeValues(OrigMatrix));
	}
	else {
		// If not, then just do the redistribution based on the the RedistMap
		
		EPETRA_CHK_ERR(RedistMatrix->PutScalar(0.0));
		// need to do this next step until we generalize the Import/Export ops for CrsMatrix
		Epetra_CrsMatrix * OrigCrsMatrix = dynamic_cast<Epetra_CrsMatrix *>(OrigMatrix); 

		if (OrigCrsMatrix==0) EPETRA_CHK_ERR(-3); // Broken for a RowMatrix at this point
		EPETRA_CHK_ERR(RedistMatrix->Export(*OrigCrsMatrix, *RedistExporter_, Add));
	}

	// Now redistribute the RHS and LHS if non-zero


	if (OrigLHS!=0) {
		EPETRA_CHK_ERR(RedistProblem_->GetLHS()->Export(*OrigLHS, *RedistExporter_, Add));
	}

	if (OrigRHS!=0) {
		EPETRA_CHK_ERR(RedistProblem_->GetRHS()->Export(*OrigRHS, *RedistExporter_, Add));
	}

  return(0);
}

//=========================================================================
// NOTE: This method should be removed and replaced with calls to Epetra_Util_ExtractHbData()
int Epetra_LinearProblemRedistor::ExtractHbData(int & M, int & N, int & nz, int * & ptr, 
																								int * & ind, double * & val, int & Nrhs, 
																								double * & rhs, int & ldrhs, 
																								double * & lhs, int & ldlhs) const {

	Epetra_CrsMatrix * RedistMatrix = dynamic_cast<Epetra_CrsMatrix *>(RedistProblem_->GetMatrix());
	
	if (RedistMatrix==0) EPETRA_CHK_ERR(-1); // This matrix is zero or not an Epetra_CrsMatrix
	if (!RedistMatrix->IndicesAreContiguous()) { // Data must be contiguous for this to work
		EPETRA_CHK_ERR(-2);
	}

	M = RedistMatrix->NumMyRows();
	N = RedistMatrix->NumMyCols();
	nz = RedistMatrix->NumMyNonzeros();
	val = (*RedistMatrix)[0];        // Dangerous, but cheap and effective way to access first element in 

	const Epetra_CrsGraph & RedistGraph = RedistMatrix->Graph();
	ind = RedistGraph[0];  // list of values and indices

	Epetra_MultiVector * LHS = RedistProblem_->GetLHS();
	Epetra_MultiVector * RHS = RedistProblem_->GetRHS();
	Nrhs = RHS->NumVectors();
	if (Nrhs>1) {
		if (!RHS->ConstantStride()) {EPETRA_CHK_ERR(-3)}; // Must have strided vectors
		if (!LHS->ConstantStride()) {EPETRA_CHK_ERR(-4)}; // Must have strided vectors
	}
	ldrhs = RHS->Stride();
	rhs = (*RHS)[0]; // Dangerous but effective (again)
	ldlhs = LHS->Stride();
	lhs = (*LHS)[0];

	// Finally build ptr vector

	if (ptr_==0) {
		ptr_ = new int[M+1];
		ptr_[0] = 0;
		for (int i=0; i<M; i++) ptr_[i+1] = ptr_[i] + RedistGraph.NumMyIndices(i);
	}
	ptr = ptr_;
	
  return(0);
}

//=========================================================================
int Epetra_LinearProblemRedistor::UpdateOriginalLHS(Epetra_MultiVector * LHS) {

	if (LHS==0) {EPETRA_CHK_ERR(-1);}
	if (!RedistProblemCreated_) EPETRA_CHK_ERR(-2);

	EPETRA_CHK_ERR(LHS->Import(*RedistProblem_->GetLHS(), *RedistExporter_, Average));

  return(0);
}

//=========================================================================
int Epetra_LinearProblemRedistor::UpdateRedistRHS(Epetra_MultiVector * RHS) {

	if (RHS==0) {EPETRA_CHK_ERR(-1);}
	if (!RedistProblemCreated_) EPETRA_CHK_ERR(-2);

	EPETRA_CHK_ERR(RedistProblem_->GetRHS()->Export(*RHS, *RedistExporter_, Insert));

  return(0);
}



