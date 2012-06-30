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

#include "Komplex_RowMatrix.hpp" 

//==========================================================================
Komplex_RowMatrix::Komplex_RowMatrix(Epetra_DataAccess CV, 
						 const Epetra_BlockMap& RowMap,
					       int NumEntriesPerRow, 
						 Komplex_KForms KForm) 
 :
{
}

//==========================================================================
Komplex_RowMatrix::Komplex_RowMatrix(Epetra_DataAccess CV, 
						 const Epetra_BlockMap& RowMap,
 						 int* NumEntriesPerRow, 
		 				 Komplex_KForms KForm)
 : 
{
}

//==========================================================================
Komplex_RowMatrix::Komplex_RowMatrix(Epetra_DataAccess CV,
						 const Epetra_BlockMap& RowMap,
 						 const Epetra_BlockMap& ColMap,
						 int NumEntriesPerRow,
						 Komplex_KForms KForm)
 :  
{
}

//==========================================================================
Komplex_RowMatrix::Komplex_RowMatrix(Epetra_DataAccess CV,
						 const Epetra_BlockMap& RowMap,
 						 const Epetra_BlockMap& ColMap,
						 int* NumEntriesPerRow,
						 Komplex_KForms KForm)
 : 
{
}

//==========================================================================
Komplex_RowMatrix::Komplex_RowMatrix(Epetra_DataAccess CV,
						 const Epetra_CrsGraph& Graph,
 						 Komplex_KForms KForm)
 : 
{
}

//==========================================================================
Komplex_RowMatrix::Komplex_RowMatrix(double cr,
						 double ci,
						 Epetra_RowMatrix* A,
 						 Komplex_KForms KForm)
 : 
{
}

//==========================================================================
Komplex_RowMatrix::Komplex_RowMatrix(double c0r,
						 double c0i,
						 Epetra_RowMatrix* A0,
 						 double c1r,
						 double c1i,
						 Epetra_RowMatrix* A1,
						 Komplex_KForms KForm)
 : 
{
}

//==========================================================================
Komplex_RowMatrix::Komplex_RowMatrix(const Komplex_RowMatrix& Matrix)
 : 
{
}

//==========================================================================
Komplex_RowMatrix::~Komplex_RowMatrix() {
}

//==========================================================================
Komplex_RowMatrix & Komplex_RowMatrix::operator=(const Komplex_RowMatrix& src) {
}

//==========================================================================
bool Komplex_RowMatrix::Filled() const {
}

//==========================================================================
bool Komplex_RowMatrix::StorageOptimized() const {
}

//==========================================================================
bool Komplex_RowMatrix::IndicesAreGlobal() const {
}

//==========================================================================
bool Komplex_RowMatrix::IndicesAreLocal() const {
}

//==========================================================================
bool Komplex_RowMatrix::IndicesAreContiguous() const {
}

//==========================================================================
bool Komplex_RowMatrix::NoDiagonal() const {
}

//==========================================================================
int Komplex_RowMatrix::IndexBase() const {
}

//==========================================================================
const Epetra_CrsGraph & Komplex_RowMatrix::Graph() const {
}

//==========================================================================
const Epetra_Map & Komplex_RowMatrix::RowMap() const {
}

//==========================================================================
const Epetra_Map & Komplex_RowMatrix::ColMap() const {
}

//==========================================================================
const Epetra_Map & Komplex_RowMatrix::DomainMap() const {
}

//==========================================================================
const Epetra_Map & Komplex_RowMatrix::RangeMap() const {
}

//==========================================================================
const Epetra_Import * Komplex_RowMatrix::Importer() const {
}

//==========================================================================
const Epetra_Export * Komplex_RowMatrix::Exporter() const {
}

//==========================================================================
const Epetra_Comm & Komplex_RowMatrix::Comm() const {
}

//==========================================================================
int Komplex_RowMatrix::PutScalar(double ScalarConstant) {
}

//==========================================================================
int Komplex_RowMatrix::Scale(double ScalarConstant) {
}

//==========================================================================
int Komplex_RowMatrix::ReplaceDiagonalValues(const Epetra_Vector& Diagonal) {
}

//==========================================================================
int Komplex_RowMatrix::FillComplete() {
}

//==========================================================================
int Komplex_RowMatrix::FillComplete(const Epetra_Map& DomainMap,
					      const Epetra_Map& RangeMap) {
}

//==========================================================================
int Komplex_RowMatrix::Multiply(bool TransA,
					  const Epetra_MultiVector& X,
 					  Epetra_MultiVector& Y) const {
return(-1);
}

//==========================================================================
int Komplex_RowMatrix::Solve(bool Upper,
				     bool Trans,
                             bool UnitDiagonal,
				     const Epetra_Vector& x,
 				     Epetra_Vector& y) const {
return(-1);
}

//==========================================================================
int Komplex_RowMatrix::Solve(bool Upper,
				     bool Trans,
				     bool UnitDiagonal,
				     const Epetra_MultiVector& X,
 				     Epetra_MultiVector& Y) const {
return(-1);
}

//==========================================================================
int Komplex_RowMatrix::InvRowSums(Epetra_Vector& x) const {
}

//==========================================================================
int Komplex_RowMatrix::LeftScale(const Epetra_Vector& x) {
}

//==========================================================================
int Komplex_RowMatrix::InvColSums(Epetra_Vector& x) const {
}

//==========================================================================
int Komplex_RowMatrix::RightScale(const Epetra_Vector& x) {
}

//==========================================================================
bool Komplex_RowMatrix::LowerTriangular() const {
}

//==========================================================================
bool Komplex_RowMatrix::UpperTriangular() const {
}

//==========================================================================
double Komplex_RowMatrix::NormInf() const {
}

//==========================================================================
double Komplex_RowMatrix::NormOne() const {
}

//==========================================================================
int Komplex_RowMatrix::NumMyNonzeros() const {
}

//==========================================================================
int Komplex_RowMatrix::NumMyCols() const {
}

//==========================================================================
int Komplex_RowMatrix::NumGlobalRows() const {
}

//==========================================================================
int Komplex_RowMatrix::NumGlobalCols() const {
}

//==========================================================================
int Komplex_RowMatrix::NumGlobalNonzeros() const {
}

//==========================================================================
int Komplex_RowMatrix::NumMyDiagonals() const {
}

//==========================================================================
int Komplex_RowMatrix::NumGlobalDiagonals() const {
}

//==========================================================================
int Komplex_RowMatrix::ExtractGlobalRowCopy(int GlobalRow, 
							  int Length,
							  int& NumEntries,
  							  double* Values,
							  int* Indices) const {
}

//==========================================================================
int Komplex_RowMatrix::ExtractGlobalRowCopy(int GlobalRow,
							  int Length,
							  int& NumEntries,
							  double* Values) const {
}

//==========================================================================
int Komplex_RowMatrix::ExtractMyRowCopy(int MyRow,
						    int Length,
						    int& NumEntries,
						    double* Values,
						    int* Indices) const {
}

//==========================================================================
int Komplex_RowMatrix::NumMyRowEntries(int MyRow, int& NumEntries) const {
}

//==========================================================================
int Komplex_RowMatrix::ExtractDiagonalCopy(Epetra_Vector& Diagonal) const {
}

//==========================================================================
int Komplex_RowMatrix::MaxNumEntries() const {
}

//==========================================================================
const Epetra_Map & Komplex_RowMatrix::RowMatrixRowMap() const {
}

//==========================================================================
const Epetra_Map & Komplex_RowMatrix::RowMatrixColMap() const {
}

//==========================================================================
const Epetra_Import * Komplex_RowMatrix::RowMatrixImporter() const {
}

//==========================================================================
void Komplex_RowMatrix::Print(ostream& os) const {
}

//==========================================================================
const char * Komplex_RowMatrix::Label() const {
}

//==========================================================================
int Komplex_RowMatrix::SetUseTranspose(bool UseTranspose) {
}

//==========================================================================
int Komplex_RowMatrix::Apply(const Epetra_MultiVector& X,
				     Epetra_MultiVector& Y) const {
}

//==========================================================================
int Komplex_RowMatrix::ApplyInverse(const Epetra_MultiVector& X,
 					      Epetra_MultiVector& Y) const {
}

//==========================================================================
bool Komplex_RowMatrix::HasNormInf() const {
}

//==========================================================================
bool Komplex_RowMatrix::UseTranspose() const {
}

//==========================================================================
const Epetra_Map & OperatorDomainMap() const {
}

//==========================================================================
const Epetra_Map & OperatorRangeMap() const {
}