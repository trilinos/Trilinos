//@HEADER
/*
************************************************************************

              Komplex: Complex Linear Solver Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
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