//@HEADER
/*
************************************************************************

              Komplex: Complex Linear Solver Package 
                Copyright (2002) Sandia Corporation

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

#include "Komplex_MultiVector.hpp"
#include "Epetra_MultiVector.hpp"
#include "Epetra_BlockMap.hpp"

//==============================================================================
Komplex_MultiVector::Komplex_MultiVector(const Epetra_BlockMap & Map, int NumVectors, bool zeroOut = true) {
  //do something
}

//==============================================================================
Komplex_MultiVector::Komplex_MultiVector(const Epetra_BlockMap & MapReal, const Epetra_BlockMap & MapImag, int NumVectors, bool zeroOut = true) {
  //do something
}

//==============================================================================
Komplex_MultiVector::Komplex_MultiVector(const Epetra_BlockMap & Map, const Epetra_MultiVector & Br, const Epetra_MultiVector & Bi) {
  //do something
}

//==============================================================================
Komplex_MultiVector::Komplex_MultiVector(const Epetra_BlockMap & MapReal, const Epetra_BlockMap & MapImag, const Epetra_MultiVector & Br, const Epetra_MultiVector & Bi) {
  //do something
}

//==============================================================================
Komplex_MultiVector::Komplex_MultiVector(const Komplex_MultiVector & Source) {
  //do something
}

//==============================================================================
Komplex_MultiVector::Komplex_MultiVector(Komplex_DataAcces CV, const Komplex_MultiVector & Source, int * Indices, int NumVectors) {
  //do something
}

//==============================================================================
Komplex_MultiVector::Komplex_MultiVector(Komplex_DataAcces CV, const Komplex_MultiVector & Source, int StartIndex, int NumVectors) {
  //do something
}

//==============================================================================
Komplex_MultiVector::~Komplex_MultiVector() {
  //do something
}

//==============================================================================
int Komplex_MultiVector::ReplaceGlobalValue(int GlobalRow, int VectorIndex, double ScalarValueReal, double ScalarValueImag) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::ReplaceGlobalValueReal(int GlobalRow, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::ReplaceGlobalValueImag(int GlobalRow, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::ReplaceGlobalValue(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValueReal, double ScalarValueImaginary) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::ReplaceGlobalValueReal(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::ReplaceGlobalValueImag(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::SumIntoGlobalValue(int GlobalRow, int VectorIndex, double ScalarValueReal, double ScalarValueImag) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::SumIntoGlobalValueReal(int GlobalRow, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::SumIntoGlobalValueImag(int GlobalRow, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::SumIntoGlobalValue(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValueReal, double ScalarValueImag) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::SumIntoGlobalValueReal(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::SumIntoGlobalValueImag(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::ReplaceMyValue(int MyRow, int VectorIndex, double ScalarValueReal, double ScalarValueImag) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::ReplaceMyValueReal(int MyRow, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::ReplaceMyValueImag(int MyRow, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::ReplaceMyValue(int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValueReal, double ScalarValueImag) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::ReplaceMyValueReal(int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::ReplaceMyValueImag(int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::SumIntoMyValue(int MyRow, int VectorIndex, double ScalarValueReal, double ScalarValueImag) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::SumIntoMyValueReal(int MyRow, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::SumIntoMyValueImag(int MyRow, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::SumIntoMyValue(int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValueReal, double ScalarValueImag) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::SumIntoMyValueReal(int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::SumIntoMyValueImag(int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::Scale(double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::Scale(double ScalarValueReal, double ScalarValueImag) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::Scale(double ScalarA, const Komplex_MultiVector & A) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::Scale(double ScalarAReal, double ScalarAImag, const Komplex_MultiVector & A) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::Norm1(double * Result) const {
  //do something
}

//==============================================================================
int Komplex_MultiVector::Norm2(double * Result) const {
  //do something
}

//==============================================================================
int Komplex_MultiVector::NormInf(double * Result) const {
  //do something
}

//==============================================================================
int Komplex_MultiVector::Multiply(char TransA, char TransB, double ScalarAB, const Komplex_MultiVector & A, const Komplex_MultiVector & B, double ScalarThis) {
  //do something
}

//==============================================================================
int Komplex_MultiVector::Multiply(double ScalarAB, const Komplex_MultiVector & A, const Komplex_MultiVector & B, double ScalarThis) {
  //do something
}

//==============================================================================
Komplex_MultiVector & Komplex_MultiVector::operator = (const Komplex_MultiVector & Source) {
  //do something
}

//==============================================================================
double *& Komplex_MultiVector::operator [] (int i) {
  //do something
}

//==============================================================================
double * const & Komplex_MultiVector::operator [] (int i) const {
  //do something
}

//==============================================================================
Komplex_Vector *& Komplex_MultiVector::operator () (int i) {
  //do something
}

//==============================================================================
const Komplex_Vector *& Komplex_MultiVector::operator() (int i) const {
  //do something
}

//==============================================================================
int Komplex_MultiVector::NumVectors() const {
  //do something
}

//==============================================================================
int Komplex_MultiVector::MyLength() const {
  //do something
}

//==============================================================================
int Komplex_MultiVector::GlobalLength() const {
  //do something
}

//==============================================================================
int Komplex_MultiVector::Stride() const {
  //do something
}

//==============================================================================
bool Komplex_MultiVector::ConstantStride() const {
  //do something
}

//==============================================================================
int Komplex_MultiVector::ReplaceMap(const Epetra_BlockMap & Map) {
  //do something
}

//==============================================================================
void Komplex_MultiVector::Print(ostream & os) const {
  //do something
}