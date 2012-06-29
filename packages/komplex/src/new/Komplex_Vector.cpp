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

#include "Komplex_Vector.hpp"
#include "Epetra_Vector.hpp"
#include "Epetra_BlockMap.hpp"
#include "Komplex_DataAccess.hpp"

//==============================================================================
Komplex_Vector::Komplex_Vector(const Epetra_BlockMap & Map, bool zeroOut = true) 
  : Komplex_MultiVector(Map, 1, zeroOut) //Vector is a special case of MultiVector
{
  SetLabel("Komplex::Vector");
}

//==============================================================================
Komplex_Vector::Komplex_Vector(const Epetra_BlockMap & Map, const Epetra_Vector & br, const Epetra_Vector & bi) {
  //do something
}

//==============================================================================
Komplex_Vector::Komplex_Vector(const Komplex_Vector & Source) 
  : Komplex_MultiVector(Source) //Vector is a special case of MultiVector
{
}

//==============================================================================
Komplex_Vector::Komplex_Vector(Komplex_DataAccess CV, const Epetra_BlockMap & Map, const Komplex_MultiVector & Source, int Index) 
  : Komplex_MultiVector(CV, Source, Index, 1) //Vector is a special case of MultiVector
{
  SetLabel("Komplex::Vector");
}

//==============================================================================
Komplex_Vector::~Komplex_Vector() {
  //do something
}

//==============================================================================
int Komplex_Vector::ReplaceGlobalValues(int NumEntries, double * Values, int * Indices) {
  //do something
}

//==============================================================================
int Komplex_Vector::ReplaceMyValues(int NumEntries, double * Values, int * Indices) {
  //do something
}

//==============================================================================
int Komplex_Vector::SumIntoGlobalValues(int NumEntries, double * Values, int * Indices) {
  //do something
}

//==============================================================================
int Komplex_Vector::SumIntoMyValues(int NumEntries, double * Values, int * Indices) {
  //do something
}

//==============================================================================
int Komplex_Vector::ReplaceGlobalValues(int NumEntries, int BlockOffset, double * Values, int * Indices) {
  //do something
}

//==============================================================================
int Komplex_Vector::ReplaceMyValues(int NumEntries, int BlockOffset, double * Values, int * Indices) {
  //do something
}

//==============================================================================
int Komplex_Vector::SumIntoGlobalValues(int NumEntries, int BlockOffset, double * Values, int * Indices) {
  //do something
}

//==============================================================================
int Komplex_Vector::SumIntoMyValues(int NumEntries, int BlockOffset, double * Values, int * Indices) {
  //do something
}

//==============================================================================
int Komplex_Vector::Scale(double ScalarValue) {
  //do something
}

//==============================================================================
int Komplex_Vector::Scale(double ScalarA, const Komplex_Vector & A) {
  //do something
}

//==============================================================================
int Komplex_Vector::Norm1(double & Result) const {
  //do something
}

//==============================================================================
int Komplex_Vector::Norm2(double & Result) const {
  //do something
}

//==============================================================================
int Komplex_Vector::NormInf(double & Result) const {
  //do something
}

//==============================================================================
Komplex_Vector & Komplex_Vector::operator = (const Komplex_Vector & Source) {
  //do something
}

//==============================================================================
double & Komplex_Vector::operator [] (int index) {
  //do something
}

//==============================================================================
const double & Komplex_Vector::operator [] (int index) const {
  //do something
}

//==============================================================================
int Komplex_Vector::Length() const {
  //do something
}

//==============================================================================
int Komplex_Vector::ReplaceMap(const Epetra_BlockMap & map) {
  //do something
}

//==============================================================================
void Komplex_Vector::Print(ostream & os) const {
  //do something
}
