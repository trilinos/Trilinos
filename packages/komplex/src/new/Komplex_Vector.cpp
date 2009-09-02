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
