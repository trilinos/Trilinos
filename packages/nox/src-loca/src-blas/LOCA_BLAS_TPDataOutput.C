// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "LOCA_BLAS_TPDataOutput.H"
#include "LOCA_Parameter_Vector.H"
#include "LOCA_Bifurcation_TPBordGroup.H"
#include "NOX_BLAS_Vector.H"

LOCA::BLAS::TPDataOutput::TPDataOutput(fstream& fs) : file(fs) {}

LOCA::BLAS::TPDataOutput::~TPDataOutput() {}

LOCA::Abstract::DataOutput&
LOCA::BLAS::TPDataOutput::operator = (const LOCA::Abstract::DataOutput& source) {
  return operator = (dynamic_cast<const LOCA::BLAS::TPDataOutput&>(source));
}

LOCA::BLAS::TPDataOutput&
LOCA::BLAS::TPDataOutput::operator = (const LOCA::BLAS::TPDataOutput& source) {
  return *this;
}

void
LOCA::BLAS::TPDataOutput::saveGroupData(const LOCA::Abstract::Group& grp) {
  saveGroupData(dynamic_cast<const LOCA::Bifurcation::TPBordGroup&>(grp));
  return;
}

void
LOCA::BLAS::TPDataOutput::saveGroupData(const LOCA::Bifurcation::TPBordGroup& grp) {
  const LOCA::Bifurcation::TPBordVector& tpx = 
    dynamic_cast<const LOCA::Bifurcation::TPBordVector&>(grp.getX());
  const NOX::BLAS::Vector& x = 
    dynamic_cast<const NOX::BLAS::Vector&>(tpx.getXVec());
  const NOX::BLAS::Vector& y = 
    dynamic_cast<const NOX::BLAS::Vector&>(tpx.getNullVec());
  const LOCA::ParameterVector& p = grp.getParams();

  for (int i=0; i<x.length(); i++)
    file << x(i) << " ";
  for (int i=0; i<y.length(); i++)
    file << y(i) << " ";
  for (int i=0; i<p.length(); i++)
    file << p[i] << " ";

  file << endl;

  return;
}
