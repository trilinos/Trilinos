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

#include "LOCA_LAPACK_AlDataOutput.H"
#include "LOCA_Parameter_Vector.H"
#include "LOCA_Bifurcation_ArcLengthGroup.H"
#include "NOX_LAPACK_Vector.H"

LOCA::LAPACK::AlDataOutput::AlDataOutput(fstream& fs) : file(fs) {}

LOCA::LAPACK::AlDataOutput::~AlDataOutput() {}

LOCA::Abstract::DataOutput&
LOCA::LAPACK::AlDataOutput::operator = (const LOCA::Abstract::DataOutput& source) {
  return operator = (dynamic_cast<const LOCA::LAPACK::AlDataOutput&>(source));
}

LOCA::LAPACK::AlDataOutput&
LOCA::LAPACK::AlDataOutput::operator = (const LOCA::LAPACK::AlDataOutput& source) {
  return *this;
}

void
LOCA::LAPACK::AlDataOutput::saveGroupData(const LOCA::Abstract::Group& grp) {
  saveGroupData(dynamic_cast<const LOCA::Bifurcation::ArcLengthGroup&>(grp));
  return;
}

void
LOCA::LAPACK::AlDataOutput::saveGroupData(const LOCA::Bifurcation::ArcLengthGroup& grp) {
  const LOCA::Bifurcation::ArcLengthVector& alx = 
    dynamic_cast<const LOCA::Bifurcation::ArcLengthVector&>(grp.getX());
  const NOX::LAPACK::Vector& x = 
    dynamic_cast<const NOX::LAPACK::Vector&>(alx.getXVec());
  const LOCA::ParameterVector& p = grp.getParams();
  double alParam = alx.getArcParam();

  for (int i=0; i<x.length(); i++)
    file << x(i) << " ";
  for (int i=0; i<p.length(); i++)
    file << p[i] << " ";
  file << alParam;

  file << endl;

  return;
}
