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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef MFLOCA_H
#define MFLOCA_H

#include <LOCA.H>
#include <LOCA_MultiContinuation_ExtendedGroup.H>
#include <LOCA_MultiContinuation_ExtendedVector.H>
typedef LOCA::MultiContinuation::ExtendedVector LMCEV; 

extern "C" {

#include <MFNSpace.h>
#include <MFNRegion.h>
#include <MFImplicitMF.h>
#include <MFNVector.h>

}

// Structure to hold info for each continuation parameter
struct ParamData {
  string name;
  int ID;
  double initialValue;
  double minValue;
  double maxValue;
  double initialStepSize;
  double minStepSize;
  double maxStepSize;
};

struct LOCAData {
  LOCAData(NOX::Solver::Generic& s, 
	   LOCA::MultiContinuation::ExtendedGroup& g, 
	   NOX::Parameter::List& par,
	   NOX::StatusTest::Generic& st,
	   list<ParamData>& conParamData);

  NOX::Solver::Generic& solver;
  LOCA::MultiContinuation::ExtendedGroup& grp;
  NOX::Parameter::List& p;
  NOX::StatusTest::Generic& status;
  list<ParamData>& paramData;
  MFNSpace space;
  int np;
  double maxNonlinearIterations;
  double aggressiveness;
  double radius;
  double maxRadius;
  double minRadius;
  double solutionMax;
};

extern "C" {

MFNSpace MFCreateLOCANSpace(LOCAData*);
MFNRegion MFNRegionCreateLOCA(LOCAData*);
MFImplicitMF MFIMFCreateLOCA(LOCAData*);
MFNVector MFCreateLOCANVectorWithData(LMCEV*);

}

#endif
