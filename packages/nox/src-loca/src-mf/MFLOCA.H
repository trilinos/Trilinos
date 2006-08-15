// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef MFLOCA_H
#define MFLOCA_H

#include <LOCA.H>
#include <LOCA_MultiContinuation_AbstractStrategy.H>
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
  LOCAData(
     const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
     const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& top_params,
     const Teuchos::RefCountPtr<NOX::Solver::Generic>& s, 
     const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractStrategy>& g, 
     const Teuchos::RefCountPtr<Teuchos::ParameterList>& par,
     const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& st,
     const Teuchos::RefCountPtr< list<ParamData> >& conParamData);

  Teuchos::RefCountPtr<LOCA::GlobalData> globalData;
  Teuchos::RefCountPtr<LOCA::Parameter::SublistParser> topParams;
  Teuchos::RefCountPtr<NOX::Solver::Generic> solver;
  Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractStrategy> grp;
  Teuchos::RefCountPtr<Teuchos::ParameterList> p;
  Teuchos::RefCountPtr<NOX::StatusTest::Generic> status;
  Teuchos::RefCountPtr< list<ParamData> > paramData;
  MFNSpace space;
  int np;
  double maxNonlinearIterations;
  double aggressiveness;
  double radius;
  double maxRadius;
  double minRadius;
  double solutionMax;
private:
  // Private to prohibit copying
  LOCAData(const LOCAData&);
  LOCAData& operator=(const LOCAData&);
};

struct LOCANVectorData {
  LOCANVectorData(const Teuchos::RefCountPtr<LMCEV>& u) : u_ptr(u) {}
  Teuchos::RefCountPtr<LMCEV> u_ptr;
};

extern "C" {

MFNSpace MFCreateLOCANSpace(LOCAData*);
MFNRegion MFNRegionCreateLOCA(LOCAData*);
MFImplicitMF MFIMFCreateLOCA(LOCAData*);
MFNVector MFCreateLOCANVectorWithData(const Teuchos::RefCountPtr<LMCEV>&);

}

#endif
