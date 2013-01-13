// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
#include <MFErrorHandler.h>

}

// Structure to hold info for each continuation parameter
struct ParamData {
  std::string name;
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
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& top_params,
     const Teuchos::RCP<NOX::Solver::Generic>& s, 
     const Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy>& g, 
     const Teuchos::RCP<Teuchos::ParameterList>& par,
     const Teuchos::RCP<NOX::StatusTest::Generic>& st,
     const Teuchos::RCP< list<ParamData> >& conParamData);

  Teuchos::RCP<LOCA::GlobalData> globalData;
  Teuchos::RCP<LOCA::Parameter::SublistParser> topParams;
  Teuchos::RCP<NOX::Solver::Generic> solver;
  Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy> grp;
  Teuchos::RCP<Teuchos::ParameterList> p;
  Teuchos::RCP<NOX::StatusTest::Generic> status;
  Teuchos::RCP< list<ParamData> > paramData;
  MFNSpace space;
  int np;
  double maxNonlinearIterations;
  double aggressiveness;
  double radius;
  double maxRadius;
  double minRadius;
  double solutionMax;
  MFErrorHandler mfErrorHandler;
private:
  // Private to prohibit copying
  LOCAData(const LOCAData&);
  LOCAData& operator=(const LOCAData&);
};

struct LOCANVectorData {
  LOCANVectorData(const Teuchos::RCP<LMCEV>& u) : u_ptr(u) {}
  Teuchos::RCP<LMCEV> u_ptr;
};

extern "C" {

MFNSpace MFCreateLOCANSpace(LOCAData*);
MFNRegion MFNRegionCreateLOCA(LOCAData*);
MFImplicitMF MFIMFCreateLOCA(LOCAData*);
MFNVector MFCreateLOCANVectorWithData(const Teuchos::RCP<LMCEV>&,MFErrorHandler);

}

#endif
