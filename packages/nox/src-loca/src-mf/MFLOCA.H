// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
     const Teuchos::RCP< std::list<ParamData> >& conParamData);

  Teuchos::RCP<LOCA::GlobalData> globalData;
  Teuchos::RCP<LOCA::Parameter::SublistParser> topParams;
  Teuchos::RCP<NOX::Solver::Generic> solver;
  Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy> grp;
  Teuchos::RCP<Teuchos::ParameterList> p;
  Teuchos::RCP<NOX::StatusTest::Generic> status;
  Teuchos::RCP< std::list<ParamData> > paramData;
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
