// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_MultiContinuation_ArcLengthGroup.H"
#include "LOCA_MultiContinuation_ArcLengthConstraint.H"
#include "LOCA_MultiContinuation_AbstractGroup.H"
#include "LOCA_MultiContinuation_ConstrainedGroup.H"
#include "LOCA_MultiPredictor_AbstractStrategy.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"

LOCA::MultiContinuation::ArcLengthGroup::ArcLengthGroup(
      const Teuchos::RCP<LOCA::GlobalData>& global_data,
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& _continuationParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const std::vector<int>& paramIDs)
  : LOCA::MultiContinuation::ExtendedGroup(global_data, topParams,
                       _continuationParams,
                       grp, pred, paramIDs),
    theta(paramIDs.size(), 1.0),
    doArcLengthScaling(true),
    gGoal(0.5),
    gMax(0.8),
    thetaMin(1.0e-3),
    isFirstRescale(true)
{
  Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface> cons
    = Teuchos::rcp(new LOCA::MultiContinuation::ArcLengthConstraint(
    globalData, Teuchos::rcp(this, false)));
  LOCA::MultiContinuation::ExtendedGroup::setConstraints(cons, false);

  double theta0 =
    continuationParams->get("Initial Scale Factor", 1.0);
  doArcLengthScaling =
    continuationParams->get("Enable Arc Length Scaling",true);
  gGoal =
    continuationParams->get("Goal Arc Length Parameter Contribution",
                     0.5);
  gMax =
    continuationParams->get("Max Arc Length Parameter Contribution",
                     0.8);
  thetaMin = continuationParams->get("Min Scale Factor", 1.0e-3);

  for (int i=0; i<numParams; i++)
    theta[i] = theta0;
}

LOCA::MultiContinuation::ArcLengthGroup::ArcLengthGroup(
             const LOCA::MultiContinuation::ArcLengthGroup& source,
             NOX::CopyType type)
  : LOCA::MultiContinuation::ExtendedGroup(source, type),
    theta(source.theta),
    doArcLengthScaling(source.doArcLengthScaling),
    gGoal(source.gGoal),
    gMax(source.gMax),
    thetaMin(source.thetaMin),
    isFirstRescale(source.isFirstRescale)
{
  Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ArcLengthConstraint>(conGroup->getConstraints())->setArcLengthGroup(Teuchos::rcp(this, false));
}


LOCA::MultiContinuation::ArcLengthGroup::~ArcLengthGroup()
{
}

NOX::Abstract::Group&
LOCA::MultiContinuation::ArcLengthGroup::operator=(
                      const NOX::Abstract::Group& source)
{
  copy(source);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::MultiContinuation::ArcLengthGroup::clone(NOX::CopyType type) const
{
  return
    Teuchos::rcp(new LOCA::MultiContinuation::ArcLengthGroup(*this, type));
}

void
LOCA::MultiContinuation::ArcLengthGroup::copy(const NOX::Abstract::Group& src)
{

  const LOCA::MultiContinuation::ArcLengthGroup& source =
    dynamic_cast<const LOCA::MultiContinuation::ArcLengthGroup&>(src);

  // Protect against A = A
  if (this != &source) {
    LOCA::MultiContinuation::ExtendedGroup::copy(source);
    theta = source.theta;
    doArcLengthScaling = source.doArcLengthScaling;
    gGoal = source.gGoal;
    gMax = source.gMax;
    thetaMin = source.thetaMin;
    isFirstRescale = source.isFirstRescale;
  }
}

void
LOCA::MultiContinuation::ArcLengthGroup::scaleTangent()
{
  double dpdsOld, dpdsNew;
  double thetaOld, thetaNew;

  scaledTangentMultiVec = tangentMultiVec;

  // Only scale the tangent if it is scalable
  if (predictor->isTangentScalable()) {

    for (int i=0; i<numParams; i++) {
      LOCA::MultiContinuation::ExtendedVector & v =
        dynamic_cast<LOCA::MultiContinuation::ExtendedVector&>(tangentMultiVec[i]);
      LOCA::MultiContinuation::ExtendedVector & sv =
        dynamic_cast<LOCA::MultiContinuation::ExtendedVector&>(scaledTangentMultiVec[i]);
      grpPtr->scaleVector(*(sv.getXVec()));
      grpPtr->scaleVector(*(sv.getXVec()));

      if (doArcLengthScaling) {

    // Estimate dpds
    thetaOld = theta[i];
    sv.getScalars()->scale(thetaOld*thetaOld);
    dpdsOld = 1.0/sqrt(sv.innerProduct(v));

    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
      globalData->locaUtils->out() <<
        std::endl << "\t" <<
        globalData->locaUtils->fill(64, '+') << std::endl << "\t" <<
        "Arc-length scaling calculation for parameter " <<
        getContinuationParameterName(i) << ": " << std::endl << "\t" <<
        "Parameter component of predictor before rescaling = " <<
        globalData->locaUtils->sciformat(dpdsOld) << std::endl << "\t" <<
        "Scale factor from previous step " <<
        globalData->locaUtils->sciformat(thetaOld) << std::endl << "\t" <<
        "Parameter contribution to arc-length equation     = " <<
        globalData->locaUtils->sciformat(thetaOld*dpdsOld) << std::endl;
    }

    // Recompute scale factor
    recalculateScaleFactor(dpdsOld, thetaOld, thetaNew);

    sv.getScalars()->scale(thetaNew*thetaNew / (thetaOld*thetaOld));

    // Calculate new dpds using new scale factor
    dpdsNew = 1.0/sqrt(sv.innerProduct(v));

    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
      globalData->locaUtils->out() << std::endl << "\t" <<
        "Parameter component of predictor after rescaling  = " <<
        globalData->locaUtils->sciformat(dpdsNew) << std::endl << "\t" <<
        "New scale factor (theta)                          = " <<
        globalData->locaUtils->sciformat(thetaNew) << std::endl << "\t" <<
        "Parameter contribution to arc-length equation     = " <<
        globalData->locaUtils->sciformat(thetaNew*dpdsNew) << std::endl <<
        "\t" << globalData->locaUtils->fill(64, '+') << std::endl;
    }

    // Rescale predictor vector
    v.scale(dpdsNew);
    sv.scale(dpdsNew);

    theta[i] = thetaNew;

    // Adjust step size scaling factor to reflect changes in
    // arc-length parameterization
    // The first time we rescale (first continuation step) we use a
    // different step size scale factor so that dpds*deltaS = step size
    // provided by user
    if (isFirstRescale) {
      stepSizeScaleFactor[i] = 1.0/dpdsNew;
    }
    else
      stepSizeScaleFactor[i] = dpdsOld/dpdsNew;
      }
    }

    if (doArcLengthScaling && isFirstRescale)
      isFirstRescale = false;

  }
}

double
LOCA::MultiContinuation::ArcLengthGroup::computeScaledDotProduct(
             const NOX::Abstract::Vector& x,
             const NOX::Abstract::Vector& y) const
{
  const LOCA::MultiContinuation::ExtendedVector& mx =
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(x);
  const LOCA::MultiContinuation::ExtendedVector& my =
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(y);

  double val = grpPtr->computeScaledDotProduct(*mx.getXVec(), *my.getXVec());
  for (int i=0; i<numParams; i++)
    val += theta[i] * theta[i] * mx.getScalar(i) * my.getScalar(i);

  return val;
}

void
LOCA::MultiContinuation::ArcLengthGroup::recalculateScaleFactor(
                                double dpds,
                                double thetaOld,
                                double& thetaNew)
{
  double g = dpds*thetaOld;

  if (g > gMax) {
    thetaNew = gGoal/dpds * sqrt( fabs(1.0 - g*g) / fabs(1.0 - gGoal*gGoal) );

    if (thetaNew < thetaMin)
      thetaNew = thetaMin;
  }
  else
    thetaNew = thetaOld;
}


