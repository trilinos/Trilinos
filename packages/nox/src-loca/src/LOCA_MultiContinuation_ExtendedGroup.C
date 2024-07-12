// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Parameter_Vector.H"
#include "LOCA_MultiContinuation_AbstractGroup.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"
#include "LOCA_MultiContinuation_ConstrainedGroup.H"
#include "LOCA_MultiPredictor_AbstractStrategy.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_MultiContinuation_CompositeConstraint.H"
#include "LOCA_MultiContinuation_CompositeConstraintMVDX.H"

LOCA::MultiContinuation::ExtendedGroup::ExtendedGroup(
             const LOCA::MultiContinuation::ExtendedGroup& source,
             NOX::CopyType type)
  : globalData(source.globalData),
    parsedParams(source.parsedParams),
    continuationParams(source.continuationParams),
    grpPtr(),
    predictor(),
    conGroup(),
    numParams(source.numParams),
    tangentMultiVec(source.tangentMultiVec, type),
    scaledTangentMultiVec(source.scaledTangentMultiVec, type),
    prevXVec(source.prevXVec, type),
    conParamIDs(source.conParamIDs),
    stepSize(source.stepSize),
    stepSizeScaleFactor(source.stepSizeScaleFactor),
    isValidPredictor(false),
    baseOnSecant(source.baseOnSecant)
{
  predictor = source.predictor->clone(type);
  conGroup = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ConstrainedGroup>(source.conGroup->clone(type));
  grpPtr = conGroup->getGroup();
  if (source.isValidPredictor && type == NOX::DeepCopy)
    isValidPredictor = true;
}


LOCA::MultiContinuation::ExtendedGroup::~ExtendedGroup()
{
}

NOX::Abstract::Group&
LOCA::MultiContinuation::ExtendedGroup::operator=(
                      const NOX::Abstract::Group& source)
{
  copy(source);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::MultiContinuation::ExtendedGroup::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new ExtendedGroup(*this, type));
}

void
LOCA::MultiContinuation::ExtendedGroup::setX(const NOX::Abstract::Vector& y)
{
  conGroup->setX(y);
}

void
LOCA::MultiContinuation::ExtendedGroup::computeX(
                          const NOX::Abstract::Group& g,
                          const NOX::Abstract::Vector& d,
                          double step)
{
  const LOCA::MultiContinuation::ExtendedGroup& mg =
    dynamic_cast<const LOCA::MultiContinuation::ExtendedGroup&>(g);

  conGroup->computeX(*(mg.conGroup), d, step);
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::computeF()
{
  return conGroup->computeF();
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::computeJacobian()
{
  return conGroup->computeJacobian();
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::computeGradient()
{
  return conGroup->computeGradient();
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::computeNewton(
                           Teuchos::ParameterList& params)
{
  return conGroup->computeNewton(params);
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobian(
                      const NOX::Abstract::Vector& input,
                      NOX::Abstract::Vector& result) const
{
  return conGroup->applyJacobian(input, result);
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianTranspose(
                      const NOX::Abstract::Vector& input,
                      NOX::Abstract::Vector& result) const
{
  return conGroup->applyJacobianTranspose(input, result);
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianInverse(
                      Teuchos::ParameterList& params,
                      const NOX::Abstract::Vector& input,
                      NOX::Abstract::Vector& result) const
{
  return conGroup->applyJacobianInverse(params, input, result);
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianMultiVector(
                     const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& result) const
{
  return conGroup->applyJacobianMultiVector(input, result);
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianTransposeMultiVector(
                     const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& result) const
{
  return conGroup->applyJacobianTransposeMultiVector(input, result);
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianInverseMultiVector(
                     Teuchos::ParameterList& params,
                     const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& result) const
{
  return conGroup->applyJacobianInverseMultiVector(params, input, result);
}

bool
LOCA::MultiContinuation::ExtendedGroup::isF() const
{
  return conGroup->isF();
}

bool
LOCA::MultiContinuation::ExtendedGroup::isJacobian() const
{
  return conGroup->isJacobian();
}

bool
LOCA::MultiContinuation::ExtendedGroup::isGradient() const
{
  return conGroup->isGradient();
}

bool
LOCA::MultiContinuation::ExtendedGroup::isNewton() const
{
  return conGroup->isNewton();
}

const NOX::Abstract::Vector&
LOCA::MultiContinuation::ExtendedGroup::getX() const
{
  return conGroup->getX();
}

const NOX::Abstract::Vector&
LOCA::MultiContinuation::ExtendedGroup::getF() const
{
  return conGroup->getF();
}

double
LOCA::MultiContinuation::ExtendedGroup::getNormF() const
{
  return conGroup->getNormF();
}

const NOX::Abstract::Vector&
LOCA::MultiContinuation::ExtendedGroup::getGradient() const
{
  return conGroup->getGradient();
}

const NOX::Abstract::Vector&
LOCA::MultiContinuation::ExtendedGroup::getNewton() const
{
  return conGroup->getNewton();
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::MultiContinuation::ExtendedGroup::getXPtr() const
{
  return conGroup->getXPtr();
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::MultiContinuation::ExtendedGroup::getFPtr() const
{
  return conGroup->getFPtr();
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::MultiContinuation::ExtendedGroup::getGradientPtr() const
{
  return conGroup->getGradientPtr();
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::MultiContinuation::ExtendedGroup::getNewtonPtr() const
{
  return conGroup->getNewtonPtr();
}

double
LOCA::MultiContinuation::ExtendedGroup::getNormNewtonSolveResidual() const
{
  return conGroup->getNormNewtonSolveResidual();
}

Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
LOCA::MultiContinuation::ExtendedGroup::getUnderlyingGroup() const
{
  return conGroup->getUnderlyingGroup();
}

Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::MultiContinuation::ExtendedGroup::getUnderlyingGroup()
{
  return conGroup->getUnderlyingGroup();
}

void
LOCA::MultiContinuation::ExtendedGroup::copy(const NOX::Abstract::Group& src)
{

  const LOCA::MultiContinuation::ExtendedGroup& source =
    dynamic_cast<const LOCA::MultiContinuation::ExtendedGroup&>(src);

  // Protect against A = A
  if (this != &source) {
    globalData = source.globalData;
    parsedParams = source.parsedParams;
    continuationParams = source.continuationParams;
    *predictor = *source.predictor;
    conGroup->copy(*source.conGroup);
    grpPtr = conGroup->getGroup();
    numParams = source.numParams;
    tangentMultiVec = source.tangentMultiVec;
    scaledTangentMultiVec = source.scaledTangentMultiVec;
    prevXVec = source.prevXVec;
    conParamIDs = source.conParamIDs;
    stepSize = source.stepSize;
    stepSizeScaleFactor = source.stepSizeScaleFactor;
    isValidPredictor = source.isValidPredictor;
    baseOnSecant = source.baseOnSecant;
  }
}

int
LOCA::MultiContinuation::ExtendedGroup::getNumParams() const
{
  return numParams;
}

void
LOCA::MultiContinuation::ExtendedGroup::preProcessContinuationStep(
                 LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  conGroup->preProcessContinuationStep(stepStatus);
}

void
LOCA::MultiContinuation::ExtendedGroup::postProcessContinuationStep(
                 LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  conGroup->postProcessContinuationStep(stepStatus);
  if (stepStatus == LOCA::Abstract::Iterator::Successful) {
    isValidPredictor = false;
    baseOnSecant = true;
  }
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::computePredictor()
{
  if (isValidPredictor)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction =
    "LOCA::MultiContinuation::ExtendedGroup::computePredictor()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute predictor
  status = predictor->compute(baseOnSecant, stepSize, *this, prevXVec,
       dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(conGroup->
                                    getX()));
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  // Fill tangent vector
  status = predictor->computeTangent(tangentMultiVec);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  scaleTangent();

  isValidPredictor = true;
  return finalStatus;
}

bool
LOCA::MultiContinuation::ExtendedGroup::isPredictor() const
{
  return isValidPredictor;
}

void
LOCA::MultiContinuation::ExtendedGroup::scaleTangent()
{
  scaledTangentMultiVec = tangentMultiVec;

  // Only scale the tangent if it is scalable
  if (predictor->isTangentScalable()) {

    for (int i=0; i<numParams; i++) {
      LOCA::MultiContinuation::ExtendedVector & v =
        dynamic_cast<LOCA::MultiContinuation::ExtendedVector&>(scaledTangentMultiVec[i]);
      grpPtr->scaleVector(*(v.getXVec()));
      grpPtr->scaleVector(*(v.getXVec()));
    }

  }
}

void
LOCA::MultiContinuation::ExtendedGroup::setPredictorTangentDirection(
             const LOCA::MultiContinuation::ExtendedVector& v,
             int i)
{
  tangentMultiVec[i] = v;
}

const LOCA::MultiContinuation::ExtendedMultiVector&
LOCA::MultiContinuation::ExtendedGroup::getPredictorTangent() const
{
  return tangentMultiVec;
}


const LOCA::MultiContinuation::ExtendedMultiVector&
LOCA::MultiContinuation::ExtendedGroup::getScaledPredictorTangent() const
{
  return scaledTangentMultiVec;
}

void
LOCA::MultiContinuation::ExtendedGroup::setPrevX(
                         const NOX::Abstract::Vector& y)
{
  prevXVec = y;
}

const LOCA::MultiContinuation::ExtendedVector&
LOCA::MultiContinuation::ExtendedGroup::getPrevX() const
{
  return prevXVec;
}

void
LOCA::MultiContinuation::ExtendedGroup::setStepSize(double deltaS, int i)
{
  stepSize[i] = deltaS;
}

double
LOCA::MultiContinuation::ExtendedGroup::getStepSize(int i) const
{
  return stepSize[i];
}

void
LOCA::MultiContinuation::ExtendedGroup::setContinuationParameter(double val,
                                 int i)
{
  conGroup->setConstraintParameter(i, val);
}

double
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameter(int i) const
{
  return conGroup->getConstraintParameter(i);
}

int
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameterID(int i) const
{
  return conParamIDs[i];
}

const std::vector<int>&
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameterIDs() const
{
  return conParamIDs;
}

std::string
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameterName(
                                  int i) const
{
  const LOCA::ParameterVector& p = grpPtr->getParams();
  return p.getLabel(conParamIDs[i]);
}

double
LOCA::MultiContinuation::ExtendedGroup::getStepSizeScaleFactor(int i) const
{
  return stepSizeScaleFactor[i];
}

void
LOCA::MultiContinuation::ExtendedGroup::printSolution() const
{
  for (int i=0; i<numParams; i++)
    grpPtr->printSolution(getContinuationParameter(i));
}

double
LOCA::MultiContinuation::ExtendedGroup::computeScaledDotProduct(
             const NOX::Abstract::Vector& x,
             const NOX::Abstract::Vector& y) const
{
  const LOCA::MultiContinuation::ExtendedVector& mx =
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(x);
  const LOCA::MultiContinuation::ExtendedVector& my =
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(y);

  double val = grpPtr->computeScaledDotProduct(*mx.getXVec(), *my.getXVec());
  for (int i=0; i<numParams; i++)
    val += mx.getScalar(i) * my.getScalar(i);

  return val;
}

int
LOCA::MultiContinuation::ExtendedGroup::projectToDrawDimension() const
{
  return numParams + grpPtr->projectToDrawDimension();
}

void
LOCA::MultiContinuation::ExtendedGroup::projectToDraw(
                const LOCA::MultiContinuation::ExtendedVector& x,
                double *px) const
{
  // first numParams components are the parameters
  for (int i=0; i<numParams; i++)
    px[i] = x.getScalar(i);

  // fill remaining solution components
  grpPtr->projectToDraw(*x.getXVec(), px+numParams);
}

int
LOCA::MultiContinuation::ExtendedGroup::getBorderedWidth() const
{
  return conGroup->getBorderedWidth();
}

Teuchos::RCP<const NOX::Abstract::Group>
LOCA::MultiContinuation::ExtendedGroup::getUnborderedGroup() const
{
  return conGroup->getUnborderedGroup();
}

bool
LOCA::MultiContinuation::ExtendedGroup::isCombinedAZero() const
{
  return conGroup->isCombinedAZero();
}

bool
LOCA::MultiContinuation::ExtendedGroup::isCombinedBZero() const
{
  return conGroup->isCombinedBZero();
}

bool
LOCA::MultiContinuation::ExtendedGroup::isCombinedCZero() const
{
  return conGroup->isCombinedCZero();
}

void
LOCA::MultiContinuation::ExtendedGroup::extractSolutionComponent(
                            const NOX::Abstract::MultiVector& v,
                                        NOX::Abstract::MultiVector& v_x) const
{
  conGroup->extractSolutionComponent(v, v_x);
}

void
LOCA::MultiContinuation::ExtendedGroup::extractParameterComponent(
               bool use_transpose,
                           const NOX::Abstract::MultiVector& v,
                           NOX::Abstract::MultiVector::DenseMatrix& v_p) const
{
  conGroup->extractParameterComponent(use_transpose, v, v_p);
}

void
LOCA::MultiContinuation::ExtendedGroup::loadNestedComponents(
               const NOX::Abstract::MultiVector& v_x,
               const NOX::Abstract::MultiVector::DenseMatrix& v_p,
               NOX::Abstract::MultiVector& v) const
{
  conGroup->loadNestedComponents(v_x, v_p, v);
}

void
LOCA::MultiContinuation::ExtendedGroup::fillA(
                                     NOX::Abstract::MultiVector& A) const
{
  conGroup->fillA(A);
}

void
LOCA::MultiContinuation::ExtendedGroup::fillB(
                                     NOX::Abstract::MultiVector& B) const
{
  conGroup->fillB(B);
}

void
LOCA::MultiContinuation::ExtendedGroup::fillC(
                         NOX::Abstract::MultiVector::DenseMatrix& C) const
{
  conGroup->fillC(C);
}

LOCA::MultiContinuation::ExtendedGroup::ExtendedGroup(
      const Teuchos::RCP<LOCA::GlobalData>& global_data,
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& conParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const std::vector<int>& paramIDs)
  : globalData(global_data),
    parsedParams(topParams),
    continuationParams(conParams),
    grpPtr(grp),
    predictor(pred),
    conGroup(),
    numParams(paramIDs.size()),
    tangentMultiVec(globalData, grp->getX(), numParams, numParams,
            NOX::ShapeCopy),
    scaledTangentMultiVec(globalData, grp->getX(), numParams, numParams,
              NOX::ShapeCopy),
    prevXVec(globalData, grp->getX(), numParams),
    conParamIDs(paramIDs),
    stepSize(numParams, 0.0),
    stepSizeScaleFactor(numParams, 1.0),
    isValidPredictor(false),
    baseOnSecant(false)
{
}

void
LOCA::MultiContinuation::ExtendedGroup::setConstraints(const Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>& constraints, bool skip_dfdp)
{
  // Form constrained group using original group and continuation constraints
  conGroup = Teuchos::rcp(new ConstrainedGroup(globalData, parsedParams,
                           continuationParams,
                           grpPtr, constraints,
                           conParamIDs,
                           skip_dfdp));
  grpPtr = conGroup->getGroup();
}
