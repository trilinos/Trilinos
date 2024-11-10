// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_MultiContinuation_AbstractGroup.H"

void
LOCA::MultiContinuation::AbstractGroup::preProcessContinuationStep(
                 LOCA::Abstract::Iterator::StepStatus /* stepStatus */)
{
}

void
LOCA::MultiContinuation::AbstractGroup::postProcessContinuationStep(
                 LOCA::Abstract::Iterator::StepStatus /* stepStatus */)
{
}

void
LOCA::MultiContinuation::AbstractGroup::projectToDraw(
                          const NOX::Abstract::Vector& x,
                          double *px) const
{
  px[0] = x.norm(NOX::Abstract::Vector::MaxNorm);
}

int
LOCA::MultiContinuation::AbstractGroup::projectToDrawDimension() const
{
  return 1;
}

double
LOCA::MultiContinuation::AbstractGroup::computeScaledDotProduct(
                     const NOX::Abstract::Vector& a,
                     const NOX::Abstract::Vector& b) const
{
  return a.innerProduct(b);
}

void
LOCA::MultiContinuation::AbstractGroup::scaleVector(NOX::Abstract::Vector& /* x */) const
{
}
