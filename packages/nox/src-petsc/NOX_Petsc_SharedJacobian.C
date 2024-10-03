// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Petsc_SharedJacobian.H" // class definition
#include "petscmat.h"

using namespace NOX::Petsc;

SharedJacobian::SharedJacobian(Mat& j) :
  owner(NULL)
{
  jacobian = &j;
}

SharedJacobian::SharedJacobian(Mat& j, Mat& p)
{
  jacobian = &j;
  prec = &p;
}

SharedJacobian::~SharedJacobian()
{
}

Mat& SharedJacobian::getJacobian(const Group* newowner)
{
  owner = newowner;
  return *jacobian;
}

const Mat& SharedJacobian::getJacobian() const
{
  return *jacobian;
}

bool SharedJacobian::isOwner(const Group* grp) const
{
  return (owner == grp);
}

Mat& SharedJacobian::getPrec(const Group* newowner)
{
  owner = newowner;
  return *prec;
}

const Mat& SharedJacobian::getPrec() const
{
  return *prec;
}


