// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef COMPRESSED_SENSING_HPP
#define COMPRESSED_SENSING_HPP

#include "ROL_LinearConstraint.hpp"
#include "ROL_NonlinearLeastSquaresObjective.hpp"
#include "random_matrix.hpp"
#include "l1_penalty_optimization_problem.hpp"

template<typename Ordinal, typename Real>
class CompressedSensingConstraint : public ROL::LinearConstraint<Real> {
public:
  CompressedSensingConstraint( Ordinal rows, Ordinal cols ) 
  : ROL::LinearConstraint<Real>(create_random_operator<Ordinal,Real>(rows,cols),
                                create_random_vector<Ordinal,Real>(rows)) {}
};

template<typename Ordinal, typename Real>
ROL::Ptr<ROL::Objective<Real>>
make_CompressedSensingObjective( Ordinal rows, Ordinal cols ) {
  using ROL::Ptr;
  using ROL::makePtr;
  Ptr<ROL::Constraint<Real>> con = makePtr<CompressedSensingConstraint<Ordinal,Real>>(rows,cols);
  auto x = ROL::TeuchosVector<Ordinal,Real>(cols);
  auto c = ROL::TeuchosVector<Ordinal,Real>(rows);
  Ptr<ROL::Objective<Real>>  obj = makePtr<ROL::NonlinearLeastSquaresObjective<Real>>(con,x,c);
  return obj;
}


#endif // COMPRESSED_SENSING_HPP

