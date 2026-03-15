// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PDEOPT_PDE_OPTBOUNDSK_HPP
#define ROL_PDEOPT_PDE_OPTBOUNDSK_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "ROL_TpetraBoundConstraint.hpp"
#include "ROL_BoundConstraint_Partitioned.hpp"

template<class Real,
         class LO=Tpetra::Map<>::local_ordinal_type, 
         class GO=Tpetra::Map<>::global_ordinal_type,
         class Node=Tpetra::Map<>::node_type>
class PDE_OptBounds : public ROL::BoundConstraint_Partitioned<Real> {
public:
  PDE_OptBounds(const ROL::Ptr<ROL::TpetraBoundConstraint<Real,LO,GO,Node>> &fieldbnd,
                const ROL::Ptr<ROL::TpetraMultiVector<Real,LO,GO,Node>> &fieldvec,
                const ROL::Ptr<ROL::BoundConstraint<Real>> &parambnd,
                const ROL::Ptr<ROL::StdVector<Real>> &paramvec)
    : ROL::BoundConstraint_Partitioned<Real>({fieldbnd,parambnd},{fieldvec,paramvec}) {}
};


#endif
