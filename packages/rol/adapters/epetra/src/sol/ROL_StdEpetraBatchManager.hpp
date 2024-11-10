// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STDEPETRABATCHMANAGER_HPP
#define ROL_STDEPETRABATCHMANAGER_HPP

#include "ROL_EpetraBatchManager.hpp"
#include "ROL_StdVector.hpp"

namespace ROL {

template<class Real> 
class StdEpetraBatchManager : public EpetraBatchManager<Real> {
public:
  StdEpetraBatchManager(ROL::Ptr<Epetra_Comm> &comm)
    : EpetraBatchManager<Real>(comm) {}

  using EpetraBatchManager<Real>::sumAll;

  void sumAll(Vector<Real> &input, Vector<Real> &output) {
    ROL::Ptr<std::vector<Real> > input_ptr
      = dynamic_cast<StdVector<Real>&>(input).getVector();
    ROL::Ptr<std::vector<Real> > output_ptr
      = dynamic_cast<StdVector<Real>&>(output).getVector();
    int dim_i = static_cast<int>(input_ptr->size());
    int dim_o = static_cast<int>(output_ptr->size());
    ROL_TEST_FOR_EXCEPTION(dim_i != dim_o, std::invalid_argument,
      ">>> (ROL::StdEpetraBatchManager::SumAll): Dimension mismatch!");
    EpetraBatchManager<Real>::sumAll(&input_ptr->front(),
                                     &output_ptr->front(),
                                     dim_i);
  }
};

}

#endif
