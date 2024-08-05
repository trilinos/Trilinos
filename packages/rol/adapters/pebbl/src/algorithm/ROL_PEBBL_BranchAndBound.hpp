// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PEBBL_BRANCHANDBOUND_HPP
#define ROL_PEBBL_BRANCHANDBOUND_HPP

#include "ROL_PEBBL_Interface.hpp"

namespace ROL {
namespace PEBBL {

template<class Real>
class BranchAndBound {
private:
  // min        obj(x)
  // subject to xl <= x <= xu
  //            econ(x) = 0         (Lagrange Multiplier: emul)
  //            cl <= icon(x) <= cu (Lagrange Multiplier: imul)
  const Ptr<IntegerProblemFactory<Real>> factory_;

  // Parameter list containing algorithmic information
  const Ptr<ParameterList> parlist_;
  
  // Application specific branching helper
  const Ptr<BranchHelper<Real>> bHelper_;

  // PEBBL Information
  Ptr<Branching<Real>> branching_;

public:

  BranchAndBound(const Ptr<Branching<Real>> &branching)
    : branching_(branching) {}

  BranchAndBound(const Ptr<IntegerProblemFactory<Real>> &factory,
                 const Ptr<ParameterList>               &parlist,
                 const Ptr<BranchHelper<Real>>          &bHelper,
                 int                                     verbosity = 0,
                 const Ptr<std::ostream>                &outStream = nullPtr)
    : factory_(factory), parlist_(parlist), bHelper_(bHelper) {
    branching_ = makePtr<Branching<Real>>(factory_,parlist_,bHelper_,verbosity,outStream);
  }

  bool solve(int &argc, char** &argv,
             std::ostream &outStream = std::cout) {
    bool flag = true;
#ifdef HAVE_MPI
    utilib::CommonIO::begin();
    utilib::CommonIO::setIOFlush(1);
  
    utilib::exception_mngr::set_stack_trace(false);
    flag = branching_->setup(argc,argv);
    if (flag) {
      utilib::exception_mngr::set_stack_trace(true);
      branching_->reset();
      branching_->printConfiguration();
      branching_->solve();
    }
    utilib::CommonIO::end();
#else
    utilib::exception_mngr::set_stack_trace(false);
    flag = branching_->setup(argc,argv);
    if (flag) {
      utilib::exception_mngr::set_stack_trace(true);
      branching_->reset();
      branching_->solve();
    }
#endif
    return flag;
  }

  const Ptr<const Vector<Real>> getSolution(void) const {
    if (branching_ != nullPtr) {
      return dynamic_cast<IntegerSolution<Real>*>(branching_->incumbent)->getVector();
    }
    else {
      return nullPtr;
    }
  }
};

} // namespace PEBBL
} // namespace ROL
#endif
