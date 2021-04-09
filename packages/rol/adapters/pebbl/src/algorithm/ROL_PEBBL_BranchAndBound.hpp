// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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
