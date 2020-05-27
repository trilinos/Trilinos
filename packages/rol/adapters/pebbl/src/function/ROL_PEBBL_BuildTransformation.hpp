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

#ifndef ROL_PEBBL_BUILDTRANSFORMATION_H
#define ROL_PEBBL_BUILDTRANSFORMATION_H

#include "ROL_AffineTransformObjective.hpp"
#include "ROL_AffineTransformConstraint.hpp"
#include "ROL_PEBBL_IntegerTransformation.hpp"

/** @ingroup func_group
    \class ROL::PEBBL::BuildTransformation
    \brief Perform integer transformation for PEBBL.

    ---
*/

namespace ROL {
namespace PEBBL {

template <class Real>
class BuildTransformation {
private:
  const Ptr<IntegerTransformation<Real>> trans_;
  const Ptr<const Vector<Real>>          x_;
  const Ptr<VectorController<Real>>      storage_;

public:
  virtual ~BuildTransformation(void) {}

  BuildTransformation(const Ptr<IntegerTransformation<Real>> &trans,
                      const Ptr<const Vector<Real>>          &x)
    : trans_(trans), x_(x),
      storage_(makePtr<VectorController<Real>>()) {}

  Ptr<Objective<Real>> transform(const Ptr<Objective<Real>> &obj) const {
    return makePtr<AffineTransformObjective<Real>>(obj,trans_,*x_,storage_);
  }

  Ptr<Constraint<Real>> transform(const Ptr<Constraint<Real>> &con) const {
    return makePtr<AffineTransformConstraint<Real>>(con,trans_,*x_,storage_);
  }

}; // class BuildTransformation

} // namespace PEBBL
} // namespace ROL

#endif
