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

#ifndef ROL_PEBBL_STDINTEGERTRANSFORMATION_H
#define ROL_PEBBL_STDINTEGERTRANSFORMATION_H

#include "ROL_PEBBL_IntegerTransformation.hpp"

/** @ingroup func_group
    \class ROL::PEBBL:IntegerStdTransformation
    \brief Defines the pebbl integer transformation operator interface for StdVectors.

    ROL's pebbl constraint interface is designed to set individual components
    of a vector to a fixed value.  The range space is the same as the domain.

    ---
*/


namespace ROL {
namespace PEBBL {

template <class Real>
class StdIntegerTransformation : public IntegerTransformation<Real> {
private:
  Ptr<std::vector<Real>> getData(Vector<Real> &x) const {
    return dynamic_cast<StdVector<Real>&>(x).getVector();
  }

 using IntegerTransformation<Real>::map_; 

public:
  StdIntegerTransformation(void)
    : IntegerTransformation<Real>() {}

  StdIntegerTransformation(const StdIntegerTransformation &T)
    : IntegerTransformation<Real>(T) {}

  void fixValues(Vector<Real> &c, bool zero = false) const {
    for (auto it=map_.begin(); it!=map_.end(); ++it) {
      (*getData(c))[it->first] = (zero ? static_cast<Real>(0) : it->second);
    }
  }

}; // class StdIntegerTransformation

} // namespace PEBBL
} // namespace ROL

#endif
