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

#ifndef ROL_TEUCHOSTRANSFORM_PEBBL_H
#define ROL_TEUCHOSTRANSFORM_PEBBL_H

#include "ROL_Transform_PEBBL.hpp"
#include "ROL_TeuchosVector.hpp"

/** @ingroup func_group
    \class ROL::TeuchosTransform_PEBBL
    \brief Defines the pebbl transform operator interface for TeuchosVectors.

    ROL's pebbl constraint interface is designed to set individual components
    of a vector to a fixed value.  The range space is the same as the domain.

    ---
*/


namespace ROL {

template <class Ordinal, class Real>
class TeuchosTransform_PEBBL : public Transform_PEBBL<Real> {
private:
  Ptr<Teuchos::SerialDenseVector<Ordinal,Real>> getData(Vector<Real> &x) const {
    return dynamic_cast<TeuchosVector<Ordinal,Real>&>(x).getVector();
  }

  Ptr<const Teuchos::SerialDenseVector<Ordinal,Real>> getConstData(const Vector<Real> &x) const {
    return dynamic_cast<const TeuchosVector<Ordinal,Real>&>(x).getVector();
  }

 using Transform_PEBBL<Real>::map_; 

public:
  TeuchosTransform_PEBBL(void)
    : Transform_PEBBL<Real>() {}

  TeuchosTransform_PEBBL(const TeuchosTransform_PEBBL &T)
    : Transform_PEBBL<Real>(T) {}

  void pruneVector(Vector<Real> &c) {
    Ptr<Teuchos::SerialDenseVector<Ordinal,Real>> cval = getData(c);
    typename std::map<int,Real>::iterator it;
    for (it=map_.begin(); it!=map_.end(); ++it) {
      (*cval)(it->first) = static_cast<Real>(0);
    }
  }

  void shiftVector(Vector<Real> &c) {
    Ptr<Teuchos::SerialDenseVector<Ordinal,Real>> cval = getData(c);
    typename std::map<int,Real>::iterator it;
    for (it=map_.begin(); it!=map_.end(); ++it) {
      (*cval)(it->first) = it->second;
    }
  }

}; // class TeuchosTransform_PEBBL

} // namespace ROL

#endif
