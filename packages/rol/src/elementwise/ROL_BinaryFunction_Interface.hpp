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

#ifndef ROL_BINARYFUNCTIONS_INTERFACE_H
#define ROL_BINARYFUNCTIONS_INTERFACE_H

namespace ROL {

namespace Elementwise {

// Forward Declarations of BinaryFunction derived types
template<typename> class Axpy;
template<typename> class Aypx;
template<typename> class Divide;
template<typename> class DivideAndInvert;
template<typename> class Greater;
template<typename> class Lesser;
template<typename> class Max;
template<typename> class Min;
template<typename> class Multiply;
template<typename> class Plus;
template<typename> class Set;


// Interface class with function of two arguments
template<class Real>
class BinaryFunction {
public:
  virtual ~BinaryFunction() {}
  virtual Real apply( const Real &x, const Real &y ) const = 0;

  struct Visitor {
    virtual void visit( const Axpy<Real>& ) {}
    virtual void visit( const Aypx<Real>& ) {}
    virtual void visit( const Divide<Real>& ) {}
    virtual void visit( const DivideAndInvert<Real>& ) {}
    virtual void visit( const Greater<Real>& ) {}
    virtual void visit( const Lesser<Real>& ) {}
    virtual void visit( const Max<Real>& ) {}
    virtual void visit( const Min<Real>& ) {}
    virtual void visit( const Multiply<Real>& ) {}
    virtual void visit( const Plus<Real>& ) {}
    virtual void visit( const Set<Real>& ) {}
  };

  virtual void accept( BinaryFunction::Visitor& ) const {}
};

} // namespace Elementwise 
} // namespace ROL

#endif // ROL_BINARYFUNCTIONS_INTERFACE_HPP

