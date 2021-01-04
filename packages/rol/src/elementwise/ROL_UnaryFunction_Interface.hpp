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

#ifndef ROL_UNARYFUNCTIONS_INTERFACE_H
#define ROL_UNARYFUNCTIONS_INTERFACE_H

namespace ROL {

namespace Elementwise {

// Forward Declarations of UnaryFunction derived types
template<typename> class AbsoluteValue;
template<typename> class Fill;
template<typename> class Heaviside;
template<typename> class Logarithm;
template<typename> class Power;
template<typename> class Reciprocal;
template<typename> class Round;
template<typename> class Scale;
template<typename> class Shift;
template<typename> class Sign;
template<typename> class SquareRoot;
template<typename> class ThresholdUpper;
template<typename> class ThresholdLower;


/** \class ROL::Elementwise::UnaryFunction
 *  \brief Interface class with function of a single argument
 * 
 */
template<class Real>
class UnaryFunction {
public:

  virtual ~UnaryFunction() {}
  virtual Real apply( const Real &x ) const = 0;

  struct Visitor {
    virtual void visit( const AbsoluteValue<Real>& ) {}
    virtual void visit( const Fill<Real>& ) {}
    virtual void visit( const Heaviside<Real>& ) {}
    virtual void visit( const Logarithm<Real>& ) {}
    virtual void visit( const Power<Real>& ) {}
    virtual void visit( const Reciprocal<Real>& ) {}
    virtual void visit( const Round<Real>& ) {}
    virtual void visit( const Scale<Real>& ) {}
    virtual void visit( const Shift<Real>& ) {}
    virtual void visit( const Sign<Real>& ) {}
    virtual void visit( const SquareRoot<Real>& ) {}
    virtual void visit( const ThresholdUpper<Real>& ) {}
    virtual void visit( const ThresholdLower<Real>& ) {}
  };

  virtual void accept( UnaryFunction::Visitor& ) const {}

};

} // namespace Elementwise 
} // namespace ROL

#endif // ROL_UNARYFUNCTIONS_INTERFACE_HPP

