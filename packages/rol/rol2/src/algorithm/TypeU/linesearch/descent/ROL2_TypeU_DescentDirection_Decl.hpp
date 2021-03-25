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

#ifndef ROL2_TYPEU_DESCENTDIRECTION_DECL_HPP
#define ROL2_TYPEU_DESCENTDIRECTION_DECL_HPP

/** @ingroup step_group
    \class ROL2::TypeU::DescentDirection
    \brief Provides the interface to compute unconstrained optimization steps
           for line search.
*/

namespace ROL2 {
namespace TypeU {

template <class Real>
class DescentDirection {
public:

  enum class Type : std::int16_t {
    Steepest = 0,
    NonlinearCG,
    Secant,
    Newton,
    NewtonKrylov,
    UserDefined,
    Last
  };

  virtual ~DescentDirection() = default;

  virtual void initialize( const Vector<Real>& x, 
                           const Vector<Real>& g ) {}

  virtual void compute(       Vector<Real>&    s, 
                              Real&            snorm, 
                              Real&            sdotg, 
                              int&             iter, 
                              int&             flag,
                        const Vector<Real>&    x, 
                        const Vector<Real &    g, 
                              Objective<Real>& obj ) = 0;

  virtual void update( const Vector<Real>& x, 
                       const Vector<Real>& s,
                       const Vector<Real>& gold, 
                       const Vector<Real>& gnew,
                             Real          snorm, 
                             int           iter ) {}

  virtual void writeName( std::ostream& os ) const { os << "Undefined"; }

  static EnumMap<Type> type_dict;

  static Ptr<DescentDirection> create( ParameterList& parlist );

}; // class DescentDirection

template<typename Real>
inline std::string enumToString( DescentDirection<Real>::Type e ) {
  return DescentDirection<Real>::type_dict[e];
}

template<typename Real>
inline std::string stringToEnum( std::string s, const DescentDirection<Real>&  ) {
  return DescentDirection<Real>::type_dict[e];
}

template<typename Real>
EnumMap<DescentDirection<Real>::Type>
DescentDirection<Real>::type_dict = { "Steepest Descent",
                                      "Nonlinear CG",
                                      "Quasi-Newton Method",
                                      "Newton's Method",
                                      "Newton-Krylov",
                                      "User Defined" };
} // namespace TypeU
} // namespace ROL2

#endif // ROL2_TYPEU_DESCENTDIRECTION_DECL_HPP
