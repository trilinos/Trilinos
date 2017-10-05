
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

#pragma once

#includes "XROL_VectorTraits.hpp"


namespace XROL {

template<class V>
class BoundConstraint {
public:

  static_assert( implements_elementwise<V>(), 
    "BoundConstraint requires either vector type V to implement elementwise "
    "functions or template specialization of BoundConstraint<V>." );
  
  BoundConstraint( const V& x, 
                   bool isLower = true,
                   magnitude_t<X> scale = 1 );

  BoundConstraint( std::unique_ptr<V> lo, 
                   std::unique_ptr<V> up,
                   magnitude_t<X> scale =1 );

  virtual ~BoundConstraint();

  virtual void update( const V& );

  void project( V& x ) const;
  
  void projectInterior( V& x ) const;

  void pruneUpperActive( V& v, 
                         const V& x, 
                         magnitude_t<V> eps = 0 ); 

  void pruneUpperActive( V& v, 
                         const dual_t<V>& g, 
                         const V& x, 
                         magnitude_t<V> eps = 0 ); 

  void pruneLowerActive( V& v, 
                         const V& x, 
                         magnitude_t<V> eps = 0 ); 

  void pruneLowerActive( V& v, 
                         const dual_t<V>& g, 
                         const V& x, 
                         magnitude_t<V> eps = 0 ); 

  const V& getUpperBound( void ) const;

  const V& getLowerBound( void ) const;
 
  bool isFeasible( const Vector<Real> &v ) const;

  void activateUpper( void );

  void activateLower( void );

  void activate( void );

  void deactivateLower( void );

  void deactivateUpper( void );

  void deactivate( void );

  bool isUpperActivated( void ) const;

  bool isLowerActivated( void ) const;

  void pruneActive( V& v, const dual_t<V>& g, const V& x, magnitude_t<V> eps = 0 ) const;

  void pruneUpperActive( V& v, const dual_t<V>& g, const V& x, magnitude_t<V> eps = 0 ) const;

  void pruneLowerActive( V& v, const dual_t<V>& g, const V& x, magnitude_t<V> eps = 0 ) const;

  void pruneInactive( V& v, const dual_t<V>& g, const V& x, magnitude_t<V> eps = 0 ) const;

  void pruneUpperInactive( V& v, const dual_t<V>& g, const V& x, magnitude_t<V> eps = 0 ) const;

  void pruneLowerInactive( V& v, const dual_t<V>& g, const V& x, magnitude_t<V> eps = 0 ) const;
  
  void computeProjectedStep( V& v, const V& x ) const;
 
  void computeProjectedGradient( dual_t<V>& g, const V& x ) const;

private:

  std::unique_ptr<V> x_lo_;
  std::unique_ptr<V> x_up_;

  magnitude_t<V> scale_;
  magnitude_t<V> min_diff_;

  bool Lactivated_;
  bool Uactivated_;

  static constexpr auto INF_ {std::numeric_limits<magnitude_t<V>>::max()};
  static constexpr auto NINF_{std::numeric_limits<magnitude_t<V>>::lowest()};

};

} // namespace XROL

