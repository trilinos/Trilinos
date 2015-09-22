// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_GROWTH_RULES
#define STOKHOS_GROWTH_RULES

namespace Stokhos {

  //! Interface for abstract growth rules
  template <typename value_type>
  class GrowthRule {
  public:

    //! Constructor
    GrowthRule() {}

    //! Destructor
    virtual ~GrowthRule() {}

    //! Evaluate growth rule
    virtual value_type operator() (const value_type& x) const = 0;
  };
  
  //! A growth rule that is the identity
  template <typename value_type>
  class IdentityGrowthRule : public GrowthRule<value_type> {
  public:
    //! Constructor
    IdentityGrowthRule() {}

    //! Destructor
    virtual ~IdentityGrowthRule() {}

    //! Evaluate growth rule
    virtual value_type operator() (const value_type& x) const { return x; }
  };

  //! A linear growth rule
  template <typename value_type>
  class LinearGrowthRule : public GrowthRule<value_type> {
  public:
    //! Constructor
    LinearGrowthRule(const value_type& a_ = value_type(1), 
		     const value_type& b_ = value_type(0)) : 
      a(a_), b(b_) {}

    //! Destructor
    virtual ~LinearGrowthRule() {}

    //! Evaluate growth rule
    virtual value_type operator() (const value_type& x) const { return a*x+b; }

  protected:

    //! Slope
    value_type a;

    //! Offset
    value_type b;
  };

  //! A growth rule that always makes the supplied order even
  /*!
   * When used in conjunction with Gaussian quadrature that generates n+1
   * points for a quadrature of order n, this always results in an odd
   * number of points, and thus includes 0.  This allows some nesting
   * in Gaussian-based sparse grids.
   */
  template <typename value_type>
  class EvenGrowthRule : public GrowthRule<value_type> {
  public:
    //! Constructor
    EvenGrowthRule() {}

    //! Destructor
    virtual ~EvenGrowthRule() {}

    //! Evaluate growth rule
    virtual value_type operator() (const value_type& x) const { 
      if (x % value_type(2) == value_type(1)) return x+value_type(1);
      return x;
    }

  };

  //! An exponential growth rule for Clenshaw-Curtis
  template <typename value_type>
  class ClenshawCurtisExponentialGrowthRule : public GrowthRule<value_type> {
  public:
    //! Constructor
    ClenshawCurtisExponentialGrowthRule() {}

    //! Destructor
    virtual ~ClenshawCurtisExponentialGrowthRule() {}

    //! Evaluate growth rule
    virtual value_type operator() (const value_type& x) const { 
      if (x == value_type(0)) return value_type(0);
      return std::pow(value_type(2),x-value_type(1));
    }

  };

  //! An exponential growth rule for Gauss-Patterson
  template <typename value_type>
  class GaussPattersonExponentialGrowthRule : public GrowthRule<value_type> {
  public:
    //! Constructor
    GaussPattersonExponentialGrowthRule() {}

    //! Destructor
    virtual ~GaussPattersonExponentialGrowthRule() {}

    //! Evaluate growth rule
    virtual value_type operator() (const value_type& x) const { 
      // Gauss-Patterson rules have precision 3*2*l-1, which is odd.
      // Since discrete orthogonality requires integrating polynomials of
      // order 2*p, setting p = 3*2*{l-1}-1 will yield the largest p such that
      // 2*p <= 3*2^l-1
      if (x == value_type(0)) return value_type(0);
      return 3*std::pow(value_type(2),x-value_type(1))-value_type(1);
    }

  };

}

#endif
