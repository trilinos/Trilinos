// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_PSEUDOSPECTRAL_ORTHOG_POLY_EXPANSION_HPP
#define STOKHOS_PSEUDOSPECTRAL_ORTHOG_POLY_EXPANSION_HPP

#include "Stokhos_OrthogPolyExpansionBase.hpp"
#include "Stokhos_PseudoSpectralOperator.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_SerialDenseVector.hpp"

namespace Stokhos {

  //! Orthogonal polynomial expansions based on numerical quadrature
  template <typename ordinal_type, typename value_type, 
	    typename point_compare_type = 
	    typename DefaultPointCompare<ordinal_type,value_type>::type,
	    typename node_type = Stokhos::StandardStorage<ordinal_type, 
							  value_type> > 
  class PseudoSpectralOrthogPolyExpansion : 
    public OrthogPolyExpansionBase<ordinal_type, value_type, node_type> {
  public:

    //! Constructor
    PseudoSpectralOrthogPolyExpansion(
    const Teuchos::RCP<const OrthogPolyBasis<ordinal_type, value_type> >& basis,
    const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk,
    const Teuchos::RCP<const PseudoSpectralOperator<ordinal_type, value_type, point_compare_type> >& ps_op,
    const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    //! Destructor
    virtual ~PseudoSpectralOrthogPolyExpansion() {}

    void timesEqual(
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& x);
    void divideEqual(
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& x);

    void timesEqual(
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const OrthogPolyApprox<ordinal_type, value_type, node_type>& x);
    void divideEqual(
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const OrthogPolyApprox<ordinal_type, value_type, node_type>& x);

    
    void times(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void times(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const value_type& a, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void times(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
               const value_type& b);
    void divide(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
                const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
                const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void divide(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
                const value_type& a, 
                const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void divide(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
                const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
                const value_type& b);

    void exp(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void log(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void log10(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void sqrt(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void pow(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void pow(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
             const value_type& a, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void pow(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
             const value_type& b);
    void cos(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void sin(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void tan(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void cosh(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void sinh(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void tanh(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void acos(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void asin(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void atan(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void atan2(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void atan2(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const value_type& a, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void atan2(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
               const value_type& b);
    void acosh(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void asinh(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void atanh(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);

    template <typename FuncT>
    void nary_op(const FuncT& func,
		 OrthogPolyApprox<ordinal_type, value_type, node_type>& c,
		 const OrthogPolyApprox<ordinal_type, value_type, node_type>** a);

    template <typename ExprT1, typename ExprT2>
    value_type compute_times_coeff(ordinal_type k, const ExprT1& a, 
				   const ExprT2& b) const;

    template <typename ExprT1, typename ExprT2>
    value_type fast_compute_times_coeff(ordinal_type k, const ExprT1& a, 
				   const ExprT2& b) const;

  private:

    // Prohibit copying
    PseudoSpectralOrthogPolyExpansion(const PseudoSpectralOrthogPolyExpansion&);

    // Prohibit Assignment
    PseudoSpectralOrthogPolyExpansion& operator=(const PseudoSpectralOrthogPolyExpansion& b);

  protected:

    //! Short-hand for Cijk
    typedef typename OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::Cijk_type Cijk_type;

    //! Pseudospectral operator
    Teuchos::RCP<const PseudoSpectralOperator<ordinal_type, value_type, point_compare_type> > ps_op;

    //! Use quadrature for times functions
    bool use_quad_for_times;

    //! Use quadrature for division functions
    bool use_quad_for_division;

    //! Expansions size
    ordinal_type sz;

    //! Number of Quad points
    ordinal_type nqp;

    //! Short-hand for SerialDenseVector
    typedef Teuchos::SerialDenseVector<ordinal_type,value_type> SDV;

    //! Temporary array for values of first argument at quad points
    SDV avals;

    //! Temporary array for values of second argument at quad points
    SDV bvals;

    //! Temporary array for values of n-ary arguments at quad points
    Teuchos::Array< Teuchos::Array< SDV > > navals;

    //! Temporary array for values of operation at quad points
    SDV fvals;

  public:

    //! Nonlinear unary function
    template <typename FuncT>
    void unary_op(
      const FuncT& func,
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);

    //! Nonlinear binary function
    template <typename FuncT>
    void binary_op(
      const FuncT& func,
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);

    //! Nonlinear binary function
    template <typename FuncT>
    void binary_op(
      const FuncT& func,
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& a, 
      const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);

    //! Nonlinear binary function
    template <typename FuncT>
    void binary_op(
      const FuncT& func,
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const value_type& b);

  protected:

    struct times_quad_func { 
      value_type operator() (const value_type& a, const value_type& b) const { 
        return a * b; 
      } 
    };

    struct div_quad_func { 
      value_type operator() (const value_type& a, const value_type& b) const { 
        return a / b; 
      } 
    };

    struct exp_quad_func { 
      value_type operator() (const value_type& a) const { 
        return std::exp(a); 
      } 
    };

    struct log_quad_func { 
      value_type operator() (const value_type& a) const { 
        return std::log(a); 
      } 
    };

    struct log10_quad_func { 
      value_type operator() (const value_type& a) const { 
        return std::log10(a); 
      } 
    };
    
    struct sqrt_quad_func { 
      value_type operator() (const value_type& a) const { 
        return std::sqrt(a); 
      } 
    };

    struct pow_quad_func { 
      value_type operator() (const value_type& a, const value_type& b) const { 
        return std::pow(a,b); 
      } 
    };

    struct cos_quad_func { 
      value_type operator() (const value_type& a) const { 
        return std::cos(a); 
      } 
    };

    struct sin_quad_func { 
      value_type operator() (const value_type& a) const { 
        return std::sin(a); 
      } 
    };

    struct tan_quad_func { 
      value_type operator() (const value_type& a) const { 
        return std::tan(a); 
      } 
    };

    struct cosh_quad_func { 
      value_type operator() (const value_type& a) const { 
        return std::cosh(a); 
      } 
    };

    struct sinh_quad_func { 
      value_type operator() (const value_type& a) const { 
        return std::sinh(a); 
      } 
    };

    struct tanh_quad_func { 
      value_type operator() (const value_type& a) const { 
        return std::tanh(a); 
      } 
    };

    struct acos_quad_func { 
      value_type operator() (const value_type& a) const { 
        return std::acos(a); 
      } 
    };

    struct asin_quad_func { 
      value_type operator() (const value_type& a) const { 
        return std::asin(a); 
      } 
    };

    struct atan_quad_func { 
      value_type operator() (const value_type& a) const { 
        return std::atan(a); 
      } 
    };

    struct atan2_quad_func { 
      value_type operator() (const value_type& a, const value_type& b) const { 
        return std::atan2(a,b); 
      } 
    };

    struct acosh_quad_func { 
      value_type operator() (const value_type & a) const { 
        return std::log(a+std::sqrt(a*a-value_type(1.0))); 
      }
    };

    struct asinh_quad_func { 
      value_type operator() (const value_type& a) const { 
        return std::log(a+std::sqrt(a*a+value_type(1.0))); 
      }
    };

    struct atanh_quad_func { 
      value_type operator() (const value_type& a) const { 
        return 0.5*std::log((value_type(1.0)+a)/(value_type(1.0)-a)); 
      } 
    };
    
  }; // class PseudoSpectralOrthogPolyExpansion

} // namespace Stokhos

#include "Stokhos_PseudoSpectralOrthogPolyExpansionImp.hpp"

#endif // STOKHOS_QUADORTHOGPOLYEXPANSION_HPP
