// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_QUADORTHOGPOLYEXPANSION_HPP
#define STOKHOS_QUADORTHOGPOLYEXPANSION_HPP

#include "Stokhos_OrthogPolyExpansionBase.hpp"
#include "Stokhos_Quadrature.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_BLAS.hpp"

namespace Stokhos {

  //! Orthogonal polynomial expansions based on numerical quadrature
  template <typename ordinal_type, typename value_type, 
	    typename node_type = Stokhos::StandardStorage<ordinal_type, 
							  value_type> > 
  class QuadOrthogPolyExpansion : 
    public OrthogPolyExpansionBase<ordinal_type, value_type, node_type> {
  public:

    //! Constructor
    QuadOrthogPolyExpansion(
    const Teuchos::RCP<const OrthogPolyBasis<ordinal_type, value_type> >& basis,
    const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk,
    const Teuchos::RCP<const Quadrature<ordinal_type, value_type> >& quad,
    const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    //! Destructor
    virtual ~QuadOrthogPolyExpansion() {}

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
    void cbrt(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
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
    QuadOrthogPolyExpansion(const QuadOrthogPolyExpansion&);

    // Prohibit Assignment
    QuadOrthogPolyExpansion& operator=(const QuadOrthogPolyExpansion& b);

  protected:

    //! Short-hand for Cijk
    typedef typename OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::Cijk_type Cijk_type;

    //! Quadrature routine
    Teuchos::RCP<const Quadrature<ordinal_type, value_type> > quad;

    //! Use quadrature for times functions
    bool use_quad_for_times;

    //! Use quadrature for division functions
    bool use_quad_for_division;

    //! Expansions size
    ordinal_type sz;
    
    //! BLAS wrappers
    Teuchos::BLAS<ordinal_type,value_type> blas;

    //! Array of Quad points
    const Teuchos::Array< Teuchos::Array<value_type> >& quad_points;

    //! Array of Quad weights
    const Teuchos::Array<value_type>& quad_weights;

    //! Values of basis at Quad points
    const Teuchos::Array< Teuchos::Array<value_type> >& quad_values;

    //! Norms of basis vectors
    const Teuchos::Array<value_type>& norms;

    //! Number of Quad points
    ordinal_type nqp;

    //! Temporary array for values of first argument at quad points
    Teuchos::Array<value_type> avals;

    //! Temporary array for values of second argument at quad points
    Teuchos::Array<value_type> bvals;

    //! Temporary array for values of n-ary arguments at quad points
    Teuchos::Array< Teuchos::Array< Teuchos::Array<value_type> > > navals;

    //! Temporary array for values of operation at quad points
    Teuchos::Array<value_type> fvals;

    //! Reshaped quad values into 1D array
    Teuchos::Array<value_type> qv;

    //! Quad values scaled by norms
    Teuchos::Array<value_type> sqv;

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

    struct cbrt_quad_func { 
      value_type operator() (const value_type& a) const { 
        return std::cbrt(a); 
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
    
  }; // class QuadOrthogPolyExpansion

} // namespace Stokhos

#include "Stokhos_QuadOrthogPolyExpansionImp.hpp"

#endif // STOKHOS_QUADORTHOGPOLYEXPANSION_HPP
