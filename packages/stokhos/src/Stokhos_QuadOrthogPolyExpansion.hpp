// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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

#ifndef STOKHOS_QUADORTHOGPOLYEXPANSION_HPP
#define STOKHOS_QUADORTHOGPOLYEXPANSION_HPP

#include "Stokhos_OrthogPolyExpansion.hpp"
#include "Stokhos_Quadrature.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_BLAS.hpp"

namespace Stokhos {

  //! Orthogonal polynomial expansions based on numerical quadrature
  template <typename ordinal_type, typename value_type> 
  class QuadOrthogPolyExpansion : 
    public OrthogPolyExpansion<ordinal_type, value_type> {
  public:

    //! Constructor
    QuadOrthogPolyExpansion(
    const Teuchos::RCP<const OrthogPolyBasis<ordinal_type, value_type> >& basis,
    const Teuchos::RCP<const Quadrature<ordinal_type, value_type> >& quad,
    bool use_quad_for_times = false);

    //! Destructor
    virtual ~QuadOrthogPolyExpansion() {}

    //! Get expansion size
    ordinal_type size() const { return sz; }

    //! Get basis
    Teuchos::RCP< const OrthogPolyBasis<ordinal_type, value_type> > 
    getBasis() const {return basis; }
 
    // Operations
    void unaryMinus(OrthogPolyApprox<ordinal_type, value_type>& c, 
                    const OrthogPolyApprox<ordinal_type, value_type>& a);

    void plusEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
		   const value_type& x);
    void minusEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
		    const value_type& x);
    void timesEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
		    const value_type& x);
    void divideEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
		     const value_type& x);

    void plusEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
                   const OrthogPolyApprox<ordinal_type, value_type>& x);
    void minusEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
                    const OrthogPolyApprox<ordinal_type, value_type>& x);
    void timesEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
                    const OrthogPolyApprox<ordinal_type, value_type>& x);
    void divideEqual(OrthogPolyApprox<ordinal_type, value_type>& c, 
                     const OrthogPolyApprox<ordinal_type, value_type>& x);

    void plus(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a, 
              const OrthogPolyApprox<ordinal_type, value_type>& b);
    void plus(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const value_type& a, 
              const OrthogPolyApprox<ordinal_type, value_type>& b);
    void plus(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a, 
              const value_type& b);
    void minus(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a,
               const OrthogPolyApprox<ordinal_type, value_type>& b);
    void minus(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const value_type& a, 
               const OrthogPolyApprox<ordinal_type, value_type>& b);
    void minus(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a, 
               const value_type& b);
    void times(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a, 
               const OrthogPolyApprox<ordinal_type, value_type>& b);
    void times(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const value_type& a, 
               const OrthogPolyApprox<ordinal_type, value_type>& b);
    void times(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a, 
               const value_type& b);
    void divide(OrthogPolyApprox<ordinal_type, value_type>& c, 
                const OrthogPolyApprox<ordinal_type, value_type>& a, 
                const OrthogPolyApprox<ordinal_type, value_type>& b);
    void divide(OrthogPolyApprox<ordinal_type, value_type>& c, 
                const value_type& a, 
                const OrthogPolyApprox<ordinal_type, value_type>& b);
    void divide(OrthogPolyApprox<ordinal_type, value_type>& c, 
                const OrthogPolyApprox<ordinal_type, value_type>& a, 
                const value_type& b);

    void exp(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a);
    void log(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a);
    void log10(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a);
    void sqrt(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void pow(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a, 
             const OrthogPolyApprox<ordinal_type, value_type>& b);
    void pow(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const value_type& a, 
             const OrthogPolyApprox<ordinal_type, value_type>& b);
    void pow(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a, 
             const value_type& b);
    void cos(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a);
    void sin(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a);
    void tan(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a);
    void cosh(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void sinh(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void tanh(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void acos(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void asin(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void atan(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void atan2(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a,
               const OrthogPolyApprox<ordinal_type, value_type>& b);
    void atan2(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const value_type& a, 
               const OrthogPolyApprox<ordinal_type, value_type>& b);
    void atan2(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a, 
               const value_type& b);
    void acosh(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a);
    void asinh(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a);
    void atanh(OrthogPolyApprox<ordinal_type, value_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type>& a);
    void abs(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a);
    void fabs(OrthogPolyApprox<ordinal_type, value_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type>& a);
    void max(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a,
             const OrthogPolyApprox<ordinal_type, value_type>& b);
    void max(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const value_type& a, 
             const OrthogPolyApprox<ordinal_type, value_type>& b);
    void max(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a, 
             const value_type& b);
    void min(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a,
             const OrthogPolyApprox<ordinal_type, value_type>& b);
    void min(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const value_type& a, 
             const OrthogPolyApprox<ordinal_type, value_type>& b);
    void min(OrthogPolyApprox<ordinal_type, value_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type>& a, 
             const value_type& b);

  private:

    // Prohibit copying
    QuadOrthogPolyExpansion(const QuadOrthogPolyExpansion&);

    // Prohibit Assignment
    QuadOrthogPolyExpansion& operator=(const QuadOrthogPolyExpansion& b);

  protected:

     //! Basis
    Teuchos::RCP<const OrthogPolyBasis<ordinal_type, value_type> > basis;

    //! Triple-product tensor
    Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk;

    //! Quadrature routine
    Teuchos::RCP<const Quadrature<ordinal_type, value_type> > quad;

    //! Use quadrature for times functions
    bool use_quad_for_times;

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

    //! Temporary array for values of operation at quad points
    Teuchos::Array<value_type> fvals;

    //! Reshaped quad values into 1D array
    Teuchos::Array<double> qv;

    //! Quad values scaled by norms
    Teuchos::Array<double> sqv;

  public:

    //! Nonlinear unary function
    template <typename FuncT>
    void unary_op(const FuncT& func,
                  OrthogPolyApprox<ordinal_type, value_type>& c, 
                  const OrthogPolyApprox<ordinal_type, value_type>& a);

    //! Nonlinear binary function
    template <typename FuncT>
    void binary_op(const FuncT& func,
                   OrthogPolyApprox<ordinal_type, value_type>& c, 
                   const OrthogPolyApprox<ordinal_type, value_type>& a, 
                   const OrthogPolyApprox<ordinal_type, value_type>& b);

    //! Nonlinear binary function
    template <typename FuncT>
    void binary_op(const FuncT& func,
                   OrthogPolyApprox<ordinal_type, value_type>& c, 
                   const value_type& a, 
                   const OrthogPolyApprox<ordinal_type, value_type>& b);

    //! Nonlinear binary function
    template <typename FuncT>
    void binary_op(const FuncT& func,
                   OrthogPolyApprox<ordinal_type, value_type>& c, 
                   const OrthogPolyApprox<ordinal_type, value_type>& a, 
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
    
  }; // class QuadOrthogPolyExpansion

} // namespace Stokhos

#include "Stokhos_QuadOrthogPolyExpansionImp.hpp"

#endif // STOKHOS_QUADORTHOGPOLYEXPANSION_HPP
