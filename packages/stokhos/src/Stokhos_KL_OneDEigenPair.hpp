// @HEADER
// ***********************************************************************
// 
//                     Stokhos Package
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_KL_ONE_D_EIGENPAIR_HPP
#define STOKHOS_KL_ONE_D_EIGENPAIR_HPP

#include <string>
#include <iostream>
#include <cmath>

#include "Teuchos_RCP.hpp"

namespace Stokhos {

  //! Namespace for analytic %KL expansions
  namespace KL {

    //! Abstract base class for one-dimensional eigenfunctions
    template <typename value_type>
    class OneDEigenFunction {
    public:

      //! Constructor
      OneDEigenFunction() {}

      //! Destructor
      virtual ~OneDEigenFunction() {}

      //! Evaluate eigenfunction
      virtual value_type evaluate(const value_type& x) const = 0;

      //! Print eigenfunction
      virtual void print(std::ostream& os) const = 0;

    }; // class OneDEigenFunction

    template <typename value_type>
    std::ostream& operator << (std::ostream& os, 
			       const OneDEigenFunction<value_type>& f) {
      f.print(os);
      return os;
    }

    //! Container for one-dimensional eigenfunction and eigenvalue
    template <typename value_type>
    struct OneDEigenPair {
      Teuchos::RCP< OneDEigenFunction<value_type> > eig_func;
      value_type eig_val;
    }; // struct OneDEigenPair

    //! One-dimensional eigenfunction for exponential covariance function
    /*!
     * Represents an eigenfunction of the form \f$A \sin(\omega (x-(b+a)/2))\f$ 
     * or  \f$A \cos(\omega (x-(b+a)/2))\f$ over the domain \f$[a,b]\f$ where
     * \f[
     *   A = \frac{1}{\sqrt{\frac{b-a}{2} \pm \frac{\sin(\omega(b-a)}{2\omega}}}
     * \f]
     * for \f$\cos\f$, \f$\sin\f$ respectively.
     */
    template <typename value_type>
    class ExponentialOneDEigenFunction : public OneDEigenFunction<value_type> {
    public:

      //! Enum identifying the type of eigenfunction
      enum TYPE { 
	SIN, ///< A*sin(omega*(x-b))
	COS  ///< A*cos(omega*(x-b))
      };

      //! Constructor
      ExponentialOneDEigenFunction(TYPE type_, const value_type& a_, 
				   const value_type& b_,
				   const value_type& omega_,
				   const std::string& dim_name_) :
	type(type_), a((b_-a_)/2.0), b((b_+a_)/2.0), omega(omega_), 
	dim_name(dim_name_) {
	if (type == SIN)
	  A = 1.0/std::sqrt(a - std::sin(2.*omega*a)/(2.*omega));
	else
	  A = 1.0/std::sqrt(a + std::sin(2.*omega*a)/(2.*omega));
      }

      //! Destructor
      ~ExponentialOneDEigenFunction() {}
      
      //! Evaluate eigenfunction
      virtual value_type evaluate(const value_type& x) const {
	if (type == SIN)
	  return A*sin(omega*(x-b));
	return A*cos(omega*(x-b));
      }
      
      //! Print eigenfunction
      virtual void print(std::ostream& os) const {
	os << A << " * ";
	if (type == SIN)
	  os << "sin(";
	else
	  os << "cos(";
	os << omega << " * (" << dim_name << " - " << b << "))";
      }

      //! Return type
      TYPE getType() const { return type; }

      //! Return frequency
      value_type getFrequency() const { return omega; }

      //! Return multiplier
      value_type getMultiplier() const { return A; }

      //! Get shift
      value_type getShift() const { return b; }

    protected:
      
      //! Type of eigenfunction (sin or cos)
      TYPE type;
      
      //! Domain length
      value_type a;

      //! Domain center
      value_type b;

      //! Multiplier for eigenfunction
      value_type A;

      //! Frequency of eigenfunction
      value_type omega;

      //! Dimesion name (e.g., x_1) for printing eigenfunction
      std::string dim_name;
    };

  } // namespace KL

} // namespace Stokhos

#endif // STOKHOS_KL_ONE_D_EIGENPAIR_HPP
