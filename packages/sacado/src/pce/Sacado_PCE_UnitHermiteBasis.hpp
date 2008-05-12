// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_PCE_UNITHERMITEBASIS_HPP
#define SACADO_PCE_UNITHERMITEBASIS_HPP

#include <vector>
#include <ostream>
#include "Sacado_PCE_StandardPoly.hpp"

namespace Sacado {

  namespace PCE {

    template <typename T>
    class UnitHermiteBasis {
    public:

      //! Typename of values
      typedef T value_type;

      //! Constructor
      UnitHermiteBasis(unsigned int degree);
      
      //! Copy constructor
      UnitHermiteBasis(const UnitHermiteBasis& b);

      //! Destructor
      ~UnitHermiteBasis();

      //! Assignment
      UnitHermiteBasis& operator=(const UnitHermiteBasis& b);

      //! Return size of basis
      unsigned int size() const;

      //! Compute norm squared of each basis element
      const std::vector<T>& norm_squared() const;

      //! Get coefficient of derivative
      T derivCoeff(unsigned int i) const;

      //! Project a polynomial into this basis
      void project(const StandardPoly<T>& poly, std::vector<T>& coeffs) const;

      //! Write polynomial in standard basis
      StandardPoly<T> toStandardBasis(const T coeffs[], 
				      unsigned int n) const;

      //! Get basis polynomial
      const StandardPoly<T>& getBasisPoly(unsigned int i) const;

      void print(std::ostream& os) const;

    protected:

      //! Degree of basis
      unsigned int d;

      //! Basis polynomials
      std::vector< StandardPoly<T> > basis;

      //! Norms
      std::vector<T> norms;

    }; // class UnitHermiteBasis

    template <typename T> 
    std::ostream& operator << (std::ostream& os, const UnitHermiteBasis<T>& b){
      b.print(os);
      return os;
    }

  } // Namespace PCE

} // Namespace Sacado

// Include template definitions
#include "Sacado_PCE_UnitHermiteBasisImp.hpp"

#endif
