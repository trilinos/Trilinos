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

#ifndef STOKHOS_ORTHOGPOLYBASISBASE_HPP
#define STOKHOS_ORTHOGPOLYBASISBASE_HPP

#include <string>
#include <vector>
#include "Stokhos_OrthogPolyBasis.hpp"

namespace Stokhos {

  template <typename T>
  class OrthogPolyBasisBase : public OrthogPolyBasis<T> {
  public:

    //! Typename of values
    typedef typename OrthogPolyBasis<T>::value_type value_type;

    //! Destructor
    ~OrthogPolyBasisBase();

    //! Return order of basis
    virtual unsigned int order() const;

    //! Return dimension of basis
    virtual unsigned int dimension() const;

    //! Return total size of basis
    virtual unsigned int size() const;

    //! Compute norm squared of each basis element
    virtual const std::vector<T>& norm_squared() const;

    //! Project product of two basis polynomials into this basis
    virtual void projectProduct(unsigned int i, unsigned int j,
				std::vector<T>& coeffs) const;

    //! Write polynomial in standard basis
    virtual Polynomial<T> toStandardBasis(const T coeffs[], 
					  unsigned int n) const;

    //! Evaluate basis polynomial at zero
    virtual T evaluateZero(unsigned int i) const;

    virtual void print(std::ostream& os) const;

  protected:

    //! Constructor
    OrthogPolyBasisBase(const std::string& name, unsigned int p);

  private:

    // Prohibit copying
    OrthogPolyBasisBase(const OrthogPolyBasisBase&);

    // Prohibit Assignment
    OrthogPolyBasisBase& operator=(const OrthogPolyBasisBase& b);
    
  protected:

    //! Name of basis
    std::string name;

    //! Order of basis
    unsigned int p;

    //! Basis polynomials
    std::vector< Polynomial<T> > basis;

    //! double-sized basis polynomials for accurately projecting products
    std::vector< Polynomial<T> > double_basis;

    //! Norms
    std::vector<T> norms;

  }; // class OrthogPolyBasisBase

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_OrthogPolyBasisBaseImp.hpp"

#endif
