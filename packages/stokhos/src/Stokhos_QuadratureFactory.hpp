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

#ifndef STOKHOS_QUADRATURE_FACTORY_HPP
#define STOKHOS_QUADRATURE_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_Quadrature.hpp"

namespace Stokhos {

  //! Factory for building multivariate quadrature strategies.
  template <typename ordinal_type, typename value_type>
  class QuadratureFactory {
  public:

    //! Constructor
    QuadratureFactory() {};

    //! Destructor
    virtual ~QuadratureFactory() {};

    //! Generate quadrature object
    static Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >
    create(Teuchos::ParameterList& sgParams);

  private:

    // Prohibit copying
    QuadratureFactory(const QuadratureFactory&);

    // Prohibit Assignment
    QuadratureFactory& operator=(const QuadratureFactory& b);

  }; // class QuadratureFactory

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_QuadratureFactoryImp.hpp"

#endif // STOKHOS_QUADRATURE_FACTORY_HPP
