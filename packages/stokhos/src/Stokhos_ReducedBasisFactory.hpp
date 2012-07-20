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

#ifndef STOKHOS_REDUCED_BASIS_FACTORY_HPP
#define STOKHOS_REDUCED_BASIS_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_ReducedPCEBasis.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_Quadrature.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

namespace Stokhos {

  /*! 
   * \brief Generate a basis from a given set of PCE expansions that is 
   * orthogonal with respect to the product measure induced by these expansions.
   */
  template <typename ordinal_type, typename value_type>
  class ReducedBasisFactory {
  public:

    //! Constructor
    /*!
     * \param params Parameter dictating choice of reduction method
     */
    ReducedBasisFactory(const Teuchos::ParameterList& params);

    //! Destructor
    virtual ~ReducedBasisFactory() {}

    //! Get reduced quadrature object
    virtual Teuchos::RCP<Stokhos::ReducedPCEBasis<ordinal_type, value_type> >
    createReducedBasis(
      ordinal_type p,
      const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
      const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
      const Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk) const;

  private:

    // Prohibit copying
    ReducedBasisFactory(const ReducedBasisFactory&);

    // Prohibit Assignment
    ReducedBasisFactory& operator=(const ReducedBasisFactory&);
    
  protected:

    //! Algorithm parameters
    Teuchos::ParameterList params;

    //! Reduction method
    std::string reduction_method;

  }; // class ReducedBasisFactory

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_ReducedBasisFactoryImp.hpp"

#endif
