// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
