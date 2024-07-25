// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_PSEUDO_SPECTRAL_OPERATOR_FACTORY_HPP
#define STOKHOS_PSEUDO_SPECTRAL_OPERATOR_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_PseudoSpectralOperator.hpp"

namespace Stokhos {

  //! Factory for building multivariate quadrature strategies.
  template <typename ordinal_type, typename value_type>
  class PseudoSpectralOperatorFactory {
  public:

    //! Constructor
    PseudoSpectralOperatorFactory() {};

    //! Destructor
    virtual ~PseudoSpectralOperatorFactory() {};

    //! Typename of generated pseudospectral operator
    typedef Stokhos::PseudoSpectralOperator<ordinal_type, value_type> psop_type;

    //! Generate pseudospectral operator object
    static Teuchos::RCP<const psop_type>
    create(Teuchos::ParameterList& sgParams);

  private:

    // Prohibit copying
    PseudoSpectralOperatorFactory(const PseudoSpectralOperatorFactory&);

    // Prohibit Assignment
    PseudoSpectralOperatorFactory& operator=(const PseudoSpectralOperatorFactory& b);

  }; // class PseudoSpectralOperatorFactory

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_PseudoSpectralOperatorFactoryImp.hpp"

#endif // STOKHOS_PSEUDO_SPECTRAL_OPERATOR_FACTORY_HPP
