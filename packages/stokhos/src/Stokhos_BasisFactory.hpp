// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_BASIS_FACTORY_HPP
#define STOKHOS_BASIS_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"

namespace Stokhos {

  //! Factory for building multivariate orthogonal polynomial bases.
  template <typename ordinal_type, typename value_type>
  class BasisFactory {
  public:

    //! Constructor
    BasisFactory() {};

    //! Destructor
    virtual ~BasisFactory() {};

    //! Generate multivariate basis
    static Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >
    create(Teuchos::ParameterList& sgParams);

  protected:

    //! Generate 1-D basis
    static Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> > create1DBasis(Teuchos::ParameterList& params);

  private:

    // Prohibit copying
    BasisFactory(const BasisFactory&);

    // Prohibit Assignment
    BasisFactory& operator=(const BasisFactory& b);

  }; // class BasisFactory

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_BasisFactoryImp.hpp"

#endif // STOKHOS_BASIS_FACTORY_HPP
