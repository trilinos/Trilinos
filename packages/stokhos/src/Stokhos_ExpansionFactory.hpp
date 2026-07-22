// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_EXPANSION_FACTORY_HPP
#define STOKHOS_EXPANSION_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_OrthogPolyExpansion.hpp"

namespace Stokhos {

  //! Factory for building multivariate expansion strategies.
  template <typename ordinal_type, typename value_type>
  class ExpansionFactory {
  public:

    //! Constructor
    ExpansionFactory() {};

    //! Destructor
    virtual ~ExpansionFactory() {};

    //! Generate multivariate expansion
    static Teuchos::RCP<Stokhos::OrthogPolyExpansion<ordinal_type, value_type> >
    create(Teuchos::ParameterList& sgParams);

  private:

    // Prohibit copying
    ExpansionFactory(const ExpansionFactory&);

    // Prohibit Assignment
    ExpansionFactory& operator=(const ExpansionFactory& b);

  }; // class ExpansionFactory

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_ExpansionFactoryImp.hpp"

#endif // STOKHOS_EXPANSION_FACTORY_HPP
