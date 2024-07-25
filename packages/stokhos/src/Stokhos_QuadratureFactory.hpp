// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
