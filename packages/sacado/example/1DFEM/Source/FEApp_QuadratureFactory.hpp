// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#ifndef FEAPP_QUADRATUREFACTORY_HPP
#define FEAPP_QUADRATUREFACTORY_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

#include "FEApp_AbstractQuadrature.hpp"

namespace FEApp {

  /*!
   * \brief A factory class to instantiate AbstractQuadrature objects
   */
  class QuadratureFactory {
  public:

    //! Default constructor
    QuadratureFactory(
	       const Teuchos::RefCountPtr<Teuchos::ParameterList>& quadParams);

    //! Destructor
    virtual ~QuadratureFactory() {}

    virtual Teuchos::RefCountPtr<FEApp::AbstractQuadrature>
    create();

  private:

    //! Private to prohibit copying
    QuadratureFactory(const QuadratureFactory&);

    //! Private to prohibit copying
    QuadratureFactory& operator=(const QuadratureFactory&);

  protected:

    //! Parameter list specifying what strategy to create
    Teuchos::RefCountPtr<Teuchos::ParameterList> quadParams;

  };

}

#endif // FEAPP_QUADRATUREFACTORY_HPP
