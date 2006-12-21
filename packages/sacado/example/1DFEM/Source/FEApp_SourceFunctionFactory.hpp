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

#ifndef FEAPP_SOURCEFUNCTIONFACTORY_HPP
#define FEAPP_SOURCEFUNCTIONFACTORY_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

#include "FEApp_AbstractSourceFunction.hpp"

namespace FEApp {

  /*!
   * \brief A factory class to instantiate AbstractSourceFunction objects
   */
  template <typename ScalarT>
  class SourceFunctionFactory {
  public:

    //! Default constructor
    SourceFunctionFactory(
	       const Teuchos::RefCountPtr<Teuchos::ParameterList>& funcParams);

    //! Destructor
    virtual ~SourceFunctionFactory() {}

    virtual Teuchos::RefCountPtr< FEApp::AbstractSourceFunction<ScalarT> >
    create();

  private:

    //! Private to prohibit copying
    SourceFunctionFactory(const SourceFunctionFactory&);

    //! Private to prohibit copying
    SourceFunctionFactory& operator=(const SourceFunctionFactory&);

  protected:

    //! Parameter list specifying what strategy to create
    Teuchos::RefCountPtr<Teuchos::ParameterList> funcParams;

  };

}

// Include implementation
#ifndef SACADO_ETI
#include "FEApp_SourceFunctionFactoryImpl.hpp"
#endif 

#endif // FEAPP_SOURCEFUNCTIONFACTORY_HPP
