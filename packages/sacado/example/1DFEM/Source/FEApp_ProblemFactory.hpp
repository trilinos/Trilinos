// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
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

#ifndef FEAPP_PROBLEMFACTORY_HPP
#define FEAPP_PROBLEMFACTORY_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

#include "FEApp_AbstractProblem.hpp"

namespace FEApp {

  /*!
   * \brief A factory class to instantiate AbstractProblem objects
   */
  class ProblemFactory {
  public:

    //! Default constructor
    ProblemFactory(
	   const Teuchos::RefCountPtr<Teuchos::ParameterList>& problemParams);

    //! Destructor
    virtual ~ProblemFactory() {}

    virtual Teuchos::RefCountPtr<FEApp::AbstractProblem>
    create();

  private:

    //! Private to prohibit copying
    ProblemFactory(const ProblemFactory&);

    //! Private to prohibit copying
    ProblemFactory& operator=(const ProblemFactory&);

  protected:

    //! Parameter list specifying what problem to create
    Teuchos::RefCountPtr<Teuchos::ParameterList> problemParams;

  };

}

#endif // FEAPP_PROBLEMFACTORY_HPP
