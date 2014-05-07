//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#ifndef RYTHMOS_RK_BUTCHER_TABLEAU_BUILDER_DECL_HPP
#define RYTHMOS_RK_BUTCHER_TABLEAU_BUILDER_DECL_HPP

#include "Rythmos_Types.hpp"

#include "Rythmos_RKButcherTableauBase.hpp"
#include "Teuchos_ObjectBuilder.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"




namespace Rythmos {

template<class Scalar>
class RKButcherTableauBuilder :
  virtual public Teuchos::ParameterListAcceptor
{
  public:
    RKButcherTableauBuilder();
    virtual ~RKButcherTableauBuilder() {}

    void setRKButcherTableauFactory(
      const RCP<const Teuchos::AbstractFactory<RKButcherTableauBase<Scalar> > > &rkbtFactory,
      const std::string &rkbtFactoryName
      );

    RCP<RKButcherTableauBase<Scalar> > create(
        const std::string &rkbt_name = ""
        ) const;

    /** \name Overridden from Teuchos::ParameterListAcceptor */
    //@{

    /** \brief . */
    void setParameterList(const RCP<Teuchos::ParameterList> & paramList);
  
    /** \brief . */
    RCP<Teuchos::ParameterList> getNonconstParameterList();
    
    /** \brief . */
    RCP<Teuchos::ParameterList> unsetParameterList();
    
    /** \brief. */
    RCP<const ParameterList> getParameterList() const;

    /** \brief. */
    RCP<const Teuchos::ParameterList> getValidParameters() const;
    
    //@}
  private:
    Teuchos::ObjectBuilder<RKButcherTableauBase<Scalar> > builder_;

    void initializeDefaults_();
};

// Nonmember constructor
template<class Scalar>
RCP<RKButcherTableauBuilder<Scalar> > rKButcherTableauBuilder();

// Nonmember helper function
template<class Scalar>
RCP<RKButcherTableauBase<Scalar> > createRKBT(const std::string& rkbt_name);


} // namespace Rythmos


#endif // RYTHMOS_RK_BUTCHER_TABLEAU_BUILDER_DECL_HPP
