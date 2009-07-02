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

#ifndef SACADO_PARAMETERREGISTRATION_HPP
#define SACADO_PARAMETERREGISTRATION_HPP

#include <sstream>
#include "Sacado_ParameterAccessor.hpp"
#include "Sacado_Traits.hpp"


namespace Sacado {
  /*!
   * @brief Parameter class for simple registration of a
   * parameter with a Parameter Library. Requires a parameter
   * name a ParameterAccessor object.
   */
  template <typename EvalType, typename EvalTypeTraits = DefaultEvalTypeTraits>
  class ParameterRegistration : 
    public Sacado::ScalarParameterEntry<EvalType, EvalTypeTraits > {

  //! Scalar type
    typedef typename EvalTypeTraits::template apply<EvalType>::type ScalarT;


  public:

    //! Constructor: Registers the parameter with the Parameter Library
    ParameterRegistration(const std::string &name_, ParameterAccessor<EvalType, EvalTypeTraits>* access_,
                            Teuchos::RCP<ParamLib> paramLib)
      : access(access_), name(name_) {

      if (paramLib != Teuchos::null) {
        if (!paramLib->isParameter(name))
          paramLib->addParameterFamily(name, true, false);
        if (!paramLib->template isParameterForType<EvalType>(name)) {
          paramLib->template addEntry<EvalType>(name, Teuchos::rcp(this,false));
        }
      }
    }

    //! Destructor
    virtual ~ParameterRegistration() {}

    //! Set real parameter value
    virtual void setRealValue(double value) { 
      setValue(ScalarT(value)); }

    //! Set parameter values using ParameterAccessor
    virtual void setValue(const ScalarT& value) { 
      access->getValue(name) =  value; }
    
    //! Get parameter value using ParameterAccessor
    virtual const ScalarT& getValue() const { 
      return access->getValue(name); 
     }
    
  protected:  
    
    //! Pointer to source function
     ParameterAccessor<EvalType, EvalTypeTraits>* access;
     const std::string name;

  };
}

#endif // SACADO_PARAMETERREGISTRATION_HPP
