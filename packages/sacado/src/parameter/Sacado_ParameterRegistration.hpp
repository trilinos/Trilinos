// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_PARAMETERREGISTRATION_HPP
#define SACADO_PARAMETERREGISTRATION_HPP

#include <sstream>
#include "Sacado_ParameterAccessor.hpp"
#include "Sacado_Traits.hpp"
#include "Sacado_ScalarParameterLibrary.hpp"

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

    typedef ScalarParameterLibrary<EvalTypeTraits> ParamLib;

    //! Constructor: Registers the parameter with the Parameter Library
    ParameterRegistration(const std::string &name_,
                          ParameterAccessor<EvalType, EvalTypeTraits>* access_,
                          ParamLib& paramLib)
      : access(access_), name(name_) {

      if (!paramLib.isParameter(name))
        paramLib.addParameterFamily(name, true, false);
      if (!paramLib.template isParameterForType<EvalType>(name))
        paramLib.template addEntry<EvalType>(name, Teuchos::rcp(this,false));
    }

    //! Constructor: Registers the parameter with the Parameter Library
    ParameterRegistration(const std::string &name_,
                          ParameterAccessor<EvalType, EvalTypeTraits>* access_,
                          const Teuchos::RCP<ParamLib>& paramLib)
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
      access->setValue(name, value);
    }

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
