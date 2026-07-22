// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_PARAMETERACCESSOR_HPP
#define SACADO_PARAMETERACCESSOR_HPP

#include <vector>
#include "Teuchos_RCP.hpp"
#include "Sacado_ScalarParameterEntry.hpp"
#include "Sacado_ScalarParameterLibrary.hpp"

namespace Sacado {

  template <typename EvalType, typename EvalTypeTraits>
  class ParameterRegistration;

  /*!
   * \brief Abstract class that provides access to a parameter
   * value in a code for the parameter library. An object of this
   * type is required to construct a ParameterRegistration object.
   */
  template<typename EvalType,
           typename EvalTypeTraits = DefaultEvalTypeTraits>
  class ParameterAccessor {
  private:
    typedef typename EvalTypeTraits::template apply<EvalType>::type ScalarT;

  public:

    typedef ScalarParameterLibrary<EvalTypeTraits> ParamLib;

    virtual ~ParameterAccessor() {};

    //! Method that returns a reference to the parameter value given the name
    //! The ParameterLibrary call this method when a parameter value changes
    virtual ScalarT& getValue(const std::string &n) = 0;

    //! Method that returns a reference to the parameter value given the name
    //! The ParameterLibrary call this method when a parameter value changes
    virtual void setValue(const std::string &n, const ScalarT& v) {
      getValue(n) = v;
    }

    void registerSacadoParameter(const std::string& name,
                                 ParamLib& paramLib);

    void registerSacadoParameter(const std::string& name,
                                 const Teuchos::RCP<ParamLib>& paramLib);

  private:
    std::vector< Teuchos::RCP< ParameterRegistration<EvalType, EvalTypeTraits> > > pr_;
  };
}

// Include implementation
#include "Sacado_ParameterAccessorImp.hpp"

#endif
