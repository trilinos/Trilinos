// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_SCALARPARAMETERFAMILY_HPP
#define SACADO_SCALARPARAMETERFAMILY_HPP

#include "Sacado_ParameterFamilyBase.hpp"
#include "Sacado_ScalarParameterEntry.hpp"

namespace Sacado {

  using namespace mpl::placeholders;

  //! Specialization of Sacado::ParameterFamilyBase for scalar parameters
  template <typename EvalTypeTraits = DefaultEvalTypeTraits>
  class ScalarParameterFamily :
    public Sacado::ParameterFamilyBase<AbstractScalarParameterEntry,
                                       ScalarParameterEntry<_,EvalTypeTraits> >
  {

    //! Typename synonym of base class
    typedef Sacado::ParameterFamilyBase<AbstractScalarParameterEntry,
                                        ScalarParameterEntry<_,EvalTypeTraits> >  BaseT;

  public:

    //! Constructor
    ScalarParameterFamily(const std::string& name_,
                          bool supports_ad_,
                          bool supports_analytic_) :
      BaseT(name_, supports_ad_, supports_analytic_) {}


    //! Destructor
    virtual ~ScalarParameterFamily() {}

    //! Set paramter value using a real number
    void setRealValueForAllTypes(double value) {
      for (typename BaseT::iterator it = this->family.begin();
           it != this->family.end(); ++it)
        (*it).second->setRealValue(value);
    }

    //! Set real parameter value
    template <class EvalType>
    void
    setRealValue(double value)
    {
      this->template getEntry<EvalType>()->setRealValue(value);
    }

    //! Set parameter to value \em value treating parameter as a constant
    template <class EvalType>
    void
    setValue(const typename EvalTypeTraits::template apply<EvalType>::type& value)
    {
      this->template getEntry<EvalType>()->setValue(value);
    }

    //! Get real parameter value
    template <class EvalType>
    double
    getRealValue() const
    {
      return this->template getEntry<EvalType>()->getRealValue();
    }

    //! Get parameter value
    template <class EvalType>
    const typename EvalTypeTraits::template apply<EvalType>::type&
    getValue() const
    {
      return this->template getEntry<EvalType>()->getValue();
    }

    //! Add a new parameter using custom entry
    /*!
     * Returns true if successful in adding entry to library, false
     * otherwise.
     */

  private:

    //! Private to prohibit copying
    ScalarParameterFamily(const ScalarParameterFamily&);

    //! Private to prohibit copying
    ScalarParameterFamily& operator = (const ScalarParameterFamily&);

  };

  /** \brief Get the value.
   *
   * \relates ScalarParameterFamily
   */
  template <typename EvalType, typename EvalTypeTraits>
  typename Sacado::ScalarParameterEntry<EvalType>::ScalarT
  getValue(const ScalarParameterFamily<EvalTypeTraits>& spf)
  {
    return spf.template getValue<EvalType>();
  }

} // namespace Sacado

#endif
