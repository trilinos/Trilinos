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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
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
