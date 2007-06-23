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

#ifndef SACADO_SCALARPARAMETERFAMILY_HPP
#define SACADO_SCALARPARAMETERFAMILY_HPP

#include "Sacado_ParameterFamilyBase.hpp"
#include "Sacado_ScalarParameterEntry.hpp"

namespace Sacado {
  
  //! Specialization of Sacado::ParameterFamilyBase for scalar parameters
  class ScalarParameterFamily : 
    public Sacado::ParameterFamilyBase<Sacado::AbstractScalarParameterEntry,
                                       Sacado::ScalarParameterEntry> 
  {

    //! Typename synonym of base class
    typedef Sacado::ParameterFamilyBase<Sacado::AbstractScalarParameterEntry,
					Sacado::ScalarParameterEntry>  BaseT;

  public:
  
    //! Constructor
    ScalarParameterFamily(const std::string& name, 
			  bool supports_ad, 
			  bool supports_analytic) : 
      BaseT(name, supports_ad, supports_analytic) {}
      

    //! Destructor
    virtual ~ScalarParameterFamily() {}

    //! Set paramter value using a real number
    void setRealValueForAllTypes(double value) {
      for (iterator it = family.begin(); it != family.end(); ++it)
	(*it).second->setRealValue(value);
    }

    //! Set parameter to value \em value treating parameter as a constant
    template <class ValueType>
    void setValueAsConstant(const ValueType& value) {
      getEntry<ValueType>()->setValueAsConstant(value);
    }

    //! Set parameter to value \em value treating parameter as an independent
    template <class ValueType>
    void setValueAsIndependent(const ValueType& value) {
      getEntry<ValueType>()->setValueAsIndependent(value);
    }

    //! Get parameter value
    template <class ValueType>
    const ValueType& getValue() const {
      return getEntry<ValueType>()->getValue();
    }

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
  template <class ValueType>
  ValueType getValue(const ScalarParameterFamily& spf)
  {
    return spf.template getValue<ValueType>();
  }

} // namespace Sacado

#endif
