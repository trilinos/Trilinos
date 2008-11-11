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

#ifndef SACADO_SCALARPARAMETERLIBRARY_HPP
#define SACADO_SCALARPARAMETERLIBRARY_HPP

#include "Sacado_ParameterLibraryBase.hpp"
#include "Sacado_ScalarParameterFamily.hpp"
#include "Sacado_ScalarParameterVector.hpp"

#include "Teuchos_TestForException.hpp"

namespace Sacado {

  /*! 
   * \brief Specialization of Sacado::ParameterLibraryBase for scalar 
   * parameters
   */
  template <typename EvalTypeTraits = DefaultEvalTypeTraits>
  class ScalarParameterLibrary : 
    public ParameterLibraryBase<ScalarParameterFamily<EvalTypeTraits>, 
                                ScalarParameterEntry<_,EvalTypeTraits> > {

  public:

    //! Typename synonym of base class
    typedef ParameterLibraryBase<ScalarParameterFamily<EvalTypeTraits>, 
                                 ScalarParameterEntry<_,EvalTypeTraits> > 
    BaseT;
  
    //! Default constructor
    ScalarParameterLibrary() {}

    //! Destructor
    virtual ~ScalarParameterLibrary() {}

    //! Set paramter value using a real number
    void setRealValueForAllTypes(const std::string& name, double value);

    //! Set real parameter to value \em value
    template <class EvalType>
    void 
    setRealValue(const std::string& name, double value);

    //! Set parameter to value \em value
    template <class EvalType>
    void 
    setValue(const std::string& name, 
             const typename EvalTypeTraits::template apply<EvalType>::type& value);

    //! Get parameter value
    template <class EvalType>
    double
    getRealValue(const std::string& name) const;

    //! Get parameter value
    template <class EvalType>
    const typename EvalTypeTraits::template apply<EvalType>::type& 
    getValue(const std::string& name) const;

    //! Returns a parameter library (singleton object).
    static ScalarParameterLibrary& getInstance() {
      static ScalarParameterLibrary instance;
      return instance;
    }

    //! Fill a vector with the supplied parameter names
    /*!
     * baseValue will be computed from each individual parameter using the
     * corresponding evaluation type EvalType
     */
    template <class EvalType>
    void
    fillVector(const Teuchos::Array<std::string>& names,
               ScalarParameterVector<EvalTypeTraits>& pv);

  private:

    //! Private to prohibit copying
    ScalarParameterLibrary(const ScalarParameterLibrary&);

    //! Private to prohibit copying
    ScalarParameterLibrary& operator = (const ScalarParameterLibrary&);

  };
  
}

template <typename EvalTypeTraits>
void
Sacado::ScalarParameterLibrary<EvalTypeTraits>::
setRealValueForAllTypes(const std::string& name, double value)
{
  typename BaseT::FamilyMap::iterator it = this->library.find(name);
  TEST_FOR_EXCEPTION(
     it == this->library.end(), 
     std::logic_error,
     std::string("Sacado::ScalararameterLibrary::setRealValueForAllTypes():  ")
     + "Invalid parameter family " + name);
  (*it).second->setRealValueForAllTypes(value);
}

template <typename EvalTypeTraits>
template <class EvalType>
void
Sacado::ScalarParameterLibrary<EvalTypeTraits>::
setRealValue(const std::string& name, double value)
{
  typename BaseT::FamilyMap::iterator it = this->library.find(name);
  TEST_FOR_EXCEPTION(
     it == this->library.end(), 
     std::logic_error,
     std::string("Sacado::ScalarParameterLibrary::setValueAsConstant():  ")
     + "Invalid parameter family " + name);
  (*it).second->setRealValue<EvalType>(value);
}

template <typename EvalTypeTraits>
template <class EvalType>
void
Sacado::ScalarParameterLibrary<EvalTypeTraits>::
setValue(
      const std::string& name, 
      const typename EvalTypeTraits::template apply<EvalType>::type& value)
{
  typename BaseT::FamilyMap::iterator it = this->library.find(name);
  TEST_FOR_EXCEPTION(
      it == this->library.end(), 
      std::logic_error,
      std::string("Sacado::ScalarParameterLibrary::setValueAsIndependent():  ")
      + "Invalid parameter family " + name);
  (*it).second->setValue<EvalType>(value);
}

template <typename EvalTypeTraits>
template <class EvalType>
double
Sacado::ScalarParameterLibrary<EvalTypeTraits>::
getRealValue(const std::string& name) const
{
  typename BaseT::FamilyMap::const_iterator it = this->library.find(name);
  TEST_FOR_EXCEPTION(
		 it == this->library.end(), 
		 std::logic_error,
		 std::string("Sacado::ScalarParameterLibrary::getValue():  ")
		 + "Invalid parameter family " + name);
  return (*it).second->getRealValue<EvalType>();
}

template <typename EvalTypeTraits>
template <class EvalType>
const typename EvalTypeTraits::template apply<EvalType>::type&
Sacado::ScalarParameterLibrary<EvalTypeTraits>::
getValue(const std::string& name) const
{
  typename BaseT::FamilyMap::const_iterator it = this->library.find(name);
  TEST_FOR_EXCEPTION(
		 it == this->library.end(), 
		 std::logic_error,
		 std::string("Sacado::ScalarParameterLibrary::getValue():  ")
		 + "Invalid parameter family " + name);
  return (*it).second->getValue<EvalType>();
}

template <typename EvalTypeTraits>
template <class EvalType>
void
Sacado::ScalarParameterLibrary<EvalTypeTraits>::
fillVector(const Teuchos::Array<std::string>& names,
           Sacado::ScalarParameterVector<EvalTypeTraits>& pv)
{
  typename BaseT::FamilyMap::iterator it;

  // Fill in parameters
  for (unsigned int i=0; i<names.size(); i++) {
    it = this->library.find(names[i]);
    TEST_FOR_EXCEPTION(
		   it == this->library.end(), 
		   std::logic_error,
		   std::string("Sacado::ParameterLibraryBase::fillVector():  ")
		   + "Invalid parameter family " + names[i]);
    pv.addParam((*it).second, 0.0);
    pv[i].baseValue = (*it).second->getRealValue<EvalType>();
  }
}


#endif
