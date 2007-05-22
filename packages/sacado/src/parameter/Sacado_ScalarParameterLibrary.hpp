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

#include "Teuchos_TestForException.hpp"

namespace Sacado {

  /*! 
   * \brief Specialization of Sacado::ParameterLibraryBase for scalar 
   * parameters
   */
  class ScalarParameterLibrary : 
    public ParameterLibraryBase<ScalarParameterFamily, ScalarParameterEntry> {

  public:
  
    //! Default constructor
    ScalarParameterLibrary() {}

    //! Destructor
    virtual ~ScalarParameterLibrary() {}

    //! Set paramter value using a real number
    void setRealValueForAllTypes(const std::string& name, double value);

    //! Set parameter to value \em value
    template <class ValueType>
    void setValue(const std::string& name, const ValueType& value);

    //! Get parameter value
    template <class ValueType>
    ValueType getValue(const std::string& name) const;

    //! Returns a parameter library (singleton object).
    static ScalarParameterLibrary& getInstance() {
      static ScalarParameterLibrary instance;
      return instance;
    }

  private:

    //! Private to prohibit copying
    ScalarParameterLibrary(const ScalarParameterLibrary&);

    //! Private to prohibit copying
    ScalarParameterLibrary& operator = (const ScalarParameterLibrary&);

  };
  
}

template <class ValueType>
void
Sacado::ScalarParameterLibrary::
setValue(const std::string& name, const ValueType& value)
{
  FamilyMap::iterator it = library.find(name);
  TEST_FOR_EXCEPTION(
		 it == library.end(), 
		 std::logic_error,
		 std::string("Sacado::ScalarParameterLibrary::setValue():  ")
		 + "Invalid parameter family " + name);
  (*it).second->setValue(value);
}

template <class ValueType>
ValueType
Sacado::ScalarParameterLibrary::
getValue(const std::string& name) const
{
  FamilyMap::const_iterator it = library.find(name);
  TEST_FOR_EXCEPTION(
		 it == library.end(), 
		 std::logic_error,
		 std::string("Sacado::ScalarParameterLibrary::getValue():  ")
		 + "Invalid parameter family " + name);
  return (*it).second->getValue<ValueType>();
}

#endif
