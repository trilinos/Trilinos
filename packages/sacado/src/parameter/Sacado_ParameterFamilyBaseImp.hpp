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

#include "Teuchos_TestForException.hpp"

template <typename EntryBase, template<typename> class EntryType>
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
ParameterFamilyBase(const std::string& name_, 
		    bool supports_ad_, 
		    bool supports_analytic_) :
  family(),
  name(name_),
  supports_ad(supports_ad_),
  supports_analytic(supports_analytic_)
{
}

template <typename EntryBase, template<typename> class EntryType>
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
~ParameterFamilyBase() 
{
}

template <typename EntryBase, template<typename> class EntryType>
std::string
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
getName() const 
{
  return name;
}

template <typename EntryBase, template<typename> class EntryType>
bool
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
supportsAD() const
{
  return supports_ad;
}

template <typename EntryBase, template<typename> class EntryType>
bool
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
supportsAnalytic() const
{
  return supports_analytic;
}

template <typename EntryBase, template<typename> class EntryType>
template <class ValueType>
bool
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
hasType() const 
{

  // Convert typename ValueType to string
  std::string valueTypeString = getTypeName<ValueType>();

  // Find entry corresponding to this ValueType
  const_iterator it = family.find(valueTypeString);
  if (it == family.end()) 
    return false;

  return true;
}

template <typename EntryBase, template<typename> class EntryType>
template <class ValueType>
bool
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
addEntry(const Teuchos::RCP< EntryType<ValueType> >& entry)
{
  // Get string representation of ValueType
  std::string valueTypeString = getTypeName<ValueType>();
    
  // Determine if entry already exists for parameter and type
  iterator it = family.find(valueTypeString);

  // If it does not, add it
  if (it == family.end()) {
    family.insert(std::pair<std::string, 
		  Teuchos::RCP<EntryBase> >(valueTypeString, entry));
  }

  else {
    return false;
  }

  return true;
}

template <typename EntryBase, template<typename> class EntryType>
template <class ValueType>
Teuchos::RCP< EntryType<ValueType> >
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
getEntry() {

  // Convert typename ValueType to string
  std::string valueTypeString = getTypeName<ValueType>();

  // Find entry corresponding to this ValueType
  iterator it = family.find(valueTypeString);
  TEST_FOR_EXCEPTION(it == family.end(), 
		     std::logic_error,
		     std::string("Sacado::ParameterFamilyBase::getEntry():  ")
		     + "Parameter entry " + name
		     + " does not have a parameter of type" 
		     + valueTypeString);

  // Cast entry to LOCA::Parameter::Entry<ValueType>
  Teuchos::RCP< EntryType<ValueType> > entry = 
    Teuchos::rcp_dynamic_cast< EntryType<ValueType> >((*it).second);
  TEST_FOR_EXCEPTION(entry == Teuchos::null, 
		     std::logic_error,
		     std::string("Sacado::ParameterFamilyBase::getEntry():  ")
		     + "Parameter entry " + name
		     + " of type" + valueTypeString
		     + " has incorrect entry type");
  
  return entry;
}

template <typename EntryBase, template<typename> class EntryType>
template <class ValueType>
Teuchos::RCP< const EntryType<ValueType> >
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
getEntry() const {

  // Convert typename ValueType to string
  std::string valueTypeString = getTypeName<ValueType>();

  // Find entry corresponding to this ValueType
  const_iterator it = family.find(valueTypeString);
  TEST_FOR_EXCEPTION(it == family.end(), 
		     std::logic_error,
		     std::string("Sacado::ParameterFamilyBase::getEntry():  ")
		     + "Parameter entry " + name
		     + " does not have a parameter of type" 
		     + valueTypeString);

  // Cast entry to LOCA::Parameter::Entry<ValueType>
  Teuchos::RCP< const EntryType<ValueType> > entry =
    Teuchos::rcp_dynamic_cast< const EntryType<ValueType> >((*it).second);
  TEST_FOR_EXCEPTION(entry == Teuchos::null, 
		     std::logic_error,
		     std::string("Sacado::ParameterFamilyBase::getEntry():  ")
		     + "Parameter entry " + name
		     + " of type" + valueTypeString
		     + " has incorrect entry type");

  return entry;
}

template <typename EntryBase, template<typename> class EntryType>
void
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
printFamily(std::ostream& os) const
{
  os << "Parameter family map:  " << name << std::endl;

  // Loop over all entries
  for (const_iterator it = family.begin(); it != family.end(); it++) {
    os << "\t" << (*it).first << std::endl;
  }

}

template <typename EntryBase, template<typename> class EntryType>
template <class ValueType>
std::string
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
getTypeName() const {
  return typeid(ValueType).name();
}
