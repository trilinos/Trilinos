// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"

template <typename FamilyType, typename EntryType>
Sacado::ParameterLibraryBase<FamilyType,EntryType>::
ParameterLibraryBase()
{
}

template <typename FamilyType, typename EntryType>
Sacado::ParameterLibraryBase<FamilyType,EntryType>::
~ParameterLibraryBase()
{
}

template <typename FamilyType, typename EntryType>
bool
Sacado::ParameterLibraryBase<FamilyType,EntryType>::
isParameter(const std::string& name) const
{
  // Get family
  typename FamilyMap::const_iterator it = library.find(name);

  return (it != library.end());
}

template <typename FamilyType, typename EntryType>
template <class EvalType>
bool
Sacado::ParameterLibraryBase<FamilyType,EntryType>::
isParameterForType(const std::string& name) const
{
  // Get family
  typename FamilyMap::const_iterator it = library.find(name);

  // First check parameter is in the library
  if (it == library.end())
    return false;

  // Determine if type is in the family
  return (*it).second->template hasType<EvalType>();
}

template <typename FamilyType, typename EntryType>
bool
Sacado::ParameterLibraryBase<FamilyType,EntryType>::
addParameterFamily(const std::string& name,
                   bool supports_ad,
                   bool supports_analytic)
{
  // Check that the parameter is not in the library
  if (isParameter(name))
    return false;

  Teuchos::RCP<FamilyType> f =
    Teuchos::rcp(new FamilyType(name, supports_ad, supports_analytic));
  library.insert(std::pair< std::string,
                 Teuchos::RCP<FamilyType> >(name, f));

  return true;
}

template <typename FamilyType, typename EntryType>
template <class EvalType>
bool
Sacado::ParameterLibraryBase<FamilyType,EntryType>::
addEntry(const std::string& name,
         const Teuchos::RCP< typename Sacado::mpl::apply<EntryType,EvalType>::type >& entry,
         const bool allow_overwrite)
{
  // Get family
  typename FamilyMap::iterator it = library.find(name);

  // First check parameter is in the library
  TEUCHOS_TEST_FOR_EXCEPTION(it == library.end(),
                     std::logic_error,
                     std::string("Sacado::ParameterLibraryBase::addEntry():  ")
                     + "Parameter family " + name
                     + " is not in the library");

  // Call family's addEntry method
  return (*it).second->template addEntry<EvalType>(entry, allow_overwrite);
}

template <typename FamilyType, typename EntryType>
template <class EvalType>
Teuchos::RCP< typename Sacado::mpl::apply<EntryType,EvalType>::type >
Sacado::ParameterLibraryBase<FamilyType,EntryType>::
getEntry(const std::string& name)
{
  // Get family
  typename FamilyMap::iterator it = library.find(name);

  // First check parameter is in the library
  TEUCHOS_TEST_FOR_EXCEPTION(it == library.end(),
                     std::logic_error,
                     std::string("Sacado::ParameterLibraryBase::getEntry():  ")
                     + "Parameter family " + name
                     + " is not in the library");

  // Call family's getEntry method
  return (*it).second->template getEntry<EvalType>();
}

template <typename FamilyType, typename EntryType>
template <class EvalType>
Teuchos::RCP< const typename Sacado::mpl::apply<EntryType,EvalType>::type >
Sacado::ParameterLibraryBase<FamilyType,EntryType>::
getEntry(const std::string& name) const
{
  // Get family
  typename FamilyMap::const_iterator it = library.find(name);

  // First check parameter is in the library
  TEUCHOS_TEST_FOR_EXCEPTION(it == library.end(),
                     std::logic_error,
                     std::string("Sacado::ParameterLibraryBase::getEntry():  ")
                     + "Parameter family " + name
                     + " is not in the library");

  // Call family's getEntry method
  return (*it).second->template getEntry<EvalType>();
}

template <typename FamilyType, typename EntryType>
template <typename BaseValueType>
void
Sacado::ParameterLibraryBase<FamilyType,EntryType>::
fillVector(const Teuchos::Array<std::string>& names,
           const Teuchos::Array<BaseValueType>& values,
           ParameterVectorBase<FamilyType,BaseValueType>& pv)
{
  typename FamilyMap::iterator it;

  // Fill in parameters
  for (unsigned int i=0; i<names.size(); i++) {
    it = library.find(names[i]);
    TEUCHOS_TEST_FOR_EXCEPTION(
                   it == library.end(),
                   std::logic_error,
                   std::string("Sacado::ParameterLibraryBase::fillVector():  ")
                   + "Invalid parameter family " + names[i]);
    pv.addParam((*it).second, values[i]);
  }
}

template <typename FamilyType, typename EntryType>
void
Sacado::ParameterLibraryBase<FamilyType,EntryType>::
print(std::ostream& os, bool print_values) const
{
  os << "Library of all registered parameters:" << std::endl;
  typename FamilyMap::const_iterator it = this->library.begin();
  for (; it != this->library.end(); ++it) {
    (*it).second->print(os, print_values);
  }
}
