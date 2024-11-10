// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"

template <typename EntryBase, typename EntryType>
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

template <typename EntryBase, typename EntryType>
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
~ParameterFamilyBase()
{
}

template <typename EntryBase, typename EntryType>
std::string
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
getName() const
{
  return name;
}

template <typename EntryBase, typename EntryType>
bool
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
supportsAD() const
{
  return supports_ad;
}

template <typename EntryBase, typename EntryType>
bool
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
supportsAnalytic() const
{
  return supports_analytic;
}

template <typename EntryBase, typename EntryType>
template <class EvalType>
bool
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
hasType() const
{

  // Convert typename EvalType to string
  std::string evalTypeString = getTypeName<EvalType>();

  // Find entry corresponding to this EvalType
  const_iterator it = family.find(evalTypeString);
  if (it == family.end())
    return false;

  return true;
}

template <typename EntryBase, typename EntryType>
template <class EvalType>
bool
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
addEntry(const Teuchos::RCP< typename Sacado::mpl::apply<EntryType,EvalType>::type >& entry,
         const bool allow_overwrite)
{
  // Get string representation of EvalType
  std::string evalTypeString = getTypeName<EvalType>();

  // Determine if entry already exists for parameter and type
  iterator it = family.find(evalTypeString);

  // If it does not, add it
  if (it == family.end()) {
    family.insert(std::pair<std::string,
                  Teuchos::RCP<EntryBase> >(evalTypeString, entry));
  }
  else if (allow_overwrite) {
    (*it).second = entry;
  }
  else {
    return false;
  }

  return true;
}

template <typename EntryBase, typename EntryType>
template <class EvalType>
Teuchos::RCP< typename Sacado::mpl::apply<EntryType,EvalType>::type >
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
getEntry() {

  // Convert typename EvalType to string
  std::string evalTypeString = getTypeName<EvalType>();

  // Find entry corresponding to this EvalType
  iterator it = family.find(evalTypeString);
  TEUCHOS_TEST_FOR_EXCEPTION(it == family.end(),
                     std::logic_error,
                     std::string("Sacado::ParameterFamilyBase::getEntry():  ")
                     + "Parameter entry " + name
                     + " does not have a parameter of type"
                     + evalTypeString);

  // Cast entry to LOCA::Parameter::Entry<EvalType>
  Teuchos::RCP<  typename Sacado::mpl::apply<EntryType,EvalType>::type > entry = Teuchos::rcp_dynamic_cast< typename Sacado::mpl::apply<EntryType,EvalType>::type >((*it).second);
  TEUCHOS_TEST_FOR_EXCEPTION(entry == Teuchos::null,
                     std::logic_error,
                     std::string("Sacado::ParameterFamilyBase::getEntry():  ")
                     + "Parameter entry " + name
                     + " of type" + evalTypeString
                     + " has incorrect entry type");

  return entry;
}

template <typename EntryBase, typename EntryType>
template <class EvalType>
Teuchos::RCP< const typename Sacado::mpl::apply<EntryType,EvalType>::type >
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
getEntry() const {

  // Convert typename EvalType to string
  std::string evalTypeString = getTypeName<EvalType>();

  // Find entry corresponding to this EvalType
  const_iterator it = family.find(evalTypeString);
  TEUCHOS_TEST_FOR_EXCEPTION(it == family.end(),
                     std::logic_error,
                     std::string("Sacado::ParameterFamilyBase::getEntry():  ")
                     + "Parameter entry " + name
                     + " does not have a parameter of type"
                     + evalTypeString);

  // Cast entry to LOCA::Parameter::Entry<EvalType>
  Teuchos::RCP< const typename Sacado::mpl::apply<EntryType,EvalType>::type > entry = Teuchos::rcp_dynamic_cast< const typename Sacado::mpl::apply<EntryType,EvalType>::type >((*it).second);
  TEUCHOS_TEST_FOR_EXCEPTION(entry == Teuchos::null,
                     std::logic_error,
                     std::string("Sacado::ParameterFamilyBase::getEntry():  ")
                     + "Parameter entry " + name
                     + " of type" + evalTypeString
                     + " has incorrect entry type");

  return entry;
}

template <typename EntryBase, typename EntryType>
void
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
print(std::ostream& os, bool print_values) const
{
  os << "\t" << name << ":  Supports AD = " << supports_ad
     << ", Supports_Analytic = " << supports_analytic << std::endl;
  if (print_values) {
    for (const_iterator it = family.begin(); it != family.end(); it++) {
      os << "\t\t" << (*it).first << " = ";
      (*it).second->print(os);
      os << std::endl;
    }
  }

}

template <typename EntryBase, typename EntryType>
template <class EvalType>
std::string
Sacado::ParameterFamilyBase<EntryBase,EntryType>::
getTypeName() const {
  return typeid(EvalType).name();
}
