// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterListModifier.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StrUtils.hpp"

namespace Teuchos {


// Constructors and/or destructors
ParameterListModifier::ParameterListModifier(const std::string &name_in)
  :name_(name_in)
{}

ParameterListModifier::~ParameterListModifier()
{}


void ParameterListModifier::printDoc(std::string const& docString, std::ostream &out) const
{
  StrUtils::printLines(out,"# ",docString);
  out << "#  Modifier Used: " << name_ << std::endl;
}


Array<std::string> ParameterListModifier::findMatchingBaseNames(const ParameterList &pl,
    const std::string &base_name, const bool &find_parameters, const bool &find_sublists) const
{
  Array<std::string> matches(0);
  ParameterList::ConstIterator itr;
  for (itr = pl.begin(); itr != pl.end(); ++itr) {
    const std::string &name = pl.name(itr);
    std::size_t found = name.find(base_name);
    if (found == 0){
      if (pl.isSublist(name)){
        if (find_sublists){
          matches.push_back(name);
        }
      } else{
        if (find_parameters){
          matches.push_back(name);
        }
      }
    }
  }
  return matches;
}


int ParameterListModifier::expandParameters(
    const std::string &param_template_name, ParameterList &pl, ParameterList &valid_pl,
    const Array<std::string> &exclude_parameters) const
{
  int replacements = 0;
  auto ignore_names = exclude_parameters;
  std::sort(ignore_names.begin(), ignore_names.end());
  ParameterList::ConstIterator itr;
  if (valid_pl.isParameter(param_template_name)){
    ParameterEntry &valid_pl_entry = valid_pl.getEntry(param_template_name);
    for (itr = pl.begin(); itr != pl.end(); ++itr) {
      const std::string &param_name = pl.name(itr);
      if (!pl.isSublist(param_name)){
        if (!std::binary_search(ignore_names.begin(), ignore_names.end(), param_name)){
          valid_pl.setEntry(param_name, valid_pl_entry);
          replacements += 1;
        }
      }
    }
    valid_pl.remove(param_template_name);
  }
  return replacements;
}


int ParameterListModifier::expandSublists(
    const std::string &sublist_template_name, ParameterList &pl, ParameterList &valid_pl,
    const Array<std::string> &exclude_sublists) const
{
  int replacements = 0;
  auto ignore_names = exclude_sublists;
  std::sort(ignore_names.begin(), ignore_names.end());
  ParameterList::ConstIterator itr;
  if (valid_pl.isSublist(sublist_template_name)){
    ParameterList &valid_pl_sublist = valid_pl.get<ParameterList>(sublist_template_name);
    for (itr = pl.begin(); itr != pl.end(); ++itr) {
      const std::string &subname = pl.name(itr);
      if (pl.isSublist(subname)){
        if (!std::binary_search(ignore_names.begin(), ignore_names.end(), subname)){
          valid_pl.set(subname, valid_pl_sublist);
          replacements += 1;
        }
      }
    }
    valid_pl.remove(sublist_template_name);
  }
  return replacements;
}


int ParameterListModifier::expandSublistsUsingBaseName(
    const std::string &base_name, ParameterList &pl, ParameterList &valid_pl,
    const bool &allow_base_name) const
{
  int replacements = 0;
  bool delete_base_name = true;
  if (valid_pl.isSublist(base_name)){
    if (pl.isSublist(base_name)){
      TEUCHOS_TEST_FOR_EXCEPTION(!allow_base_name, std::logic_error,
          "Sublist can't have the same name as the parameter template name "
          "without `allow_base_name=true`.");
      delete_base_name = false;
    }
    Array<std::string> matches = findMatchingBaseNames(pl, base_name, false, true);
    replacements = matches.length();
    for (const std::string &match_name : matches){
      valid_pl.set(match_name, valid_pl.get<ParameterList>(base_name));
    }
    if (delete_base_name){
      valid_pl.remove(base_name);
    }
  }
  return replacements;
}


int ParameterListModifier::setDefaultsInSublists(const std::string &param_name,
    ParameterList &pl, const Array<std::string> &sublist_names,
    const bool remove_param) const
{
  int num_defaults = 0;
  if (pl.isParameter(param_name)){
    for (const std::string &sublist_name : sublist_names){
      if (pl.isSublist(sublist_name)){
        ParameterList &sublist = pl.sublist(sublist_name);
        if (!sublist.isParameter(param_name)){
          ParameterEntry &pl_entry = pl.getEntry(param_name);
          sublist.setEntry(param_name, pl_entry);
        }
      }
    }
    if (remove_param){
      pl.remove(param_name);
    }
  }
  return num_defaults;
}


} // namespace Teuchos
