// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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


void ParameterListModifier::printDoc(std::string const& docString, std::ostream &out) const{
  StrUtils::printLines(out,"# ",docString);
  out << "#  Modifier Used: " << name_ << std::endl;
}


Array<std::string> ParameterListModifier::findMatchingBaseNames(const ParameterList &pl,
    const std::string &base_name, const bool &find_parameters, const bool &find_sublists) const{
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
      if (!allow_base_name){
        throw std::logic_error("Sublist can't have the same name as the parameter template name without `allow_base_name=true`.");
      } else{
        delete_base_name = false;
      }
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


int ParameterListModifier::setDefaultsInSublists(const std::string &param_name, ParameterList &pl,
    const Array<std::string> &sublist_names) const{
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
    pl.remove(param_name);
  }
  return num_defaults;
}


template<typename T>
bool ParameterListModifier::replaceScalarParameterWithArray(const std::string &param_name, const std::string &new_name,
    ParameterList &pl, const bool &throw_if_new_name_exists) const{
  bool param_exists = false;
  if (pl.isParameter(param_name)){
    param_exists = true;
    if (!pl.isType<T>(param_name)){
      throw std::logic_error("The parameter " + new_name + " is not of type " + typeid(T).name());
    }
    if ((pl.isParameter(new_name) || pl.isSublist(new_name)) && throw_if_new_name_exists){
      throw std::logic_error("The parameter " + new_name + " already exists in this parameter list.");
    }
    Array<T> params = tuple<T>(pl.get<T>(param_name));
    pl.set(new_name, params);
    pl.remove(param_name);
  }
  return param_exists;
}

// Instantiation of template member function for `replaceScalarParameterWithArray`
// We usually only convert scalar strings, integers, and floats to arrays
template bool ParameterListModifier::replaceScalarParameterWithArray<std::string>(const std::string&, const std::string&,
    ParameterList&, const bool&) const;
template bool ParameterListModifier::replaceScalarParameterWithArray<int>(const std::string&, const std::string&,
    ParameterList&, const bool&) const;
template bool ParameterListModifier::replaceScalarParameterWithArray<double>(const std::string&, const std::string&,
    ParameterList&, const bool&) const;


} // namespace Teuchos
