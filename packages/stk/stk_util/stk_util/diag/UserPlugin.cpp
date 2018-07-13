/*
// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
 */

#include <stdexcept>
#include <sstream>
#include <string>

#include <stk_util/diag/UserPlugin.hpp>
#include <stk_util/diag/SlibDiagWriter.hpp>

#ifdef SIERRA_DLOPEN_ENABLED
#include <dlfcn.h>
#endif

namespace sierra {
namespace Plugin {

std::string
derived_id_name(
  int			derived_id)
{
  std::ostringstream derived_name;
  derived_name << "enum id " << derived_id;
  return derived_name.str();
}


Registry::RegistryMap &
Registry::getRegistryMap()
{
  static RegistryMap s_registryMap;

  return s_registryMap;
}


Registry &
Registry::rootInstance()
{
  static Registry registry;

  return registry;
}



void
Registry::registerIt(
  const NamePair &	name_pair,
  void *		func_ptr)
{
  slibout.m(Slib::LOG_PLUGIN) << "Registering " << name_pair.second
			      << " of type " << demangle(name_pair.first->name())
			      << " at " << func_ptr << stk::diag::dendl;

  RegistryMap::const_iterator registry_entry = getRegistryMap().find(name_pair);
  if (registry_entry != getRegistryMap().end() && (*registry_entry).second != func_ptr) {
    std::ostringstream strout;
    strout << "Function with signature " << demangle((*registry_entry).first.first->name())
	   << " and derived name '" << (*registry_entry).first.second
	   << "' already registered to create function at address " << (*registry_entry).second;
    throw std::invalid_argument(strout.str());
  }
  getRegistryMap()[name_pair] = func_ptr;
}


void *
Registry::getPluginPtr(
  const NamePair &	name_pair) const
{
  void *creator_function = getFuncPtr(name_pair);
  if (creator_function)
    return creator_function;
  else {
    std::ostringstream strout;

    strout << "User plugin creator function with base class '" << demangle(name_pair.first->name())
	   << "' and derived class name '" << name_pair.second
	   << "' not found in registry";
    throw std::invalid_argument(strout.str());
  }
}


void *
Registry::getFunctionPtr(
  const NamePair &	name_pair) const
{
  void *creator_function = getFuncPtr(name_pair);
  if (creator_function)
    return creator_function;
  else {
    std::ostringstream strout;
    strout << "User subroutine " << name_pair.second << "\n"
	   << " with signature " << demangle(name_pair.first->name()) << "\n"
           << " not found in registry";
    throw std::invalid_argument(strout.str());
  }
}


Registry *
Registry::getFactoryPtr(
  const NamePair &	name_pair) const
{
  Registry *creator_function = reinterpret_cast<Registry *>(getFuncPtr(name_pair));
  if (creator_function)
    return creator_function;
  else {
    std::ostringstream strout;
    strout << "Registry does not contain function with signature " << demangle(name_pair.first->name())
	   << " and derived name '" << name_pair.second << "'";
    throw std::invalid_argument(strout.str());
  }
}


void *
Registry::getFuncPtr(
  const NamePair &	name_pair) const
{
  RegistryMap::const_iterator registry_entry = getRegistryMap().find(name_pair);
  return registry_entry == getRegistryMap().end() ? nullptr : (*registry_entry).second;
}


std::vector<std::string>
Registry::getDerivedNames(
  const std::type_info &	type) const
{
  std::vector<std::string> derived_names;

  for (RegistryMap::const_iterator it = getRegistryMap().begin(); it != getRegistryMap().end(); ++it)
    if (*(*it).first.first == type)
      derived_names.push_back((*it).first.second);

  return derived_names;
}


typedef void (*dl_register_t)();

void
Registry::registerDL(
  const char *		so_path,
  const char *		function_name)
{
#ifdef SIERRA_DLOPEN_ENABLED
  slibout.m(Slib::LOG_PLUGIN) << "Loading dynamic library " << so_path << stk::diag::dendl;
  void *dl = dlopen(so_path, RTLD_NOW);
  if (!dl){
    throw std::runtime_error(dlerror());
  }

  if (function_name) {
    std::string s = std::strlen(function_name) ? function_name : "dl_register";

    dl_register_t f = reinterpret_cast<dl_register_t>(dlsym(dl, s.c_str()));
    if (!f) {
      s = s + SIERRA_FORTRAN_SUFFIX;

      f = reinterpret_cast<dl_register_t>(dlsym(dl, s.c_str()));
    }

    if (f) {
      slibout.m(Slib::LOG_PLUGIN) << "Executing dynamic library " << so_path << " function " << s << "()" << stk::diag::dendl;
      (*f)();
    }
    else {
      if (std::strlen(function_name)) {
        std::ostringstream str;
        str << "Registration function " << function_name << " not found in " << so_path;
        throw std::runtime_error(str.str().c_str());
      }
    }
  }

#else
  throw std::runtime_error("Dynamic linkage is not supported on this platform");
#endif
}


template <>
void *
Registry::getsym<void *>(
  const char *  sym)
{
#ifdef SIERRA_DLOPEN_ENABLED
  void *s = nullptr;
  void *dl = dlopen(nullptr, RTLD_LAZY);
  if (dl) {
    s = dlsym(dl, sym);
    dlclose(dl);
  }

  return s;
#else
  return nullptr;
#endif
}


std::ostream &
Registry::verbose_print(
  std::ostream &		os) const
{
  for (RegistryMap::const_iterator it = getRegistryMap().begin(); it != getRegistryMap().end(); ++it)
    os << (*it).first.second << " of type " << demangle((*it).first.first->name()) << " at " << (*it).second << std::endl;
  return os;
}


stk::diag::Writer &
Registry::verbose_print(
  stk::diag::Writer &		dout) const
{
  if (dout.shouldPrint()) {
    dout << "Registry, size " << getRegistryMap().size() << stk::diag::push << stk::diag::dendl;

    for (RegistryMap::const_iterator it = getRegistryMap().begin(); it != getRegistryMap().end(); ++it)
      dout << (*it).first.second << " of type " << demangle((*it).first.first->name()) << " at " << (*it).second << stk::diag::dendl;
    dout << stk::diag::pop;
  }
  return dout;
}


extern "C" {

void SIERRA_FORTRAN(register_user_subroutine)(
  type_info_func *		type_id,
  void *			user_subroutine,
  const char *			name,
  int				name_len)
{
  sierra::Plugin::Registry::rootInstance().registerIt(std::make_pair(type_id(), std::string(name, name_len)), user_subroutine);
}

} // extern "C"

} // namespace Plugin
} // namespace sierra
