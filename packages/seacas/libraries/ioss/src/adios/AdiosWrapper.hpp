// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
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
//     * Neither the name of NTESS nor the names of its
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

#ifndef IOSS_Ioad_AdiosWrapper_hpp
#define IOSS_Ioad_AdiosWrapper_hpp

#include "adios/AdiosWrapper.h"
#include <Ioss_Utils.h> // for Utils, IOSS_ERROR, etc

namespace Ioad {

  template <typename T>
  void AdiosWrapper::DefineVariable(const std::string &name, const adios2::Dims &shape,
                                    const adios2::Dims &start, const adios2::Dims &count,
                                    const bool constantDims)
  {
    adios2::Variable<T> var = this->IO::InquireVariable<T>(name);
    if (!var) {
      this->IO::DefineVariable<T>(name, shape, start, count);
    }
  }

  template <typename T>
  T AdiosWrapper::GetAttribute(const std::string &attribute_name, bool ignore_missing,
                               T default_value)
  {
    adios2::Attribute<T> attribute = this->IO::InquireAttribute<T>(attribute_name);
    if (!attribute) {
      if (ignore_missing) {
        return default_value;
      }
      else {
        std::ostringstream errmsg;
        errmsg << "ERROR: " << attribute_name << " not found.\n";
        IOSS_ERROR(errmsg);
      }
    }
    return attribute.Data()[0];
  }

  template <typename T> void AdiosWrapper::GetSync(adios2::Variable<T> var, T *data)
  {
    this->Engine::Get(var, data, adios2::Mode::Sync);
  }

  template <typename T> void AdiosWrapper::GetSync(adios2::Variable<T> var, T &data)
  {
    this->Engine::Get(var, data, adios2::Mode::Sync);
  }

  template <typename T> void AdiosWrapper::GetSync(std::string var_name, T *data)
  {
    this->Engine::Get(var_name, data, adios2::Mode::Sync);
  }

  template <typename T> void AdiosWrapper::GetSync(std::string var_name, T &data)
  {
    this->Engine::Get(var_name, data, adios2::Mode::Sync);
  }

  template <typename T> void AdiosWrapper::InquireAndPut(const std::string &name, const T *value)
  {
    adios2::Variable<T> var = this->IO::InquireVariable<T>(name);
    if (var) {
      this->Engine::Put<T>(name, value,
                           adios2::Mode::Sync); // If not Sync, variables are not saved correctly.
    }
    else {
      std::ostringstream errmsg;
      errmsg << "ERROR: Could not find variable '" << name << "'\n";
      IOSS_ERROR(errmsg);
    }
  }

  template <typename T> void AdiosWrapper::DefineAttribute(const std::string &name, const T &value)
  {
    adios2::Attribute<T> attr = this->IO::InquireAttribute<T>(name);
    if (!attr) {
      this->IO::DefineAttribute<T>(name, value);
    }
  }

  template <typename T>
  void AdiosWrapper::DefineMetaVariable(const std::string &meta_name,
                                        const std::string &variable_name)
  {
    std::string         encoded_name = EncodeMetaVariable(meta_name, variable_name);
    adios2::Variable<T> var          = this->IO::InquireVariable<T>(encoded_name);
    if (!var) {
      this->IO::DefineVariable<T>(encoded_name);
    }
  }

  template <typename T>
  void AdiosWrapper::PutMetaVariable(const std::string &meta_name, T value,
                                     const std::string &variable_name)
  {
    std::string         name = EncodeMetaVariable(meta_name, variable_name);
    adios2::Variable<T> var  = this->IO::InquireVariable<T>(name);
    if (var) {
      this->Engine::Put<T>(var, &value,
                           adios2::Mode::Sync); // If not Sync, variables are not saved correctly.
    }
    else {
      std::ostringstream errmsg;
      errmsg << "ERROR: " << name << " variable not defined.\n";
      IOSS_ERROR(errmsg);
    }
  }

  template <typename T>
  T AdiosWrapper::GetMetaVariable(const std::string &meta_name, const std::string &variable_name)
  {
    T variable;
    this->Engine::Get<T>(EncodeMetaVariable(meta_name, variable_name), variable,
                         adios2::Mode::Sync); // If not Sync, variables are not saved correctly.
    return variable;
  }

} // end of namespace Ioad

#endif
