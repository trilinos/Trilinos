/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#ifndef _fei_Param_hpp_
#define _fei_Param_hpp_

#include <fei_macros.hpp>
#include <string>

namespace fei {
/** Simple container that pairs a name with a parameter that can be one
    of several different primitive types. This container is used as the
    value-type of the class fei::ParameterSet.

    Instances of fei::Param are fully defined at construction, and
    can not be altered later, except by assignment.
*/
class Param {
  public:
  /** enumeration for parameter-types */
  enum ParamType {
    STRING = 0,
    DOUBLE = 1,
    INT = 2,
    BOOL = 3,
    VOID = 4,
    BAD_TYPE = 5
  };

  /** Constructor */
  Param(const char* name, const char* value);
  /** Constructor */
  Param(const char* name, double value);
  /** Constructor */
  Param(const char* name, int value);
  /** Constructor */
  Param(const char* name, const void* value);
  /** Constructor */
  Param(const char* name, bool value);
  /** Copy Constructor */
  Param(const Param& src);

  /** Destructor */
  virtual ~Param();

  /** Assignment operator */
  Param& operator=(const Param& src);

  /** Query for the name of the parameter */
  const std::string& getName() const;

  /** Query for the type of the parameter */
  ParamType getType() const;

  /** Query for string value. Returned string is empty if
      getType() doesn't return Param::STRING */
  const std::string& getStringValue() const;

  /** Query for double value. Returned double is meaningless if
      getType() doesn't return Param::DOUBLE */
  double getDoubleValue() const;

  /** Query for int value. Returned int is meaningless if
      getType() doesn't return Param::INT */
  int getIntValue() const;

  /** Query for bool value. Returned bool is meaningless if
      getType() doesn't return Param::BOOL */
  bool getBoolValue() const;

  /** Query for void-pointer value. Returned void-pointer is
      meaningless if getType() doesn't return Param::VOID */
  const void* getVoidValue() const;

  private:
  ParamType type_;
  std::string name_;
  std::string string_value_;
  double double_value_;
  int int_value_;
  bool bool_value_;
  const void* void_value_;
};
}//namespace fei

inline
const std::string& fei::Param::getName() const
{
  return name_;
}

inline
fei::Param::ParamType fei::Param::getType() const
{
  return type_;
}

inline
const std::string& fei::Param::getStringValue() const
{
  return string_value_;
}

inline
double fei::Param::getDoubleValue() const
{
  return double_value_;
}

inline
int fei::Param::getIntValue() const
{
  return int_value_;
}

inline
bool fei::Param::getBoolValue() const
{
  return bool_value_;
}

inline
const void* fei::Param::getVoidValue() const
{
  return void_value_;
}
#endif
