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



#include <fei_macros.hpp>

#include <fei_Param.hpp>

fei::Param::Param(const char* name,
		  const char* value)
  : type_(STRING),
    name_(),
    string_value_(),
    double_value_(0.0),
    int_value_(0),
    bool_value_(false),
    void_value_(NULL)
{
  if (name != 0) name_ = name;
  if (value != 0) string_value_ = value;
}

fei::Param::Param(const char* name,
		  double value)
  : type_(DOUBLE),
    name_(),
    string_value_(),
    double_value_(value),
    int_value_(0),
    bool_value_(false),
    void_value_(NULL)
{
  if (name != 0) name_ = name;
}

fei::Param::Param(const char* name,
		  int value)
  : type_(INT),
    name_(),
    string_value_(),
    double_value_(0.0),
    int_value_(value),
    bool_value_(false),
    void_value_(NULL)
{
  if (name != 0) name_ = name;
}

fei::Param::Param(const char* name,
		  bool value)
  : type_(BOOL),
    name_(),
    string_value_(),
    double_value_(0.0),
    int_value_(0),
    bool_value_(value),
    void_value_(NULL)
{
  if (name != 0) name_ = name;
}

fei::Param::Param(const char* name,
		  const void* value)
  : type_(VOID),
    name_(),
    string_value_(),
    double_value_(0.0),
    int_value_(0),
    bool_value_(false),
    void_value_(value)
{
  if (name != 0) name_ = name;
}

fei::Param::Param(const fei::Param& src)
  : type_(src.type_),
    name_(),
    string_value_(),
    double_value_(0.0),
    int_value_(0),
    void_value_(NULL)
{
  *this = src;
}

fei::Param::~Param()
{
}

fei::Param& fei::Param::operator=(const fei::Param& src)
{
  name_ = src.name_;
  switch(type_) {
  case STRING:
    string_value_ = src.string_value_; break;
  case DOUBLE:
    double_value_ = src.double_value_; break;
  case INT:
    int_value_ = src.int_value_; break;
  case VOID:
    void_value_ = src.void_value_; break;
  case BOOL:
    bool_value_ = src.bool_value_; break;
  case BAD_TYPE:
    break;
  }

  return(*this);
}
