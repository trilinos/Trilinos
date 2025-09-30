// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#ifndef STK_STK_SEARCH_STK_SEARCH_OBJECTOUTSIDEDOMAINPOLICY_HPP_
#define STK_STK_SEARCH_STK_SEARCH_OBJECTOUTSIDEDOMAINPOLICY_HPP_

#include <string>

namespace stk {
namespace search {

//BEGINObjectOutsideDomainPolicy
enum class ObjectOutsideDomainPolicy { IGNORE, EXTRAPOLATE, TRUNCATE, PROJECT, ABORT, UNDEFINED_OBJFLAG = 0xff };
//ENDObjectOutsideDomainPolicy

inline ObjectOutsideDomainPolicy get_object_outside_domain_policy(const std::string& id)
{
  if(id == "IGNORE") return ObjectOutsideDomainPolicy::IGNORE;
  if(id == "EXTRAPOLATE") return ObjectOutsideDomainPolicy::EXTRAPOLATE;
  if(id == "TRUNCATE") return ObjectOutsideDomainPolicy::TRUNCATE;
  if(id == "PROJECT") return ObjectOutsideDomainPolicy::PROJECT;
  if(id == "ABORT") return ObjectOutsideDomainPolicy::ABORT;

  return ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG;
}

inline std::string get_object_outside_domain_policy(const ObjectOutsideDomainPolicy id)
{
  switch(id) {
  case ObjectOutsideDomainPolicy::IGNORE:
    return "IGNORE";
  case ObjectOutsideDomainPolicy::EXTRAPOLATE:
    return "EXTRAPOLATE";
  case ObjectOutsideDomainPolicy::TRUNCATE:
    return "TRUNCATE";
  case ObjectOutsideDomainPolicy::PROJECT:
    return "PROJECT";
  case ObjectOutsideDomainPolicy::ABORT:
    return "ABORT";
  case ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG:
    return "UNDEFINED";
  }

  return std::string("");
}

} // namespace search
} // namespace stk


#endif /* STK_STK_SEARCH_STK_SEARCH_OBJECTOUTSIDEDOMAINPOLICY_HPP_ */
