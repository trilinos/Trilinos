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

#ifndef BOXIDENT_HPP
#define BOXIDENT_HPP

#include "Kokkos_Core.hpp"
namespace stk::search {

template <typename BoxType, typename IdentType>
struct BoxIdent
{
  using box_type = BoxType;
  using ident_type = IdentType;
  using first_type = BoxType;
  using second_type = IdentType;

  BoxType box;
  IdentType ident;
};

template <typename BoxType, typename IdentProcType>
struct BoxIdentProc
{
  using box_type = BoxType;
  using ident_proc_type = IdentProcType;
  using first_type = BoxType;
  using second_type = IdentProcType;

  BoxType box;
  IdentProcType identProc;
};

template <typename DomainIdentType, typename RangeIdentType>
struct IdentIntersection
{
  using domain_ident_type = DomainIdentType;
  using range_ident_type = RangeIdentType;

  DomainIdentType domainIdent;
  RangeIdentType rangeIdent;

  KOKKOS_FORCEINLINE_FUNCTION bool operator==(IdentIntersection<DomainIdentType, RangeIdentType> const& rhs) const
  {
    return domainIdent == rhs.domainIdent && rangeIdent == rhs.rangeIdent;
  }

  KOKKOS_FORCEINLINE_FUNCTION bool operator<(IdentIntersection<DomainIdentType, RangeIdentType> const& rhs) const
  {
    return domainIdent < rhs.domainIdent ||
          (!(rhs.domainIdent < domainIdent) && rangeIdent < rhs.rangeIdent);
  }

};

template <typename DomainIdentProcType, typename RangeIdentProcType>
struct IdentProcIntersection
{
  DomainIdentProcType domainIdentProc;
  RangeIdentProcType rangeIdentProc;

  KOKKOS_FORCEINLINE_FUNCTION bool operator==(IdentProcIntersection<DomainIdentProcType, RangeIdentProcType> const& rhs) const
  {
    return domainIdentProc == rhs.domainIdentProc && rangeIdentProc == rhs.rangeIdentProc;
  }

  KOKKOS_FORCEINLINE_FUNCTION bool operator<(IdentProcIntersection<DomainIdentProcType, RangeIdentProcType> const& rhs) const
  {
    return domainIdentProc < rhs.domainIdentProc ||
          (!(rhs.domainIdentProc < domainIdentProc) && rangeIdentProc < rhs.rangeIdentProc);
  }

};

template<class T>
struct Comparator {

KOKKOS_FUNCTION
bool operator()(T a, T b) const {
  return a < b;
}

};

}

#endif // BOXIDENT_HPP
