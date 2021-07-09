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

#ifndef stk_search_CoarseSearch_hpp
#define stk_search_CoarseSearch_hpp

#include <stk_util/stk_config.h>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/CoarseSearchKdTree.hpp>
#include <stk_search/SearchMethod.hpp>
#include <vector>
#include <utility>


namespace stk { namespace search {

inline
std::ostream& operator<<(std::ostream &out, SearchMethod method)
{
  switch( method )   {
  case KDTREE:                 out << "KDTREE"; break;
  case MORTON_LINEARIZED_BVH:  out << "MORTON_LINEARIZED_BVH"; break;
  }
  return out;
}

// THIS MIGHT BE WHAT WE ACTUALLY WANT FOR THE INTERFACE.
template <typename DomainBox, typename DomainIdent, typename RangeBox, typename RangeIdent>
void coarse_search(
                    std::vector<std::pair<DomainBox,DomainIdent> > const& domain,
                    std::vector<std::pair<RangeBox, RangeIdent> >  const& range,
                    SearchMethod                                          method,
                    stk::ParallelMachine                                  comm,
                    std::vector< std::pair< IdentProc<DomainIdent, unsigned int>,
                                            IdentProc<RangeIdent, unsigned int> > > &  intersections
                  )
{
  std::cerr << "Future version of coarse_search called" << std::endl;
  abort();
}

// intersections will be those of distributed domain boxes associated with this
// processor rank via get_proc<DomainIdent>(.) that intersect distributed range
// boxes.  Optionally, also include intersections of distributed domain boxes
// with distributed range boxes associated with this processor rank via
// get_proc<RangeIdent>(.).
// Note that the search results that are placed in the intersections result argument
// are not sorted. If the caller needs this vector to be sorted, you must
// call std::sort(intersections.begin(), intersections.end()) or similar.
template <typename DomainBox, typename DomainIdent, typename RangeBox, typename RangeIdent>
void coarse_search( std::vector<std::pair<DomainBox,DomainIdent> > const& domain,
               std::vector<std::pair<RangeBox,RangeIdent> >   const& range,
               SearchMethod                                          method,
               stk::ParallelMachine                                  comm,
               std::vector< std::pair< DomainIdent, RangeIdent> > &  intersections,
               bool communicateRangeBoxInfo=true
             )
{
  switch( method )
  {
  case KDTREE:
    coarse_search_kdtree_driver(domain,range,comm,intersections,communicateRangeBoxInfo);
    break;
  default:
    std::cerr << "coarse_search(..) interface used does not support SearchMethod " << method << std::endl;
    abort();
    break;
  }
}

}} // namespace stk::search

#endif // stk_search_CoarseSearch_hpp
