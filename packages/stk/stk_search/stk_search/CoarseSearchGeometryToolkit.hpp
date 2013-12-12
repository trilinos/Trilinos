/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_SEARCH_COARSE_SEARCH_GEOMETRY_TOOLKIT_HPP
#define STK_SEARCH_COARSE_SEARCH_GEOMETRY_TOOLKIT_HPP

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/environment/ReportHandler.hpp>

namespace stk { namespace search {

template <typename DomainBox, typename DomainIdent, typename RangeBox, typename RangeIdent>
void coarse_search_geometry_toolkit( std::vector< std::pair<DomainBox,DomainIdent> > const& local_domain,
                                std::vector< std::pair<RangeBox,RangeIdent> > const& local_range,
                                stk::ParallelMachine comm,
                                std::vector<std::pair<DomainIdent, RangeIdent> >& output
                              )
{

}


} }

#endif
