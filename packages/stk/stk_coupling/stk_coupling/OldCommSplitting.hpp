/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_COUPLING_OLD_COMM_SPLITTING_HPP
#define STK_COUPLING_OLD_COMM_SPLITTING_HPP

#include <stk_util/parallel/Parallel.hpp>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include "stk_util/stk_config.h"

#ifndef STK_HIDE_DEPRECATED_CODE  // delete October 2022
namespace stk
{
namespace coupling
{

STK_DEPRECATED
std::pair<int, int> calc_my_root_and_other_root_ranks(MPI_Comm global, MPI_Comm local);

STK_DEPRECATED_MSG("prefer stk::couping::are_comms_unequal")
bool has_split_comm(MPI_Comm global, MPI_Comm local);

STK_DEPRECATED
MPI_Comm split_comm(MPI_Comm parentCommunicator, int color);

}
}

#endif
#endif /* STK_COUPLING_OLD_COMM_SPLITTING_HPP */
