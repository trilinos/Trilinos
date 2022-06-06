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

namespace stk
{
namespace coupling
{

std::pair<int, int> calc_my_root_and_other_root_ranks(MPI_Comm global, MPI_Comm local);

bool has_split_comm(MPI_Comm global, MPI_Comm local);

MPI_Comm split_comm(MPI_Comm parentCommunicator, int color);

}
}

#endif /* STK_COUPLING_OLD_COMM_SPLITTING_HPP */
