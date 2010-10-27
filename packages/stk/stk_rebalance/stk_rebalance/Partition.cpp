/*--------------------------------------------------------------------*/
/*    Copyright 2001, 2002 Sandia Corporation.                        */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// Copyright 2001, 2002 Sandia Corporation, Albuquerque, NM.

#include <stdexcept>

#include <stk_rebalance/Partition.hpp>

using namespace stk;
namespace stk {
  using namespace rebalance;
}

Partition::Partition(ParallelMachine comm) :
  comm_(comm)
{
}

Partition::Partition(const Partition & p) :
  comm_(p.comm_)
{
}

Partition::~Partition()
{
}



