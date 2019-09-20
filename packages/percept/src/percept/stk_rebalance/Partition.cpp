// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.




#include <stdexcept>

#include <percept/stk_rebalance/Partition.hpp>

using namespace stk;
namespace stk {
  using namespace rebalance;
}

Partition::Partition(ParallelMachine comm) : comm_(comm) { }
Partition::~Partition() { }



