/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include "exodusII.h"

/*!
\deprecated Use ex_get_set_dist_fact()(exoid, EX_NODE_SET, node_set_id,
node_set_dist_fact)

The function ex_get_node_set_dist_fact() returns the node set
distribution factors for a single node set. Memory must be allocated
for the list of distribution factors(num_dist_in_set in length)
before this function is invoked.

Because the distribution factors are floating point values, the
application code must declare the array passed to be the appropriate
type (float or double) to match the compute word size passed in
ex_create() or ex_open().

\return In case of an error, ex_get_node_set_dist_fact() returns a
negative number; a warning will return a positive number. Possible
causes of errors include:
  -  a warning value is returned if no distribution factors were stored.

\param[in]  exoid               exodus file ID returned from a previous call to
ex_create() or ex_open().
\param[in]  node_set_id         The node set ID.
\param[out] node_set_dist_fact  Returned array containing the distribution
factors in the node set.

Refer to the description of ex_get_node_set_param() for a sample code
segment to read a node set's distribution factors.
*/

int ex_get_node_set_dist_fact(int exoid, ex_entity_id node_set_id, void *node_set_dist_fact)
{
  return ex_get_set_dist_fact(exoid, EX_NODE_SET, node_set_id, node_set_dist_fact);
}
