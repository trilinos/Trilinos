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
\deprecated Use ex_get_set_param()(exoid, EX_NODE_SET, node_set_id,
num_nodes_in_set, num_df_in_set)

The function ex_get_node_set_param() reads the number of nodes which
describe a single node set and the number of distribution factors for
the node set.

\return In case of an error, ex_get_node_set_param() returns a
negative number; a warning will return a positive number. Possible
causes of errors include:
  -  data file not properly opened with call to ex_create() or ex_open()
  -  a warning value is returned if no node sets are stored in the file.
  -  incorrect node set ID.

\param[in]  exoid            exodus file ID returned from a previous call to
ex_create() or ex_open().
\param[in]  node_set_id      The node set ID.
\param[out] num_nodes_in_set Returned number of nodes in the node set.
\param[out] num_df_in_set    Returned number of distribution factors in the node
set.

The following code segment will read a node set from an open exodus
file :
~~~{.c}
int error, exoid, id, num_nodes_in_set, num_df_in_set, *node_list;
float *dist_fact;

\comment{read node set parameters}
id = 100;

error = ex_get_node_set_param(exoid, id, &num_nodes_in_set,
                              &num_df_in_set);

\comment{read node set node list}
node_list = (int *) calloc(num_nodes_in_set, sizeof(int));
error = ex_get_node_set(exoid, id, node_list);

\comment{read node set distribution factors}
if (num_df_in_set > 0) {
   dist_fact = (float *) calloc(num_nodes_in_set, sizeof(float));
   error = ex_get_node_set_dist_fact(exoid, id, dist_fact);
}

\comment{Same result using non-deprecated functions}
error = ex_get_set_param(exoid, EX_NODE_SET, id, &num_nodes_in_set,
&num_df_in_set);
error = ex_get_set(exoid, EX_NODE_SET, id, node_list);
if (num_df_in_set > 0) {
   error = ex_get_set_dist_fact(exoid, EX_NODE_SET, id, dist_fact);
}

~~~

 */

int ex_get_node_set_param(int exoid, ex_entity_id node_set_id, void_int *num_nodes_in_set,
                          void_int *num_df_in_set)
{
  return ex_get_set_param(exoid, EX_NODE_SET, node_set_id, num_nodes_in_set, num_df_in_set);
}
