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

#include "exodusII.h" // for ex_get_var_time, etc
#include <stdint.h>   // for int64_t

/*!
\deprecated Use ex_get_var_time()(exoid, EX_NODAL, nodal_var_index,
node_number, beg_time_step, end_time_step, node_var_vals)
The function ex_get_nodal_var_time() reads the values of a nodal
variable for a single node through a specified number of time
steps. Memory must be allocated for the nodal variable values array
before this function is invoked.

Because nodal variables are floating point values, the application
code must declare the array passed to be the appropriate type
(float or ``double'') to match the compute word size passed in
ex_create() or ex_open().

\return In case of an error, ex_get_nodal_var_time() returns a
negative number; a warning will return a positive number. Possible
causes of errors include:
  -  specified nodal variable does not exist.
  -  a warning value is returned if no nodal variables are stored in the file.

\param[in] exoid             exodus file ID returned from a previous call to
ex_create() or
                             ex_open().

\param[in] nodal_var_index   The index of the desired nodal variable. The first
variable has an
                             index of 1.

\param[in] node_number       The internal ID (see  Section LocalNodeIds) of the
desired
                             node. The first node is 1.

\param[in] beg_time_step     The beginning time step for which a nodal variable
value
                             is desired. This is not a time value but rather a
time step number,
                             as described under ex_put_time(). The first time
step is 1.

\param[in] end_time_step     The last time step for which a nodal variable value
is desired. If
                             negative, the last time step in the database will
be used. The first
                             time step is 1.

\param[out]  nodal_var_vals  Returned array of(end_time_step {-}
beg_time_step +1) values
                             of the node_number-th node for the
nodal_var_index-th nodal
                             variable.


For example, the following code segment will read the values of the
first nodal variable for node number one for all time steps stored in
the data file:

~~~{.c}
int num_time_steps, var_index, node_num, beg_time, end_time, error, exoid;
float *var_values;

\comment{determine how many time steps are stored}
num_time_steps = ex_inquire_int(exoid, EX_INQ_TIME);

\comment{read a nodal variable through time}
var_values = (float *) calloc (num_time_steps, sizeof(float));

var_index = 1; node_num = 1; beg_time = 1; end_time = -1;
error = ex_get_var_time(exoid, EX_NODAL, var_index, node_num, beg_time,
                              end_time, var_values);

~~~
*/

int ex_get_nodal_var_time(int exoid, int nodal_var_index, int64_t node_number, int beg_time_step,
                          int end_time_step, void *nodal_var_vals)
{
  return ex_get_var_time(exoid, EX_NODAL, nodal_var_index, node_number, beg_time_step,
                         end_time_step, nodal_var_vals);
}
