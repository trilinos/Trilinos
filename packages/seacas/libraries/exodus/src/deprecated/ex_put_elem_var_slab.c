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
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function(s) contained in this file:
 *
 *      ex_put_elem_var_slab()
 *
 *****************************************************************************
 *
 *  Variable Index:
 *
 *      exoid               - The NetCDF ID of an already open NemesisI file.
 *      time_step          - The time step to write this data to.
 *      elem_var_index     - The index of this elemental variable.
 *      elem_blk_id        - The ID of the element block being written to.
 *      start_pos          - The start point for outputting data. The
 *                           first value is 0.
 *      num_vals           - The number of values to be output.
 *      elem_var_vals      - Pointer to the vector of values to be output.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#include "exodusII.h" // for ex_put_partial_var, etc
#include <stdint.h>   // for int64_t

/*
 * writes the values of a single element variable for one element block,
 * starting at start_pos, at one time step to the database; assume the
 * first time step and element variable index are 1
 */

/*!
 * \deprecated Use ex_put_partial_var()(exoid, time_step, EX_ELEM_BLOCK, elem_var_index,
 elem_blk_id, start_pos, num_vals, elem_var_vals)
 */
int ex_put_elem_var_slab(int exoid, int time_step, int elem_var_index, ex_entity_id elem_blk_id,
                         int64_t start_pos, int64_t num_vals, void *elem_var_vals)
{
  return ex_put_partial_var(exoid, time_step, EX_ELEM_BLOCK, elem_var_index, elem_blk_id, start_pos,
                            num_vals, elem_var_vals);
}
