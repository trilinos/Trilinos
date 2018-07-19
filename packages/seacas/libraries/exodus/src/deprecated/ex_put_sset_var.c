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
/*****************************************************************************
 *
 * expev - ex_put_sset_var
 *
 * entry conditions -
 *   input parameters:
 *       int     exoid                   exodus file id
 *       int     time_step               time step number
 *       int     sset_var_index          sideset variable index
 *       int     sset_id                 sideset id
 *       int     num_faces_this_sset     number of faces in this sideset
 *
 * exit conditions -
 *
 *
 * exit conditions -
 *
 * revision history -
 *
 *
 *****************************************************************************/

#include "exodusII.h" // for ex_put_var, ex_entity_id, etc
#include <stdint.h>   // for int64_t

/*!
 * writes the values of a single sideset variable for one sideset at
 * one time step to the database; assume the first time step and
 * sideset variable index are 1
 * \param      exoid                   exodus file id
 * \param      time_step               time step number
 * \param      sset_var_index          sideset variable index
 * \param      sset_id                 sideset id
 * \param      num_faces_this_sset     number of faces in this sideset
 * \param      sset_var_vals           the variable values to be written
 * \deprecated Use ex_put_var()(exoid, time_step, EX_SIDE_SET, sset_var_index,
 * sset_id, num_faces_this_sset, sset_var_vals)
 */

int ex_put_sset_var(int exoid, int time_step, int sset_var_index, ex_entity_id sset_id,
                    int64_t num_faces_this_sset, const void *sset_var_vals)
{
  return ex_put_var(exoid, time_step, EX_SIDE_SET, sset_var_index, sset_id, num_faces_this_sset,
                    sset_var_vals);
}
