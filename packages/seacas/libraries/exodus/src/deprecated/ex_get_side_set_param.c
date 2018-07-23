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
 * exgsp - ex_get_side_set_param
 *
 * entry conditions -
 *   input parameters:
 *       int     exoid                   exodus file id
 *       int     side_set_id             side set id
 *
 * exit conditions -
 *       int*    num_side_in_set         number of sides in the side set
 *       int*    num_dist_fact_in_set    number of distribution factors in the
 *                                       side set
 *
 * revision history -
 *
 *
 *****************************************************************************/

#include "exodusII.h"

/*!
 * reads the number of sides and the number of distribution factors which
 * describe a single side set
 * \param      exoid                   exodus file id
 * \param      side_set_id             side set id
 * \param[out] num_side_in_set         number of sides in the side set
 * \param[out] num_dist_fact_in_set    number of distribution factors in the
 * \deprecated Use ex_get_set_param()(exoid, EX_SIDE_SET, side_set_id,
 * num_side_in_set, num_dist_fact_in_set)
 */

int ex_get_side_set_param(int exoid, ex_entity_id side_set_id, void_int *num_side_in_set,
                          void_int *num_dist_fact_in_set)
{
  return ex_get_set_param(exoid, EX_SIDE_SET, side_set_id, num_side_in_set, num_dist_fact_in_set);
}
