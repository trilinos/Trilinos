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
* expss - ex_put_side_set
*
* entry conditions -
*   input parameters:
*       int     exoid                   exodus file id
*       int     side_set_id             side set id
*       int*    side_set_elem_list      array of elements in side set
*       int*    side_set_side_list      array of sides in side set

* exit conditions -
*
* revision history -
*
*
*****************************************************************************/

#include "exodusII.h" // for ex_put_set, void_int, etc

/*!
 * writes the side set element list and side set side list for a single side set
 * \param   exoid                   exodus file id
 * \param   side_set_id             side set id
 * \param  *side_set_elem_list      array of elements in side set
 * \param  *side_set_side_list      array of sides in side set
 * \deprecated  Use ex_put_set()(exoid, EX_SIDE_SET, side_set_id,
 * side_set_elem_list, side_set_side_list)
 */

int ex_put_side_set(int exoid, ex_entity_id side_set_id, const void_int *side_set_elem_list,
                    const void_int *side_set_side_list)
{
  return ex_put_set(exoid, EX_SIDE_SET, side_set_id, side_set_elem_list, side_set_side_list);
}
