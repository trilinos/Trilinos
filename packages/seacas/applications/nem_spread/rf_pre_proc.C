/*
 * Copyright (C) 2009-2017 National Technology & Engineering Solutions of
 * Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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

#include "rf_allo.h"    // for array_alloc
#include <cstddef>      // for size_t
#include <nem_spread.h> // for NemSpread

template void NemSpread<double, int>::create_elem_types(void);
template void NemSpread<float, int>::create_elem_types(void);
template void NemSpread<double, int64_t>::create_elem_types(void);
template void NemSpread<float, int64_t>::create_elem_types(void);

template <typename T, typename INT> void NemSpread<T, INT>::create_elem_types()
/*
 *      Function which creates a vector of element types for each element.
 */
{
  globals.Elem_Type = (int **)array_alloc(__FILE__, __LINE__, 1, Proc_Info[2], sizeof(int *));

  for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {

    globals.Elem_Type[iproc] = (int *)array_alloc(
        __FILE__, __LINE__, 1, globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc],
        sizeof(int));

    /*
     *     Loop through all the element blocks on the processor,
     *     setting the element type for each element in each block
     */
    size_t ielem_count = 0;
    for (int i = 0; i < globals.Proc_Num_Elem_Blk[iproc]; i++) {
      int ielem_type = globals.Proc_Elem_Blk_Types[iproc][i];
      for (int j = 0; j < globals.Proc_Num_Elem_In_Blk[iproc][i]; j++) {
        globals.Elem_Type[iproc][ielem_count++] = ielem_type;
      }
    }
  }
}
