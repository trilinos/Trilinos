/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
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
