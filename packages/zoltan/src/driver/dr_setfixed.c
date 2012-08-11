
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "dr_const.h"
#include "dr_loadbal_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*****************************************************************************/
void setup_fixed_obj(MESH_INFO_PTR mesh, int Num_Global_Parts)
{
int i;
int proc, nprocs;

  MPI_Comm_rank(MPI_COMM_WORLD, &proc);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /* Initialize fixed elements in mesh */
  for (i = 0; i < mesh->num_elems; i++) mesh->elements[i].fixed_part = -1;

  switch (Test.Fixed_Objects) {
  case 1:
      /* Fix 100% of objects */
      for (i = 0; i < mesh->num_elems; i++) {
        mesh->elements[i].fixed_part = i % Num_Global_Parts;
      }
      break;
  case 2:
      /* Fix 50% of objects */
      for (i = 1; i < mesh->num_elems; i+=2) {
        mesh->elements[i].fixed_part = i % Num_Global_Parts;
      }
      break;
  case 3:
      /* Fix 10% of objects */
      for (i = 0; i < mesh->num_elems; i+=10) {
        mesh->elements[i].fixed_part = i % Num_Global_Parts;
      }
      break;
  case 4:
      /* Fix 100% of objects to partition 1 */
      for (i = 0; i < mesh->num_elems; i++) {
        mesh->elements[i].fixed_part = 1;
      }
      break;
  case 5:
      /* Fix 50% of objects to partition 1 */
      for (i = mesh->num_elems/2; i < mesh->num_elems; i++) {
        mesh->elements[i].fixed_part = 1;
      }
      break;
  case 6:
      /* Fix 10% of objects to partition 1 */
      for (i = 0; i < mesh->num_elems/10; i++) {
        mesh->elements[i].fixed_part = 1;
      }
      break;
  case 7:
      /* Fix one object on each proc. */
      for (i = 0; i < MIN(1, mesh->num_elems); i++) {
        mesh->elements[i].fixed_part = (proc%Num_Global_Parts);
      }
      break;
  case 8:
      /* Fix half the objects on half the procs, 25% overall */
      if (proc&1){
        for (i = 0; i < mesh->num_elems/2; i++) {
          mesh->elements[i].fixed_part = (proc%Num_Global_Parts);
        }
      }
      break;
  case 9:
      /* Fix first two objects on proc 0 only. */
      if (proc == 0){
        for (i = 0; i < MIN(2, mesh->num_elems); i++) {
          mesh->elements[i].fixed_part = i;
        }
      }
      break;
  case 10:
      /* Fix 0% of objects */
      break;
  }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
