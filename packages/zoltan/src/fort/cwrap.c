// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/*****************************************************************************
 * KDD:  2/7/11
 * KDD:  This interface uses pointer aliasing to translate C pointers to 
 * KDD:  something that F90 can store and return.  It worked fine until 
 * KDD:  gcc 4.5, where compiler optimizations removed assignment 
 * KDD:  statements in Zfw_Set_Fn (assignments to the Fortran callbacks), 
 * KDD:  causing segmentation faults when Zoltan attempted to call those
 * KDD:  callbacks.  This problem was fixed by declaring all Zoltan_Structs
 * KDD:  using this aliasing to be "volatile".  A better fix would use the
 * KDD:  F90 "C_PTR" type in fwrap.f90.  We can do that fix as time permits.
 * KDD:  See Bugzilla bug 5077 for more info.
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_const.h"
#include "all_allo_const.h"
#include "cwrap_fmangle.h"

/* reconstruct the lb pointer from the nbyte 1-byte integers in addr_lb */
#define ADDR_TO_LB(addr_lb, lb) { \
   unsigned char *p; \
   int i; \
   p = (unsigned char *) &(lb); \
   for (i=0; i<(int)(sizeof(struct Zoltan_Struct *)); i++) \
     {*p = (unsigned char)(addr_lb)[i]; p++;} }

/* construct the nbyte 1-byte integers from the lb pointer */
/* Always assuming 64-bit pointers in F90 interface. */
/* If 32-bit pointers, pad remaining fields with 0.  */
#define LB_TO_ADDR(lb, addr_lb, nbytes) { \
   unsigned char *p; \
   int i; \
   p = (unsigned char *) &(lb); \
   for (i = 0; i < (int)(sizeof(struct Zoltan_Struct *)); i++) { \
     (addr_lb)[i] = (int) *p;  \
     p++; \
   } \
   for (i = (int)(sizeof(struct Zoltan_Struct *)); i < (nbytes); i++) \
     (addr_lb)[i] = 0; \
 }

/*--------------------------------------------------------------------*/
/* Variables                                                          */

static struct Zoltan_Struct *Zoltan_Current;
void Zoltan_Reftree_Get_Child_Order(struct Zoltan_Struct *, int *, int *);

/*--------------------------------------------------------------------*/
/* Utilities                                                          */

/* some MPI implementations may require conversion between a Fortran
   communicator and a C communicator.  This routine is used to perform the
   conversion.  It may need different forms for different MPI libraries. */

MPI_Comm Zoltan_comm_f2c(int *f_comm)
{
#ifndef NO_MPI2
/* MPI 2 provides a standard way of doing this */
   return MPI_Comm_f2c((MPI_Fint)(*f_comm));
#else
/* will probably need some special cases here */
/* when in doubt, just return the input */
   return (MPI_Comm)(*f_comm);
#endif
}

/*****************************************************************************/
/* These routines get the address of an array allocated by fortran and
   return it */
void Zfw_Get_Address_int(int *addr,
                         long *ret_addr)
{
   /* Assuming sizeof(long) == pointer size.  True on most linux systems. */
   /* May not be true in 64-bit Windows, but does anyone use F90 on Windoze? */
   if (sizeof(long) != sizeof(int *)) {
     ZOLTAN_PRINT_ERROR(-1, "Zfw_Get_Address_int", 
       "sizeof(long) != sizeof(int *); F90 allocation will not work properly.\n Contact Zoltan developers for help.");
   }
   *ret_addr = (long)addr;
}

void Zfw_Get_Address_struct(int *addr,
			    long *ret_addr)
{
   /* Assuming sizeof(long) == pointer size.  True on most linux systems. */
   /* May not be true in 64-bit Windows, but does anyone use F90 on Windoze? */
   if (sizeof(long) != sizeof(int *)) {
     ZOLTAN_PRINT_ERROR(-1, "Zfw_Get_Address_struct", 
       "sizeof(long) != sizeof(int *); F90 allocation will not work properly.\n Contact Zoltan developers for help.");
   }
   *ret_addr = (long)addr;
}

/*****************************************************************************/
int Zfw_Get_Wgt_Dim(int *addr_lb, int *nbytes)
{
   volatile struct Zoltan_Struct *lb;
   ADDR_TO_LB(addr_lb, lb);
   return lb->Obj_Weight_Dim;
}

/*****************************************************************************/
int Zfw_Get_Comm_Dim(int *addr_lb, int *nbytes)
{
   volatile struct Zoltan_Struct *lb;
   ADDR_TO_LB(addr_lb, lb);
   return lb->Edge_Weight_Dim;
}

/*--------------------------------------------------------------------*/
/* Reverse wrappers for callbacks                                     */
/*--------------------------------------------------------------------*/

void Zoltan_Part_Multi_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, int num_obj,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *parts,
  int *ierr)
{
   Zoltan_Current->Get_Part_Multi_Fort(data,
                       &num_gid_entries, &num_lid_entries, &num_obj,
                       global_id, local_id, parts, ierr);
}

/*****************************************************************************/

int Zoltan_Part_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  int *ierr)
{
   return Zoltan_Current->Get_Part_Fort(data,
                                            &num_gid_entries, &num_lid_entries,
                                            global_id, local_id, ierr);
}

/*****************************************************************************/
int Zoltan_Num_Edges_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  int *ierr)
{
   return Zoltan_Current->Get_Num_Edges_Fort(data,
                                            &num_gid_entries, &num_lid_entries,
                                            global_id, local_id, ierr);
}

/*****************************************************************************/
void Zoltan_Num_Edges_Multi_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, int num_obj,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *num_edges,
  int *ierr)
{
   Zoltan_Current->Get_Num_Edges_Multi_Fort(data,
                                            &num_gid_entries, &num_lid_entries,
                                            &num_obj, global_id, local_id, 
                                            num_edges, ierr);
}

/*****************************************************************************/
void Zoltan_Edge_List_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  ZOLTAN_ID_PTR nbor_global_id, int *nbor_procs,
  int wdim, float *nbor_ewgts, int *ierr)
{
   Zoltan_Current->Get_Edge_List_Fort(data, &num_gid_entries, &num_lid_entries,
                                     global_id, local_id,
                                     nbor_global_id, nbor_procs, &wdim,
                                     nbor_ewgts, ierr);
}

/*****************************************************************************/
void Zoltan_Edge_List_Multi_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, int num_obj,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *num_edges,
  ZOLTAN_ID_PTR nbor_global_id, int *nbor_procs,
  int wdim, float *nbor_ewgts, int *ierr)
{
   Zoltan_Current->Get_Edge_List_Multi_Fort(data,
                                            &num_gid_entries, &num_lid_entries,
                                            &num_obj, global_id, local_id,
                                            num_edges,
                                            nbor_global_id, nbor_procs, &wdim,
                                            nbor_ewgts, ierr);
}

/*****************************************************************************/
int Zoltan_Num_Geom_Fort_Wrapper(void *data, int *ierr)
{
   return Zoltan_Current->Get_Num_Geom_Fort(data,ierr);
}

/*****************************************************************************/
void Zoltan_Geom_Multi_Fort_Wrapper(
  void *data, int num_gid_entries, int num_lid_entries, int num_obj,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int num_dim,
  double *geom_vec, int *ierr)
{
   Zoltan_Current->Get_Geom_Multi_Fort(data, &num_gid_entries, &num_lid_entries,
                                       &num_obj, global_id, local_id,
                                       &num_dim, geom_vec, ierr);
}
/*****************************************************************************/
void Zoltan_Geom_Fort_Wrapper(
  void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  double *geom_vec, int *ierr)
{
   Zoltan_Current->Get_Geom_Fort(data, &num_gid_entries, &num_lid_entries,
                                global_id, local_id, geom_vec, ierr);
}

/*****************************************************************************/
int Zoltan_Num_Obj_Fort_Wrapper(void *data, int *ierr)
{
   return Zoltan_Current->Get_Num_Obj_Fort(data, ierr);
}

/*****************************************************************************/
void Zoltan_Obj_List_Fort_Wrapper(void *data,
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
  int wdim, float *objwgts, int *ierr)
{
   Zoltan_Current->Get_Obj_List_Fort(data, &num_gid_entries, &num_lid_entries,
                                    global_ids, local_ids, &wdim,
                                    objwgts, ierr);
}

/*****************************************************************************/
int Zoltan_First_Obj_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, 
  ZOLTAN_ID_PTR first_global_id,
  ZOLTAN_ID_PTR first_local_id,
  int wdim, float *first_obj_wgt, int *ierr)
{
   return Zoltan_Current->Get_First_Obj_Fort(data, 
                                            &num_gid_entries, &num_lid_entries,
                                            first_global_id,
                                            first_local_id, &wdim,
                                            first_obj_wgt, ierr);
}

/*****************************************************************************/
int Zoltan_Next_Obj_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, 
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  ZOLTAN_ID_PTR next_global_id, ZOLTAN_ID_PTR next_local_id,
  int wdim, float *next_obj_wgt, int *ierr)
{
   return Zoltan_Current->Get_Next_Obj_Fort(data, 
                                           &num_gid_entries, &num_lid_entries, 
                                           global_id, local_id,
                                           next_global_id, next_local_id,
                                           &wdim, next_obj_wgt, ierr);
}

/*****************************************************************************/
int Zoltan_Num_Border_Obj_Fort_Wrapper(void *data, int nbor_proc, int *ierr)
{
   return Zoltan_Current->Get_Num_Border_Obj_Fort(data, &nbor_proc, ierr);
}

/*****************************************************************************/
void Zoltan_Border_Obj_List_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, 
  int nbor_proc,
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
  int wdim, float *objwgts, int *ierr)
{
   Zoltan_Current->Get_Border_Obj_List_Fort(data, 
                                           &num_gid_entries, &num_lid_entries, 
                                           &nbor_proc, global_ids,
                                           local_ids, &wdim, objwgts, ierr);
}

/*****************************************************************************/
int Zoltan_First_Border_Obj_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  int nbor_proc,
  ZOLTAN_ID_PTR first_global_id,
  ZOLTAN_ID_PTR first_local_id,
  int wdim, float *first_obj_wgt,
  int *ierr)
{
   return Zoltan_Current->Get_First_Border_Obj_Fort(data, 
                                                   &num_gid_entries, 
                                                   &num_lid_entries, 
                                                   &nbor_proc,
                                                   first_global_id,
                                                   first_local_id, &wdim,
                                                   first_obj_wgt, ierr);
}

/*****************************************************************************/
int Zoltan_Next_Border_Obj_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id,
  ZOLTAN_ID_PTR local_id, int nbor_proc,
  ZOLTAN_ID_PTR next_global_id,
  ZOLTAN_ID_PTR next_local_id,
  int wdim, float *next_obj_wgt,
  int *ierr)
{
   return Zoltan_Current->Get_Next_Border_Obj_Fort(data, 
                                                  &num_gid_entries,
                                                  &num_lid_entries,
                                                  global_id, local_id,
                                                  &nbor_proc, next_global_id,
                                                  next_local_id, &wdim,
                                                  next_obj_wgt, ierr);
}

/*****************************************************************************/
int Zoltan_Obj_Size_Fort_Wrapper(void *data, int num_gid_entries,
  int num_lid_entries, ZOLTAN_ID_PTR global_id, 
  ZOLTAN_ID_PTR local_id, int *ierr)
{
   return Zoltan_Current->Get_Obj_Size_Fort(data,
             &num_gid_entries, &num_lid_entries,
             global_id, local_id, ierr);
}

/*****************************************************************************/
void Zoltan_Obj_Size_Multi_Fort_Wrapper(
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  int num_ids,
  ZOLTAN_ID_PTR global_ids,
  ZOLTAN_ID_PTR local_ids,
  int *num_bytes,
  int *ierr)
{
   Zoltan_Current->Get_Obj_Size_Multi_Fort(data,
             &num_gid_entries, &num_lid_entries, &num_ids,
             global_ids, local_ids, num_bytes, ierr);
}

/*****************************************************************************/
void Zoltan_Pre_Migrate_PP_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids, int *import_procs, int *import_to_proc,
  int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_proc,
  int *ierr)
{
   Zoltan_Current->Migrate.Pre_Migrate_PP_Fort(data, 
                   &num_gid_entries, &num_lid_entries,
                   &num_import, import_global_ids, import_local_ids, 
                   import_procs, import_to_proc,
                   &num_export, export_global_ids, export_local_ids,
                   export_procs, export_to_proc, ierr);
}

/*****************************************************************************/
void Zoltan_Mid_Migrate_PP_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids, int *import_procs, int *import_to_proc,
  int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_proc,
  int *ierr)
{
   Zoltan_Current->Migrate.Mid_Migrate_PP_Fort(data,
                   &num_gid_entries, &num_lid_entries,
                   &num_import, import_global_ids, import_local_ids,
                   import_procs, import_to_proc,
                   &num_export, export_global_ids, export_local_ids,
                   export_procs, export_to_proc, ierr);
}

/*****************************************************************************/
void Zoltan_Post_Migrate_PP_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids, int *import_procs, int *import_to_proc,
  int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_proc,
  int *ierr)
{
   Zoltan_Current->Migrate.Post_Migrate_PP_Fort(data,
                   &num_gid_entries, &num_lid_entries,
                   &num_import, import_global_ids, import_local_ids,
                   import_procs, import_to_proc,
                   &num_export, export_global_ids, export_local_ids,
                   export_procs, export_to_proc, ierr);

}

/*****************************************************************************/
void Zoltan_Pre_Migrate_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids, int *import_procs,
  int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs,
  int *ierr)
{
   Zoltan_Current->Migrate.Pre_Migrate_Fort(data, 
                                           &num_gid_entries,
                                           &num_lid_entries,
                                           &num_import,
                                           import_global_ids,
                                           import_local_ids, import_procs,
                                           &num_export, export_global_ids,
                                           export_local_ids, export_procs,
                                           ierr);
}

/*****************************************************************************/
void Zoltan_Mid_Migrate_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids, int *import_procs,
  int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs,
  int *ierr)
{
   Zoltan_Current->Migrate.Mid_Migrate_Fort(data,
                                           &num_gid_entries,
                                           &num_lid_entries,
                                           &num_import,
                                           import_global_ids,
                                           import_local_ids, import_procs,
                                           &num_export, export_global_ids,
                                           export_local_ids, export_procs,
                                           ierr);
}

/*****************************************************************************/
void Zoltan_Post_Migrate_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids, int *import_procs,
  int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs,
  int *ierr)
{
   Zoltan_Current->Migrate.Post_Migrate_Fort(data,
                                            &num_gid_entries, &num_lid_entries, 
                                            &num_import,
                                            import_global_ids,
                                            import_local_ids, import_procs,
                                            &num_export, export_global_ids,
                                            export_local_ids, export_procs,
                                            ierr);
}

/*****************************************************************************/
void Zoltan_Pack_Obj_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  int dest_proc, int size, char *buf, int *ierr)
{
   Zoltan_Current->Pack_Obj_Fort(data, 
                                        &num_gid_entries, &num_lid_entries, 
                                        global_id, local_id,
                                        &dest_proc, &size, buf, ierr);
}

/*****************************************************************************/
void Zoltan_Pack_Obj_Multi_Fort_Wrapper(
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  int num_ids,
  ZOLTAN_ID_PTR global_ids,
  ZOLTAN_ID_PTR local_ids,
  int *dest_proc,
  int *size,
  int *index,
  char *buffer,
  int *ierr)
{
  int factor = sizeof(int) / sizeof(char);
  int i;
 
  /* Convert index array from indices into char * to indices into int *. */
  /* Add 1 for F90 one-based indexing. */
  for (i = 0; i < num_ids; i++) {
    /* Sanity check */
    if (index[i] % factor != 0) {
      ZOLTAN_PRINT_ERROR(-1, "Zoltan_Pack_Obj_Multi_Fort_Wrapper", 
                         "Alignment problem in index array.");
      
      *ierr = ZOLTAN_FATAL;
      return;
    }
    index[i] = index[i]/factor + 1;
  }
     
  Zoltan_Current->Pack_Obj_Multi_Fort(data, &num_gid_entries, &num_lid_entries,
                                      &num_ids, global_ids, local_ids,
                                      dest_proc, size, index, buffer, ierr);

  /* Restore index array to original condition. */
  for (i = 0; i < num_ids; i++) 
    index[i] = (index[i] - 1) * factor;
}

/*****************************************************************************/
void Zoltan_Unpack_Obj_Fort_Wrapper(void *data, int num_gid_entries,
                                ZOLTAN_ID_PTR global_id, int size,
                                char *buf, int *ierr)
{
   Zoltan_Current->Unpack_Obj_Fort(data, &num_gid_entries, 
                                          global_id, &size, buf, ierr);
}

/*****************************************************************************/
void Zoltan_Unpack_Obj_Multi_Fort_Wrapper(
  void *data,
  int num_gid_entries,
  int num_ids,
  ZOLTAN_ID_PTR global_ids,
  int *size,
  int *index,
  char *buffer,
  int *ierr)
{
  int factor = sizeof(int) / sizeof(char);
  int i;
 
  /* Convert index array from indices into char * to indices into int *. */
  /* Add 1 for F90 one-based indexing. */
  for (i = 0; i < num_ids; i++) {
    /* Sanity check */
    if (index[i] % factor != 0) {
      ZOLTAN_PRINT_ERROR(-1, "Zoltan_Pack_Obj_Multi_Fort_Wrapper", 
                         "Alignment problem in index array.");
      
      *ierr = ZOLTAN_FATAL;
      return;
    }
    index[i] = index[i]/factor + 1;
  }
     
  Zoltan_Current->Unpack_Obj_Multi_Fort(data, &num_gid_entries, &num_ids,
                                        global_ids, size, index, buffer, ierr);

  /* Restore index array to original condition. */
  for (i = 0; i < num_ids; i++) 
    index[i] = (index[i] - 1) * factor;
}


/*****************************************************************************/
int Zoltan_Num_Coarse_Obj_Fort_Wrapper(void *data, int *ierr)
{
   return Zoltan_Current->Get_Num_Coarse_Obj_Fort(data, ierr);
}

/*****************************************************************************/
void Zoltan_Coarse_Obj_List_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_ids,
  ZOLTAN_ID_PTR local_ids, int *assigned, int *num_vert,
  ZOLTAN_ID_PTR vertices, int *in_order, ZOLTAN_ID_PTR in_vertex,
  ZOLTAN_ID_PTR out_vertex, int *ierr)
{
   Zoltan_Current->Get_Coarse_Obj_List_Fort(data, 
                                           &num_gid_entries, &num_lid_entries,
                                           global_ids, local_ids,
                                           assigned, num_vert, vertices,
                                           in_order, in_vertex, out_vertex,
                                           ierr);
}

/*****************************************************************************/
void Zoltan_HG_Size_CS_Fort_Wrapper(void *data, 
  int *num_lists, int *num_pins, int *format, int *ierr)
{
   Zoltan_Current->Get_HG_Size_CS_Fort(data, num_lists, num_pins, format,
                                           ierr);
}
/*****************************************************************************/
void Zoltan_HG_Size_Edge_Wts_Fort_Wrapper(void *data, int *num_edges, int *ierr)
{
   Zoltan_Current->Get_HG_Size_Edge_Wts_Fort(data, num_edges, ierr);
}
/*****************************************************************************/
void Zoltan_HG_CS_Fort_Wrapper(void *data, 
  int num_gid_entries, int nrowcol, int npins, int format,
  ZOLTAN_ID_PTR rowcol_GID, int *rowcol_ptr, ZOLTAN_ID_PTR pin_GID,
  int *ierr)
{
   Zoltan_Current->Get_HG_CS_Fort(data,
                       &num_gid_entries, &nrowcol, &npins, &format,
                       rowcol_GID, rowcol_ptr, pin_GID, ierr);
}
/*****************************************************************************/
void Zoltan_HG_Edge_Wts_Fort_Wrapper(void *data,
  int num_gid_entries, int num_lid_entries, int nedges, int edge_weight_dim,
  ZOLTAN_ID_PTR edge_GID, ZOLTAN_ID_PTR edge_LID, float *edge_weight,
  int *ierr)

{
   Zoltan_Current->Get_HG_Edge_Wts_Fort(data,
       &num_gid_entries, &num_lid_entries, &nedges, &edge_weight_dim,
        edge_GID, edge_LID, edge_weight, ierr);
}
/*****************************************************************************/
int Zoltan_Num_Fixed_Obj_Fort_Wrapper(void *data, int *ierr)
{
   return Zoltan_Current->Get_Num_Fixed_Obj_Fort(data, ierr);
}
/*****************************************************************************/
void Zoltan_Fixed_Obj_List_Fort_Wrapper(void *data,
  int num_fixed_obj, int num_gid_entries, 
  ZOLTAN_ID_PTR fixed_gids, int *fixed_part, int *ierr)
{
   Zoltan_Current->Get_Fixed_Obj_List_Fort(data,
              &num_fixed_obj, &num_gid_entries, 
              fixed_gids, fixed_part, ierr);
}
/*****************************************************************************/
int Zoltan_First_Coarse_Obj_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, 
  ZOLTAN_ID_PTR global_id,
  ZOLTAN_ID_PTR local_id, int *assigned,
  int *num_vert, ZOLTAN_ID_PTR vertices,
  int *in_order, ZOLTAN_ID_PTR in_vertex,
  ZOLTAN_ID_PTR out_vertex, int *ierr)
{
   return Zoltan_Current->Get_First_Coarse_Obj_Fort(data, 
                                                   &num_gid_entries, 
                                                   &num_lid_entries,
                                                   global_id, local_id,
                                                   assigned, num_vert, vertices,
                                                   in_order, in_vertex,
                                                   out_vertex, ierr);
}

/*****************************************************************************/
int Zoltan_Next_Coarse_Obj_Fort_Wrapper(void *data, int num_gid_entries, 
  int num_lid_entries, ZOLTAN_ID_PTR global_id,
  ZOLTAN_ID_PTR local_id, 
  ZOLTAN_ID_PTR next_global_id, 
  ZOLTAN_ID_PTR next_local_id,
  int *assigned,
  int *num_vert, ZOLTAN_ID_PTR vertices,
  ZOLTAN_ID_PTR in_vertex, ZOLTAN_ID_PTR out_vertex, int *ierr)
{
   return Zoltan_Current->Get_Next_Coarse_Obj_Fort(data, &num_gid_entries,
                                                  &num_lid_entries,
                                                  global_id, local_id,
                                                  next_global_id, next_local_id,
                                                  assigned, num_vert, vertices,
                                                  in_vertex, out_vertex, ierr);
}

/*****************************************************************************/
int Zoltan_Num_Child_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, 
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  int *ierr)
{
   return Zoltan_Current->Get_Num_Child_Fort(data, 
                                            &num_gid_entries, &num_lid_entries,
                                            global_id, local_id, ierr);
}

/*****************************************************************************/
void Zoltan_Child_List_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, 
  ZOLTAN_ID_PTR parent_gid,
  ZOLTAN_ID_PTR parent_lid, ZOLTAN_ID_PTR child_gids,
  ZOLTAN_ID_PTR child_lids, int *assigned,
  int *num_vert, ZOLTAN_ID_PTR vertices,
  ZOLTAN_REF_TYPE *ref_type, ZOLTAN_ID_PTR in_vertex,
  ZOLTAN_ID_PTR out_vertex, int *ierr)
{
   Zoltan_Current->Get_Child_List_Fort(data, &num_gid_entries, &num_lid_entries,
                                      parent_gid, parent_lid,
                                      child_gids, child_lids, assigned,
                                      num_vert, vertices,
                                      ref_type, in_vertex, out_vertex, ierr);
}

/*****************************************************************************/
void Zoltan_Child_Weight_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  int wgt_dim, float *obj_wgt, int *ierr)
{
   Zoltan_Current->Get_Child_Weight_Fort(data, 
                                        &num_gid_entries, &num_lid_entries,
                                        global_id, local_id, &wgt_dim,
                                        obj_wgt, ierr);
}

/*****************************************************************************/
int Zoltan_Hier_Num_Levels_Fort_Wrapper(void *data, int *ierr)
{
  return Zoltan_Current->Get_Hier_Num_Levels_Fort(data, ierr);
}

/*****************************************************************************/
int Zoltan_Hier_Part_Fort_Wrapper(void *data, int level, int *ierr)
{
  return Zoltan_Current->Get_Hier_Part_Fort(data, &level, ierr);
}

/*****************************************************************************/
void Zoltan_Hier_Method_Fort_Wrapper(void *data, int level, 
				     struct Zoltan_Struct *zz, int *ierr)
{
  int *fort_zz; /* maybe this should be void *? */
  extern ZOLTAN_FORT_MALLOC_SET_STRUCT_FN fort_malloc_set_struct;
  int zz_addr_bytes[8];
  int nbytes = 8;
  LB_TO_ADDR(zz, zz_addr_bytes, nbytes);

  /* create a Fortran Zoltan_Struct for zz */
  Zoltan_Special_Fort_Malloc_Set_Struct(zz_addr_bytes,&fort_zz);

  /* call the callback */
  Zoltan_Current->Get_Hier_Method_Fort(data, &level, fort_zz, ierr);

}

/*****************************************************************************/
/*--------------------------------------------------------------------*/
/* C wrapper functions                                                */
/*--------------------------------------------------------------------*/

/*****************************************************************************/
int Zfw_Initialize(float *ver)
{
   int myArgc;
   char **myArgv;
   int result;
   myArgc = 1;
   myArgv = (char **) ZOLTAN_MALLOC((myArgc+1)*sizeof(char *));
   myArgv[0] = "unknown";
   myArgv[1] = NULL;
   result = Zoltan_Initialize(myArgc,myArgv,ver);
   ZOLTAN_FREE(&myArgv);
   return result;
}

/*****************************************************************************/
int Zfw_Initialize1(int *argc, int *argv, int *starts, float *ver)
{
   int i, j, result;
   char **myArgv;
   myArgv = (char **) ZOLTAN_MALLOC(((*argc)+1)*sizeof(char *));
   for (i=0; i<(*argc); i++) {
      myArgv[i] = (char *) ZOLTAN_MALLOC((starts[i+1]-starts[i]+1)*sizeof(char));
      for (j=0; j<starts[i+1]-starts[i]; j++) {
         myArgv[i][j] = (char) argv[starts[i]+j-1];
      }
      myArgv[i][starts[i+1]-starts[i]] = '\0';
   }
   myArgv[*argc] = NULL;
   result = Zoltan_Initialize(*argc,myArgv,ver);
   for (i=0; i<(*argc); i++) 
     ZOLTAN_FREE(&(myArgv[i]));
   ZOLTAN_FREE(&myArgv);
   return result;
}

/*****************************************************************************/
void Zfw_Create(int *f_communicator, int *addr_lb, int *nbytes)
{
   volatile struct Zoltan_Struct *lb;
   MPI_Comm c_communicator;
   c_communicator = Zoltan_comm_f2c(f_communicator);
   lb = Zoltan_Create(c_communicator);
   lb->Fortran = 1;
   LB_TO_ADDR(lb, addr_lb, *nbytes);
}

/*****************************************************************************/
void Zfw_Copy(int *addr_lb1, int *addr_lb2, int *nbytes)
{
   volatile struct Zoltan_Struct *in, *out;

   ADDR_TO_LB(addr_lb1, in);

   out = Zoltan_Copy(in);
   out->Fortran = 1;

   LB_TO_ADDR(out, addr_lb2, *nbytes);
}

/*****************************************************************************/
int Zfw_Copy_To(int *addr_lb1, int *addr_lb2, int *nbytes)
{
   volatile struct Zoltan_Struct *to, *from;
   ADDR_TO_LB(addr_lb1, to);
   ADDR_TO_LB(addr_lb2, from);
   return Zoltan_Copy_To(to, from);
}

/*****************************************************************************/
void Zfw_Destroy(int *addr_lb, int *nbytes)
{
   volatile struct Zoltan_Struct *lb;
   ADDR_TO_LB(addr_lb, lb);
   Zoltan_Destroy(&lb);
}

/*****************************************************************************/
int Zfw_Align(int *size)
{
   return Zoltan_Align(*size);
}

/*****************************************************************************/
void Zfw_Memory_Stats()
{
   Zoltan_Memory_Stats();
}

/*****************************************************************************/
int Zfw_Set_Fn(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
               void *data)
{
   volatile struct Zoltan_Struct *lb;
   ADDR_TO_LB(addr_lb, lb);
   switch(*type) {
   case ZOLTAN_PART_MULTI_FN_TYPE:
      lb->Get_Part_Multi_Fort = (ZOLTAN_PART_MULTI_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Part_Multi_Fort_Wrapper, data);
      break;
   case ZOLTAN_PART_FN_TYPE:
      lb->Get_Part_Fort = (ZOLTAN_PART_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Part_Fort_Wrapper, data);
      break;
   case ZOLTAN_NUM_EDGES_MULTI_FN_TYPE:
      lb->Get_Num_Edges_Multi_Fort = (ZOLTAN_NUM_EDGES_MULTI_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Num_Edges_Multi_Fort_Wrapper, data);
      break;
   case ZOLTAN_EDGE_LIST_MULTI_FN_TYPE:
      lb->Get_Edge_List_Multi_Fort = (ZOLTAN_EDGE_LIST_MULTI_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Edge_List_Multi_Fort_Wrapper, data);
      break;
   case ZOLTAN_NUM_EDGES_FN_TYPE:
      lb->Get_Num_Edges_Fort = (ZOLTAN_NUM_EDGES_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Num_Edges_Fort_Wrapper, data);
      break;
   case ZOLTAN_EDGE_LIST_FN_TYPE:
      lb->Get_Edge_List_Fort = (ZOLTAN_EDGE_LIST_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Edge_List_Fort_Wrapper, data);
      break;
   case ZOLTAN_NUM_GEOM_FN_TYPE:
      lb->Get_Num_Geom_Fort = (ZOLTAN_NUM_GEOM_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Num_Geom_Fort_Wrapper, data);
      break;
   case ZOLTAN_GEOM_MULTI_FN_TYPE:
      lb->Get_Geom_Multi_Fort = (ZOLTAN_GEOM_MULTI_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Geom_Multi_Fort_Wrapper, data);
      break;
   case ZOLTAN_GEOM_FN_TYPE:
      lb->Get_Geom_Fort = (ZOLTAN_GEOM_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Geom_Fort_Wrapper, data);
      break;
   case ZOLTAN_NUM_OBJ_FN_TYPE:
      lb->Get_Num_Obj_Fort = (ZOLTAN_NUM_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Num_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_OBJ_LIST_FN_TYPE:
      lb->Get_Obj_List_Fort = (ZOLTAN_OBJ_LIST_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Obj_List_Fort_Wrapper, data);
      break;
   case ZOLTAN_FIRST_OBJ_FN_TYPE:
      lb->Get_First_Obj_Fort = (ZOLTAN_FIRST_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_First_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_NEXT_OBJ_FN_TYPE:
      lb->Get_Next_Obj_Fort = (ZOLTAN_NEXT_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Next_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_NUM_BORDER_OBJ_FN_TYPE:
      lb->Get_Num_Border_Obj_Fort = (ZOLTAN_NUM_BORDER_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Num_Border_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_BORDER_OBJ_LIST_FN_TYPE:
      lb->Get_Border_Obj_List_Fort = (ZOLTAN_BORDER_OBJ_LIST_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Border_Obj_List_Fort_Wrapper, data);
      break;
   case ZOLTAN_FIRST_BORDER_OBJ_FN_TYPE:
      lb->Get_First_Border_Obj_Fort = (ZOLTAN_FIRST_BORDER_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_First_Border_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_NEXT_BORDER_OBJ_FN_TYPE:
      lb->Get_Next_Border_Obj_Fort = (ZOLTAN_NEXT_BORDER_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Next_Border_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_PRE_MIGRATE_PP_FN_TYPE:
      lb->Migrate.Pre_Migrate_PP_Fort = (ZOLTAN_PRE_MIGRATE_PP_FORT_FN *)fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Pre_Migrate_PP_Fort_Wrapper, data);
      break;
   case ZOLTAN_MID_MIGRATE_PP_FN_TYPE:
      lb->Migrate.Mid_Migrate_PP_Fort = (ZOLTAN_MID_MIGRATE_PP_FORT_FN *)fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Mid_Migrate_PP_Fort_Wrapper, data);
      break;
   case ZOLTAN_POST_MIGRATE_PP_FN_TYPE:
      lb->Migrate.Post_Migrate_PP_Fort =(ZOLTAN_POST_MIGRATE_PP_FORT_FN*)fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Post_Migrate_PP_Fort_Wrapper, data);
      break;
   case ZOLTAN_PRE_MIGRATE_FN_TYPE:
      lb->Migrate.Pre_Migrate_Fort = (ZOLTAN_PRE_MIGRATE_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Pre_Migrate_Fort_Wrapper, data);
      break;
   case ZOLTAN_MID_MIGRATE_FN_TYPE:
      lb->Migrate.Mid_Migrate_Fort = (ZOLTAN_MID_MIGRATE_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Mid_Migrate_Fort_Wrapper, data);
      break;
   case ZOLTAN_POST_MIGRATE_FN_TYPE:
      lb->Migrate.Post_Migrate_Fort = (ZOLTAN_POST_MIGRATE_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Post_Migrate_Fort_Wrapper, data);
      break;
   case ZOLTAN_OBJ_SIZE_FN_TYPE:
      lb->Get_Obj_Size_Fort = (ZOLTAN_OBJ_SIZE_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Obj_Size_Fort_Wrapper, data);
      break;
   case ZOLTAN_PACK_OBJ_FN_TYPE:
      lb->Pack_Obj_Fort = (ZOLTAN_PACK_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Pack_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_UNPACK_OBJ_FN_TYPE:
      lb->Unpack_Obj_Fort = (ZOLTAN_UNPACK_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Unpack_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE:
      lb->Get_Obj_Size_Multi_Fort = (ZOLTAN_OBJ_SIZE_MULTI_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Obj_Size_Multi_Fort_Wrapper, data);
      break;
   case ZOLTAN_PACK_OBJ_MULTI_FN_TYPE:
      lb->Pack_Obj_Multi_Fort = (ZOLTAN_PACK_OBJ_MULTI_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Pack_Obj_Multi_Fort_Wrapper, data);
      break;
   case ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE:
      lb->Unpack_Obj_Multi_Fort = (ZOLTAN_UNPACK_OBJ_MULTI_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Unpack_Obj_Multi_Fort_Wrapper, data);
      break;
   case ZOLTAN_NUM_COARSE_OBJ_FN_TYPE:
      lb->Get_Num_Coarse_Obj_Fort = (ZOLTAN_NUM_COARSE_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Num_Coarse_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_COARSE_OBJ_LIST_FN_TYPE:
      lb->Get_Coarse_Obj_List_Fort = (ZOLTAN_COARSE_OBJ_LIST_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Coarse_Obj_List_Fort_Wrapper, data);
      break;
   case ZOLTAN_FIRST_COARSE_OBJ_FN_TYPE:
      lb->Get_First_Coarse_Obj_Fort = (ZOLTAN_FIRST_COARSE_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_First_Coarse_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_NEXT_COARSE_OBJ_FN_TYPE:
      lb->Get_Next_Coarse_Obj_Fort = (ZOLTAN_NEXT_COARSE_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Next_Coarse_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_NUM_CHILD_FN_TYPE:
      lb->Get_Num_Child_Fort = (ZOLTAN_NUM_CHILD_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Num_Child_Fort_Wrapper, data);
      break;
   case ZOLTAN_CHILD_LIST_FN_TYPE:
      lb->Get_Child_List_Fort = (ZOLTAN_CHILD_LIST_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Child_List_Fort_Wrapper, data);
      break;
   case ZOLTAN_CHILD_WEIGHT_FN_TYPE:
      lb->Get_Child_Weight_Fort = (ZOLTAN_CHILD_WEIGHT_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Child_Weight_Fort_Wrapper, data);
      break;
   case ZOLTAN_HG_SIZE_CS_FN_TYPE:
      lb->Get_HG_Size_CS_Fort = (ZOLTAN_HG_SIZE_CS_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_HG_Size_CS_Fort_Wrapper, data);
      break;
   case ZOLTAN_HG_SIZE_EDGE_WTS_FN_TYPE:
      lb->Get_HG_Size_Edge_Wts_Fort = (ZOLTAN_HG_SIZE_EDGE_WTS_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_HG_Size_Edge_Wts_Fort_Wrapper, data);
      break;
   case ZOLTAN_HG_CS_FN_TYPE:
      lb->Get_HG_CS_Fort = (ZOLTAN_HG_CS_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_HG_CS_Fort_Wrapper, data);
      break;
   case ZOLTAN_HG_EDGE_WTS_FN_TYPE:
      lb->Get_HG_Edge_Wts_Fort = (ZOLTAN_HG_EDGE_WTS_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_HG_Edge_Wts_Fort_Wrapper, data);
      break;
   case ZOLTAN_NUM_FIXED_OBJ_FN_TYPE:
      lb->Get_Num_Fixed_Obj_Fort = (ZOLTAN_NUM_FIXED_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type,
               (void (*)())Zoltan_Num_Fixed_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_FIXED_OBJ_LIST_FN_TYPE:
      lb->Get_Fixed_Obj_List_Fort = (ZOLTAN_FIXED_OBJ_LIST_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type,
               (void (*)())Zoltan_Fixed_Obj_List_Fort_Wrapper, data);
      break;
   case ZOLTAN_HIER_NUM_LEVELS_FN_TYPE:
      lb->Get_Hier_Num_Levels_Fort = (ZOLTAN_HIER_NUM_LEVELS_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Hier_Num_Levels_Fort_Wrapper, data);
      break;
   case ZOLTAN_HIER_PART_FN_TYPE:
      lb->Get_Hier_Part_Fort = (ZOLTAN_HIER_PART_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Hier_Part_Fort_Wrapper, data);
      break;
   case ZOLTAN_HIER_METHOD_FN_TYPE:
      lb->Get_Hier_Method_Fort = (ZOLTAN_HIER_METHOD_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Hier_Method_Fort_Wrapper, data);
      break;
   default:
      return Zoltan_Set_Fn(lb, *type, (void (*)())NULL, data);
      break;
   }
}

/*****************************************************************************/
int Zfw_Set_Fn0f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)())
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)NULL);
}

/*****************************************************************************/
int Zfw_Set_Fn1f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                  int *data)
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn2f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 float *data)
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn3f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 double *data)
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn4f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_1) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn5f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_2) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn6f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_3) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn7f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_4) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn8f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_1) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn9f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_2) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_FnAf(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_3) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_FnBf(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_4) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn0s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)())
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)NULL);
}

/*****************************************************************************/
int Zfw_Set_Fn1s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                  int *data)
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn2s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                  float *data)
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn3s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                  double *data)
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn4s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_1) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn5s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_2) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn6s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_3) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn7s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_4) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn8s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_1) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn9s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_2) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_FnAs(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_3) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_FnBs(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_4) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Param(int *addr_lb, int *nbytes, int *int_param_name,
                   int *param_name_len, int *int_new_value, int *new_value_len)
{
   volatile struct Zoltan_Struct *lb;
   char *param_name, *new_value;
   int i, result;
   param_name = (char *)ZOLTAN_MALLOC(*param_name_len+1);
   new_value = (char *)ZOLTAN_MALLOC(*new_value_len+1);
   ADDR_TO_LB(addr_lb, lb);
   for (i=0; i<(*param_name_len); i++) param_name[i] = (char)int_param_name[i];
   param_name[*param_name_len] = '\0';
   for (i=0; i<(*new_value_len); i++) new_value[i] = (char)int_new_value[i];
   new_value[*new_value_len] = '\0';
   result = Zoltan_Set_Param(lb, param_name, new_value);
   ZOLTAN_FREE(&param_name);
   ZOLTAN_FREE(&new_value);
   return result;
}

/*****************************************************************************/
int Zfw_Set_Param_Vec(int *addr_lb, int *nbytes, int *int_param_name,
                   int *param_name_len, int *int_new_value, int *new_value_len,
                   int index)
{
   volatile struct Zoltan_Struct *lb;
   char *param_name, *new_value;
   int i, result;
   param_name = (char *)ZOLTAN_MALLOC(*param_name_len+1);
   new_value = (char *)ZOLTAN_MALLOC(*new_value_len+1);
   ADDR_TO_LB(addr_lb, lb);
   for (i=0; i<(*param_name_len); i++) param_name[i] = (char)int_param_name[i];
   param_name[*param_name_len] = '\0';
   for (i=0; i<(*new_value_len); i++) new_value[i] = (char)int_new_value[i];
   new_value[*new_value_len] = '\0';
   result = Zoltan_Set_Param_Vec(lb, param_name, new_value, index);
   ZOLTAN_FREE(&param_name);
   ZOLTAN_FREE(&new_value);
   return result;
}

/*****************************************************************************/
int Zfw_LB_Partition(int *addr_lb, int *nbytes, int *changes, 
  int *num_gid_entries, int *num_lid_entries,
  int *num_import,
  ZOLTAN_ID_PTR *import_global_ids, ZOLTAN_ID_PTR *import_local_ids,
  int **import_procs, int **import_to_part, int *num_export,
  ZOLTAN_ID_PTR *export_global_ids, ZOLTAN_ID_PTR *export_local_ids,
  int **export_procs, int **export_to_part
#ifdef PGI
/* PGI uses hidden arguments when it passes pointers */
   ,int *imp_gid_hide, int *imp_lid_hide, int *imp_proc_hide,
    int *imp_to_part_hide,
    int *exp_gid_hide, int *exp_lid_hide, int *exp_proc_hide,
    int *exp_to_part_hide
#endif
#ifdef FUJITSU
/* Fujitsu and Lahey use a hidden argument for every argument */
/* TEMP need to verify this with Fujitsu or Lahey */
   ,int *addr_lb_hide, int *nbytes_hide, int *changes_hide,
    int *num_gid_entries_hide, int *num_lid_entries_hide,
    int *num_import_hide, int *imp_gid_hide, int *imp_lid_hide,
    int *imp_proc_hide, int *imp_to_part_hide,
    int *num_export_hide, int *exp_gid_hide,
    int *exp_lid_hide, int *exp_proc_hide, int *exp_to_part_hide
#endif
)
{
   volatile struct Zoltan_Struct *lb;
#if defined (PGI) || defined (FUJITSU)
#define F90LB_TEMP 3
#else
#define F90LB_TEMP 2
#endif
   ZOLTAN_ID_PTR temp_imp_gid[F90LB_TEMP], temp_exp_gid[F90LB_TEMP];
   ZOLTAN_ID_PTR temp_imp_lid[F90LB_TEMP], temp_exp_lid[F90LB_TEMP];
   int *temp_imp_proc[F90LB_TEMP], *temp_exp_proc[F90LB_TEMP];
   int *temp_imp_to_part[F90LB_TEMP], *temp_exp_to_part[F90LB_TEMP];
#undef F90LB_TEMP

/* reconstruct the lb pointer from the nbyte 1-byte integers in addr_lb */

   ADDR_TO_LB(addr_lb, lb);
   Zoltan_Current = lb;

/* put the address of the Fortran pointer into temp_*[1] to be passed to
   Fortran for allocation.  The address of the allocated space will be
   in temp_*[0] so it can be used by C without messing up the Fortran pointer*/

   temp_imp_gid[1] = (ZOLTAN_ID_PTR)import_global_ids;
   temp_imp_lid[1] = (ZOLTAN_ID_PTR)import_local_ids;
   temp_imp_proc[1] = (int *)import_procs;
   temp_imp_to_part[1] = (int *)import_to_part;
   temp_exp_gid[1] = (ZOLTAN_ID_PTR)export_global_ids;
   temp_exp_lid[1] = (ZOLTAN_ID_PTR)export_local_ids;
   temp_exp_proc[1] = (int *)export_procs;
   temp_exp_to_part[1] = (int *)export_to_part;

/* for PGI and FUJITSU, put the hidden argument in temp_*[2] */

#if defined (PGI) || defined (FUJITSU)
   temp_imp_gid[2] = (ZOLTAN_ID_PTR)imp_gid_hide;
   temp_imp_lid[2] = (ZOLTAN_ID_PTR)imp_lid_hide;
   temp_imp_proc[2] = (int *)imp_proc_hide;
   temp_imp_to_part[2] = (int *)imp_to_part_hide;
   temp_exp_gid[2] = (ZOLTAN_ID_PTR)exp_gid_hide;
   temp_exp_lid[2] = (ZOLTAN_ID_PTR)exp_lid_hide;
   temp_exp_proc[2] = (int *)exp_proc_hide;
   temp_exp_to_part[2] = (int *)exp_to_part_hide;
#endif

/* call Zoltan_LB_Partition */

   return Zoltan_LB_Partition(lb, changes, num_gid_entries, num_lid_entries, 
                     num_import, temp_imp_gid, temp_imp_lid,
                     temp_imp_proc, temp_imp_to_part,
                     num_export, temp_exp_gid, temp_exp_lid,
                     temp_exp_proc, temp_exp_to_part);
}

/*****************************************************************************/
int Zfw_LB_Eval(int *addr_lb, int *nbytes, int *print_stats)
{
/* KDD Note that we are not yet ready to support the return arguments
 * KDD of Zoltan_LB_Eval in the F90 interface yet.  For now, we'll only
 * KDD print stats.
 */
   volatile struct Zoltan_Struct *lb;
   ADDR_TO_LB(addr_lb, lb);
   Zoltan_Current = lb;

   return  Zoltan_LB_Eval(lb, *print_stats, NULL, NULL, NULL);
}

/*****************************************************************************/
int Zfw_LB_Set_Part_Sizes(int *addr_lb, int *nbytes, int *global_part, int *len,
                          int *partids, int *wgtidx, float *partsizes)
{
   volatile struct Zoltan_Struct *lb;
   ADDR_TO_LB(addr_lb, lb);

   return Zoltan_LB_Set_Part_Sizes(lb, *global_part, *len, partids, wgtidx,
                                   partsizes);
}

/*****************************************************************************/
int Zfw_LB_Point_Assign(int *addr_lb, int *nbytes, double *coords, int *proc)
{
   volatile struct Zoltan_Struct *lb;
   ADDR_TO_LB(addr_lb, lb);

   return Zoltan_LB_Point_Assign(lb, coords, proc);
}

/*****************************************************************************/
int Zfw_LB_Point_PP_Assign(int *addr_lb, int *nbytes, double *coords, int *proc,
                           int *part)
{
   volatile struct Zoltan_Struct *lb;
   ADDR_TO_LB(addr_lb, lb);

   return Zoltan_LB_Point_PP_Assign(lb, coords, proc, part);
}

/*****************************************************************************/
int Zfw_LB_Box_Assign(int *addr_lb, int *nbytes, double *xmin, double *ymin,
                     double *zmin, double *xmax, double *ymax, double *zmax,
                     int *procs, int *numprocs)
{
   volatile struct Zoltan_Struct *lb;
   ADDR_TO_LB(addr_lb, lb);

   return Zoltan_LB_Box_Assign(lb, *xmin, *ymin, *zmin, *xmax, *ymax, *zmax, 
                               procs, numprocs);
}

/*****************************************************************************/
int Zfw_LB_Box_PP_Assign(int *addr_lb, int *nbytes, double *xmin, double *ymin,
                     double *zmin, double *xmax, double *ymax, double *zmax,
                     int *procs, int *numprocs, int *parts, int *numparts)
{
   volatile struct Zoltan_Struct *lb;
   ADDR_TO_LB(addr_lb, lb);

   return Zoltan_LB_Box_PP_Assign(lb, *xmin, *ymin, *zmin, *xmax, *ymax, *zmax,
                                  procs, numprocs, parts, numparts);
}

/*****************************************************************************/
int Zfw_Invert_Lists(int *addr_lb, int *nbytes, 
  int *num_input,
  ZOLTAN_ID_PTR input_global_ids, ZOLTAN_ID_PTR input_local_ids,
  int *input_procs, int *input_to_part, int *num_output,
  ZOLTAN_ID_PTR *output_global_ids, ZOLTAN_ID_PTR *output_local_ids,
  int **output_procs, int **output_to_part
#ifdef PGI
  ,int *output_gid_hide, int *output_lid_hide, int *output_proc_hide, 
   int *output_to_part_hide
#endif
#ifdef FUJITSU
 ,int *addr_lb_hide, int *nbytes_hide,
  int *num_input_hide,
  int *input_global_ids_hide, int *input_local_ids_hide,
  int *input_procs_hide, int *input_to_part_hide,
  int *num_output_hide,
  int *output_gid_hide, int *output_lid_hide, 
  int *output_proc_hide, int *output_to_part_hide
#endif
)
{
   volatile struct Zoltan_Struct *lb;
#if defined (PGI) || defined(FUJITSU)
#define F90LB_TEMP 3
#else
#define F90LB_TEMP 2
#endif
   ZOLTAN_ID_PTR temp_output_gid[F90LB_TEMP];
   ZOLTAN_ID_PTR temp_output_lid[F90LB_TEMP];
   int *temp_output_proc[F90LB_TEMP];
   int *temp_output_to_part[F90LB_TEMP];
#undef F90LB_TEMP

/* reconstruct the lb pointer from the nbyte 1-byte integers in addr_lb */
   ADDR_TO_LB(addr_lb, lb);

/* put the address of the Fortran pointer into temp_*[1] to be passed to
   Fortran for allocation.  The address of the allocated space will be
   in temp_*[0] so it can be used by C without messing up the Fortran pointer*/

   temp_output_gid[1] = (ZOLTAN_ID_PTR)output_global_ids;
   temp_output_lid[1] = (ZOLTAN_ID_PTR)output_local_ids;
   temp_output_proc[1] = (int *)output_procs;
   temp_output_to_part[1] = (int *)output_to_part;

/* for PGI and FUJITSU, put the hidden argument in temp_*[2] */

#if defined (PGI) || defined(FUJITSU)
   temp_output_gid[2] = (ZOLTAN_ID_PTR)output_gid_hide;
   temp_output_lid[2] = (ZOLTAN_ID_PTR)output_lid_hide;
   temp_output_proc[2] = (int *)output_proc_hide;
   temp_output_to_part[2] = (int *)output_to_part_hide;
#endif

/* call Zoltan_Invert_Lists */

   return Zoltan_Invert_Lists(lb, 
                     *num_input, input_global_ids,
                     input_local_ids, input_procs, input_to_part,
                     num_output, temp_output_gid, temp_output_lid,
                     temp_output_proc, temp_output_to_part);
}

/*****************************************************************************/
int Zfw_Compute_Destinations(int *addr_lb, int *nbytes, 
  int *num_input,
  ZOLTAN_ID_PTR input_global_ids, ZOLTAN_ID_PTR input_local_ids,
  int *input_procs, int *num_output,
  ZOLTAN_ID_PTR *output_global_ids, ZOLTAN_ID_PTR *output_local_ids,
  int **output_procs
#ifdef PGI
  ,int *output_gid_hide, int *output_lid_hide, int *output_proc_hide
#endif
#ifdef FUJITSU
 ,int *addr_lb_hide, int *nbytes_hide,
  int *num_input_hide,
  int *input_global_ids_hide, int *input_local_ids_hide,
  int *input_procs_hide, int *num_output_hide,
  int *output_gid_hide, int *output_lid_hide, int *output_proc_hide
#endif
)
{
   volatile struct Zoltan_Struct *lb;
#if defined (PGI) || defined(FUJITSU)
#define F90LB_TEMP 3
#else
#define F90LB_TEMP 2
#endif
   ZOLTAN_ID_PTR temp_output_gid[F90LB_TEMP];
   ZOLTAN_ID_PTR temp_output_lid[F90LB_TEMP];
   int *temp_output_proc[F90LB_TEMP];
#undef F90LB_TEMP

/* reconstruct the lb pointer from the nbyte 1-byte integers in addr_lb */
   ADDR_TO_LB(addr_lb, lb);

/* put the address of the Fortran pointer into temp_*[1] to be passed to
   Fortran for allocation.  The address of the allocated space will be
   in temp_*[0] so it can be used by C without messing up the Fortran pointer*/

   temp_output_gid[1] = (ZOLTAN_ID_PTR)output_global_ids;
   temp_output_lid[1] = (ZOLTAN_ID_PTR)output_local_ids;
   temp_output_proc[1] = (int *)output_procs;

/* for PGI and FUJITSU, put the hidden argument in temp_*[2] */

#if defined (PGI) || defined(FUJITSU)
   temp_output_gid[2] = (ZOLTAN_ID_PTR)output_gid_hide;
   temp_output_lid[2] = (ZOLTAN_ID_PTR)output_lid_hide;
   temp_output_proc[2] = (int *)output_proc_hide;
#endif

/* call Zoltan_Compute_Destinations */

   return Zoltan_Compute_Destinations(lb, 
                     *num_input, input_global_ids,
                     input_local_ids, input_procs, 
                     num_output, temp_output_gid, temp_output_lid,
                     temp_output_proc);
}


/*****************************************************************************/
int Zfw_Migrate(int *addr_lb, int *nbytes, 
 int *num_import,
 ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
 int *import_procs, int *import_to_part, int *num_export,
 ZOLTAN_ID_PTR export_global_ids, ZOLTAN_ID_PTR export_local_ids,
 int *export_procs, int *export_to_part)
{
   volatile struct Zoltan_Struct *lb;
   ADDR_TO_LB(addr_lb, lb);

   Zoltan_Current = lb;
   return Zoltan_Migrate(lb,
                         *num_import,import_global_ids,import_local_ids,
                         import_procs,import_to_part,
                         *num_export,export_global_ids,
                         export_local_ids,export_procs,export_to_part);
}

/*****************************************************************************/
int Zfw_Help_Migrate(int *addr_lb, int *nbytes, 
 int *num_import,
 ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
 int *import_procs, int *num_export,
 ZOLTAN_ID_PTR export_global_ids, ZOLTAN_ID_PTR export_local_ids,
 int *export_procs)
{
   volatile struct Zoltan_Struct *lb;
   ADDR_TO_LB(addr_lb, lb);

   Zoltan_Current = lb;
   return Zoltan_Help_Migrate(lb,
                          *num_import,import_global_ids,import_local_ids,
                          import_procs,*num_export,export_global_ids,
                          export_local_ids,export_procs);
}

/*****************************************************************************/
int Zfw_Order(
 int *addr_lb, int *nbytes,
 int *num_gid_entries,
 int *num_obj,
 ZOLTAN_ID_PTR gids, 
 ZOLTAN_ID_PTR perm)
{
   volatile struct Zoltan_Struct *lb;
   int ierr;
   ADDR_TO_LB(addr_lb, lb);

   Zoltan_Current = lb;
   ierr = Zoltan_Order(lb,*num_gid_entries,*num_obj,
                       gids, perm);
   return ierr;
}

/*****************************************************************************/
int Zfw_Color(
 int *addr_lb, int *nbytes,
 int *num_gid_entries, 
 int *num_obj,
 ZOLTAN_ID_PTR gids,
 int *color_exp)
{
   volatile struct Zoltan_Struct *lb;
   int ierr;
   ADDR_TO_LB(addr_lb, lb);

   Zoltan_Current = lb;
   ierr = Zoltan_Color(lb,*num_gid_entries,*num_obj,
                       gids, color_exp);
   return ierr;
}

/*****************************************************************************/
int Zfw_Color_Test(
 int *addr_lb, int *nbytes,
 int *num_gid_entries, int *num_lid_entries,
 int *num_obj,
 ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
 int *color_exp)
{
   volatile struct Zoltan_Struct *lb;
   int ierr;
   ADDR_TO_LB(addr_lb, lb);

   Zoltan_Current = lb;
   ierr = Zoltan_Color_Test(lb,num_gid_entries,num_lid_entries,*num_obj,
                       gids, lids, color_exp);
   return ierr;
}

/*****************************************************************************/
int Zfw_Generate_Files(int *addr_lb, int *nbytes, int *int_filename,
                   int *filename_len, int *base_index, int *gen_geom,
                   int *gen_graph, int *gen_hg)
{
   volatile struct Zoltan_Struct *lb;
   char *filename;
   int i, result;
   filename = (char *)ZOLTAN_MALLOC(*filename_len+1);
   ADDR_TO_LB(addr_lb, lb);

   Zoltan_Current = lb;
   for (i=0; i<(*filename_len); i++) filename[i] = (char)int_filename[i];
   filename[*filename_len] = '\0';
   result = Zoltan_Generate_Files(lb, filename, *base_index, *gen_geom,
                                  *gen_graph, *gen_hg);
   ZOLTAN_FREE(&filename);
   return result;
}

/*****************************************************************************/
int Zfw_RCB_Box(int *addr_lb, int *nbytes, int *part, int *ndim,
                double *xmin, double *ymin, double *zmin, 
                double *xmax, double *ymax, double *zmax)
{
   volatile struct Zoltan_Struct *lb;
   ADDR_TO_LB(addr_lb, lb);

   return Zoltan_RCB_Box(lb, *part, ndim, xmin, ymin, zmin, xmax, ymax, zmax);
}
/*****************************************************************************/
void Zfw_Register_Fort_Malloc(ZOLTAN_FORT_MALLOC_INT_FN *fort_malloc_int,
			      ZOLTAN_FORT_FREE_INT_FN *fort_free_int,
			      ZOLTAN_FORT_MALLOC_SET_STRUCT_FN *fort_malloc_set_struct)
{
   Zoltan_Register_Fort_Malloc(fort_malloc_int,fort_free_int,
			       fort_malloc_set_struct);
}

/*****************************************************************************/
void Zfw_Reftree_Get_Child_Order(
  int *addr_lb, 
  int *nbytes, 
  int *order, 
  int *ierr)
{
   volatile struct Zoltan_Struct *lb;
   ADDR_TO_LB(addr_lb, lb);

   Zoltan_Reftree_Get_Child_Order(lb,order,ierr);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
