/*
//@HEADER
// *****************************************************************************
//                        Adelus
//
// Copyright 2020 NTESS and the Adelus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER
*/

#ifndef __ADELUS_MACROS_H__
#define __ADELUS_MACROS_H__

//#define grey_c(P)     ((P)^((P)>>1))

#define lrow_to_grow(R) ( (mesh_row(me) + nprocs_col*(R))  )

#define grow_to_lrow(R) ( (R/nprocs_col)  )

//// #define col_owner(C)  (((C)%nprocs_row) + (me - me%nprocs_row))
#define col_owner(C)  ( proc_num(mesh_row(me) , (C)%nprocs_row) )

//// #define row_owner(R)  ((((R)%nprocs_col)*nprocs_row) + (me%nprocs_row))
#define row_owner(R)  ( proc_num((R)%nprocs_col , mesh_col(me)) )

//#define owner(R, C)   ((((R)%nprocs_col)*nprocs_row) + ((C)%nprocs_row))

#define mesh_col(P)   ((P)%nprocs_row)

#define mesh_row(P)   ((P)/nprocs_row)

#define proc_num(R,C) ((R)*nprocs_row + (C))

//#define mac_send_msg(D,B,S,T)  MPI_Send(B,S,MPI_CHAR,D,T,MPI_COMM_WORLD)

#endif
