// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

int Zoltan_LocalHSFC_Order(
	ZZ *zz,               /* Zoltan structure */
	int num_obj,          /* Number of (local) objects to order. */
	ZOLTAN_ID_PTR gids,   /* List of global ids (local to this proc) */
	                      /* The application must allocate enough space */
        ZOLTAN_ID_PTR lids,   /* List of local ids (local to this proc) */
	   		      /* The application must allocate enough space */
	ZOLTAN_ID_PTR rank,      /* rank[i] is the rank of gids[i] */
	int *iperm,
	ZOOS *order_opt       /* Ordering options, parsed by Zoltan_Order */
			   );
