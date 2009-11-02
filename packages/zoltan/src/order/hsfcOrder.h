/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 *****************************************************************************/

int Zoltan_LocalHSFC_Order(
	ZZ *zz,               /* Zoltan structure */
	int num_obj,          /* Number of (local) objects to order. */
	ZOLTAN_ID_PTR gids,   /* List of global ids (local to this proc) */
	                      /* The application must allocate enough space */
        ZOLTAN_ID_PTR lids,   /* List of local ids (local to this proc) */
	   		      /* The application must allocate enough space */
	int *rank,            /* rank[i] is the rank of gids[i] */
	int *iperm,
	ZOOS *order_opt       /* Ordering options, parsed by Zoltan_Order */
			   );

