/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
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
 *     * Neither the name of Sandia Corporation nor the names of its
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
/* mdcgi - multiple simultaneous devices routines for cgi  */
#include <stdio.h>
#include "stdtyp.h"
#include "fortyp.h"
#include "cgi.h"
#include "ifdefx.h"
#include "mdcgi.h"
/******************************************************************************/
/*									      */
/*	Global variables						      */
/*									      */
/******************************************************************************/

/* these are shared with sdcgi.c, and are defined there */

extern device_struct	devices[MAX_DEVICES];	/* list of currently active */
extern anything	*in_params [MAX_IN_PARAMS];	/* params sent to driver */
extern short	num_devices;		/* how many items in devices*/
extern anything	*sol_surf;		/* current solicitation surface, */

/* these aren't shared with anybody */
static int	num_oldest = 0;	/* # surfs in oldest_surfs */
static anything	*oldest_surfs [MAX_SURFACES];	/* states of oldest surfaces */

/******************************************************************************/
/*									      */
/*	Added functions to allow for multiple simultaneous devices	      */
/*									      */
/******************************************************************************/



/******************************************************************************/
/*									      */
/*	xcoon - turn output on for a cgi display surface.		      */
/*   Note: it is not an error to turn output on for a surface which is	      */
/*	already on.							      */
/*   Note: it is an error to turn output on for a surface which is not	      */
/*	initialized. 							      */
/*									      */
/******************************************************************************/
void xcoon_ (anything **surface_id) /* which surface to turn output on for*/
{
   /* does surface_id point to a valid surface? */
      /* if not, do error action */
   /* find out which device this surface is on */
   /* rearrange the surface list for that device */

   short	dev;			/* which device to look at now */
   short	dev_found = 0;		/* which device was it found on */
   short	surf;			/* which surface on device to look at */
   short	surf_found = 0;		/* which active_surface was found */
   anything	*temp;			/* for swap */

   /* search active devices for this surface */
   dev_found = -1;
   for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
      for (surf = 0; surf < devices [dev]. num_active_surfaces; ++surf) {
         if (*surface_id == devices[dev]. statelist [surf]) {
            dev_found = dev; surf_found = surf; break;
         } /* end if found on list */
      } /* end for */
   } /* end for all devices */

   /* does surface_id point to a valid surface? */
   if (dev_found < 0) {
      fprintf (stderr, "xcoon: surface id not a valid surface\n");
      return;
   } /* end if surface not found */

   /* if surface was off */
   if (surf_found >= devices [dev_found]. num_on_surfaces) { 
      ++devices [dev_found]. num_on_surfaces;
      /* swap target surface with first off surface if there is one */
      if (devices [dev_found]. num_on_surfaces !=
          devices [dev_found]. num_active_surfaces) {
         temp = devices [dev_found]. statelist [surf_found];
         devices [dev_found]. statelist [surf_found] =
                      devices [dev_found].
                       statelist [devices [dev_found].num_on_surfaces-1];
         devices [dev_found].
                      statelist [devices [dev_found].num_on_surfaces-1] =
                      temp;
      } /* end if there is an off surface */
   } /* end if surface was off */
} /* end xcoon */



/******************************************************************************/
/*									      */
/*	xcact - initialize and activate a cgi display surface		      */
/*									      */
/******************************************************************************/
void xcact_ (void (*device_fn)(), anything **p_surface_id)
{
   short	arg1;			/* arg passed to device routine */
   short	i;
   anything	*temp_surface[1];	/* for calling driver */
   short	which_device;		/* index of this device in devices */
   short	which_surface;		/* index of surface in devices*/


   /* is device already initialized? */
   which_device = -1;
   for (i = 0; i < num_devices; ++i) {
      if (device_fn == devices [i]. device_fn) {
         which_device = i; break;
      } /* end if */
   } /* end for */

   if (which_device < 0) {			/* if device not initialized */
      if (num_devices >= MAX_DEVICES) {	/* if no room */
         fprintf (stderr,
                  "xcact: can't activate: too many initialized devices\n");
         *p_surface_id = NULL;
         return;
      } /* end if no room for device */
   } /* end if device not initialized */

   /* call the device driver with ACTIVATE, so it can allocate a state list */
   arg1 = ACTIVATE_FN;
   in_params [0] = (anything *)&arg1;
   (*device_fn) (in_params, 1, temp_surface);

   if (temp_surface [0] == NULL) {	/* error */
      fprintf (stderr, "xcact: can't activate surface: driver error\n");
      *p_surface_id = NULL;
      return;
   } /* end if error on activate */

   if (which_device < 0) {		/* if device not initialized */
      /* then initialize it */
      which_device = num_devices;
      ++num_devices;
      devices [which_device]. num_active_surfaces = 0;
      devices [which_device]. num_on_surfaces = 0;
      devices [which_device]. device_fn = device_fn;
   } /* end if */

   /* add new surface to device */
   which_surface = devices [which_device]. num_active_surfaces;
   ++devices [which_device]. num_active_surfaces;
   devices [which_device]. statelist [which_surface] = temp_surface [0];

   /* if new surface is the oldest surface, make it the solicitation surface */
   if (num_oldest == 0) sol_surf = temp_surface [0];

   /* put this surface as newest surf on oldest_surfs list */
   oldest_surfs [num_oldest] = temp_surface [0];
   ++num_oldest;

   /* return id of new surface */
   *p_surface_id = temp_surface [0];

} /* end xcact */



/******************************************************************************/
/*									      */
/*	xcooff - turn output off for a display surface			      */
/*   Note: it is not an error to turn off a surface which is		      */
/*	already off.							      */
/*   Note: it is an error to turn off a surface which is not active.	      */
/*									      */
/******************************************************************************/
void xcooff_ (anything **surface_id)
{
   /* does surface_id point to a valid surface? */
      /* if not, do error action */
   /* find out which device this surface is on */
   /* rearrange the surface list for that device */

   short	dev;			/* which device to look at now */
   short	dev_found;		/* which device was it found on */
   short	surf;			/* which surface on device to look at */
   short	surf_found;		/* which active_surface was found */
   anything	*temp;			/* used for swap */


   /* search devices for this surface */
   dev_found = -1;
   for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
      for (surf = 0; surf < devices [dev]. num_active_surfaces; ++surf) {
         if (*surface_id == devices[dev]. statelist [surf]) {
            dev_found = dev; surf_found = surf; break;
         } /* end if found on list */
      } /* end for */
   } /* end for all devices */

   /* does surface_id point to a valid surface? */
   if (dev_found < 0) {
      fprintf (stderr, "xcooff: surface id not a valid surface\n");
      return;
   } /* end if surface not found */

   /* if output for surface was on */
   if (surf_found < devices [dev_found]. num_on_surfaces) {	
      --devices [dev_found]. num_on_surfaces;
      /* swap target with last on surface.  Since num_on_surfaces has */
      /* been decremented, this slot is now first off surface */
      temp = devices [dev_found]. statelist [surf_found];
      devices [dev_found]. statelist [surf_found] =
         devices [dev_found]. statelist
                             [devices[dev_found]. num_on_surfaces];
      devices [dev_found]. statelist
                             [devices[dev_found]. num_on_surfaces] = temp;
   } /* end if surface was on */
} /* end xcooff */




/******************************************************************************/
/*									      */
/*	xcdact - deactivate a cgi display surface			      */
/*	Note: it is an error to deactivate a surface which is not active.     */
/*									      */
/******************************************************************************/
void xcdact_ (anything **surface_id)
{
   /* does surface_id point to a valid surface? */
      /* if not, do error action */
   /* call the device driver to allow it to deallocate the state list */
   /* remove the surface from the device list */
   /* add the surface to the free surface list */
   /* if there are no more surfaces for this device, remove the device */

   short	arg1;			/* arg passed to device routine */
   short	dev;			/* which device to look at now */
   short	dev_found;		/* which device was it found on */
   short	hole;			/* index in device of hole to fill */
   short	i;
   short	old_ptr=0;		/* where in oldest_surfs is this one */
   short	surf;			/* which surface on device to look at */
   short	surf_found=0;		/* which surface was it */
   anything	*temp_surface[1];	/* for calling driver */


   /* search devices for this surface */
   dev_found = -1;
   for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
      for (surf = 0; surf < devices [dev]. num_active_surfaces; ++surf) {
         if (*surface_id == devices[dev]. statelist [surf]) {
            dev_found = dev; surf_found = surf;
            break;
         } /* end if found on list */
      } /* end for */
   } /* end for all devices */

   /* does surface_id point to a valid surface? */
   if (dev_found < 0) {
      fprintf (stderr, "xcdact: surface id not a valid surface\n");
      return;
   } /* end if surface not found */

   /* call the device driver to allow it to deallocate the state list */
   arg1 = DEACTIVATE_FN;
   in_params [0] = (anything *)&arg1;
   temp_surface [0] = devices[dev_found]. statelist [surf_found];
   (*devices [dev_found]. device_fn) (in_params, 1, temp_surface);

   /* remove the surface from the device list */

   /* if more than 1 active surface on this device */
   if (devices [dev_found]. num_active_surfaces > 1) {

      /* if surf_found is not off or last on surface for this device */
         /* move last on surface to hole left by indexed surface */
      if (surf_found < (devices [dev_found]. num_on_surfaces - 1)) {
         devices [dev_found]. statelist [surf_found] =
                   devices [dev_found]. statelist
                                [devices [dev_found]. num_on_surfaces - 1];
      } /* end if surf_found not off or last on surface */

      /* move last off surface to hole left by either last active */
      /* surface or this surface if surface was off.  */
      /* if no off surfaces, this moves last on surface onto itself */
      hole = (surf_found < devices [dev_found]. num_on_surfaces) ?
                 devices [dev_found]. num_on_surfaces-1 : surf_found;
      devices [dev_found]. statelist [hole] =
         devices [dev_found].
              statelist [devices [dev_found]. num_active_surfaces - 1];

   } /* end if more than 1 active surface */

   /* if surface was on, now there is one less on */
   if (surf_found < devices [dev_found]. num_on_surfaces) {
      --devices [dev_found]. num_on_surfaces;
   } /* end if surface was on */

   /* clean up so things look nice while debugging (no necessary function) */
   devices [dev_found].
       statelist [devices [dev_found]. num_active_surfaces - 1] = NULL;

   /* no list manipulating necessary, just change count */
   --devices [dev_found]. num_active_surfaces;

   /* if no more surfaces on this device, remove the device */
   if (devices [dev_found]. num_active_surfaces == 0) {
      /* move last device to hole left by this device */
      /* if no more devices, this just moves last device onto itself */
      devices [dev_found] = devices [num_devices - 1];
      --num_devices;
   } /* end if no more surfaces */

   /* remove this surface from the list of oldest surfaces */
   /* first, find the right surface index in oldest_surfs */
   for (i = 0; i < num_oldest; ++i) {
      if (oldest_surfs [i] == *surface_id) {
         old_ptr = i; break;
      } /* end if */
   } /* end for */
   /* then compress the list starting at that point */
   for (i = old_ptr; i < num_oldest-1; ++i) {
      oldest_surfs [i] = oldest_surfs [i+1];
   } /* end if */
   --num_oldest;
   oldest_surfs [num_oldest] = NULL;

   /* if surface was the solicitation surface, set the solicitation surface */
   /* to the oldest surface still initialized */
   if (*surface_id == sol_surf) {
      sol_surf = oldest_surfs [0];
   } /* end if */

} /* end xcdact */



/******************************************************************************/
/*									      */
/*	xcsol - set solicitation surface				      */
/*									      */
/******************************************************************************/
void xcsol_ (anything **surface_id)
{
   /* does surface_id point to a valid surface? */
      /* if not, do error action */
   /* set the solicitation surface to surface_id */

   short	dev;			/* which device to look at now */
   short	dev_found;		/* which device was it found on */
   short	surf;			/* which surface on device to look at */
   short	surf_found;		/* which active_surface was found */


   /* search devices for this surface */
   dev_found = -1;
   for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
      for (surf = 0; surf < devices [dev]. num_active_surfaces; ++surf) {
         if (*surface_id == devices[dev]. statelist [surf]) {
            dev_found = dev; surf_found = surf; break;
         } /* end if found on list */
      } /* end for */
   } /* end for all devices */

   /* does surface_id point to a valid surface? */
   if (dev_found < 0) {
      fprintf (stderr, "xcsol: surface id not a valid surface\n");
      return;
   } /* end if surface not found */

   sol_surf = *surface_id;

} /* end xcsol */

