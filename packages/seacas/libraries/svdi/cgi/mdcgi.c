/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/* mdcgi - multiple simultaneous devices routines for cgi  */
#include "mdcgi.h"
#include "stdtyp.h"
#include <stdio.h> // for fprintf, stderr, NULL
/******************************************************************************/
/*                                                                            */
/*      Global variables                                                      */
/*                                                                            */
/******************************************************************************/

/* these are shared with sdcgi.c, and are defined there */

extern device_struct devices[MAX_DEVICES];     /* list of currently active */
extern anything *    in_params[MAX_IN_PARAMS]; /* params sent to driver */
extern short         num_devices;              /* how many items in devices*/
extern anything *    sol_surf;                 /* current solicitation surface, */

/* these aren't shared with anybody */
static int       num_oldest = 0;             /* # surfs in oldest_surfs */
static anything *oldest_surfs[MAX_SURFACES]; /* states of oldest surfaces */

/******************************************************************************/
/*                                                                            */
/*      Added functions to allow for multiple simultaneous devices            */
/*                                                                            */
/******************************************************************************/

/******************************************************************************/
/*                                                                            */
/*      xcoon - turn output on for a cgi display surface.                     */
/*   Note: it is not an error to turn output on for a surface which is        */
/*      already on.                                                           */
/*   Note: it is an error to turn output on for a surface which is not        */
/*      initialized.                                                          */
/*                                                                            */
/******************************************************************************/
#if defined(ADDC_)
void xcoon_(anything **surface_id) /* which surface to turn output on for*/
#else
void xcoon(anything **surface_id) /* which surface to turn output on for*/
#endif
{
  /* does surface_id point to a valid surface? */
  /* if not, do error action */
  /* find out which device this surface is on */
  /* rearrange the surface list for that device */

  short     dev;            /* which device to look at now */
  short     dev_found = 0;  /* which device was it found on */
  short     surf;           /* which surface on device to look at */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp;           /* for swap */

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (*surface_id == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  /* does surface_id point to a valid surface? */
  if (dev_found < 0) {
    fprintf(stderr, "xcoon: surface id not a valid surface\n");
    return;
  } /* end if surface not found */

  /* if surface was off */
  if (surf_found >= devices[dev_found].num_on_surfaces) {
    ++devices[dev_found].num_on_surfaces;
    /* swap target surface with first off surface if there is one */
    if (devices[dev_found].num_on_surfaces != devices[dev_found].num_active_surfaces) {
      temp = devices[dev_found].statelist[surf_found];
      devices[dev_found].statelist[surf_found] =
          devices[dev_found].statelist[devices[dev_found].num_on_surfaces - 1];
      devices[dev_found].statelist[devices[dev_found].num_on_surfaces - 1] = temp;
    } /* end if there is an off surface */
  }   /* end if surface was off */
} /* end xcoon */

/******************************************************************************/
/*                                                                            */
/*      xcact - initialize and activate a cgi display surface                 */
/*                                                                            */
/******************************************************************************/
#if defined(ADDC_)
void xcact_(void (*device_fn)(), anything **p_surface_id)
#else
void xcact(void (*device_fn)(), anything **p_surface_id)
#endif
{
  short     arg1; /* arg passed to device routine */
  short     i;
  anything *temp_surface[1]; /* for calling driver */
  short     which_device;    /* index of this device in devices */
  short     which_surface;   /* index of surface in devices*/

  /* is device already initialized? */
  which_device = -1;
  for (i = 0; i < num_devices; ++i) {
    if (device_fn == devices[i].device_fn) {
      which_device = i;
      break;
    } /* end if */
  }   /* end for */

  if (which_device < 0) {             /* if device not initialized */
    if (num_devices >= MAX_DEVICES) { /* if no room */
      fprintf(stderr, "xcact: can't activate: too many initialized devices\n");
      *p_surface_id = NULL;
      return;
    } /* end if no room for device */
  }   /* end if device not initialized */

  /* call the device driver with ACTIVATE, so it can allocate a state list */
  arg1         = ACTIVATE_FN;
  in_params[0] = (anything *)&arg1;
  (*device_fn)(in_params, 1, temp_surface);

  if (temp_surface[0] == NULL) { /* error */
    fprintf(stderr, "xcact: can't activate surface: driver error\n");
    *p_surface_id = NULL;
    return;
  } /* end if error on activate */

  if (which_device < 0) { /* if device not initialized */
    /* then initialize it */
    which_device = num_devices;
    ++num_devices;
    devices[which_device].num_active_surfaces = 0;
    devices[which_device].num_on_surfaces     = 0;
    devices[which_device].device_fn           = device_fn;
  } /* end if */

  /* add new surface to device */
  which_surface = devices[which_device].num_active_surfaces;
  ++devices[which_device].num_active_surfaces;
  devices[which_device].statelist[which_surface] = temp_surface[0];

  /* if new surface is the oldest surface, make it the solicitation surface */
  if (num_oldest == 0) {
    sol_surf = temp_surface[0];
  }

  /* put this surface as newest surf on oldest_surfs list */
  oldest_surfs[num_oldest] = temp_surface[0];
  ++num_oldest;

  /* return id of new surface */
  *p_surface_id = temp_surface[0];

} /* end xcact */

/******************************************************************************/
/*                                                                            */
/*      xcooff - turn output off for a display surface                        */
/*   Note: it is not an error to turn off a surface which is                  */
/*      already off.                                                          */
/*   Note: it is an error to turn off a surface which is not active.          */
/*                                                                            */
/******************************************************************************/

/******************************************************************************/
/*                                                                            */
/*      xcdact - deactivate a cgi display surface                             */
/*      Note: it is an error to deactivate a surface which is not active.     */
/*                                                                            */
/******************************************************************************/

/******************************************************************************/
/*                                                                            */
/*      xcsol - set solicitation surface                                      */
/*                                                                            */
/******************************************************************************/
#if defined(ADDC_)
void xcsol_(anything **surface_id)
#else
void xcsol(anything **surface_id)
#endif
{
  /* does surface_id point to a valid surface? */
  /* if not, do error action */
  /* set the solicitation surface to surface_id */

  short dev;       /* which device to look at now */
  short dev_found; /* which device was it found on */
  short surf;      /* which surface on device to look at */

  /* search devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (*surface_id == devices[dev].statelist[surf]) {
        dev_found = dev;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  /* does surface_id point to a valid surface? */
  if (dev_found < 0) {
    fprintf(stderr, "xcsol: surface id not a valid surface\n");
    return;
  } /* end if surface not found */

  sol_surf = *surface_id;

} /* end xcsol */
