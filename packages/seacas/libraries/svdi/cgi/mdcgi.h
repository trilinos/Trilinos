/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "cgi.h"                    // for device_struct, MAX_DEVICES, etc
void xcoon_(anything **surface_id); /* which surface to turn output on for*/
void xcact_(void (*device_fn)(), anything **p_surface_id);
