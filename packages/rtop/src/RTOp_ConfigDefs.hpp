// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef _RTOP_CONFIGDEFS_H_
#define _RTOP_CONFIGDEFS_H_

/* Let Teuchos' configure process do all of the work! */
#include <Teuchos_ConfigDefs.hpp>

#include <RTOp_Config.h>

#ifdef HAVE_RTOP_DEBUG
#  define RTOP_DEBUG
#endif

#endif /*_RTOP_CONFIGDEFS_H_*/
