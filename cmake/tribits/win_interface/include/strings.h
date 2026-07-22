// @HEADER
// *****************************************************************************
//            TriBITS: Tribal Build, Integrate, and Test System
//
// Copyright 2013-2016 NTESS and the TriBITS contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// On windows stricmp and strnicmp
// are strcasecmp and strncasecmp, and are
// include in other header files
#define strcasecmp stricmp
#define strncasecmp strnicmp
