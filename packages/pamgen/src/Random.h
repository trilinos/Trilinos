// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// the system defined random number generator
#ifdef _WIN32
# define SRANDOM(i) srand(i)
# define RANDOM() rand()
#else
# define SRANDOM(i) srandom(i)
# define RANDOM() random()
#endif
