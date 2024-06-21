// @HEADER
// *****************************************************************************
//            TriBITS: Tribal Build, Integrate, and Test System
//
// Copyright 2013-2016 NTESS and the TriBITS contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef _MSC_VER
# define NOMINMAX
# include <Winsock2.h>
# include <process.h>
# define getpid _getpid
inline void sleep(int sec)
{
  Sleep(sec * 1000);
}
#pragma comment(lib, "Ws2_32.lib")
#endif
