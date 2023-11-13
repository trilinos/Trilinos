// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <sys_info.h>

#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64)
#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#undef IN
#undef OUT
#include <fmt/ostream.h>
#include <sstream>
#else
#include <sys/utsname.h>
#endif
#include <fmt/format.h>

std::string sys_info(const std::string &codename)
{
  // Return 'uname' output.  This is used as information data records
  // in output exodus files to help in tracking when/where/... the
  // file was created

#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64)
  char  machine_name[MAX_COMPUTERNAME_LENGTH + 1] = {0};
  DWORD buf_len                                   = MAX_COMPUTERNAME_LENGTH + 1;
  ::GetComputerName(machine_name, &buf_len);

  std::string info = codename + ": ";
  info += machine_name;
  info += ", OS: ";

  std::string   os = "Microsoft Windows";
  OSVERSIONINFO osvi;

  ZeroMemory(&osvi, sizeof(OSVERSIONINFO));
  osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);

  if (GetVersionEx(&osvi)) {
    DWORD             build = osvi.dwBuildNumber & 0xFFFF;
    std::stringstream str;
    fmt::print(info, " {}.{} {} (Build {})", osvi.dwMajorVersion, osvi.dwMinorVersion,
               osvi.szCSDVersion, build);
    os += str.str();
  }
  info += os;
#else
  struct utsname sys_info
  {
  };
  uname(&sys_info);

  std::string info =
      fmt::format("{}: {}, OS: {} {}, {}, Machine: {}", codename, sys_info.nodename,
                  sys_info.sysname, sys_info.release, sys_info.version, sys_info.machine);
#endif
  return info;
}
