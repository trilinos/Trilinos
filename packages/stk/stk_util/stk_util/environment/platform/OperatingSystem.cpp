/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/environment/platform/OperatingSystem.hpp>

#include <pwd.h>
#include <unistd.h>

#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <stdexcept>
#include <iomanip>
#include <algorithm>
#include <locale>

#include <fcntl.h>

#if defined(__GNUC__)
#include <fstream>
#include <malloc.h>
#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>
#if __GNUC__ == 3 || __GNUC__ == 4
#include <cxxabi.h>
#endif

#elif defined(__PGI)
#include <fstream>
#include <malloc.h>
#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>

#elif defined(__sun)
#include <fstream>
#include <procfs.h>
#include <sys/resource.h>
#include <sys/systeminfo.h>
#include <sys/utsname.h>
#include <sys/time.h>

#elif defined(__SUNPRO_CC)
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/utsname.h>
#include <netdb.h>
#endif

#if defined(__PUMAGON__)
extern "C" {
#include <util.h>
#include <sys/param.h>
}

#elif defined(__sgi)
#include <sys/time.h>
#include <sys/resource.h>

#elif defined(REDS)
#include <sys/param.h>
#include <sys/utsname.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <catamount/catmalloc.h>
extern void *__heap_start;		///< Magic address of end of instruction/static data
extern void *__heap_end;		///< Magic address of end of heap

#elif defined(__JVN)
#include <sys/param.h>
#include <sys/utsname.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

#else
#include <sys/utsname.h>
#include <sys/time.h>
#include <netdb.h>
#endif

namespace stk {

std::string
hostname()
{
  char buf[255];
  ::gethostname(buf, sizeof(buf));
  return std::string(buf);
}


std::string
domainname()
{
#if defined(__PUMAGON__) || defined(REDS)
  return std::string(".sandia.gov");

#elif defined(__sun)
  std::string domain(".");
  char buf[255];

  ::sysinfo(SI_SRPC_DOMAIN, buf, sizeof(buf));
  if (std::strlen(buf)) {
    domain += buf;
  }
  return domain;

#else //#elif defined(__linux) || defined(__sgi) || defined(_AIX)
  std::string domain(".");
  char buf[255];

  ::getdomainname(buf, sizeof(buf));
  if (::strlen(buf)) {
    domain += buf;
  }
  return domain;

#endif
}


std::string
username()
{
#if defined(REDS)
  return "unknown";
#else
  struct passwd *user_info = ::getpwuid(::geteuid());

  return (user_info ? user_info->pw_name : "unknown");
#endif
}


std::string
hardware()
{
#ifndef __sgi
  struct utsname	uts_name;

  uname(&uts_name);

  return uts_name.machine;
#else
  std::string s;
  return s;
#endif
}


std::string
osname()
{
#ifndef __sgi
  struct utsname	uts_name;

  uname(&uts_name);

  return uts_name.sysname;
#else
  std::string s;
  return s;
#endif
}


std::string
osversion()
{
#ifndef __sgi
  struct utsname	uts_name;

  uname(&uts_name);

  return uts_name.release;
#else
  std::string s;
  return s;
#endif
}


int
pid()
{
  return ::getpid();
}


int
pgrp()
{
#if defined(__PUMAGON__) || defined(REDS)
  return 0;
#else
  return ::getpgrp();
#endif
}


} // namespace stk
