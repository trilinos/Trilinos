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
bool
path_access(
  const std::string &	name,
  int			mode)

{
  return !name.empty() && ::access(name.c_str(), mode) == 0;
}


bool
path_exists(
  const std::string &	name)
{
  return path_access(name, F_OK);
}


bool
path_read_access(
  const std::string &	name)
{
  return path_access(name, R_OK);
}


bool
path_write_access(
  const std::string &	name)
{
  return path_access(name, W_OK);
}


namespace {

struct flock *
file_lock(
  short	type,
  short	whence)
{
//  /* %TRACE[SPEC]% */ Tracespec trace__("sierra::Fmwk::<empty>::file_lock( short type, short whence)"); /* %TRACE% */
  static struct flock ret;
  ret.l_type = type;
  ret.l_start = 0;
  ret.l_whence = whence;
  ret.l_len = 0;
  ret.l_pid = 0; //getpid();
  return &ret;
}

} // namespace <empty>

bool
write_lock(
  int		fd)
{
  int i =::fcntl(fd, F_SETLK, file_lock(F_WRLCK, SEEK_SET));
//   if (i == -1)
//     fmwkout << "Write lock failed " << errno << dendl;

  return i != -1;
}


bool
release_lock(
  int		fd)
{
  int i =::fcntl(fd, F_SETLK, file_lock(F_UNLCK, SEEK_SET));
//   if (i == -1)
//     fmwkout << "Release lock failed " << errno << dendl;

  return i != -1;
}


bool
read_lock(
  int		fd)
{
  return ::fcntl(fd, F_SETLK, file_lock(F_RDLCK, SEEK_SET)) != -1;
}


bool
append_lock(
  int		fd)
{
  return ::fcntl(fd, F_SETLK, file_lock(F_WRLCK, SEEK_END)) != -1;
}

} // namespace stk
