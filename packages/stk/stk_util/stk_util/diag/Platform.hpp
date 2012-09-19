/*--------------------------------------------------------------------*/
/*    Copyright 2003, 2008 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_UTIL_DIAG_Platform_h
#define STK_UTIL_DIAG_Platform_h

#include <iosfwd>
#include <mpi.h>

#include <stk_util/util/FeatureTest.hpp>
#include <stk_util/diag/Env.hpp>
#include <stk_util/diag/String.hpp>
#include <stk_util/diag/Option.hpp>
#include <stk_util/environment/RegisterProduct.hpp>

namespace sierra {
namespace Env {

///
/// @addtogroup EnvDetail
/// @{
///

/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>startup_preparallel_platform</b> ...
 *
 */
void startup_preparallel_platform();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>hostname</b> returns the hostname of the host running the
 * application.
 *
 * @return			a <b>String</b> value of the host name obtained from
 *				the operating system.
 */
std::string hostname();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>domainname</b> returns the domainname of the domain running the
 * application.
 *
 * @return			a <b>String</b> value of the domain name obtained from
 *				the operating system.
 */
std::string domainname();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>username</b> returns the username of the user running the
 * application.
 *
 * @return			a <b>String</b> value of the username obtained from
 *				the operating system.
 */
std::string username();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>hardware</b> returns the hardware type of the host running the
 * application.
 *
 * @return			a <b>String</b> value of the <b>machine</b>
 *				field of the <b>uname</b> system call or equivalent
 *				obtained from the operating system.
 */
std::string hardware();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>osname</b> returns the operating system nameof the host running the
 * application.
 *
 * @return			a <b>String</b> value of the <b>sysname</b>
 *				field of the <b>uname</b> system call or equivalent
 *				obtained from the operating system.
 */
std::string osname();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>osversion</b> returns the hardware type of the host running the
 * application.
 *
 * @return			a <b>String</b> value of the <b>release</b>
 *				field of the <b>uname</b> system call or equivalent
 *				obtained from the operating system.
 */
std::string osversion();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>pid</b> returns the process id of the process running the
 * application.
 *
 * @return			a <b>int</b> value of the process id obtained from
 *				the operating system.
 */
int pid();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>pgrp</b> returns the process group id of the process running
 * the application.
 *
 * @return			a <b>int</b> value of the process group id obtained from
 *				the operating system.
 */
int pgrp();


/**
 * @brief Member function <b>get_heap_info</b> returns the amount of heap
 * memory used in bytes and the largest free block of memory in bytes.
 *
 * @param heap_size		a <b>size_t</b> returns the amount of heap
 * memory used in bytes.
 *
 * @param largest_free		a <b>size_t</b> returns the largest free block
 * of memory.
 *
 */
void get_heap_info(size_t &heap_size, size_t &largest_free);


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>get_heap_usage</b> returns the number of bytes used by the heap.
 *
 * @return			a <b>size_t</b> value of the number of bytes used by
 *				the heap.
 */
inline size_t get_heap_usage() {
  size_t heap_size;
  size_t largest_free;
  get_heap_info(heap_size, largest_free);

  return heap_size;
}

/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>get_available_memory</b> returns an estimation of the amount of memory available to the process.
 *
 * @return			a <b>size_t</b> value of the number of bytes available to the process.
 */
size_t get_available_memory();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>get_memory_info</b> returns the total memory usage of the
 * process and the number of page faults accumulated by the process.
 *
 * @param memory_usage		a <b>size_t</b> reference to receive the number of
 *				bytes currently used by the process.
 *
 * @param faults		a <b>size_t</b> reference to treceive the number of
 *				page faults incurred by the process.
 *
 */
void get_memory_info(size_t &memory_usage, size_t &faults);


/**
 * @ingroup EnvRuntimeInformationDetail EnvOutput
 * @brief Function <b>path_exists</b> returns true if the path exists.
 *
 * @param path			a <b>String</b> const reference to the path to have
 *				existence tested.
 *
 * @return			a <b>bool</b> value of true if the path exists.
 */
bool path_exists(const std::string &path);


/**
 * @ingroup EnvRuntimeInformation EnvOutputDetail
 * @brief Function <b>path_access</b> returns true if the process has permission to
 * access path with the specified mode.
 *
 * @param path			a <b>String</b> const reference to the path to check
 *				for <b>mode</b> access.
 *
 * @param mode			an <b>int</b> value of the mode to test.
 *
 * @return			a <b>bool</b> value of true of the process has
 *				permission to access the path with <b>mode</b>
 *				access.
 */
bool path_access(const std::string &path, int mode);


/**
 * @ingroup EnvRuntimeInformationDetail EnvOutputDetail
 * @brief Function <b>path_read_access</b> returns true if the process has read
 * access to the path.
 *
 * @param path			a <b>String</b> const reference to the path to check
 *				for read access.
 *
 * @return			a <b>bool</b> value of true of the process has
 *				permission to access the path with read access.
 */
bool path_read_access(const std::string &path);

/**
 * @ingroup EnvRuntimeInformationDetail EnvOutputDetail
 * @brief Function <b>path_write_access</b> returns true if the process has write
 * access to the path.
 *
 * @param path			a <b>String</b> const reference to the path to check
 *				for write access.
 *
 * @return			a <b>bool</b> value of true of the process has
 *				permission to access the path with write access.
 */
bool path_write_access(const std::string &path);


/**
 * @ingroup EnvRuntimeInformationDetail EnvOutputDetail
 * @brief Function <b>read_lock</b> returns true if the process was able to place
 * a shared lock on the specified file descriptor.
 *
 * @param fd			an <b>int</b> value of the file description to
 *				attempt to lock.
 *
 * @return			a <b>bool</b> value of true of the lock succeeded.
 */
bool read_lock(int fd);


/**
 * @ingroup EnvRuntimeInformationDetail EnvOutputDetail
 * @brief Function <b>write_lock</b> returns true if the process was able to place
 * an exclusive lock on the specified file descriptor.
 *
 * @param fd			an <b>int</b> value of the file description to
 *				attempt to lock.
 *
 * @return			a <b>bool</b> value of true of the lock succeeded.
 */
bool write_lock(int fd);


/**
 * @ingroup EnvRuntimeInformationDetail EnvOutputDetail
 * @brief Function <b>append_lock</b> returns true if the process was able to place
 * an exclusive lock on the end of the specified file descriptor.  Existing records may
 * still be accessed.
 *
 * @param fd			an <b>int</b> value of the file description to
 *				attempt to lock.
 *
 * @return			a <b>bool</b> value of true of the lock succeeded.
 */
bool append_lock(int fd);


/**
 * @ingroup EnvRuntimeInformationDetail EnvOutputDetail
 * @brief Function <b>release_lock</b> returns true if the process was able to
 * release a lock previously palced on the specified file descriptor.
 *
 * @param fd			an <b>int</b> value of the file description to have
 *				the lock released.
 *
 * @return			a <b>bool</b> value of true of the lock release
 *				succeeded.
 */
bool release_lock(int fd);

///
/// @}
///

} // namespace Env
} // namespace sierra

#endif // STK_UTIL_DIAG_Platform_h
