/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_environment_platform_OperatingSystem_hpp
#define stk_util_environment_platform_OperatingSystem_hpp

#include <string>

namespace stk {

///
/// @addtogroup EnvDetail
/// @{
///

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

///
/// @}
///

} // namespace stk

#endif // stk_util_environment_platform_OperatingSystem_hpp
