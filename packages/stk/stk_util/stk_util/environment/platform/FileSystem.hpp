#ifndef stk_util_environment_platform_FileSystem_hpp
#define stk_util_environment_platform_FileSystem_hpp

namespace sierra stk {

///
/// @addtogroup EnvDetail
/// @{
///

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

} // namespace stk

#endif // stk_util_environment_platform_FileSystem_hpp
