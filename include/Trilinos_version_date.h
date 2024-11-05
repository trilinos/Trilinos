/* TRILINOS_VERSION_DATE
 *
 * This macro gives the version date of the Trilinos git repository at
 * the time it was configured.  It gives the date YYYY-MM-DD and the 2-digit
 * (24) hour hh of the commit as a 10-digit integer:
 *
 *   YYYYYMMDDhh
 *
 * As long as the branch for the git repo is not hard reset, the first-parent
 * history should give a monotonically increasing 10-digit integer.  This
 * 10-digit date/time integer YYYYMMDDhh will fit in a signed 32-bit integer
 * with a maximum value of 2^32 / 2 - 1 = 2147483647.  Therefore, the maximum
 * date that can be handled is the year 2147 with the max date/time of 2147 12
 * 31 23 = 2147123123.  Modern C/C++ compilers (and other implementations of
 * the C preprocessor) should be able to handle 32-bit integers.
 * 
 * This header file is meant to be included by downstream codes to determine
 * the version of Trilinos being used and allows
 * TRILINOS_VERSION_DATE to be used in C/C++ ifdefs like:
 *
 * #if defined(TRILINOS_VERSION_DATE) && TRILINOS_VERSION_DATE >= 2019032704
 *   // The version is newer than 2019-03-27 04:00:00 UTC
 *   ...
 * #else
 *   // The version is older than 2019-03-27 04:00:00 UTC
 *     ...
 *  #endif
 *
 * This allows downstream codes to know the fine-grained version of
 * Trilinos at configure and build time to adjust for the addition of
 * new features, deprecation of code, or breaks in backward compatibility
 * (which occur in specific commits with unique commit dates).
 *
 */

#define TRILINOS_VERSION_DATE 2024110204
