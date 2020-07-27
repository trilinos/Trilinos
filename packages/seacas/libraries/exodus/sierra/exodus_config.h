/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/* NOTE: Definitions in this file only control building via the Sierra
 *       bjam system For SEACAS and Trilinos builds, these are all controlled
 *       in the cmake configure step
 */

/* Enable this block to build a thread-safe version of the exodus library */
#if 0
#define EXODUS_THREADSAFE
#endif

/* Deprecated Code Handling Options:
 * 1. Ignore -- treat deprecated functions as normal non-deprecated functions (default)
 * 2. Delete -- the deprecated functions are not defined or compiled (SEACAS_HIDE_DEPRECATED_CODE is
 * defined)
 * 3. Warn   -- if used in client code, issue a warning. (SEACAS_WARN_DEPRECATED_CODE is defined)
 *
 * The symbols SEACAS_HIDE_DEPRECATED_CODE and SEACAS_DEPRECATED are defined in exodus_config.h
 * In a TriBITs-based system, this include file is generated from cmake-variable definitions.
 * In other build systems, the exodus_config.h file is hard-wired.
 */

/* Enable this block to eliminate the building of deprecated functions in library */
#if 0
#define SEACAS_HIDE_DEPRECATED_CODE
#endif

/* Enable this block to enable warnings if any deprecated functions are used */
#if 0
#ifndef SEACAS_DEPRECATED
#if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#define SEACAS_DEPRECATED __attribute__((__deprecated__))
#else
#define SEACAS_DEPRECATED
#endif
#endif
#endif
