/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
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
