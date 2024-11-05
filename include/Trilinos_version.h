/*

            Trilinos: An Object-Oriented Solver Framework
                    Copyright (c) 2001-2024 NTESS

Copyright 2001-2024 National Technology & Engineering Solutions of Sandia,
LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the
U.S. Government retains certain rights in this software.

Copyright the Trilinos contributors.


*/

#ifndef TRILINOS_VERSION_H
#define TRILINOS_VERSION_H


/* Trilinos version numbering convention.
 *
 * Trilinos version numbers take the form X.Y.Z where:
 *
 *   X: The major version number that defines a window of (perfect) backward
 *   compatibility (see below).
 *
 *   Y: The release version number within a backward-compatible set of
 *   versions X.  Even numbers (0, 2, 4, ...) are used for releases and odd
 *   numbers (1, 3, 5, ...) are used for development versions in-between
 *   releases.
 *
 *   Z: The minor release version number for minor releases taken off off a
 *   release branch X.Y.  Even numbers (0, 2, 4, ...) are used for customer
 *   releases and odd numbers (1, 3, 5, ...) are used for the code on the
 *   release X.Y branch in-between minor releases.
 *
 * All Trilinos releases (i.e. X.Y where Y is even) are taken off of the
 * development branch (i.e. the dev version X-1.R or X.Y-1) and are given a
 * name containing the version number X.Y.  The initial releases in a backward
 * compatible set are then given the release numbers:
 *
 *   X.0.0, X.2.0, X.4.0, ...
 *
 * The intermediate development versions are given the release numbers:
 *
 *   X.1.0, X.3.0, X.5.0, ....
 *
 * For development versions, the minor release version number Z is always 0.
 *
 * The minor releases for a given release branch X.Y are given the version
 * numbers:
 *
 *   X.Y.0, X.Y.2, X.Y.4, ...
 *
 * The version numbers given to the code in the release branch X.Y in-between
 * minor releases (which are not branched, only tagged) are:
 *
 *   X.Y.1, X.Y.3, X.Y.5, ...
 *
 * In this way, client code can just examine the version number in this file
 * and know exactly what version of Trilinos they are working with with no
 * ambiguity no matter what.
 */


/* The major version number xx (allows up 99 major Trilinos release version
 * numbers).
 *
 * The major Trilinos version number defines a window of backward
 * compatibility.
 */
#define TRILINOS_MAJOR_VERSION 16

/* The major, release, and minor release version numbers (i.e. xx.yy.zz).
*
* NOTE: When numbers are less than 10, it is padded with a 0.  For example,
* development version 10.1 of Trilinos is designated 100100 and the release
* version 10.2.4 is designated 100204.  This preserves the comparability of
* these version numbers with simple comparison operators used in #ifdef tests.
*/
#define TRILINOS_MAJOR_MINOR_VERSION 160100

/* NOTE: These macros are given long int values to allow comparisons in
 * preprocessor #if statements.  For example, you can do comparisons with ==,
 * <, <=, >, and >=.
 *
 * NOTE: The C++ standard for the C preprocessor requires that the arguments
 * for #if must be convertible into a long int.  Expressions that convert to 1
 * are true and expressions that convert to 0 are false.
 */

/* \brief Version string for Trilinos.
 *
 * NOTE: This string is to be used for outputting, not for comparison logic.
 */
#define TRILINOS_VERSION_STRING "16.1.0-dev"

#endif /* TRILINOS_VERSION_H */
