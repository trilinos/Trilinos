// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_CONFIGDEFS_HPP
#define TPETRA_CONFIGDEFS_HPP

#ifndef __cplusplus
#define __cplusplus
#endif // ifndef __cplusplus

/* this section undefines all the things autotools defines for us that we wish it didn't. */

#ifdef PACKAGE
#undef PACKAGE
#endif // ifdef PACKAGE

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif // ifdef PACKAGE_NAME

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif // ifdef PACKAGE_BUGREPORT

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif // ifdef PACKAGE_STRING

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif // ifdef PACKAGE_TARNAME

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif // ifdef PACKAGE_VERSION

#ifdef VERSION
#undef VERSION
#endif // ifdef VERSION

// end of undoing autoconf's work section

#include <Tpetra_config.h>
#include <Teuchos_ConfigDefs.hpp>

// Used in all Tpetra code that explicitly must a type (like a loop index)
// that is used with the Teuchos::Array[View,RCP] classes.
namespace Tpetra { typedef Teuchos_Ordinal Array_size_type; }

// these make some of the macros in Tpetra_Util.hpp much easier to describe
#ifdef HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS
  #define TPETRA_THROWS_EFFICIENCY_WARNINGS 1
#else
  #define TPETRA_THROWS_EFFICIENCY_WARNINGS 0
#endif

#ifdef HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS
  #define TPETRA_PRINTS_EFFICIENCY_WARNINGS 1
#else
  #define TPETRA_PRINTS_EFFICIENCY_WARNINGS 0
#endif

#ifdef HAVE_TPETRA_THROW_ABUSE_WARNINGS
  #define TPETRA_THROWS_ABUSE_WARNINGS 1
#else
  #define TPETRA_THROWS_ABUSE_WARNINGS 0
#endif

#ifdef HAVE_TPETRA_PRINT_ABUSE_WARNINGS
  #define TPETRA_PRINTS_ABUSE_WARNINGS 1
#else
  #define TPETRA_PRINTS_ABUSE_WARNINGS 0
#endif

namespace Tpetra {
  typedef size_t global_size_t;

  /*! Local versus global allocation of Map elements */
  enum LocalGlobal {
    LocallyReplicated,  /*!< Indicates that map elements are locally replicated across all nodes */
    GloballyDistributed /*!< Indicates that map elements are globally distributed across all nodes */
  };

  /*! Return status of Map lookup */
  enum LookupStatus {
    AllIDsPresent, /*!< Indicates that all queried IDs were present in the Map */
    IDNotPresent   /*!< Indicates that at least one of the specified IDs was not present in the Map */
  };

  /*! Allocation profile for matrix/graph entries */
  enum ProfileType {
    StaticProfile,  /*!< Indicates a single, static allocation (strict and more efficient) */
    DynamicProfile  /*!< Indicates multiple, dynamic allocations (flexibile, but less efficient) */
  };

  /*! Optimize storage option */
  enum OptimizeOption {
    DoOptimizeStorage,   /*!< Indicates that storage should be optimized */
    DoNotOptimizeStorage /*!< Indicates that storage should not be optimized */
  };

  /*! Combine mode */
	enum CombineMode {
		ADD,    /*!< Existing values will be summed with new values. */
		INSERT, /*!< Insert new values that don't currently exist. */
		REPLACE /*!< Existing values will be replaced with new values. */
	};

}

#endif // TPETRA_CONFIGDEFS_HPP
