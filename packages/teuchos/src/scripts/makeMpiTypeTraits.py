# @HEADER
# ***********************************************************************
#
#                    Teuchos: Common Tools Package
#                 Copyright (2004) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ***********************************************************************
# @HEADER

'''Generate C++ code for converting from a C++ type to its MPI_Datatype.

Author: Mark Hoemmen <mhoemme@sandia.gov>
Date: Aug, Sep 2012

Introduction
============

This module generates C++ code for MpiTypeTraits, a traits class
which, for a given C++ type Packet, returns the corresponding
MPI_Datatype.  The MPI_Datatype tells MPI (Message Passing Interface)
how to send and receive objects of type Packet.  The traits class also
tells users whether or not the MPI_Datatype must be freed (via
MPI_Type_free()) after use.

Usage
=====

If you run this code as a script with no arguments, it will write two
files to the current working directory:

Teuchos_MpiTypeTraits.hpp
Teuchos_MpiTypeTraits.cpp

The first file contains the generic declaration and documentation of
MpiTypeTraits, declarations of full specializations of MpiTypeTraits,
and declarations and definitions of some partial specializations.
The second file contains definitions of the full specializations of
MpiTypeTraits.

After creating these C++ files, one should copy them manually into the
parent directory (teuchos/src).  They should also be listed in
teuchos/src/CMakeLists.txt, conditionally if MPI is enabled.  Look for
"IF (${PACKAGE_NAME}_ENABLE_MPI)", and append Teuchos_MpiTypeTraits.hpp
to the HEADERS set and Teuchos_MpiTypeTraits.cpp to the SOURCES set.'''

# We actually reuse templates, which is why we use string.Template
# instead of "new"-style format statements.  This is hopefully more
# efficient.
from string import Template
# Redirect stdout to a file.
import sys

def emitCopyrightNotice ():
    '''Print the Teuchos package's copyright notice.
    This belongs at the top of every Teuchos header or source file.'''
    
    print '''// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
'''

def makeTypeMap ():
    '''Return the map from C++ Packet type to the corresponding MPI_Datatype.

Each key is a C++ type t for which we are generating a specialization
of MpiTypeTraits<Packet=t>.  Its corresponding value is the
MPI_Datatype which tells MPI how to send and receive data of that
type.  The map only includes keys for C++ types for which the
corresponding MPI_Datatype is a "basic" (not derived) type.  A derived
MPI_Datatype must be freed after use via MPI_Type_free(); a "basic"
MPI_Datatype need not and must not be freed after use.'''

    # Start with an empty dictionary; fill it below.
    d = {}

    # Many built-in C++ types have a simple conversion to their
    # standard MPI_Datatype values: Capitalize, replace each space
    # with an underscore, and prepend MPI_.
    simpleNames = ['char', 'signed char', 'unsigned char',
		   'short', 'unsigned short',
		   'int', 'unsigned int',
		   'long', 'unsigned long',
		   'long long', 'unsigned long long', \
		   'size_t', 'ptrdiff_t',
		   'float', 'double', 'long double']
    d.update (dict ((name, 'MPI_' + name.upper ().replace (' ', '_')) for name in simpleNames))

    # Other C++ types don't have a straightforward conversion.
    # MPI_CXX_COMPLEX is in the MPI 3.0 standard (thanks in small part
    # to Jeff Hammond and Jed Brown thinking of me!) but
    # MPI_C_COMPLEX, MPI_C_FLOAT_COMPLEX, MPI_C_DOUBLE_COMPLEX, and
    # MPI_C_LONG_DOUBLE_COMPLEX work perfectly well, in case your MPI
    # implementation does not support the 3.0 standard.  We deal with
    # bool elsewhere as a separate case, since there is no basic MPI
    # datatype which standardly supports bool.
    d.update ({'std::complex<float>': 'MPI_C_FLOAT_COMPLEX', \
	       'std::complex<double>': 'MPI_C_DOUBLE_COMPLEX', \
	       'std::complex<long double>': 'MPI_C_LONG_DOUBLE_COMPLEX'})
    return d

def emit (s, indent=0):
    '''Print code to stdout, with the proper indent level.

    s (string): the string to print.  It may contain one or more lines.
    indent (integer): the number of spaces to indent each line.'''
    for line in s.split('\n'):
	print ' '*indent + line

def emitGenericDecl (d, indent=0):
    '''Emit the declaration of MpiTypeTraits, with documentation.

    d (dictionary): return value of makeTypeMap().
    indent (integer): the number of spaces to indent.'''
    
    emit ('''/// \class MpiTypeTraits
/// \\brief Traits class mapping a given Packet type to its MPI_Datatype.
/// \\tparam Packet The type of data being received and sent.
///
/// \section Summary
///
/// For a given C++ type Packet, this traits class gives you the
/// corresponding MPI_Datatype.  The latter tells MPI how to send and
/// receive objects of type Packet.  MpiTypeTraits also tells you
/// whether MPI_Datatype is a basic or derived type.  Derived types
/// must be freed after use via MPI_Type_free().  For more details,
/// please refer to the MPI standard.
///
/// This class works by specialization.  If there is no specialization
/// of MpiTypeTraits for a type Packet, then expressions involving
/// MpiTypeTraits<Packet> will not compile.  You are welcome to
/// specialize this class for your own Packet types.  Be sure to check
/// this header file for the list of specializations that we provide.
/// Duplicating an existing specialization will probably result in a
/// compile error.  Also, be sure that the values of needFree and
/// constantSize (see below) in your specialization are compile-time
/// constants whose values are visible in the header file.
///
/// \\note We assume that on all MPI processes, for all built-in
///   C++ types T, sizeof(T) has the same value.  It is technically
///   possible to run on clusters with heterogeneous nodes, such that
///   sizeof(T) for types such as bool or size_t may have different
///   values on different nodes.  However, Teuchos and its client
///   packages do not currently work for this use case, so we do not
///   support it here.  We make no effort to check whether different
///   processes have different values of sizeof(T) for all built-in
///   C++ types T.
///
/// \section Design discussion
///
/// This section is likely only of interest to implementers of
/// specializations of MpiTypeTraits, or to those who have questions
/// about the design choices made by the authors of MpiTypeTraits.
///
/// \subsection Notes on derived datatypes
///
/// Please refer to the <a href="http://www.mpi-forum.org/docs/mpi22-report/node78.htm#Node78">relevant part of the MPI 2.2 standard</a>.
///
/// MPI_TYPE_FREE "[m]arks the datatype object associated with datatype for
/// deallocation and sets datatype to MPI_DATATYPE_NULL. Any communication
/// that is currently using this datatype will complete normally.  Freeing
/// a datatype does not affect any other datatype that was built from the
/// freed datatype. The system behaves as if input datatype arguments to
/// derived datatype constructors are passed by value."
///
/// I'll unpack these statements:
///
/// 1. We may call MPI_Type_free on any MPI_Datatype, even if a
///    nonblocking operation (e.g., an MPI_Irecv) with that datatype
///    hasn't yet completed.
///
/// This lets us use RCP<OpaqueWrapper<MPI_Datatype> > to handle
/// calling MPI_Type_free automatically for any datatype.  We don't
/// have to try to track MPI_Irecv or other nonblocking operations
/// (such as the new nonblocking collectives in MPI 3.0).
///
/// 2. We may safely create composite datatypes (e.g., std::pair<X,Y>)
///    from previously created derived datatypes.
///
/// The composite datatype will work correctly even if its component
/// datatypes are freed later.  Furthermore, we may free the composite
/// datatype either before or after its component datatypes are freed.
/// This lets <tt>MpiTypeTraits<T>::makeType()</tt> return a "raw" 
/// MPI_Datatype rather than an <tt>RCP<OpaqueWrapper<MPI_Datatype> ></tt>.
///
/// \subsection Why does makeType() return a "raw" MPI_Datatype?
///
/// Users might wonder why makeType() returns a "raw" MPI_Datatype
/// rather than a "wrapped" <tt>RCP<OpaqueWrapper<MPI_Datatype> ></tt>
/// with a suitable deallocator.  That would be helpful for creating 
/// composite datatypes; for example, it would make it easier to write 
/// exception-safe code.  We chose to return the "raw" MPI_Datatype first, 
/// because it adds overhead for the overwhelmingly common case of basic 
/// (not derived) datatypes, and second, to give clients of this class a
/// choice about their deallocation mechanism.
template<class Packet>
class MpiTypeTraits {
 public:
  //! The type of the data being received and sent.
  typedef Packet packet_type;
  
  /// Whether the MPI_Datatype returned by makeType() is derived or basic.
  /// Derived types must be freed after use with MPI_Type_free(); basic types
  /// must _not_ be freed.  See the
  /// <a href="http://www.mcs.anl.gov/research/projects/mpi/www/www3/MPI_Type_free.html">documentation of MPI_Type_free()</a>.
  ///
  /// \\note The class datum needFree is a compile-time constant so that
  ///   you may use it in template metaprogramming (as a bool template
  ///   parameter).  If its value in your specialization of MpiTypeTraits
  ///   depends on some system parameters, you must resolve its value at
  ///   configure or compile time.
  static const bool needFree = false;

  /// \\brief Whether all Packet objects on all MPI processes
  ///   have the same size throughout the lifetime of a program.
  ///
  /// If this is true, then throughout the lifetime of a program:
  ///
  /// 1. For every MPI process, the same MPI_Datatype (which is the value
  ///    returned by makeType()) is valid for sending or receiving one
  ///    object of type Packet.
  ///
  /// 2. On every MPI process, sending an object x of type Packet using 
  ///    that MPI_Datatype results in the receiver receiving an object y of
  ///    type Packet such that x and y are equal (using whatever definition 
  ///    of "equal" makes sense for objects of type Packet).
  ///
  /// It is possible for constantSize to be false, even if sizeof(Packet) is
  /// a compile-time constant.  For example, Packet may contain a pointer to a
  /// dynamically allocated array, whose length may differ for different
  /// Packet objects.  If that dynamically array has the same length on all
  /// MPI processes at all times in a single execution of a program, then
  /// constantSize would be true.
  ///
  /// If constantSize is false, we make no guarantees of correctness when
  /// using the MPI_Datatype returned from makeTyope().  Any of the following
  /// may have different MPI_Datatypes:
  ///
  /// 1. Different Packet instances on the same process
  /// 2. Different Packet instances on different processes
  /// 3. The same Packet instance on the same process at different times
  ///
  /// \\note To Trilinos developers: If implementing a cache for MPI datatypes,
  ///   it would be wise not to cache datatypes corresponding to Packet types
  ///   for which constantSize is false.  This will take care of the common case
  ///   where all processes agree on the size of Packet in a bulk synchronous
  ///   phase of computation, even though Packet's interface allows its size to
  ///   change at any time.
  static const bool constantSize = false;

  /// \\brief Return the MPI_Datatype and length corresponding to the Packet type.
  ///
  /// If constantSize is true, then you may use the returned MPI_Datatype
  /// to send and receive objects of type Packet, throughout the lifetime
  /// of the program.  Otherwise, we make no guarantees that this is
  /// possible; for example, it may be necessary to agree on the amount
  /// of data to send before sending or receiving it.
  ///
  /// This function requires a Packet instance.  In general, for maximum
  /// portability, one may need this in order to compute the MPI_Datatype.
  /// For example, one might compute the MPI_Datatype of an std::pair<X,Y> 
  /// instance \c x using MPI_Type_struct(), by measuring address offsets
  /// of x.first and x.second from the address of x itself.  C++ does not
  /// provide a way to perform introspection on the layout of a class or
  /// struct, without actually making an instance of it.
  ///
  /// If this function did not take a Packet instance, then Packet would 
  /// need to be default constructible (so that we could make an actual 
  /// Packet instance without knowing what the type Packet is).
  static std::pair<MPI_Datatype, size_t> makeType (const Packet& example) {
    // Raise a compile-time error in case no specialization of this
    // traits class has been defined for the given Packet type.
    Packet::noSpecializationForThisType();

    // Return a valid value, so the compiler doesn't complain.
    return std::make_pair (MPI_TYPE_NULL, static_cast<size_t> (1));
  }
};
''', indent)
	
def emitSpecializationDecls (d, indent=0):
    '''Emit declarations of the specializations of MpiTypeTraits.

    d (dictionary): return value of makeTypeMap().
    indent (integer): the number of spaces to indent.'''
    
    # It's a template template.  Whee!!!
    # Fill it in for each C++ type key of the dictionary.
    tmpl = Template('''// Specialization of MpiTypeTraits<T> for T=${key}.
template<>
class MpiTypeTraits< ${key} > {
 public:
  typedef ${key} packet_type;
  static const bool needFree = false;
  static const bool constantSize = true;

  static std::pair<MPI_Datatype, size_t> makeType (const ${key}& example);
};
''')
    for k,v in d.iteritems():
        if k != 'long long' and k != 'unsigned long long' and k != 'size_t' and k != 'ptrdiff_t': # we handle these specially
	    emit (tmpl.substitute (key=k), indent)
    emit ('''#if defined(HAVE_TEUCHOS_LONG_LONG_INT)

// Partial specializations for long long and unsigned long long.
// On platforms with sizeof(size_t) <= sizeof(unsigned long long), 
// this should take care of the size_t and ptrdiff_t specializations
// as well, since we've covered all built-in unsigned integer types
// above with size <= sizeof(unsigned long long).
''', indent)
    emit (tmpl.substitute (key='long long', value=d['long long']), indent)
    emit (tmpl.substitute (key='unsigned long long', value=d['unsigned long long']), indent)
    emit ('''// The C preprocessor does not allow "sizeof(T)" expressions in #if
// statements, even if T is a built-in type.  Otherwise, we could test
// for 'sizeof(size_t) > sizeof(unsigned long int)'.  The constants
// below are defined in the <cstdint> header file.
#elif SIZE_MAX > ULONG_MAX

// We already have an unsigned long int specialization above.  If
// Teuchos support for "long long" is enabled, then we've taken care
// of all possible lengths of size_t: unsigned (char, short, int,
// long, long long).  If "long long" is _not_ enabled, we need to
// check if sizeof(size_t) > sizeof(unsigned long).  If so, then we
// need a specialization for size_t.  Ditto for ptrdiff_t (which is a
// signed type of the same length as size_t).
''', indent)
    emit (tmpl.substitute (key='ptrdiff_t', value=d['ptrdiff_t']), indent)
    emit (tmpl.substitute (key='size_t', value=d['size_t']), indent)
    emit ('#endif // HAVE_TEUCHOS_LONG_LONG_INT\n', indent)

def emitSpecializationDefs (d, indent=0):
    '''Emit definitions of the specializations of MpiTypeTraits.

    d (dictionary): return value of makeTypeMap().
    indent (integer): the number of spaces to indent.'''
    
    # It's a template template.  Whee!!!
    # Fill it in for each (C++ type, MPI_Datatype) pair in the dictionary.
    tmpl = Template('''// Specialization of MpiTypeTraits<T> for T=${key}.
template<>
MPI_Datatype MpiTypeTraits< ${key} >::makeType () {
  return std::make_pair (${value}, static_cast<size_t> (1));
}
''')
    for k,v in d.iteritems():
        if k != 'long long' and k != 'unsigned long long' and k != 'size_t' and k != 'ptrdiff_t': # we handle these specially
	    emit (tmpl.substitute (key=k, value=v), indent)
    emit ('''#if defined(HAVE_TEUCHOS_LONG_LONG_INT)

// Partial specializations for long long and unsigned long long.
// On platforms with sizeof(size_t) <= sizeof(unsigned long long), 
// this should take care of the size_t and ptrdiff_t specializations
// as well, since we've covered all built-in unsigned integer types
// above with size <= sizeof(unsigned long long).
''', indent)
    emit (tmpl.substitute (key='long long', value=d['long long']), indent)
    emit (tmpl.substitute (key='unsigned long long', value=d['unsigned long long']), indent)
    emit ('''// The C preprocessor does not allow "sizeof(T)" expressions in #if
// statements, even if T is a built-in type.  Otherwise, we could test
// for 'sizeof(size_t) > sizeof(unsigned long int)'.  The constants
// below are defined in the <cstdint> header file.
#elif SIZE_MAX > ULONG_MAX

// We already have an unsigned long int specialization above.  If
// Teuchos support for "long long" is enabled, then we've taken care
// of all possible lengths of size_t: unsigned (char, short, int,
// long, long long).  If "long long" is _not_ enabled, we need to
// check if sizeof(size_t) > sizeof(unsigned long).  If so, then we
// need a specialization for size_t.  Ditto for ptrdiff_t (which is a
// signed type of the same length as size_t).
''', indent)
    emit (tmpl.substitute (key='ptrdiff_t', value=d['ptrdiff_t']), indent)
    emit (tmpl.substitute (key='size_t', value=d['size_t']), indent)
    emit ('#endif // HAVE_TEUCHOS_LONG_LONG_INT\n', indent)

def headerizeFilename (filename):
    '''Convert filename to a suitable string for a header file's #ifndef ... #endif.

    filename (string): Name of the header file.'''
    return filename.upper ().replace ('.', '_')

def emitStdPairDeclAndDef (indent=0):
    '''Emit decl. and def. of MpiTypeTraits<std::pair<X,Y> > partial specialization.

    indent (integer): the number of spaces to indent.'''
    emit ('''//! Partial specialization for an std::pair of two types X and Y.
template<class X, class Y>
class MpiTypeTraits<std::pair<X, Y> > {
 public:
  typedef std::pair<X, Y> packet_type;
  static const bool needFree = true;
  static const bool constantSize = MpiTypeTraits<X>::constantSize && MpiTypeTraits<Y>::constantSize;

  static std::pair<MPI_Datatype, size_t> makeType (const std::pair<X, Y>& example) {
    int err = 0;
    MPI_Datatype t; // the return value
    int blkLens[2] = {1, 1}; // one each of X and Y
    MPI_Aint disps[2] = {&(example.first) - &example, &(example.second) - &example};

    // NOTE (mfh 30 Aug 2012) If we call makeType() for X and Y, we
    // might be wasting some effort if X and Y are derived (instead of
    // basic) types.  However, this behavior is still correct; it does
    // not cause problems if X's or Y's type gets freed but the type
    // of their pair is not.  This comes from the MPI 2.2 standard,
    // in particular
    // <a href="http://www.mpi-forum.org/docs/mpi22-report/node78.htm#Node78">this section</a>.
    //
    // "Freeing a datatype does not affect any other datatype that was
    // built from the freed datatype.  The system behaves as if input
    // datatype arguments to derived datatype constructors are passed by
    // value."
    //
    // It's unlikely that our users will be constructing any composite
    // datatypes from datatypes that themselves are expensive to
    // compute.  If that does become the case, however, we could
    // refactor type creation to use a datatype cache that returns
    // RCP<OpaqueWrapper<MPI_Datatype> > for the datatype
    // corresponding to a given type T.  This could even be optional,
    // so that makeType() could take in a datatype cache if one is
    // available, and otherwise construct datatypes from scratch.
    MPI_Datatype x, y;

    MPI_Datatype x = MpiTypeTraits<X>::makeType ();
    // If creating y doesn't succeed, makeType will throw.  If x is a
    // derived type, this will leak memory.  Thus, if creating y
    // throws, we first free x and then rethrow.  Another way to deal
    // with this would be for makeType() to return
    // RCP<OpaqueWrapper<MPI_Datatype> > with a suitable deallocator,
    // but we chose not to do that because it adds overhead for the
    // overwhelmingly common case of basic datatypes.
    try { 
      MPI_Datatype y = MpiTypeTraits<Y>::makeType ();
    } catch (...) {
      if (MpiTypeTraits<X>::needFree) {
	MPI_Type_free (&x);
      }
      throw;
    }
    MPI_Datatype types[2] = { x, y };

    err = MPI_Type_create_struct (2, blkLens, disps, types, &t);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, 
      "MpiTypeTraits<std::pair<" << TypeNameTraits<X>::name() << ", " 
      << TypeNameTraits<Y>::name() << ">: MPI_Type_create_struct failed.  "
      "This may mean either Teuchos programmer error or something wrong "
      "with the current state of MPI.");

    err = MPI_Type_commit (&t);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, 
      "MpiTypeTraits<std::pair<" << TypeNameTraits<X>::name() << ", " 
      << TypeNameTraits<Y>::name() << ">: MPI_Type_commit failed.  "
      "This may mean either Teuchos programmer error or something wrong "
      "with the current state of MPI.");
    return std::make_pair (t, static_cast<size_t> (1));
  }
};
''', indent)


def emitBoolDecl (indent=0):
    '''Emit declaration of the specialization of MpiTypeTraits for bool.

    indent (integer): the number of spaces to indent.'''

    emit ('''/// \\brief Specialization of MpiTypeTraits for Packet type \c bool.
///
/// \\note For various reasons, we cannot promise that makeType() for
///   the Packet=bool specialization returns a basic datatype.  Even
///   if it does, we cannot promise that the MPI_Op reduction
///   functions for typical Boolean operations are not user-defined.
///   This means that bool should not in general be used in one-sided
///   communication operations, since we cannot guarantee this is
///   possible on all platforms.  In general, you should prefer using
///   C++ built-in integer types to encode Boolean values.  In most
///   cases, \\c int should suffice.  For large arrays of Boolean 
///   values, you might consider using \\c char to save space.
///
/// \warning std::vector<bool> is not an array of bool; it is a bit
///   set with a different representation than std::vector<T> for
///   other types T.  Users who attempt to pass the address of the
///   first entry of an std::vector<bool> into an MPI function will
///   experience undefined and likely undesirable behavior.
///
/// \warning We assume that on all MPI processes, sizeof(bool) has
///   the same value.  It is technically possible to run on clusters
///   with heterogeneous nodes, such that sizeof(bool) may have
///   different values on different nodes.  However, we do not support
///   this (currently uncommon) use case.
template<>
class MpiTypeTraits<bool> {
 public:
  typedef bool packet_type;
  // The right-hand side of this expression is a compile-time constant.
  static const bool needFree = sizeof(bool) != sizeof(char) &&
                               sizeof(bool) != sizeof(short) &&
                               sizeof(bool) != sizeof(int) &&
                               sizeof(bool) != sizeof(long);
  // We assume that every bool on every MPI process has the same size.
  static const bool constantSize = true;

  static std::pair<MPI_Datatype, size_t> makeType (const bool& x);
};
''', indent)

def emitDeclBody (d, indent=0):
    emitGenericDecl (d, indent)
    emitSpecializationDecls (d, indent)
    emitBoolDecl (indent)
    emitStdPairDeclAndDef (indent)

def emitBoolDef (indent=0):
    '''Emit definition of the specialization of MpiTypeTraits for bool.

    indent (integer): the number of spaces to indent.'''

    emit ('''// C++ bool is problematic to MPI because MPI's C++ bindings are
// deprecated (as of MPI 2.2).  This means that we can't rely on a
// native standard MPI_Datatype (in this case, MPI::BOOL) for bool.
// It also means that we don't get standard reduction operators for
// bool; we have to define our own custom operators, or hope that one
// of the standard operators for an integer type of the same size
// works correctly.  Finally, std::vector<bool> is not at all an array
// of bool; it's a bit set.  Users who attempt to pass the address of
// the first entry of an std::vector<bool> into an MPI function will
// experience undefined and likely undesirable behavior.
//
// Note that MPI_C_BOOL is for C99 '_Bool', which is not guaranteed to
// have the same representation (or even the same sizeof value) as C++
// 'bool'.  The 2003 version of the C++ standard (Section 5.3.3/1)
// says that sizeof(bool) is implementation defined, whereas
// sizeof(char) == sizeof(signed char) == sizeof(unsigned char) == 1.
// (There are legitimate performance reasons to want bool to be as big
// as a word; some CPUs are inefficient at extracting and using data
// smaller than a word.) 
//
// We combine two approaches in order to ensure that MpiTypeTraits 
// is correct on all platforms. 
//
// 1. Find the integer type of the same length as bool.
//
// It's reasonable to assume that bool has the same length as some
// kind of standard integer type.  If it doesn't, we fall back to the
// second approach:
//
// 2. Create a derived datatype, using sizeof(bool) contiguous MPI_CHAR.
//
// This works because bool carries information (so sizeof(bool) > 0),
// and sizeof(T) for any type T must be a nonnegative integer.  Thus,
// sizeof(bool) must be a multiple of sizeof(char) (which == 1 by the
// 2003 version of the C++ standard).

namespace { // anonymous
  //! Fall-back for returning a derived MPI_Datatype for Packet=bool.
  static MPI_Datatype makeDerivedTypeForBool () {
    int err = 0;
    MPI_Datatype t;
    // This is always guaranteed to work.  sizeof(T) for any type T
    // must be a nonnegative integer, and thus must be a multiple of
    // sizeof(char) (which == 1 by the C++03 standard).
    err = MPI_Type_contiguous (sizeof(bool), MPI_CHAR, &t);
    if (err != MPI_SUCCESS) {
      TEUCHOS_TEST_FOR_EXCEPTION(err == MPI_ERR_INTERN, std::runtime_error, 
        "MpiTypeTraits<bool>: MPI_Type_contiguous returned MPI_ERR_INTERN.  "
        "This means that MPI failed to acquire memory for allocating the "
        "custom MPI_Datatype.");
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::logic_error, 
        "MpiTypeTraits<bool>: MPI_Type_contiguous failed, "
        "probably due to Teuchos programmer error.  "
	"Please report this to the Teuchos developers.");
    } 
    err = MPI_Type_commit (&t);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, 
      "MpiTypeTraits<bool>: MPI_Type_commit failed.  This may mean either "
      "Teuchos programmer error or something wrong with the current state of "
      "MPI.");
    return t;
  }
} // namespace (anonymous)

template<>
std::pair<MPI_Datatype, size_t> MpiTypeTraits<bool>::makeType () {
  // The compiler should be able to optimize away the false branches.
  if (sizeof(bool) == sizeof(char)) {
    return std::make_pair (MPI_CHAR, static_cast<size_t> (1));
  } else if (sizeof(bool) == sizeof(short)) {
    return std::make_pair (MPI_SHORT, static_cast<size_t> (1));
  } else if (sizeof(bool) == sizeof(int)) {
    return std::make_pair (MPI_INT, static_cast<size_t> (1));
  } else if (sizeof(bool) == sizeof(long)) {
    return std::make_pair (MPI_LONG, static_cast<size_t> (1));
  } else { // Fall-back if bool isn't the same size as a common integer type.
    return std::make_pair (makeDerivedTypeForBool (), static_cast<size_t> (1));
  }
}
''', indent)

def emitDefBody (d, indent=0):
    emitSpecializationDefs (d, indent)
    emitBoolDef (indent)

def emitDeclHeader (d, filename, indent=0):
    '''Emit the entire file of MpiTypeTraits declarations (a .hpp file).

    d (dictionary): return value of makeTypeMap().
    filename (string): name of the header file.
    indent (integer): the number of spaces to indent.'''

    hf = headerizeFilename (filename)
    
    emitCopyrightNotice ()
    topMatter = '''#ifndef {0}
#define {0}

#include <mpi.h>
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <utility> // std::pair

/// \\file {1}
/// \\brief Declaration of MpiTypeTraits and its specializations
/// \\warning This file is automatically generated from a script.
///   Any desired changes should be made to the script, not to this
///   file.  Changes to this file will be lost the next time the
///   script is run.
'''.format (hf, filename) + '''namespace Teuchos {
namespace details {
'''
    emit (topMatter, indent)
    emitDeclBody (d, indent)
    emit ('} // namespace details\n} // namespace Teuchos\n', indent)
    emit ('#endif // #ifndef {0}'.format (hf), indent)

def emitDefSource (d, filename, indent=0):
    '''Emit the entire file of MpiTypeTraits definitions (a .cpp file).

    d (dictionary): return value of makeTypeMap().
    filename (string): name of the header file.
    indent (integer): the number of spaces to indent.'''

    emitCopyrightNotice ()
    topMatter = '''/// \\file {0}
/// \\brief Definition of specializations of MpiTypeTraits
/// \\warning This file is automatically generated from a script.
///   Any desired changes should be made to the script, not to this
///   file.  Changes to this file will be lost the next time the
///   script is run.
'''.format (filename) + '''namespace Teuchos {
namespace details {
'''
    emit (topMatter, indent)
    emitDefBody (d, indent)
    emit ('} // namespace details\n} // namespace Teuchos\n', indent)

def run ():
    '''Generate the two header files mentioned in the module's documentation.

    This writes the header file of MpiTypeTraits declarations
    'Teuchos_MpiTypeTraits.hpp', and the source file of definitions
    'Teuchos_MpiTypeTraits.cpp', for all specializations of
    MpiTypeTraits that this module knows how to generate.  Both files
    are written to the current working directory.'''

    d = makeTypeMap ()
    declsFilename = 'Teuchos_MpiTypeTraits.hpp'
    defsFilename = 'Teuchos_MpiTypeTraits.cpp'

    with open(declsFilename, 'w') as outfile:
        try:
	    sys.stdout = outfile # Rebind stdout to the output file.
	    emitDeclHeader (d, declsFilename)
	finally:
	    # If something went wrong, be sure to restore the original
	    # binding of stdout.
	    sys.stdout = sys.__stdout__

    with open(defsFilename, 'w') as outfile:
	try:
	    sys.stdout = outfile
	    emitDefSource (d, defsFilename)
	finally:
	    sys.stdout = sys.__stdout__


# If running the module as an executable script, just call run().
if __name__ == "__main__":
    if len (sys.argv) > 1:
        raise ValueError ('This script does not currently take any command-line arguments.')
    else:
        run ()
