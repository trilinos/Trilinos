//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//@HEADER
#ifndef EPETRAEXT_PACKTRAITS_H
#define EPETRAEXT_PACKTRAITS_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

// ---------- Standard Includes ----------

#include <string>
#include <vector>

// ----------   Xyce Includes   ----------

// ----------  Other Includes   ----------

// ----------  Fwd Declarations ----------

namespace EpetraExt {

  /** Traits for packing and unpacking of data into char buffers for
      communication.  Should be replaced by something more generic in
      Teuchos.
   */
template <typename T>
struct PackTraits
{
  /** Returns size in byte necessary to pack datatype

      @param object Input, object to be sized for packing.

      @return Size in bytes needed for packed object.
   */
  static size_t size( T const & object )
  { return object.packedByteCount(); }

  /** Packs object into char buffer

      @param object data to be packed.
      @param buf buffer to be used for packed data.
      @param size total size of buffer (for overrun check).
      @param pos current position in buffer for packing.
   */
  static void pack( T const & object, char * buf, size_t size, int & pos )
  { object.pack( buf, size, pos ); }

  /** Unpacks object from char buffer

      @param object data to be unpacked.
      @param buf buffer to be used for unpacking data.
      @param size total size of buffer (for overrun check).
      @param pos current position in buffer for unpacking.
   */
  static void unpack( T & object, char * buf, size_t size, int & pos )
  { object.unpack( buf, size, pos ); }
};

///
/** Full specialization of <tt>PackTraits</tt> for std::string
 */
template <>
struct PackTraits<std::string>
{
  static size_t size( std::string const & object )
  { return object.length() + sizeof(size_t); }

  static void pack( std::string const & object, char * buf, size_t size, int & pos )
  {
    size_t len = object.length();
    std::memcpy( buf+pos, &len, sizeof(size_t) );
    pos += sizeof(size_t);
    std::memcpy( buf+pos, object.c_str(), len );
    pos += len;
  }

  static void unpack( std::string & object, char * buf, size_t size, int & pos )
  {
    size_t len;
    std::memcpy( &len, buf+pos, sizeof(size_t) );
    pos += sizeof(size_t);
    object = std::string( buf+pos, len );
    pos += len;
  }
};

///
/** Partial specialization of <tt>PackTraits</tt> for std::vector<>
    containing a primitive type
 */
template <typename T>
struct PackTraits< std::vector<T> >
{
  static size_t size( std::vector<T> const & object )
  { return object.size() * sizeof(T) + sizeof(size_t); }

  static void pack( std::vector<T> const & object, char * buf, size_t size, int & os )
  {
    size_t len = object.size();
    std::memcpy( buf+pos, &len, sizeof(size_t) );
    pos += sizeof(size_t);
    std::memcpy( buf+pos, &object[0], len*sizeof(T) );
    pos += len*sizeof(T);
  }

  static void unpack( std::vector<T> & object, char * buf, int size, int & pos )
  {
    size_t len;
    std::memcpy( &len, buf+pos, sizeof(size_t) );
    pos += sizeof(size_t);
    object.resize(len);
    std::memcpy( &object[0], buf+pos, len*sizeof(T) );
    pos += len*sizeof(T);
  }
};

} //namespace EpetraExt

#endif
