// @HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
#ifndef EPETRAEXT_PACKTRAITS_H
#define EPETRAEXT_PACKTRAITS_H

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
