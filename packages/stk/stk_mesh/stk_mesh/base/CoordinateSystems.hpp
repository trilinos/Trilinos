// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef stk_mesh_fem_CoordinateSystems_hpp
#define stk_mesh_fem_CoordinateSystems_hpp

//----------------------------------------------------------------------

#include <stk_mesh/base/Types.hpp>
#include <Shards_Array.hpp>             // for ArrayDimTag::size_type, etc
#include <string>                       // for string

#define DEPRECATED_SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( ADT ) \
  class ADT : public shards::ArrayDimTag { \
  public: \
    const char * name() const { static const char n[] = # ADT; return n; } \
    const ADT & tag() { static const ADT self ; return self ; } \
  private: \
  ~ADT() {} \
  ADT() { \
    std::cerr << "Warning: The stk::mesh::" #ADT " type is deprecated and will soon be removed." << std::endl; \
  } \
    ADT( const ADT & ); \
    ADT & operator = ( const ADT & ); \
  };

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_field_dimension_tags
 *  \{
 *
 * ArrayDimTags are required for multidimensional Field types; they specify
 * the dimensions of the field and the intent of each dimension. Note that
 * scalar Field types do not involve ArrayDimTags. This file defines a number
 * of ArrayDimTags that we believe will be widely useful for STK users. Clients
 * have the freedom to define their own ArrayDimTags as well.
 *
 * Example use of Cartesian ArrayDimTag to create a field type:
 *   stk::mesh::Field<double, stk::mesh::Cartesian>
 */

DEPRECATED_SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( SimpleArrayTag )

/**
 *   \brief Implement an shards::ArrayDimTag for Cartesian coordinate dimensions.
 *
 *   A Cartesian coordinate has up to three dimensions in X, Y, Z order.
 */
struct Cartesian3d : public shards::ArrayDimTag {

  enum { Size = 3 };                    ///< default size

  enum { X = 0 , Y = 1 , Z = 2 };       ///< Identifiers for each dimension

  const char * name() const {
    static const char n[] = "Cartesian3d";
    return n;
  }

  static const Cartesian3d & tag() {
    static const Cartesian3d self;
    return self;
  }

private:
  Cartesian3d() {
    std::cerr << "Warning: The stk::mesh::Cartesian3d type is deprecated and will soon be removed." << std::endl;
  }
  Cartesian3d( const Cartesian3d & );
  Cartesian3d & operator = ( const Cartesian3d & );
};

/**
 *   \brief Implement an shards::ArrayDimTag for Cartesian 2d coordinate dimensions.
 *
 *   A Cartesian coordinate has up to two dimensions in X, Y order.
 */
struct Cartesian2d: public shards::ArrayDimTag {

  enum { Size = 2 };                    ///< default size

  enum { X = 0 , Y = 1 };       ///< Identifiers for each dimension

  const char * name() const {
    static const char n[] = "Cartesian2d";
    return n;
  }

  static const Cartesian2d & tag() {
    static const Cartesian2d self;
    return self;
  }

private:
  Cartesian2d() {
    std::cerr << "Warning: The stk::mesh::Cartesian2d type is deprecated and will soon be removed." << std::endl;
  }
  Cartesian2d( const Cartesian2d & );
  Cartesian2d & operator = ( const Cartesian2d & );
};

typedef Cartesian3d Cartesian;
/**
 *   \brief Implement an shards::ArrayDimTag for Cylindrical coordinate dimensions.
 *
 *   A Cylindral coordinate has up to three dimensions in
 *   radius, angle, and longitudinal-distance order.
 */
struct Cylindrical : public shards::ArrayDimTag {

  enum { Radius = 0 , R = 0 ,           ///< Identifiers for each dimension
         Angle = 1 ,  A = 1 ,
         Z = 2 };

  const char * name() const {
    static const char n[] = "Cylindrical";
    return n;
  }

  static const Cylindrical & tag() {
    static const Cylindrical self;
    return self;
  }

private:
  Cylindrical() {
    std::cerr << "Warning: The stk::mesh::Cylindrical type is deprecated and will soon be removed." << std::endl;
  }
  Cylindrical( const Cylindrical & );
  Cylindrical & operator = ( const Cylindrical & );
};

/**
 *  \brief Implement an shards::ArrayDimTag for FullTensor.
 *
 * \todo REFACTOR  Where should FullTensor live, in the application,
 *                 in the toolkit or a common application header?
 */
struct FullTensor36 : public shards::ArrayDimTag {

  enum { Size = 9 };

/*
 * Note on Ordering:  This is the ordering as used in the old 
 * Sierra Framework and is somewhat standard in that a symmetric
 * tensor is the first six values of a full tensor and a diagonal
 * only tensor is the first three values of that.  
 *
 * I think this is actually in ERROR in that (XZ,YX,ZY) SHOULD
 * be (6,7,8) NOT (8,6,7).  But backwards compatibility is useful.
 *
 * \todo Look at the proper ordering of a full second order tensor.
 */
  enum { XX = 0 , XY = 3 , XZ = 8 ,
         YX = 6 , YY = 1 , YZ = 4 ,
         ZX = 5 , ZY = 7 , ZZ = 2 };

  const char * name() const {
    static const char n[] = "FullTensor36";
    return n;
  }

  static const FullTensor36 & tag() {
    static const FullTensor36 self;
    return self;
  }

private:
  FullTensor36() {
    std::cerr << "Warning: The stk::mesh::FullTensor36 type is deprecated and will soon be removed." << std::endl;
  }
  FullTensor36( const FullTensor36 & );
  FullTensor36 & operator = ( const FullTensor36 & );
};

typedef FullTensor36 FullTensor;

/**
 *  \brief Implement an shards::ArrayDimTag for FullTensor.
 */
struct FullTensor22 : public shards::ArrayDimTag {

  enum { Size = 4 };

  enum { XX = 0 , XY = 2 , 
         YX = 3 , YY = 1};

  const char * name() const {
    static const char n[] = "FullTensor22";
    return n;
  }

  static const FullTensor22 & tag() {
    static const FullTensor22 self;
    return self;
  }

private:
  FullTensor22() {
    std::cerr << "Warning: The stk::mesh::FullTensor22 type is deprecated and will soon be removed." << std::endl;
  }
  FullTensor22( const FullTensor22 & );
  FullTensor22 & operator = ( const FullTensor22 & );
};

//----------------------------------------------------------------------

/**
 *  \brief Implement an shards::ArrayDimTag for SymmetricTensor.
 *
 * \todo REFACTOR  Where should SymmetricTensor live, in the application,
 *                 in the toolkit or a common application header?
 */
struct SymmetricTensor33 : public shards::ArrayDimTag {

  enum { Size = 6 };

  enum { XX = 0 , XY = 3,  XZ = 5,
         YX = 3 , YY = 1,  YZ = 4, 
         ZX = 5 , ZY = 4,  ZZ = 2};

  const char * name() const {
    static const char n[] = "SymmetricTensor33";
    return n;
  }

  static const SymmetricTensor33 & tag() {
    static const SymmetricTensor33 self;
    return self;
  }

private:
  SymmetricTensor33() {
    std::cerr << "Warning: The stk::mesh::SymmetricTensor33 type is deprecated and will soon be removed." << std::endl;
  }
  SymmetricTensor33( const SymmetricTensor33 & );
  SymmetricTensor33 & operator = ( const SymmetricTensor33 & );
};

typedef SymmetricTensor33 SymmetricTensor;

/**
 *  \brief Implement an shards::ArrayDimTag for SymmetricTensor.
 *    
 *  SymmetricTensor31 is an axisymmetric tensor in 3D.  It 
 *  has the radius and height of the cylindrical coordinate
 *  system but with no theta coordinate.
 */
struct SymmetricTensor31 : public shards::ArrayDimTag {

  enum { Size = 4 };

  enum { rr = 0 , rz = 2 , 
         zr = 3 , zz = 1};

  const char * name() const {
    static const char n[] = "SymmetricTensor31";
    return n;
  }

  static const SymmetricTensor31 & tag() {
    static const SymmetricTensor31 self;
    return self;
  }

private:
  SymmetricTensor31() {
    std::cerr << "Warning: The stk::mesh::SymmetricTensor31 type is deprecated and will soon be removed." << std::endl;
  }
  SymmetricTensor31( const SymmetricTensor31 & );
  SymmetricTensor31 & operator = ( const SymmetricTensor31 & );
};

/**
 *  \brief Implement an shards::ArrayDimTag for SymmetricTensor.
 */
struct SymmetricTensor21 : public shards::ArrayDimTag {

  enum { Size = 3 };

  enum { XX = 0 , XY = 2 , 
         YX = 2 , YY = 1 };

  const char * name() const {
    static const char n[] = "SymmetricTensor21";
    return n;
  }

  static const SymmetricTensor21 & tag() {
    static const SymmetricTensor21 self;
    return self;
  }

private:
  SymmetricTensor21() {
    std::cerr << "Warning: The stk::mesh::SymmetricTensor21 type is deprecated and will soon be removed." << std::endl;
  }
  SymmetricTensor21( const SymmetricTensor21 & );
  SymmetricTensor21 & operator = ( const SymmetricTensor21 & );
};

/**
 *  \brief Implement an shards::ArrayDimTag for AsymmetricTensor.
 *
 * Note: I think by Axymmetric is ment Skew-symmetric. 
 * Asymmetric would be any non-symmetric tensor while skew-symmetric
 * means it is equal to the negative of it's transpose. This 
 * forces the diagonals to be zero and only the three off-diagonal
 * elements are useful.
 */
struct AsymmetricTensor03 : public shards::ArrayDimTag {

  enum { Size = 3 };

  enum {  /* XX = 0 */  XY = 0 ,   XZ = 2 ,
             YX = 0 ,/* YY = 0 */  YZ = 1 , 
             ZX = 2 ,   ZY = 1  /* ZZ=0 */ };

  const char * name() const {
    static const char n[] = "AsymmetricTensor03";
    return n;
  }

  static const AsymmetricTensor03 & tag() {
    static const AsymmetricTensor03 self;
    return self;
  }

private:
  AsymmetricTensor03() {
    std::cerr << "Warning: The stk::mesh::AsymmetricTensor03 type is deprecated and will soon be removed." << std::endl;
  }
  AsymmetricTensor03( const AsymmetricTensor03 & );
  AsymmetricTensor03 & operator = ( const AsymmetricTensor03 & );
};

typedef AsymmetricTensor03 AsymmetricTensor;

/**
 *  \brief Implement an shards::ArrayDimTag for Matrix.
 */
struct Matrix22 : public shards::ArrayDimTag {

  enum { Size = 4 };

  enum { XX = 0 , XY = 2 , 
         YX = 1,  YY = 3 };

  const char * name() const {
    static const char n[] = "Matrix22";
    return n;
  }

  static const Matrix22 & tag() {
    static const Matrix22 self;
    return self;
  }

private:
  Matrix22() {
    std::cerr << "Warning: The stk::mesh::Matrix22 type is deprecated and will soon be removed." << std::endl;
  }
  Matrix22( const Matrix22 & );
  Matrix22 & operator = ( const Matrix22 & );
};

/**
 *  \brief Implement an shards::ArrayDimTag for Matrix.
 */
struct Matrix33 : public shards::ArrayDimTag {

  enum { Size = 9 };

  enum { XX = 0 , XY = 3 , XZ = 6 ,
         YX = 1 , YY = 4 , YZ = 7 ,
         ZX = 2 , ZY = 5 , ZZ = 8 };

  const char * name() const {
    static const char n[] = "Matrix33";
    return n;
  }

  static const Matrix33 & tag() {
    static const Matrix33 self;
    return self;
  }

private:
  Matrix33() {
    std::cerr << "Warning: The stk::mesh::Matrix33 type is deprecated and will soon be removed." << std::endl;
  }
  Matrix33( const Matrix33 & );
  Matrix33 & operator = ( const Matrix33 & );
};

typedef Matrix33 Matrix;

//----------------------------------------------------------------------

/** \} */

} //namespace mesh
} //namespace stk

#endif //stk_mesh_fem_CoordinateSystems_hpp
