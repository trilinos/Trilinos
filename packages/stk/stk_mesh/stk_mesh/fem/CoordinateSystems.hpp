/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_fem_CoordinateSystems_hpp
#define stk_mesh_fem_CoordinateSystems_hpp

//----------------------------------------------------------------------

#include <Shards_Array.hpp>

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_field_dimension_tags
 *  \{
 */

/**
 *   \brief Implement an shards::ArrayDimTag for Cartesian coordinate dimensions.
 *
 *   A Cartesian coordinate has up to three dimensions in X, Y, Z order.
 */
struct Cartesian3d : public shards::ArrayDimTag {

  enum { Size = 3 };                    ///< default size

  enum { X = 0 , Y = 1 , Z = 2 };       ///< Identifiers for each dimension

  const char * name() const ;
  std::string to_string( size_type size , size_type index ) const ;
  size_type    to_index(  size_type size , const std::string & ) const ;
  static const Cartesian3d & tag();       ///< Singleton

private:
  Cartesian3d() {}
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

  const char * name() const ;
  std::string to_string( size_type size , size_type index ) const ;
  size_type    to_index(  size_type size , const std::string & ) const ;
  static const Cartesian2d & tag();       ///< Singleton

private:
  Cartesian2d() {}
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

  const char * name() const ;
  std::string to_string( size_type size , size_type index ) const ;
  size_type    to_index(  size_type size , const std::string & ) const ;
  static const Cylindrical & tag(); ///< Singleton

private:
  Cylindrical() {}
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

  const char * name() const ;
  std::string to_string( size_type, size_type) const  ;
  size_type    to_index(  size_type, const std::string & ) const  ;
  static const FullTensor36 & tag(); ///< Singleton

private:
  FullTensor36() {}
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

  const char * name() const ;
  std::string to_string( size_type, size_type) const  ;
  size_type    to_index(  size_type, const std::string & ) const  ;
  static const FullTensor22 & tag(); ///< Singleton

private:
  FullTensor22() {}
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

  const char * name() const  ;
  std::string to_string( size_type, size_type) const  ;
  size_type    to_index(  size_type , const std::string & ) const ;
  static const SymmetricTensor33 & tag(); ///< Singleton

private:
  SymmetricTensor33() {}
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

  const char * name() const  ;
  std::string to_string( size_type, size_type) const  ;
  size_type    to_index(  size_type , const std::string & ) const ;
  static const SymmetricTensor31 & tag(); ///< Singleton

private:
  SymmetricTensor31() {}
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

  const char * name() const  ;
  std::string to_string( size_type, size_type) const  ;
  size_type    to_index(  size_type , const std::string & ) const ;
  static const SymmetricTensor21 & tag(); ///< Singleton

private:
  SymmetricTensor21() {}
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

  const char * name() const  ;
  std::string to_string( size_type, size_type) const  ;
  size_type    to_index(  size_type , const std::string & ) const ;
  static const AsymmetricTensor03 & tag(); ///< Singleton

private:
  AsymmetricTensor03() {}
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

  const char * name() const  ;
  std::string to_string( size_type, size_type) const  ;
  size_type    to_index(  size_type , const std::string & ) const ;
  static const Matrix22 & tag(); ///< Singleton

private:
  Matrix22() {}
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

  const char * name() const  ;
  std::string to_string( size_type, size_type) const  ;
  size_type    to_index(  size_type , const std::string & ) const ;
  static const Matrix33 & tag(); ///< Singleton

private:
  Matrix33() {}
  Matrix33( const Matrix33 & );
  Matrix33 & operator = ( const Matrix33 & );
};

//----------------------------------------------------------------------

/** \} */

} //namespace mesh
} //namespace stk

#endif //stk_mesh_fem_CoordinateSystems_hpp
