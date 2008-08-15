/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

#ifndef util_DimTagCommon_hpp
#define util_DimTagCommon_hpp

#include <Dimension.hpp>

namespace phdmesh {

/** \brief Cartesian coordinates multi-index map ordinate tag.
 */
class Cartesian : public DimTag {
public:

  enum { X = 0 , Y = 1 , Z = 2 };

  const char * name() const ;

  std::string to_string( unsigned size , unsigned index ) const ;

  unsigned to_index( unsigned size , const std::string & ) const ;

  static const Cartesian & descriptor();
};

/** \brief Cylindrical coordinates multi-index map ordinate tag.
 */
class Cylindrical : public DimTag {
public:

  enum { Radius = 0 , R = 0 ,
         Angle  = 1 , A = 1 ,
                      Z = 2 };

  const char * name() const ;

  std::string to_string( unsigned size , unsigned index ) const ;

  unsigned to_index( unsigned size , const std::string & ) const ;

  static const Cylindrical & descriptor();
};

//----------------------------------------------------------------------

/** \brief Diagonal Cartesian tensor multi-index map ordinate tag.
 */
class DiagonalTensor : public DimTag {
public:

  enum { XX = 0 , YY = 1 , ZZ = 2 };

  const char * name() const ;

  std::string to_string( unsigned size , unsigned index ) const ;

  unsigned to_index( unsigned size , const std::string & ) const ;

  static const Cylindrical & descriptor();
};

/** \brief Symmetric Cartesian tensor multi-index map ordinate tag.
 */
template< unsigned > class SymmetricTensor ;

template<>
class SymmetricTensor<3> : public DimTag {
public:

  enum { XX = 0 ,
         XY = 3 , YY = 1 ,
         XZ = 4 , YZ = 5 , ZZ = 2 };

  const char * name() const ;

  std::string to_string( unsigned size , unsigned index ) const ;

  unsigned to_index( unsigned size , const std::string & ) const ;

  static const SymmetricTensor<3> & descriptor();
};

template<>
class SymmetricTensor<2> : public DimTag {
public:

  enum { XX = 0 ,
         XY = 2 , YY = 1 };

  const char * name() const ;

  std::string to_string( unsigned size , unsigned index ) const ;

  unsigned to_index( unsigned size , const std::string & ) const ;

  static const SymmetricTensor<3> & descriptor();
};

//----------------------------------------------------------------------

}

#endif

