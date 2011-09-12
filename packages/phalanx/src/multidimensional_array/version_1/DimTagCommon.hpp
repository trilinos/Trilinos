/*
// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

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

