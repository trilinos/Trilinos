// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_CubatureCompositeTetDef.hpp
    \brief  Definition file for the Intrepid::CubatureCompositeTet class.
    \author Created by P. Bochev, D. Ridzal, D. Day, and J. Ostien.
*/

namespace Intrepid {

template <class Scalar, class ArrayPoint, class ArrayWeight>
CubatureCompositeTet<Scalar,ArrayPoint,ArrayWeight>::CubatureCompositeTet(const int degree) {
  this->degree_    = degree;
  this->dimension_ = 3;
  TEST_FOR_EXCEPTION(degree != 3,
                     std::out_of_range,
                     ">>> ERROR (CubatureDirectTetDefault): No direct cubature rule implemented for the desired polynomial degree.");
} // end constructor



template <class Scalar, class ArrayPoint, class ArrayWeight>
const CubatureTemplate *  CubatureCompositeTet<Scalar,ArrayPoint,ArrayWeight>::exposeCubatureData() const {
  return cubature_data_;
}



template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureCompositeTet<Scalar,ArrayPoint,ArrayWeight>::getMaxAccuracy() const {
  return INTREPID_CUBATURE_COMPOSITE_TET_MAX_ENUM;
}



template <class Scalar, class ArrayPoint, class ArrayWeight>
const char* CubatureCompositeTet<Scalar,ArrayPoint,ArrayWeight>::getName() const {
  return cubature_name_;
} // end getName



template <class Scalar, class ArrayPoint, class ArrayWeight>
const char* CubatureCompositeTet<Scalar,ArrayPoint,ArrayWeight>::cubature_name_ = "INTREPID_CUBATURE_COMPOSITE_TET";


//-------------------------------------------------------------------------------------//
//                          Definition of cubature templates                           //
//-------------------------------------------------------------------------------------//

/*
   Cubature templates for lines are defined the reference Tetrahedron:

   Tetrahedron    -> (0,0,0), (1,0,0), (0,1,0), (0,0,1)
*/

/*
   This static const member contains templates for default tetrahedron rules.
*/

template <class Scalar, class ArrayPoint, class ArrayWeight>
const CubatureTemplate CubatureCompositeTet<Scalar,ArrayPoint,ArrayWeight>::cubature_data_[INTREPID_CUBATURE_COMPOSITE_TET_MAX_ENUM+1] =
{
  // Cubature templates for the reference tet {(0,0,0), (1,0,0), (0,1,0), (0,0,1)}
  //
  {
    1,
    {{1./4., 1./4., 1./4.}},
    {1./6.}
  },
  {
    1,
    {{1./4., 1./4., 1./4.}},
    {1./6.}
  },
  {
    4,
    {{0.1381966011250105151795413165634361882280, 0.1381966011250105151795413165634361882280, 0.1381966011250105151795413165634361882280},
     {0.5854101966249684544613760503096914353161, 0.1381966011250105151795413165634361882280, 0.1381966011250105151795413165634361882280},
     {0.1381966011250105151795413165634361882280, 0.5854101966249684544613760503096914353161, 0.1381966011250105151795413165634361882280},
     {0.1381966011250105151795413165634361882280, 0.1381966011250105151795413165634361882280, 0.5854101966249684544613760503096914353161}},
    {1./24.,
     1./24.,
     1./24.,
     1./24.}
  },
  {
    5,
    {{1./4., 1./4., 1./4.},
     {1./6., 1./6., 1./6.},
     {1./6., 1./6., 1./2.},
     {1./6., 1./2., 1./6.},
     {1./2., 1./6., 1./6.}},
    { 1./16.,
      5./192.,
      5./192.,
      5./192.,
      5./192.}
  }

};
    
} // end namespace Intrepid
