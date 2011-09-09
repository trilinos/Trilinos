// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
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
