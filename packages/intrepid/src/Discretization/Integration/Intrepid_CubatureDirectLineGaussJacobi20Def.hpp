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

/** \file   Intrepid_CubatureDirectLineGaussJacobi20Def.hpp
    \brief  Definition file for the Intrepid::CubatureDirectLineGaussJacobi20 class.
    \author Created by P. Bochev, D. Ridzal, and D. Day.
*/

namespace Intrepid {

template <class Scalar, class ArrayPoint, class ArrayWeight>
CubatureDirectLineGaussJacobi20<Scalar,ArrayPoint,ArrayWeight>::CubatureDirectLineGaussJacobi20(const int degree) {
  this->degree_    = degree;
  this->dimension_ = 1;
  TEUCHOS_TEST_FOR_EXCEPTION((degree < 0) || (degree > INTREPID_CUBATURE_LINE_GAUSSJACOBI20_MAX_ENUM),
                     std::out_of_range,
                     ">>> ERROR (CubatureDirectLineGaussJacobi20): No cubature rule implemented for the desired polynomial degree.");
} // end constructor



template <class Scalar, class ArrayPoint, class ArrayWeight>
const CubatureTemplate *  CubatureDirectLineGaussJacobi20<Scalar,ArrayPoint,ArrayWeight>::exposeCubatureData() const {
  return cubature_data_;
}



template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureDirectLineGaussJacobi20<Scalar,ArrayPoint,ArrayWeight>::getMaxAccuracy() const {
  return INTREPID_CUBATURE_LINE_GAUSSJACOBI20_MAX_ENUM;
}



template <class Scalar, class ArrayPoint, class ArrayWeight>
const char* CubatureDirectLineGaussJacobi20<Scalar,ArrayPoint,ArrayWeight>::getName() const {
  return cubature_name_;
} // end getName



template <class Scalar, class ArrayPoint, class ArrayWeight>
const char* CubatureDirectLineGaussJacobi20<Scalar,ArrayPoint,ArrayWeight>::cubature_name_ = "INTREPID_CUBATURE_LINE_GAUSSJACOBI20";


//-------------------------------------------------------------------------------------//
//                          Definition of cubature templates                           //
//-------------------------------------------------------------------------------------//

/*
   Cubature templates for lines are defined the reference cell:

   Line -> (-1,0,0),(1,0,0)
*/

/*
   This static const member contains templates for GaussJacobi20(-Legendre) rules.
*/

template <class Scalar, class ArrayPoint, class ArrayWeight>
const CubatureTemplate CubatureDirectLineGaussJacobi20<Scalar,ArrayPoint,ArrayWeight>::cubature_data_[INTREPID_CUBATURE_LINE_GAUSSJACOBI20_MAX_ENUM+1] =
{

  // Collection of GaussJacobi20 rules on [-1,1]
  // The rule with array index i is exact for polynomials up to order i
  {
    1,
    {{-0.5, 0.0, 0.0}},
    {2.66666666666666666666666666}
  },
  {
    1,
    {{-0.5, 0.0, 0.0}},
    {2.66666666666666666666666666}
  },
  {
    2,
    {{-7.549703546891172e-1, 0.0, 0.0},
     {8.830368802245062e-2, 0.0, 0.0}},
    {1.860379610028064,
     8.062870566386037e-01}
  },
  {
    2,
    {{-7.549703546891172e-1, 0.0, 0.0},
     {8.830368802245062e-2, 0.0, 0.0}},
    {1.860379610028064,
     8.062870566386037e-01}
  },
  {
    3,
    {{-8.540119518537008e-01, 0.0, 0.0},
     {-3.059924679232963e-01, 0.0, 0.0},
     { 4.100044197769969e-01, 0.0, 0.0}},
    {1.257090888519093e+00,
     1.169970154078928e+00,
     2.396056240686456e-01}
  },
  {
    3,
    {{-8.540119518537008e-01, 0.0, 0.0},
     {-3.059924679232963e-01, 0.0, 0.0},
     { 4.100044197769969e-01, 0.0, 0.0}},
    {1.257090888519093e+00,
     1.169970154078928e+00,
     2.396056240686456e-01}
  },
  {
    4,
    {{-9.029989011060054e-01, 0.0, 0.0},
     {-5.227985248962754e-01, 0.0, 0.0},
     {3.409459020873505e-02, 0.0, 0.0},
     {5.917028357935457e-01, 0.0, 0.0}},
    {8.871073248902235e-01,
     1.147670318393715e+00,
     5.490710973833849e-01,
     8.281792599934450e-02}
  },
  {
    4,
    {{-9.029989011060054e-01, 0.0, 0.0},
     {-5.227985248962754e-01, 0.0, 0.0},
     {3.409459020873505e-02, 0.0, 0.0},
     {5.917028357935457e-01, 0.0, 0.0}},
    {8.871073248902235e-01,
     1.147670318393715e+00,
     5.490710973833849e-01,
     8.281792599934450e-02}
  },
  {
    5,
    {{-9.308421201635699e-01, 0.0, 0.0},
     {-6.530393584566085e-01, 0.0, 0.0},
     {-2.202272258689614e-01, 0.0, 0.0},
     {2.686669452617736e-01, 0.0, 0.0},
     {7.021084258940329e-01, 0.0, 0.0}},
    {6.541182742861678e-01,
     1.009591695199292e+00,
     7.136012897727201e-01,
     2.564448057836956e-01,
     3.291060162479211e-02}
  },
  {
    5,
    {{-9.308421201635699e-01, 0.0, 0.0},
     {-6.530393584566085e-01, 0.0, 0.0},
     {-2.202272258689614e-01, 0.0, 0.0},
     {2.686669452617736e-01, 0.0, 0.0},
     {7.021084258940329e-01, 0.0, 0.0}},
    {6.541182742861678e-01,
     1.009591695199292e+00,
     7.136012897727201e-01,
     2.564448057836956e-01,
     3.291060162479211e-02}
  },
  {
    6,
    {{-9.481908898126656e-01, 0.0, 0.0},
     {-7.368721166840297e-01, 0.0, 0.0},
     {-3.951261639542174e-01, 0.0, 0.0},
     {1.807282632950432e-02, 0.0, 0.0},
     {4.313622546234276e-01, 0.0, 0.0},
     {7.736112323551237e-01, 0.0, 0.0}},
    {5.003096218126469e-01,
     8.590119978942462e-01,
     7.566174939883307e-01,
     4.103165690369299e-01,
     1.257623774795603e-01,
     1.464860645495425e-02}
  },
  {
    6,
    {{-9.481908898126656e-01, 0.0, 0.0},
     {-7.368721166840297e-01, 0.0, 0.0},
     {-3.951261639542174e-01, 0.0, 0.0},
     {1.807282632950432e-02, 0.0, 0.0},
     {4.313622546234276e-01, 0.0, 0.0},
     {7.736112323551237e-01, 0.0, 0.0}},
    {5.003096218126469e-01,
     8.590119978942462e-01,
     7.566174939883307e-01,
     4.103165690369299e-01,
     1.257623774795603e-01,
     1.464860645495425e-02}
  } // end GaussJacobi20

};
    
} // end namespace Intrepid

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

