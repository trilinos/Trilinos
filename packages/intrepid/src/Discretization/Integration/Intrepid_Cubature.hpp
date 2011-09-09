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

/** \file   Intrepid_Cubature.hpp
    \brief  Header file for the Intrepid::Cubature class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CUBATURE_HPP
#define INTREPID_CUBATURE_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_Types.hpp"

namespace Intrepid {

/** \class Intrepid::Cubature
    \brief Defines the base class for cubature (integration) rules in Intrepid.

    Cubature template (rule) consists of cubature points and cubature weights.
    Intrepid provides a small collection of frequently used cubature rule templates
    for FEM reconstructions on simplices (edge, tri, tet) and the pyramid cell,
    defined in the derived classes of CubatureDirect.

    For quad, hex, and triprism cells cubature templates are tensor products of CubatureDirect
    templates. The tensor-product cubatures are defined in the derived class CubatureTensor.
*/
template<class Scalar, class ArrayPoint = FieldContainer<Scalar>, class ArrayWeight = ArrayPoint>
class Cubature {
  private:

  public:

  Cubature() {}

  virtual ~Cubature() {}


  /** \brief Returns cubature points and weights
             (return arrays must be pre-sized/pre-allocated).

      \param cubPoints       [out]     - Array containing the cubature points.
      \param cubWeights      [out]     - Array of corresponding cubature weights.
  */
  virtual void getCubature(ArrayPoint  & cubPoints,
                           ArrayWeight & cubWeights) const = 0;


  /** \brief Returns the number of cubature points.
  */
  virtual int getNumPoints() const = 0;


  /** \brief Returns dimension of the integration domain.
  */
  virtual int getDimension() const = 0;


  /** \brief Returns algebraic accuracy (e.g. max. degree of polynomial
             that is integrated exactly). For tensor-product or sparse
             rules, algebraic accuracy for each coordinate direction
             is returned.

             Since the size of the return argument need not be known
             ahead of time, other return options are possible, depending
             on the type of the cubature rule. 
  */
  virtual void getAccuracy(std::vector<int> & accuracy) const = 0;

}; // class Cubature 

}// end namespace Intrepid

#endif
