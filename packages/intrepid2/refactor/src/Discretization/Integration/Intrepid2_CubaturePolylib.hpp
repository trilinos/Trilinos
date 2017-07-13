// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_CubaturePolylib.hpp
    \brief  Header file for the Intrepid2::CubaturePolylib class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_POLYLIB_HPP__
#define __INTREPID2_CUBATURE_POLYLIB_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_CubatureDirect.hpp"

#include "Intrepid2_Polylib.hpp"


namespace Intrepid2 {

  /** \class Intrepid2::CubaturePolylib
      \brief Utilizes cubature (integration) rules contained in the library Polylib
      (Spencer Sherwin, Aeronautics, Imperial College London) within Intrepid.

      They are based on zeros of Jacobi polynomials,
      e.g. Legendre (alpha=beta=0, default), Chebyshev (alpha=beta=-0.5), etc.
      They are given on the interval [-1,1] and are optimal with respect to
      the following requirements, yielding 4 subclasses:
      \li Gauss             - no restrictions (all integration points are contained in the open interval (-1,1))
      \li Gauss-Radau-Left  - left-most integration point is fixed at -1
      \li Gauss-Radau-Right - right-most integration point is fixed at +1
      \li Gauss-Lobatto     - left-most and right-most integration points are fixed at -1 and +1, respectively
  */
  template<typename ExecSpaceType = void,
           typename pointValueType = double,
           typename weightValueType = double>
  class CubaturePolylib
    : public CubatureDirect<ExecSpaceType,pointValueType,weightValueType> {
  public:
    typedef typename CubatureDirect<ExecSpaceType,pointValueType,weightValueType>::pointViewType  pointViewType; 
    typedef typename CubatureDirect<ExecSpaceType,pointValueType,weightValueType>::weightViewType weightViewType;

    CubaturePolylib( const ordinal_type degree,
                     const EPolyType polytype = POLYTYPE_GAUSS,
                     const double alpha = 0.0,
                     const double beta = 0.0 );

    /** \brief Returns cubature name.
     */
    virtual
    const char*
    getName() const {
      return "CubaturePolylib";
    }


  };

}

#include "Intrepid2_CubaturePolylibDef.hpp"

#endif
