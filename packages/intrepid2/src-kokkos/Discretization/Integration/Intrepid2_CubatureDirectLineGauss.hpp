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

/** \file   Intrepid_CubatureDirectLineGauss.hpp
    \brief  Header file for the Intrepid2::CubatureDirectLineGauss class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_DIRECT_LINE_GAUSS_HPP__
#define __INTREPID2_CUBATURE_DIRECT_LINE_GAUSS_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_CubatureDirect.hpp"

namespace Intrepid2 {

  /** \class Intrepid2::CubatureDirectLineGauss
      \brief Defines Gauss integration rules on a line.
  */
  template<typename ExecSpaceType>
  class CubatureDirectLineGauss : public CubatureDirect<ExecSpaceType> {
  private:

    /** \brief Complete set of data defining line Gauss(-Legendre) rules.
     */
    static const CubatureDataStatic cubatureDataStatic_[Parameters::MaxCubatureDegreeEdge+1]; // initialized once
    /**/         CubatureData       cubatureData_      [Parameters::MaxCubatureDegreeEdge+1]; // copied to memory space
    
  public:
    
    /** \brief Constructor.
        \param degree           [in]     - The degree of polynomials that are integrated
                                           exactly by this cubature rule. Default: 0.
    */
    CubatureDirectLineGauss(const ordinal_type degree = 0);
    ~CubatureDirectLineGauss() = default;
    
    /** \brief Exposes cubature data.
     */
    CubatureData getCubatureData() const;
    
    /** \brief Returns maximum cubature accuracy.
     */
    ordinal_type getMaxAccuracy() const;
    
    /** \brief Returns cubature name.
     */
    const char* getName() const;
  };
  
} // end namespace Intrepid2

// include templated definitions
#include <Intrepid2_CubatureDirectLineGaussDef.hpp>

#endif
