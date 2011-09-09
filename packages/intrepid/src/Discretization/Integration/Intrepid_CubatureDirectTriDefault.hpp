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

/** \file   Intrepid_CubatureDirectTriDefault.hpp
    \brief  Header file for the Intrepid::CubatureDirectTriDefault class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CUBATURE_DIRECT_TRI_DEFAULT_HPP
#define INTREPID_CUBATURE_DIRECT_TRI_DEFAULT_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_CubatureDirect.hpp"
#include "Teuchos_TestForException.hpp"

/** \def   INTREPID_CUBATURE_TRI_DEFAULT_MAX
  \brief The maximum degree of the polynomial that can be integrated exactly by
         a direct triangle rule of the default type.
*/
// srkenno@sandia.gov 6/21/10:
// see below comment for the enum
#define INTREPID_CUBATURE_TRI_DEFAULT_MAX 20


namespace Intrepid {

/** \class Intrepid::CubatureDirectTriDefault
    \brief Defines direct integration rules on a triangle.
*/
template<class Scalar, class ArrayPoint = FieldContainer<Scalar>, class ArrayWeight = ArrayPoint >
class CubatureDirectTriDefault : public Intrepid::CubatureDirect<Scalar,ArrayPoint,ArrayWeight> {
  public:
  // srkenno@sandia.gov 6/21/10:
  // This indirection is to workaround a compiler bug on the sun platform, 5.7 toolset, SunOS 10.
  enum {INTREPID_CUBATURE_TRI_DEFAULT_MAX_ENUM = INTREPID_CUBATURE_TRI_DEFAULT_MAX};
  private:

  /** \brief Complete set of data defining default cubature rules on a triangle.
  */
  static const CubatureTemplate cubature_data_[INTREPID_CUBATURE_TRI_DEFAULT_MAX_ENUM+1];

  /** \brief Names of templates for frequently used direct cubature rules.
  */
  static const char *cubature_name_;


  public:

  ~CubatureDirectTriDefault() {}

  /** \brief Constructor.

      \param degree           [in]     - The degree of polynomials that are integrated
                                         exactly by this cubature rule. Default: 0.
  */
  CubatureDirectTriDefault(const int degree = 0);

  /** \brief Returns cubature name.
  */
  const char* getName() const;

  /** \brief Exposes cubature data.
  */
  const CubatureTemplate * exposeCubatureData() const;

  /** \brief Returns maximum cubature accuracy.
  */
  int getMaxAccuracy() const;

  /** \brief Exposes cubature data, accessible without construction.
  */
  static const CubatureTemplate (& exposeCubatureDataStatic())[INTREPID_CUBATURE_TRI_DEFAULT_MAX_ENUM+1];

}; // end class CubatureDirect 

template<class Scalar, class ArrayPoint, class ArrayWeight>
inline const CubatureTemplate (& CubatureDirectTriDefault<Scalar,ArrayPoint,ArrayWeight>::exposeCubatureDataStatic())[INTREPID_CUBATURE_TRI_DEFAULT_MAX_ENUM+1] {
  return cubature_data_;
}

} // end namespace Intrepid


// include templated definitions
#include <Intrepid_CubatureDirectTriDefaultDef.hpp>

#endif
