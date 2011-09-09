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

/** \file   Intrepid_DefaultCubatureFactory.hpp
\brief  Header file for the abstract base class Intrepid::DefaultCubatureFactory.
\author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_DEFAULT_CUBATURE_FACTORY_HPP
#define INTREPID_DEFAULT_CUBATURE_FACTORY_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Cubature.hpp"
#include "Shards_CellTopology.hpp"
#include "Teuchos_RCP.hpp"

/////   list of default cubature includes   /////

#include "Intrepid_CubatureDirectLineGauss.hpp"
#include "Intrepid_CubatureDirectTriDefault.hpp"
#include "Intrepid_CubatureDirectTetDefault.hpp"
#include "Intrepid_CubatureCompositeTet.hpp"
#include "Intrepid_CubatureTensor.hpp"

///// end of list of default cubature includes /////


namespace Intrepid {
  
/** \class Intrepid::DefaultCubatureFactory
    \brief A factory class that generates specific instances of cubatures.
*/
template<class Scalar, class ArrayPoint=FieldContainer<Scalar>, class ArrayWeight=ArrayPoint >
class DefaultCubatureFactory {
  private:

  public:
    
  /** \brief Default constructor.
  */
  DefaultCubatureFactory() {};

  /** \brief Destructor.
  */
  virtual ~DefaultCubatureFactory() {};

  /** \brief Factory method.

      \param cell        [in]    - Cell topology.
      \param degree      [in]    - Array of polynomial degrees, one for each component cubature.

      \return
              - RCP to cubature with given specifications.
  */
  Teuchos::RCP<Cubature<Scalar,ArrayPoint,ArrayWeight> > create(const shards::CellTopology & cellTopology,
                                                                const std::vector<int> & degree);

  /** \brief Factory method.

      \param cell        [in]    - Cell topology.
      \param degree      [in]    - A single polynomial degree, used for all component cubatures.

      \return
              - RCP to cubature with given specifications.
  */
  Teuchos::RCP<Cubature<Scalar,ArrayPoint,ArrayWeight> > create(const shards::CellTopology & cellTopology,
                                                                int   degree);
    
};
  
}// namespace Intrepid

#include "Intrepid_DefaultCubatureFactoryDef.hpp"

#endif
