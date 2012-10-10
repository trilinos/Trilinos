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

/** \file   Intrepid_Basis.hpp
    \brief  Header file for the abstract base class Intrepid::Basis.
    \author Created by P. Bochev and D. Ridzal.
 */

#ifndef INTREPID_TENSORBASIS_HPP
#define INTREPID_TENSORBASIS_HPP
#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Intrepid_Utils.hpp"
#include "Shards_CellTopology.hpp"
#include "Intrepid_Basis.hpp"

using Teuchos::Array;
using Teuchos::RCP;

namespace Intrepid {
    
/** \class  Intrepid::TensorBasis
    \brief  An abstract base class that defines interface for bases
            that are tensor products of simpler bases.
*/
template<class Scalar, class ArrayScalar>
class TensorBasis: public Basis<Scalar,ArrayScalar> 
{
private:

  
protected:  
  Array< Array< RCP< Basis< Scalar , ArrayScalar > > > > bases_; 

  void setBases( Array< Array< RCP< Basis< Scalar , ArrayScalar > > > > & bases )
  {
    bases_.resize( bases.size() );
    for (int i=0;i<bases.size();i++)
      {
	bases_[i].resize( bases[i].size() );
	for (int j=0;j<bases[i].size();j++)
	  {
	    bases_[i][j] = bases[i][j];
	  }
      }
  }

public:
  Array< Array< RCP< Basis< Scalar , ArrayScalar > > > > &getBases() 
  { return bases_; }

  /** \brief  Destructor
   */
  virtual ~TensorBasis() {}

}; // class Tensor Basis

}
#endif
