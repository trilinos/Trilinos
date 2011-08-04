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

using Teuchos::Array;
using Teuchos::RCP;

namespace Intrepid {
    
/** \class  Intrepid::TensorBasis
    \brief  An abstract base class that defines interface for bases
            that are tensor products of simpler bases.
*/
template<class Scalar, class ArrayScalar>
class TensorBasis: public Basis<Scalar,ArrayScalar> {
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
