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

/** \file   Intrepid_CubatureTensorPyrDef.hpp
    \brief  Definition file for the Intrepid::CubatureTensorPyr class.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

template <class Scalar, class ArrayPoint, class ArrayWeight>
CubatureTensorPyr<Scalar,ArrayPoint,ArrayWeight>::CubatureTensorPyr(std::vector< Teuchos::RCP<Cubature<Scalar,ArrayPoint,ArrayWeight> > > cubatures) :
                                                            CubatureTensor<Scalar,ArrayPoint,ArrayWeight>(cubatures) {}


template <class Scalar, class ArrayPoint, class ArrayWeight>
CubatureTensorPyr<Scalar,ArrayPoint,ArrayWeight>::CubatureTensorPyr(Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature1,
                                                            Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature2) :
                                                            CubatureTensor<Scalar,ArrayPoint,ArrayWeight>(cubature1,cubature2) {}


template <class Scalar, class ArrayPoint, class ArrayWeight>
CubatureTensorPyr<Scalar,ArrayPoint,ArrayWeight>::CubatureTensorPyr(Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature1,
                                                            Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature2,
                                                            Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature3) :
                                                              	  CubatureTensor<Scalar,ArrayPoint,ArrayWeight>(cubature1,cubature2, cubature3) {}


template <class Scalar, class ArrayPoint, class ArrayWeight>
CubatureTensorPyr<Scalar,ArrayPoint,ArrayWeight>::CubatureTensorPyr(Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature, int n) :
CubatureTensor<Scalar,ArrayPoint,ArrayWeight>(cubature, n) {}


template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureTensorPyr<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint  & cubPoints,
                                                                  ArrayWeight & cubWeights) const {
  CubatureTensor<Scalar,ArrayPoint,ArrayWeight>::getCubature(cubPoints, cubWeights);
  
  int numCubPoints = getNumPoints();

  for(int i=0; i<numCubPoints; i++)
  {
	  Scalar zTerm = 0.5*(1.0 - cubPoints(i,2));
	  cubPoints(i,0) = zTerm*cubPoints(i,0); //   zTerm*(cubPoints(i,0)+1.0 )-1.0;
	  cubPoints(i,1) = zTerm*cubPoints(i,1); //zTerm*(cubPoints(i,1)+1.0 )-1.0;
	  cubPoints(i,2) = 1.0 - zTerm;   // 0.5 + 0.5 z
	  cubWeights(i) /= 8.;
  }

} // end getCubature

template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureTensorPyr<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint& cubPoints,
                                                                   ArrayWeight& cubWeights,
                                                                   ArrayPoint& cellCoords) const
{
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (CubatureTensorPyr): Cubature defined in reference space calling method for physical space cubature.");
}

template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureTensorPyr<Scalar,ArrayPoint,ArrayWeight>::getNumPoints() const {
	return CubatureTensor<Scalar,ArrayPoint,ArrayWeight>::getNumPoints();
} // end getNumPoints


template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureTensorPyr<Scalar,ArrayPoint,ArrayWeight>::getDimension() const {
	return CubatureTensor<Scalar,ArrayPoint,ArrayWeight>::getDimension();
} // end dimension


template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureTensorPyr<Scalar,ArrayPoint,ArrayWeight>::getAccuracy(std::vector<int> & degree) const {
	CubatureTensor<Scalar,ArrayPoint,ArrayWeight>::getAccuracy(degree);
} // end getAccuracy

} // end namespace Intrepid

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

