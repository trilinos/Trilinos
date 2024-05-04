#ifndef INTREPID_HGRAD_TET_COMP12_FEMDEF_HPP
#define INTREPID_HGRAD_TET_COMP12_FEMDEF_HPP
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

/** \file   Intrepid_HGRAD_TET_COMP12_FEMDef.hpp
    \brief  Definition file for the composite H(grad)-compatible FEM basis of degree 1
            on Tetrahedron cell with 12 sub-tetrahedrons.
    \author Created by P. Bochev, J. Ostien, K. Peterson and D. Ridzal.
*/

namespace Intrepid {

  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_TET_COMP12_FEM<Scalar, ArrayScalar>::Basis_HGRAD_TET_COMP12_FEM()
  {
    this -> basisCardinality_  = 10;
    this -> basisDegree_       = 1;
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<11> >() );
    this -> basisType_         = BASIS_FEM_DEFAULT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;
  }


template<class Scalar, class ArrayScalar>
void Basis_HGRAD_TET_COMP12_FEM<Scalar, ArrayScalar>::initializeTags() {

  // Basis-dependent intializations
  int tagSize  = 4;        // size of DoF tag
  int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
  int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
  int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

  // An array with local DoF tags assigned to basis functions, in the order of their local enumeration
  int tags[]  = { 0, 0, 0, 1,
                  0, 1, 0, 1,
                  0, 2, 0, 1,
                  0, 3, 0, 1,
                  1, 0, 0, 1,
                  1, 1, 0, 1,
                  1, 2, 0, 1,
                  1, 3, 0, 1,
                  1, 4, 0, 1,
                  1, 5, 0, 1,
  };

  // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
  Intrepid::setOrdinalTagData(this -> tagToOrdinal_,
                              this -> ordinalToTag_,
                              tags,
                              this -> basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);
}



template<class Scalar, class ArrayScalar>
void Basis_HGRAD_TET_COMP12_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
                                                                const ArrayScalar &  inputPoints,
                                                                const EOperator      operatorType) const {

  // Verify arguments
#ifdef HAVE_INTREPID_DEBUG
  Intrepid::getValues_HGRAD_Args<Scalar, ArrayScalar>(outputValues,
                                                      inputPoints,
                                                      operatorType,
                                                      this -> getBaseCellTopology(),
                                                      this -> getCardinality() );
#endif

  // Number of evaluation points = dim 0 of inputPoints
  int dim0 = inputPoints.dimension(0);

  // Temporaries: (x,y,z) coordinates of the evaluation point
  Scalar r = 0.0;
  Scalar s = 0.0;
  Scalar t = 0.0;

  // Temporary for the auriliary node basis function
  Scalar aux = 0.0;

  // Array to store all the subtets containing the given pt
  Teuchos::Array<int> pt_tets;


  switch (operatorType) {

    case OPERATOR_VALUE:
      for (int i0 = 0; i0 < dim0; i0++) {
        r = inputPoints(i0, 0);
        s = inputPoints(i0, 1);
        t = inputPoints(i0, 2);

	pt_tets = getLocalSubTetrahedra(r,s,t);

	// idependent verification shows that a given point will produce
	// the same shape functions for each tet that contains it, so
	// we only need to use the first one returned.
	//for (int pt = 0; pt < pt_tets.size(); ++pt)
        if (pt_tets[0] != -1) {
          int subtet = pt_tets[0];
          aux = 0.0;
          // outputValues is a rank-2 array with dimensions (basisCardinality_, dim0)
          switch (subtet) {
          case 0:
            outputValues(0, i0) = 1. - 2. * (r + s + t);
            outputValues(4, i0) = 2. * r;
            outputValues(6, i0) = 2. * s;
            outputValues(7, i0) = 2. * t;
            break;
          case 1:
            outputValues(1, i0) = 2. * r - 1.;
            outputValues(4, i0) = 2. - 2. * (r + s + t);
            outputValues(5, i0) = 2. * s;
            outputValues(8, i0) = 2. * t;
            break;
          case 2:
            outputValues(2, i0) = 2. * s - 1.;
            outputValues(5, i0) = 2. * r;
            outputValues(6, i0) = 2. - 2. * (r + s + t);
            outputValues(9, i0) = 2. * t;
            break;
          case 3:
            outputValues(3, i0) = 2. * t - 1.;
            outputValues(7, i0) = 2. - 2. * (r + s + t);
            outputValues(8, i0) = 2. * r;
            outputValues(9, i0) = 2. * s;
            break;
          case 4:
            outputValues(4, i0) = 1. - 2. * (s + t);
            outputValues(5, i0) = 2. * (r + s) - 1.;
            outputValues(8, i0) = 2. * (r + t) - 1.;
            aux = 2. - 4. * r;
            break;
          case 5:
            outputValues(5, i0) = 2. * (r + s) - 1.;
            outputValues(8, i0) = 2. * (r + t) - 1.;
            outputValues(9, i0) = 2. * (s + t) - 1.;
            aux = 4. - 4. * (r + s + t);
            break;
          case 6:
            outputValues(7, i0) = 1. - 2. * (r + s);
            outputValues(8, i0) = 2. * (r + t) - 1.;
            outputValues(9, i0) = 2. * (s + t) - 1.;
            aux = 2. - 4. * t;
            break;
          case 7:
            outputValues(4, i0) = 1. - 2. * (s + t);
            outputValues(7, i0) = 1. - 2. * (r + s);
            outputValues(8, i0) = 2. * (r + t) - 1.;
            aux = 4. * s;
            break;
          case 8:
            outputValues(4, i0) = 1. - 2. * (s + t);
            outputValues(5, i0) = 2. * (r + s) - 1.;
            outputValues(6, i0) = 1. - 2. * (r + t);
            aux = 4. * t;
            break;
          case 9:
            outputValues(5, i0) = 2. * (r + s) - 1.;
            outputValues(6, i0) = 1. - 2. * (r + t);
            outputValues(9, i0) = 2. * (s + t) - 1.;
            aux = 2. - 4. * s;
            break;
          case 10:
            outputValues(6, i0) = 1. - 2. * (r + t);
            outputValues(7, i0) = 1. - 2. * (r + s);
            outputValues(9, i0) = 2. * (s + t) - 1.;
            aux = 4. * r;
            break;
          case 11:
            outputValues(4, i0) = 1. - 2. * (s + t);
            outputValues(6, i0) = 1. - 2. * (r + t);
            outputValues(7, i0) = 1. - 2. * (r + s);
            aux = 4. * (r + s + t) - 2.;
            break;
          }
          outputValues(4, i0) += aux/6.0;
          outputValues(5, i0) += aux/6.0;
          outputValues(6, i0) += aux/6.0;
          outputValues(7, i0) += aux/6.0;
          outputValues(8, i0) += aux/6.0;
          outputValues(9, i0) += aux/6.0;
        }
      }
      break;

  case OPERATOR_GRAD:
  case OPERATOR_D1:
    {
      // initialize to 0.0 since we will be accumulating
      outputValues.initialize(0.0);

      FieldContainer<Scalar> Lopt(10,3);
      for (int pt=0; pt < dim0; ++pt) {

        r = inputPoints(pt, 0);
        s = inputPoints(pt, 1);
        t = inputPoints(pt, 2);

        Lopt(0,0) = (-17 + 20*r + 20*s + 20*t)/8.;
        Lopt(0,1) = (-17 + 20*r + 20*s + 20*t)/8.;
        Lopt(0,2) = (-17 + 20*r + 20*s + 20*t)/8.;
        Lopt(1,0) = -0.375 + (5*r)/2.;
        Lopt(1,1) = 0.;
        Lopt(1,2) = 0.;
        Lopt(2,0) = 0.;
        Lopt(2,1) = -0.375 + (5*s)/2.;
        Lopt(2,2) = 0.;
        Lopt(3,0) = 0.;
        Lopt(3,1) = 0.;
        Lopt(3,2) = -0.375 + (5*t)/2.;
        Lopt(4,0) = (-35*(-1 + 2*r + s + t))/12.;
        Lopt(4,1) = (-4 - 35*r + 5*s + 10*t)/12.;
        Lopt(4,2) = (-4 - 35*r + 10*s + 5*t)/12.;
        Lopt(5,0) = (-1 + 5*r + 40*s - 5*t)/12.;
        Lopt(5,1) = (-1 + 40*r + 5*s - 5*t)/12.;
        Lopt(5,2) = (-5*(-1 + r + s + 2*t))/12.;
        Lopt(6,0) = (-4 + 5*r - 35*s + 10*t)/12.;
        Lopt(6,1) = (-35*(-1 + r + 2*s + t))/12.;
        Lopt(6,2) = (-4 + 10*r - 35*s + 5*t)/12.;
        Lopt(7,0) = (-4 + 5*r + 10*s - 35*t)/12.;
        Lopt(7,1) = (-4 + 10*r + 5*s - 35*t)/12.;
        Lopt(7,2) = (-35*(-1 + r + s + 2*t))/12.;
        Lopt(8,0) = (-1 + 5*r - 5*s + 40*t)/12.;
        Lopt(8,1) = (-5*(-1 + r + 2*s + t))/12.;
        Lopt(8,2) = (-1 + 40*r - 5*s + 5*t)/12.;
        Lopt(9,0) = (-5*(-1 + 2*r + s + t))/12.;
        Lopt(9,1) = (-1 - 5*r + 5*s + 40*t)/12.;
        Lopt(9,2) = (-1 - 5*r + 40*s + 5*t)/12.;

        for (int node=0; node < 10; ++node) {
          for (int dim=0; dim < 3; ++dim) {
            outputValues(node,pt,dim) = Lopt(node,dim);
          }
        }
      }
    }
    break;

    case OPERATOR_CURL:
      TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_CURL), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_TET_COMP12_FEM): CURL is invalid operator for rank-0 (scalar) functions in 3D");
      break;

    case OPERATOR_DIV:
      TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_TET_COMP12_FEM): DIV is invalid operator for rank-0 (scalar) functions in 3D");
      break;

    case OPERATOR_D2:
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
      {
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, DkCardinality)
        int DkCardinality = Intrepid::getDkCardinality(operatorType,
                                                       this -> basisCellTopology_.getDimension() );
        for(int dofOrd = 0; dofOrd < this -> basisCardinality_; dofOrd++) {
          for (int i0 = 0; i0 < dim0; i0++) {
            for(int dkOrd = 0; dkOrd < DkCardinality; dkOrd++){
              outputValues(dofOrd, i0, dkOrd) = 0.0;
            }
          }
        }
      }
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( !( Intrepid::isValidOperator(operatorType) ), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_TET_COMP12_FEM): Invalid operator type");
  }
}



template<class Scalar, class ArrayScalar>
void Basis_HGRAD_TET_COMP12_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                             const ArrayScalar &    inputPoints,
                                                             const ArrayScalar &    cellVertices,
                                                             const EOperator        operatorType) const {
  TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (Basis_HGRAD_TET_COMP12_FEM): FEM Basis calling an FVD member function");

}

template<class Scalar, class ArrayScalar>
void Basis_HGRAD_TET_COMP12_FEM<Scalar, ArrayScalar>::getDofCoords(ArrayScalar & DofCoords) const {
#ifdef HAVE_INTREPID_DEBUG
  // Verify rank of output array.
  TEUCHOS_TEST_FOR_EXCEPTION( !(DofCoords.rank() == 2), std::invalid_argument,
                      ">>> ERROR: (Intrepid::Basis_HGRAD_TET_COMP12_FEM::getDofCoords) rank = 2 required for DofCoords array");
  // Verify 0th dimension of output array.
  TEUCHOS_TEST_FOR_EXCEPTION( !( DofCoords.dimension(0) == this -> basisCardinality_ ), std::invalid_argument,
                      ">>> ERROR: (Intrepid::Basis_HGRAD_TET_COMP12_FEM::getDofCoords) mismatch in number of DoF and 0th dimension of DofCoords array");
  // Verify 1st dimension of output array.
  TEUCHOS_TEST_FOR_EXCEPTION( !( DofCoords.dimension(1) == (int)(this -> basisCellTopology_.getDimension()) ), std::invalid_argument,
                      ">>> ERROR: (Intrepid::Basis_HGRAD_TET_COMP12_FEM::getDofCoords) incorrect reference cell (1st) dimension in DofCoords array");
#endif

  DofCoords(0,0) = 0.0;   DofCoords(0,1) = 0.0; DofCoords(0,2) = 0.0;
  DofCoords(1,0) = 1.0;   DofCoords(1,1) = 0.0; DofCoords(1,2) = 0.0;
  DofCoords(2,0) = 0.0;   DofCoords(2,1) = 1.0; DofCoords(2,2) = 0.0;
  DofCoords(3,0) = 0.0;   DofCoords(3,1) = 0.0; DofCoords(3,2) = 1.0;
  DofCoords(4,0) = 0.5;   DofCoords(4,1) = 0.0; DofCoords(4,2) = 0.0;
  DofCoords(5,0) = 0.5;   DofCoords(5,1) = 0.5; DofCoords(5,2) = 0.0;
  DofCoords(6,0) = 0.0;   DofCoords(6,1) = 0.5; DofCoords(6,2) = 0.0;
  DofCoords(7,0) = 0.0;   DofCoords(7,1) = 0.0; DofCoords(7,2) = 0.5;
  DofCoords(8,0) = 0.5;   DofCoords(8,1) = 0.0; DofCoords(8,2) = 0.5;
  DofCoords(9,0) = 0.0;   DofCoords(9,1) = 0.5; DofCoords(9,2) = 0.5;
}

template<class Scalar, class ArrayScalar>
Teuchos::Array<int>
Basis_HGRAD_TET_COMP12_FEM<Scalar, ArrayScalar>::getLocalSubTetrahedra(Scalar x, Scalar y, Scalar z) const
{

  Teuchos::Array<int> subTets;
  int count(0);

  // local coords
  Scalar xyz = x + y + z;
  Scalar xy = x + y;
  Scalar xz = x + z;
  Scalar yz = y + z;

  // cycle through each subdomain and push back if the point lies within

  // subtet #0 E0 := 0.0 <= r + s + t <= 0.5 && 0.0 <= r <= 0.5 && 0.0 <= s <= 0.5 && 0.0 <= t <= 0.5
  if ( (0.0 <= xyz && xyz <= 0.5) && (0.0 <= x && x <= 0.5) && (0.0 <= y && y <= 0.5) && (0.0 <= z && z <= 0.5) ) {
    count++;
    subTets.push_back(0);
  }

  // subtet #1 E1 := 0.5 <= r + s + t <= 1.0 && 0.5 <= r <= 1.0 && 0.0 <= s <= 0.5 && 0.0 <= t <= 0.5
  if ( (0.5 <= xyz && xyz <= 1.0) && (0.5 <= x && x <= 1.0) && (0.0 <= y && y <= 0.5) && (0.0 <= z && z <= 0.5) ) {
    count++;
    subTets.push_back(1);
  }

  // subtet #2 E2 := 0.5 <= r + s + t <= 1.0 && 0.0 <= r <= 0.5 && 0.5 <= s <= 1.0 && 0.0 <= t <= 0.5
  if ( (0.5 <= xyz && xyz <= 1.0) && (0.0 <= x && x <= 0.5) && (0.5 <= y && y <= 1.0) && (0.0 <= z && z<= 0.5) ) {
    count++;
    subTets.push_back(2);
  }

  // subtet #3 E3 := 0.5 <= r + s + t <= 1.0 && 0.0 <= r <= 0.5 && 0.0 <= s <= 0.5 && 0.5 <= t <= 1.0
  if ( (0.5 <= xyz && xyz <= 1.0) && (0.0 <= x && x <= 0.5) && (0.0 <= y && y <= 0.5) && (0.5 <= z && z <= 1.0) ) {
    count++;
    subTets.push_back(3);
  }

  // subtet #4 E4 := 0.0 <= s + t <= 0.5 && 0.5 <= r + s <= 1.0 && 0.5 <= r + t <= 1.0 && 0.0 <= r <= 0.5
  if ( (0.0 <= yz && yz <= 0.5) && (0.5 <= xy && xy <= 1.0) && (0.5 <= xz && xz <= 1.0) && (0.0 <= x && x <= 0.5) ) {
    count++;
    subTets.push_back(4);
  }

  // subtet #5 E5 := 0.5 <= r + s <= 1.0 && 0.5 <= s + t <= 1.0 && 0.5 <= r + t <= 1.0 && 0.75 <= r + s + t <= 1.0
  if ( (0.5 <= xy && xy <= 1.0) && (0.5 <= yz && yz <= 1.0) && (0.5 <= xz && xz <= 1.0) && (0.75 <= xyz && xyz <= 1.0) ) {
    count++;
    subTets.push_back(5);
  }

  // subtet #6 E6 := 0.5 <= s + t <= 1.0 && 0.0 <= r + s <= 0.5 && 0.5 <= r + t <= 1.0 && 0.0 <= t <= 0.5
  if ( (0.5 <= yz && yz <= 1.0) && (0.0 <= xy && xy <= 0.5) && (0.5 <= xz && xz <= 1.0) && (0.0 <= z && z <= 0.5) ) {
    count++;
    subTets.push_back(6);
  }

  // subtet #7 E7 := 0.0 <= s + t <= 0.5 && 0.0 <= r + s <= 0.5 && 0.5 <= r + t <= 1.0 && 0.0 <= s <= 0.25
  if ( (0.0 <= yz && yz <= 0.5) && (0.0 <= xy && xy <= 0.5) && (0.5 <= xz && xz <= 1.0) && (0.0 <= y && y <= 0.25) ) {
    count++;
    subTets.push_back(7);
  }

  // subtet #8 E8 := 0.0 <= r + t <= 0.5 && 0.0 <= s + t <= 0.5 &&  0.5 <= r + s <= 1.0 && 0.0 <= t <= 0.25
  if ( (0.0 <= xz && xz <= 0.5) && (0.0 <= yz && yz <= 0.5) && (0.5 <= xy && xy <= 1.0) && (0.0 <= z && z <= 0.25) ) {
    count++;
    subTets.push_back(8);
  }

  // subtet #9 E9 := 0.0 <= r + t <= 0.5 && 0.5 <= r + s <= 1.0 &&  0.5 <= s + t <= 1.0 && 0.0 <= s <= 0.5
  if ( (0.0 <= xz && xz <= 0.5) && (0.5 <= xy && xy <= 1.0) && (0.5 <= yz && yz <= 1.0) && (0.0 <= y && y <= 0.5) ) {
    count++;
    subTets.push_back(9);
  }

  // subtet #10 E10 := 0.0 <= r + t <= 0.5 && 0.5 <= s + t <= 1.0 && 0.0 <= r + s <= 0.5 && 0.0 <= r <= 0.25
  if ( (0.0 <= xz && xz <= 0.5) && (0.5 <= yz && yz <= 1.0) && (0.0 <= xy && xy <= 0.5) && (0.0 <= x && x <= 0.25) ) {
    count++;
    subTets.push_back(10);
  }

  // subtet #11 E11 := 0.5 <= r + s + t <= 0.75 && 0.0 <= r + t <= 0.5 && 0.0 <= s + t <= 0.5 && 0.0 <= r + s <= 0.5
  if ( (0.5 <= xyz && xyz <= 0.75) && (0.0 <= xz && xz <= 0.5) && (0.0 <= yz && yz <= 0.5) && (0.0 <= xy && xy <= 0.5) ) {
    count++;
    subTets.push_back(11);
  }

  // if the point doesn't lie in the parent domain return -1
  if (count == 0) {
    subTets.push_back(-1);
  }

  return subTets;
}


template<class Scalar, class ArrayScalar>
Intrepid::FieldContainer<Scalar>
Basis_HGRAD_TET_COMP12_FEM<Scalar, ArrayScalar>::getWeights(const ArrayScalar & inPts) const
{
  int numPoints = inPts.dimension(0);
  Intrepid::FieldContainer<Scalar> w(numPoints, 12);
  w.initialize(0.0);
  Teuchos::Array< Teuchos::Array<int> > pt_tets;

  for (int pt = 0; pt < numPoints; ++pt)
    pt_tets.push_back(getLocalSubTetrahedra(inPts(pt,0), inPts(pt,1), inPts(pt,2)));

  Teuchos::Array<int> flat;
  FieldContainer<int> count(12);

  for (int pt = 0; pt < numPoints; ++pt)
    for (int i = 0; i < pt_tets[pt].size(); ++i)
      flat.push_back(pt_tets[pt][i]);

  for (int i = 0; i < flat.size(); ++i)
    count(flat[i])++;

  for (int pt = 0; pt < numPoints; ++pt)
    for (int i = 0; i < pt_tets[pt].size(); ++i)
      w(pt, pt_tets[pt][i]) = 1.0/count(pt_tets[pt][i]);

  return w;
}



template<class Scalar, class ArrayScalar>
Intrepid::FieldContainer<Scalar>
Basis_HGRAD_TET_COMP12_FEM<Scalar, ArrayScalar>::getBarycentricCoords(const ArrayScalar & inPts) const
{
  int numPoints = inPts.dimension(0);
  Intrepid::FieldContainer<Scalar> lambda(numPoints, 4);

  for (int pt = 0; pt < numPoints; ++pt)
  {
    lambda(pt,0) = 1. - inPts(pt,0) - inPts(pt,1) - inPts(pt,2);
    lambda(pt,1) = inPts(pt,0);
    lambda(pt,2) = inPts(pt,1);
    lambda(pt,3) = inPts(pt,2);
  }

  return lambda;
}

template<class Scalar, class ArrayScalar>
Scalar
Basis_HGRAD_TET_COMP12_FEM<Scalar, ArrayScalar>::det44(const Intrepid::FieldContainer<Scalar> a) const
{
  Scalar det = a(0,3) * a(1,2) * a(2,1) * a(3,0)
    - a(0,2) * a(1,3) * a(2,1) * a(3,0)
    - a(0,3) * a(1,1) * a(2,2) * a(3,0)
    + a(0,1) * a(1,3) * a(2,2) * a(3,0)
    + a(0,2) * a(1,1) * a(2,3) * a(3,0)
    - a(0,1) * a(1,2) * a(2,3) * a(3,0)
    - a(0,3) * a(1,2) * a(2,0) * a(3,1)
    + a(0,2) * a(1,3) * a(2,0) * a(3,1)
    + a(0,3) * a(1,0) * a(2,2) * a(3,1)
    - a(0,0) * a(1,3) * a(2,2) * a(3,1)
    - a(0,2) * a(1,0) * a(2,3) * a(3,1)
    + a(0,0) * a(1,2) * a(2,3) * a(3,1)
    + a(0,3) * a(1,1) * a(2,0) * a(3,2)
    - a(0,1) * a(1,3) * a(2,0) * a(3,2)
    - a(0,3) * a(1,0) * a(2,1) * a(3,2)
    + a(0,0) * a(1,3) * a(2,1) * a(3,2)
    + a(0,1) * a(1,0) * a(2,3) * a(3,2)
    - a(0,0) * a(1,1) * a(2,3) * a(3,2)
    - a(0,2) * a(1,1) * a(2,0) * a(3,3)
    + a(0,1) * a(1,2) * a(2,0) * a(3,3)
    + a(0,2) * a(1,0) * a(2,1) * a(3,3)
    - a(0,0) * a(1,2) * a(2,1) * a(3,3)
    - a(0,1) * a(1,0) * a(2,2) * a(3,3)
    + a(0,0) * a(1,1) * a(2,2) * a(3,3);

  return det;
}

template<class Scalar, class ArrayScalar>
Intrepid::FieldContainer<Scalar>
Basis_HGRAD_TET_COMP12_FEM<Scalar, ArrayScalar>::inverse44(const Intrepid::FieldContainer<Scalar> a) const
{
  Intrepid::FieldContainer<Scalar> ai(4,4);
  Scalar xj = det44(a);

  ai(0,0) = (1/xj) * (-a(1,3) * a(2,2) * a(3,1) + a(1,2) * a(2,3) * a(3,1) + a(1,3) * a(2,1) * a(3,2) - a(1,1) * a(2,3) * a(3,2) - a(1,2) * a(2,1) * a(3,3) + a(1,1) * a(2,2) * a(3,3));
  ai(0,1) = (1/xj) * ( a(0,3) * a(2,2) * a(3,1) - a(0,2) * a(2,3) * a(3,1) - a(0,3) * a(2,1) * a(3,2) + a(0,1) * a(2,3) * a(3,2) + a(0,2) * a(2,1) * a(3,3) - a(0,1) * a(2,2) * a(3,3));
  ai(0,2) = (1/xj) * (-a(0,3) * a(1,2) * a(3,1) + a(0,2) * a(1,3) * a(3,1) + a(0,3) * a(1,1) * a(3,2) - a(0,1) * a(1,3) * a(3,2) - a(0,2) * a(1,1) * a(3,3) + a(0,1) * a(1,2) * a(3,3));
  ai(0,3) = (1/xj) * ( a(0,3) * a(1,2) * a(2,1) - a(0,2) * a(1,3) * a(2,1) - a(0,3) * a(1,1) * a(2,2) + a(0,1) * a(1,3) * a(2,2) + a(0,2) * a(1,1) * a(2,3) - a(0,1) * a(1,2) * a(2,3));

  ai(1,0) = (1/xj) * ( a(1,3) * a(2,2) * a(3,0) - a(1,2) * a(2,3) * a(3,0) - a(1,3) * a(2,0) * a(3,2) + a(1,0) * a(2,3) * a(3,2) + a(1,2) * a(2,0) * a(3,3) - a(1,0) * a(2,2) * a(3,3));
  ai(1,1) = (1/xj) * (-a(0,3) * a(2,2) * a(3,0) + a(0,2) * a(2,3) * a(3,0) + a(0,3) * a(2,0) * a(3,2) - a(0,0) * a(2,3) * a(3,2) - a(0,2) * a(2,0) * a(3,3) + a(0,0) * a(2,2) * a(3,3));
  ai(1,2) = (1/xj) * ( a(0,3) * a(1,2) * a(3,0) - a(0,2) * a(1,3) * a(3,0) - a(0,3) * a(1,0) * a(3,2) + a(0,0) * a(1,3) * a(3,2) + a(0,2) * a(1,0) * a(3,3) - a(0,0) * a(1,2) * a(3,3));
  ai(1,3) = (1/xj) * (-a(0,3) * a(1,2) * a(2,0) + a(0,2) * a(1,3) * a(2,0) + a(0,3) * a(1,0) * a(2,2) - a(0,0) * a(1,3) * a(2,2) - a(0,2) * a(1,0) * a(2,3) + a(0,0) * a(1,2) * a(2,3));

  ai(2,0) = (1/xj) * (-a(1,3) * a(2,1) * a(3,0) + a(1,1) * a(2,3) * a(3,0) + a(1,3) * a(2,0) * a(3,1) - a(1,0) * a(2,3) * a(3,1) - a(1,1) * a(2,0) * a(3,3) + a(1,0) * a(2,1) * a(3,3));
  ai(2,1) = (1/xj) * ( a(0,3) * a(2,1) * a(3,0) - a(0,1) * a(2,3) * a(3,0) - a(0,3) * a(2,0) * a(3,1) + a(0,0) * a(2,3) * a(3,1) + a(0,1) * a(2,0) * a(3,3) - a(0,0) * a(2,1) * a(3,3));
  ai(2,2) = (1/xj) * (-a(0,3) * a(1,1) * a(3,0) + a(0,1) * a(1,3) * a(3,0) + a(0,3) * a(1,0) * a(3,1) - a(0,0) * a(1,3) * a(3,1) - a(0,1) * a(1,0) * a(3,3) + a(0,0) * a(1,1) * a(3,3));
  ai(2,3) = (1/xj) * ( a(0,3) * a(1,1) * a(2,0) - a(0,1) * a(1,3) * a(2,0) - a(0,3) * a(1,0) * a(2,1) + a(0,0) * a(1,3) * a(2,1) + a(0,1) * a(1,0) * a(2,3) - a(0,0) * a(1,1) * a(2,3));

  ai(3,0) = (1/xj) * ( a(1,2) * a(2,1) * a(3,0) - a(1,1) * a(2,2) * a(3,0) - a(1,2) * a(2,0) * a(3,1) + a(1,0) * a(2,2) * a(3,1) + a(1,1) * a(2,0) * a(3,2) - a(1,0) * a(2,1) * a(3,2));
  ai(3,1) = (1/xj) * (-a(0,2) * a(2,1) * a(3,0) + a(0,1) * a(2,2) * a(3,0) + a(0,2) * a(2,0) * a(3,1) - a(0,0) * a(2,2) * a(3,1) - a(0,1) * a(2,0) * a(3,2) + a(0,0) * a(2,1) * a(3,2));
  ai(3,2) = (1/xj) * ( a(0,2) * a(1,1) * a(3,0) - a(0,1) * a(1,2) * a(3,0) - a(0,2) * a(1,0) * a(3,1) + a(0,0) * a(1,2) * a(3,1) + a(0,1) * a(1,0) * a(3,2) - a(0,0) * a(1,1) * a(3,2));
  ai(3,3) = (1/xj) * (-a(0,2) * a(1,1) * a(2,0) + a(0,1) * a(1,2) * a(2,0) + a(0,2) * a(1,0) * a(2,1) - a(0,0) * a(1,2) * a(2,1) - a(0,1) * a(1,0) * a(2,2) + a(0,0) * a(1,1) * a(2,2));

  return ai;
}

template<class Scalar, class ArrayScalar>
Intrepid::FieldContainer<Scalar>
Basis_HGRAD_TET_COMP12_FEM<Scalar, ArrayScalar>::getSubTetGrads() const
{
  Intrepid::FieldContainer<Scalar> dx(3,11,12);
  dx.initialize(0.0);
  // fill in dx
  dx(0,0,0)   = -2.0;
  dx(0,1,1)   =  2.0;
  dx(0,4,0)   =  2.0;
  dx(0,4,1)   = -2.0;
  dx(0,5,2)   =  2.0;
  dx(0,5,4)   =  2.0;
  dx(0,5,5)   =  2.0;
  dx(0,5,8)   =  2.0;
  dx(0,5,9)   =  2.0;
  dx(0,6,2)   = -2.0;
  dx(0,6,8)   = -2.0;
  dx(0,6,9)   = -2.0;
  dx(0,6,10)  = -2.0;
  dx(0,6,11)  = -2.0;
  dx(0,7,3)   = -2.0;
  dx(0,7,6)   = -2.0;
  dx(0,7,7)   = -2.0;
  dx(0,7,10)  = -2.0;
  dx(0,7,11)  = -2.0;
  dx(0,8,3)   =  2.0;
  dx(0,8,4)   =  2.0;
  dx(0,8,5)   =  2.0;
  dx(0,8,6)   =  2.0;
  dx(0,8,7)   =  2.0;
  dx(0,10,4)  = -4.0;
  dx(0,10,5)  = -4.0;
  dx(0,10,10) =  4.0;
  dx(0,10,11) =  4.0;

  // fill in dy
  dx(1,0,0)   = -2.0;
  dx(1,2,2)   =  2.0;
  dx(1,4,1)   = -2.0;
  dx(1,4,4)   = -2.0;
  dx(1,4,7)   = -2.0;
  dx(1,4,8)   = -2.0;
  dx(1,4,11)  = -2.0;
  dx(1,5,1)   =  2.0;
  dx(1,5,4)   =  2.0;
  dx(1,5,5)   =  2.0;
  dx(1,5,8)   =  2.0;
  dx(1,5,9)   =  2.0;
  dx(1,6,0)   =  2.0;
  dx(1,6,2)   = -2.0;
  dx(1,7,3)   = -2.0;
  dx(1,7,6)   = -2.0;
  dx(1,7,7)   = -2.0;
  dx(1,7,10)  = -2.0;
  dx(1,7,11)  = -2.0;
  dx(1,9,3)   =  2.0;
  dx(1,9,5)   =  2.0;
  dx(1,9,6)   =  2.0;
  dx(1,9,9)   =  2.0;
  dx(1,9,10)  =  2.0;
  dx(1,10,5)  = -4.0;
  dx(1,10,7)  =  4.0;
  dx(1,10,9)  = -4.0;
  dx(1,10,11) =  4.0;

  // fill in dz
  dx(2,0,0)   = -2.0;
  dx(2,3,3)   =  2.0;
  dx(2,4,1)   = -2.0;
  dx(2,4,4)   = -2.0;
  dx(2,4,7)   = -2.0;
  dx(2,4,8)   = -2.0;
  dx(2,4,11)  = -2.0;
  dx(2,6,2)   = -2.0;
  dx(2,6,8)   = -2.0;
  dx(2,6,9)   = -2.0;
  dx(2,6,10)  = -2.0;
  dx(2,6,11)  = -2.0;
  dx(2,7,0)   =  2.0;
  dx(2,7,3)   = -2.0;
  dx(2,8,1)   =  2.0;
  dx(2,8,4)   =  2.0;
  dx(2,8,5)   =  2.0;
  dx(2,8,6)   =  2.0;
  dx(2,8,7)   =  2.0;
  dx(2,9,2)   =  2.0;
  dx(2,9,5)   =  2.0;
  dx(2,9,6)   =  2.0;
  dx(2,9,9)   =  2.0;
  dx(2,9,10)  =  2.0;
  dx(2,10,5)  = -4.0;
  dx(2,10,6)  = -4.0;
  dx(2,10,8)  =  4.0;
  dx(2,10,11) =  4.0;

  return dx;
}

template<class Scalar, class ArrayScalar>
Intrepid::FieldContainer<Scalar>
Basis_HGRAD_TET_COMP12_FEM<Scalar, ArrayScalar>::getSubTetDetF() const
{
  Intrepid::FieldContainer<Scalar> xJ(12);
  // set sub-elem jacobians
  xJ(0) = 1./48.; xJ(1) = 1./48.; xJ(2) = 1./48.; xJ(3) = 1./48.;
  xJ(4) = 1./96.; xJ(5) = 1./96.; xJ(6) = 1./96.; xJ(7) = 1./96.;
  xJ(8) = 1./96.; xJ(9) = 1./96.; xJ(10) = 1./96.; xJ(11) = 1./96.;

  return xJ;
}


}// namespace Intrepid
#endif

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

