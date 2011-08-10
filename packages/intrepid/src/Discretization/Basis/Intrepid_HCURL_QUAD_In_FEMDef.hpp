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
//                    Denis Ridzal (dridzal@sandia.gov)
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HCURL_QUAD_In_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(div) functions on HEX cells.
    \author Created by R. Kirby, P. Bochev, D. Ridzal and K. Peterson.
*/

namespace Intrepid {

  template<class Scalar, class ArrayScalar>
  Basis_HCURL_QUAD_In_FEM<Scalar,ArrayScalar>::Basis_HCURL_QUAD_In_FEM( int order ,
								      const ArrayScalar & ptsClosed ,
								      const ArrayScalar & ptsOpen):
    closedBasis_( order , ptsClosed ),
    openBasis_( order-1 , ptsOpen ) ,
    closedPts_( ptsClosed ),
    openPts_( ptsOpen )
  {
    this -> basisDegree_       = order;
    this -> basisCardinality_  = 2 * closedBasis_.getCardinality() * openBasis_.getCardinality(); 
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
    this -> basisType_         = BASIS_FEM_FIAT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;
  }

  template<class Scalar, class ArrayScalar>
  Basis_HCURL_QUAD_In_FEM<Scalar,ArrayScalar>::Basis_HCURL_QUAD_In_FEM( int order , const EPointType &pointType ):
    closedBasis_( order , pointType==POINTTYPE_SPECTRAL?POINTTYPE_SPECTRAL:POINTTYPE_EQUISPACED ),
    openBasis_( order-1 , pointType==POINTTYPE_SPECTRAL?POINTTYPE_SPECTRAL_OPEN:POINTTYPE_EQUISPACED ) ,
    closedPts_( order+1 , 1 ),
    openPts_( order , 1 )
  {
    this -> basisDegree_       = order;
    this -> basisCardinality_  = 2 * closedBasis_.getCardinality() * openBasis_.getCardinality(); 
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
    this -> basisType_         = BASIS_FEM_FIAT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;

    PointTools::getLattice<Scalar,FieldContainer<Scalar> >( closedPts_ ,
                                                            shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >()) ,
                                                            order ,
                                                            0 ,
                                                            pointType==POINTTYPE_SPECTRAL?POINTTYPE_WARPBLEND:POINTTYPE_EQUISPACED );

    if (pointType == POINTTYPE_SPECTRAL)
      {
	PointTools::getGaussPoints<Scalar,FieldContainer<Scalar> >( openPts_ ,
								    order - 1 );
      }
    else
      {
	PointTools::getLattice<Scalar,FieldContainer<Scalar> >( openPts_ ,
								shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >()) ,
								order - 1,
								0 ,
								POINTTYPE_EQUISPACED );
      }

  }
  
  template<class Scalar, class ArrayScalar>
  void Basis_HCURL_QUAD_In_FEM<Scalar, ArrayScalar>::initializeTags() {
    
    // Basis-dependent intializations
    int tagSize  = 4;        // size of DoF tag
    int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
    int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
    int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
    
    std::vector<int> tags( tagSize * this->getCardinality() );
    
    const std::vector<std::vector<int> >& closedDofTags = closedBasis_.getAllDofTags();
    const std::vector<std::vector<int> >& openDofTags = openBasis_.getAllDofTags();

    std::map<int,std::map<int,int> > total_dof_per_entity;
    std::map<int,std::map<int,int> > current_dof_per_entity;

    for (int i=0;i<4;i++) {
      total_dof_per_entity[0][i] = 0;
      current_dof_per_entity[0][i] = 0;
    }
    for (int i=0;i<4;i++) {
      total_dof_per_entity[1][i] = 0;
      current_dof_per_entity[1][i] = 0;
    }
    total_dof_per_entity[2][0] = 0;
    current_dof_per_entity[2][0] = 0;

    // tally dof on each facet.  none on vertex
    for (int i=0;i<4;i++) {
      total_dof_per_entity[1][i] = openBasis_.getCardinality();
    }

    total_dof_per_entity[2][0] = this->getCardinality() - 4 * openBasis_.getCardinality();

    int tagcur = 0;
    // loop over the x-face tangent basis functions, which are (0,phi(x)psi(y))
    // for psi in the closed basis and phi in the open
    for (int j=0;j<openBasis_.getCardinality();j++) {
      const int odim = openDofTags[j][0];
      const int oent = openDofTags[j][1];
      for (int i=0;i<closedBasis_.getCardinality();i++) {
	const int cdim = closedDofTags[i][0];
	const int cent = closedDofTags[i][1];
	int dofdim;
	int dofent;
	ProductTopology::lineProduct2d(cdim,cent,odim,oent,dofdim,dofent);
	tags[4*tagcur] = dofdim;
	tags[4*tagcur+1] = dofent;
	tags[4*tagcur+2] = current_dof_per_entity[dofdim][dofent];
	current_dof_per_entity[dofdim][dofent]++;
	tags[4*tagcur+3] = total_dof_per_entity[dofdim][dofent];
	tagcur++;
      }
    }
    // now we have to do it for the other basis functions, which are
    // (psi(x)phi(y)) for psi in the closed basis and phi in the open
    for (int j=0;j<closedBasis_.getCardinality();j++) {
      const int cdim = closedDofTags[j][0];
      const int cent = closedDofTags[j][1];
      for (int i=0;i<openBasis_.getCardinality();i++) {
	const int odim = openDofTags[i][0];
	const int oent = openDofTags[i][1];
	int dofdim;
	int dofent;
	ProductTopology::lineProduct2d(odim,oent,cdim,cent,dofdim,dofent);
	tags[4*tagcur] = dofdim;
	tags[4*tagcur+1] = dofent;
	tags[4*tagcur+2] = current_dof_per_entity[dofdim][dofent];
	current_dof_per_entity[dofdim][dofent]++;
	tags[4*tagcur+3] = total_dof_per_entity[dofdim][dofent];
	tagcur++;
      }
    }

    // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
    Intrepid::setOrdinalTagData(this -> tagToOrdinal_,
				this -> ordinalToTag_,
				&(tags[0]),
				this -> basisCardinality_,
				tagSize,
				posScDim,
				posScOrd,
				posDfOrd);
  }


  template<class Scalar, class ArrayScalar>
  void Basis_HCURL_QUAD_In_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
							      const ArrayScalar &  inputPoints,
							      const EOperator      operatorType) const {
    
    // Verify arguments
#ifdef HAVE_INTREPID_DEBUG
    Intrepid::getValues_HCURL_Args<Scalar, ArrayScalar>(outputValues,
						       inputPoints,
						       operatorType,
						       this -> getBaseCellTopology(),
						       this -> getCardinality() );
#endif
    
    // Number of evaluation points = dim 0 of inputPoints
    int dim0 = inputPoints.dimension(0);
    
    // separate out points
    FieldContainer<Scalar> xPoints(dim0,1);
    FieldContainer<Scalar> yPoints(dim0,1);
    
    for (int i=0;i<dim0;i++) {
      xPoints(i,0) = inputPoints(i,0);
      yPoints(i,0) = inputPoints(i,1);
    }
    
    switch (operatorType) {
    case OPERATOR_VALUE:
      {
	FieldContainer<Scalar> closedBasisValsXPts( closedBasis_.getCardinality() , dim0 );
	FieldContainer<Scalar> closedBasisValsYPts( closedBasis_.getCardinality() , dim0 );
	FieldContainer<Scalar> openBasisValsXPts( openBasis_.getCardinality() , dim0 );
	FieldContainer<Scalar> openBasisValsYPts( openBasis_.getCardinality() , dim0 );
	
	closedBasis_.getValues( closedBasisValsXPts , xPoints , OPERATOR_VALUE );
	closedBasis_.getValues( closedBasisValsYPts , yPoints , OPERATOR_VALUE );
	openBasis_.getValues( openBasisValsXPts , xPoints , OPERATOR_VALUE );
	openBasis_.getValues( openBasisValsYPts , yPoints , OPERATOR_VALUE );
	
	int bfcur = 0;
	// x component bfs are (closed(x) open(y),0)
	for (int j=0;j<openBasis_.getCardinality();j++) {
	  for (int i=0;i<closedBasis_.getCardinality();i++) {
	    for (int l=0;l<dim0;l++) {
	      outputValues(bfcur,l,0) = 0.0;
	      outputValues(bfcur,l,1) = closedBasisValsXPts(i,l) * openBasisValsYPts(j,l);
	    }
	    bfcur++;
	  }
	}
	
	// y component bfs are (0,open(x) closed(y))
	for (int j=0;j<closedBasis_.getCardinality();j++) {
	  for (int i=0;i<openBasis_.getCardinality();i++) {
	    for (int l=0;l<dim0;l++) {
	      outputValues(bfcur,l,0) = openBasisValsXPts(i,l) * closedBasisValsYPts(j,l);
	      outputValues(bfcur,l,1) = 0.0;
	    }
	    bfcur++;
	  }
	}
      }
      break;
    case OPERATOR_CURL:
      {
	FieldContainer<Scalar> closedBasisDerivsXPts( closedBasis_.getCardinality() , dim0 , 1 );
	FieldContainer<Scalar> closedBasisDerivsYPts( closedBasis_.getCardinality() , dim0 , 1 );
	FieldContainer<Scalar> openBasisValsXPts( openBasis_.getCardinality() , dim0 );
	FieldContainer<Scalar> openBasisValsYPts( openBasis_.getCardinality() , dim0 );
	
	closedBasis_.getValues( closedBasisDerivsXPts , xPoints , OPERATOR_D1 );
	closedBasis_.getValues( closedBasisDerivsYPts , yPoints , OPERATOR_D1 );
	openBasis_.getValues( openBasisValsXPts , xPoints , OPERATOR_VALUE );
	openBasis_.getValues( openBasisValsYPts , yPoints , OPERATOR_VALUE );
	
	int bfcur = 0;
	
	// x component basis functions first
	for (int j=0;j<openBasis_.getCardinality();j++) {
	  for (int i=0;i<closedBasis_.getCardinality();i++) {
	    for (int l=0;l<dim0;l++) {
	      outputValues(bfcur,l) = closedBasisDerivsXPts(i,l,0) * openBasisValsYPts(j,l);
	    }
	    bfcur++;
	  }
	}
	
	// now y component basis functions
	for (int j=0;j<closedBasis_.getCardinality();j++) {
	  for (int i=0;i<openBasis_.getCardinality();i++) {
	    for (int l=0;l<dim0;l++) {
	      outputValues(bfcur,l) = openBasisValsXPts(i,l) * closedBasisDerivsYPts(j,l,0);
	    }
	    bfcur++;
	  }
	}
      }
      break;
    case OPERATOR_DIV:
      TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
			  ">>> ERROR (Basis_HCURL_QUAD_In_FEM): DIV is invalid operator for HCURL Basis Functions");
      break;
      
    case OPERATOR_GRAD:
      TEST_FOR_EXCEPTION( (operatorType == OPERATOR_GRAD), std::invalid_argument,
			  ">>> ERROR (Basis_HCURL_QUAD_In_FEM): GRAD is invalid operator for HCURL Basis Functions");
      break;
      
    case OPERATOR_D1:
    case OPERATOR_D2:
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
      TEST_FOR_EXCEPTION( ( (operatorType == OPERATOR_D1)    ||
			    (operatorType == OPERATOR_D2)    ||
			    (operatorType == OPERATOR_D3)    ||
			    (operatorType == OPERATOR_D4)    ||
			    (operatorType == OPERATOR_D5)    ||
			    (operatorType == OPERATOR_D6)    ||
			    (operatorType == OPERATOR_D7)    ||
			    (operatorType == OPERATOR_D8)    ||
			    (operatorType == OPERATOR_D9)    ||
			    (operatorType == OPERATOR_D10) ),
			  std::invalid_argument,
			  ">>> ERROR (Basis_HCURL_QUAD_In_FEM): Invalid operator type");
      break;
      
    default:
      TEST_FOR_EXCEPTION( ( (operatorType != OPERATOR_VALUE) &&
			    (operatorType != OPERATOR_GRAD)  &&
			    (operatorType != OPERATOR_CURL)  &&
			    (operatorType != OPERATOR_CURL)   &&
			    (operatorType != OPERATOR_D1)    &&
			    (operatorType != OPERATOR_D2)    &&
			    (operatorType != OPERATOR_D3)    &&
			    (operatorType != OPERATOR_D4)    &&
			    (operatorType != OPERATOR_D5)    &&
			    (operatorType != OPERATOR_D6)    &&
			    (operatorType != OPERATOR_D7)    &&
			    (operatorType != OPERATOR_D8)    &&
			    (operatorType != OPERATOR_D9)    &&
			    (operatorType != OPERATOR_D10) ),
			  std::invalid_argument,
			  ">>> ERROR (Basis_HCURL_QUAD_In_FEM): Invalid operator type");
    }
  }
  
  
  
  template<class Scalar, class ArrayScalar>
  void Basis_HCURL_QUAD_In_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
							      const ArrayScalar &    inputPoints,
							      const ArrayScalar &    cellVertices,
							      const EOperator        operatorType) const {
    TEST_FOR_EXCEPTION( (true), std::logic_error,
			">>> ERROR (Basis_HCURL_QUAD_In_FEM): FEM Basis calling an FVD member function");
  }

  template<class Scalar, class ArrayScalar>
  void Basis_HCURL_QUAD_In_FEM<Scalar,ArrayScalar>::getDofCoords(ArrayScalar & DofCoords) const
  {
    // x-component basis functions
    int cur = 0;

    for (int j=0;j<closedPts_.dimension(0);j++)
      {
        for (int i=0;i<openPts_.dimension(0);i++)
          {
            DofCoords(cur,0) = openPts_(i,0);
            DofCoords(cur,1) = closedPts_(j,0);
            cur++;
          }
      }

    // y-component basis functions
    for (int j=0;j<openPts_.dimension(0);j++)
      {
        for (int i=0;i<closedPts_.dimension(0);i++)
          {
            DofCoords(cur,0) = closedPts_(i,0);
            DofCoords(cur,1) = openPts_(j,0);
            cur++;
          }
      }

    return;
  }


  
}// namespace Intrepid
