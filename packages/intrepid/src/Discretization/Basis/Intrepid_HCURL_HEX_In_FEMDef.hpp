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

/** \file   Intrepid_HCURL_HEX_In_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(div) functions on HEX cells.
    \author Created by R. Kirby, P. Bochev, D. Ridzal and K. Peterson.
*/

namespace Intrepid {

  template<class Scalar, class ArrayScalar>
  Basis_HCURL_HEX_In_FEM<Scalar,ArrayScalar>::Basis_HCURL_HEX_In_FEM( int order ,
								      const ArrayScalar & ptsClosed ,
								      const ArrayScalar & ptsOpen):
    closedBasis_( order , ptsClosed ),
    openBasis_( order-1 , ptsOpen ) ,
    closedPts_( ptsClosed ),
    openPts_( ptsOpen )    
  {
    this -> basisDegree_       = order;
    this -> basisCardinality_  = 3 * closedBasis_.getCardinality() 
      * closedBasis_.getCardinality() * openBasis_.getCardinality(); 
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
    this -> basisType_         = BASIS_FEM_FIAT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;

    Array<Array<RCP<Basis<Scalar,ArrayScalar > > > > bases(3);
    bases[0].resize(3); bases[1].resize(3); bases[2].resize(3);
    bases[0][0] = rcp( &openBasis_ , false );
    bases[0][1] = rcp( &closedBasis_ , false );
    bases[0][2] = rcp( &closedBasis_ , false );
    bases[1][0] = rcp( &closedBasis_ , false );
    bases[1][1] = rcp( &openBasis_ , false );
    bases[1][2] = rcp( &closedBasis_ , false );
    bases[2][0] = rcp( &closedBasis_ , false );
    bases[2][1] = rcp( &closedBasis_ , false );
    bases[2][2] = rcp( &openBasis_ , false );
    this->setBases( bases );

  }

  template<class Scalar, class ArrayScalar>
  Basis_HCURL_HEX_In_FEM<Scalar,ArrayScalar>::Basis_HCURL_HEX_In_FEM( int order , const EPointType &pointType ):
    closedBasis_( order , pointType==POINTTYPE_SPECTRAL?POINTTYPE_SPECTRAL:POINTTYPE_EQUISPACED ),
    openBasis_( order-1 , pointType==POINTTYPE_SPECTRAL?POINTTYPE_SPECTRAL_OPEN:POINTTYPE_EQUISPACED ),
    closedPts_( order+1 , 1 ),
    openPts_( order , 1 )
  {
    this -> basisDegree_       = order;
    this -> basisCardinality_  = 3 * closedBasis_.getCardinality() 
      * closedBasis_.getCardinality() * openBasis_.getCardinality(); 
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
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

    Array<Array<RCP<Basis<Scalar,ArrayScalar > > > > bases(3);
    bases[0].resize(3); bases[1].resize(3); bases[2].resize(3);
    bases[0][0] = rcp( &openBasis_ , false );
    bases[0][1] = rcp( &closedBasis_ , false );
    bases[0][2] = rcp( &closedBasis_ , false );
    bases[1][0] = rcp( &closedBasis_ , false );
    bases[1][1] = rcp( &openBasis_ , false );
    bases[1][2] = rcp( &closedBasis_ , false );
    bases[2][0] = rcp( &closedBasis_ , false );
    bases[2][1] = rcp( &closedBasis_ , false );
    bases[2][2] = rcp( &openBasis_ , false );
    this->setBases( bases );

  }
  
  template<class Scalar, class ArrayScalar>
  void Basis_HCURL_HEX_In_FEM<Scalar, ArrayScalar>::initializeTags() {
    
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

    const int order = this->basisDegree_;

    // vertices
    for (int i=0;i<4;i++) {
      total_dof_per_entity[0][i] = 0;
      current_dof_per_entity[0][i] = 0;
    }
    // edges
    for (int i=0;i<12;i++) {
      total_dof_per_entity[1][i] = 0;
      current_dof_per_entity[1][i] = 0;
    }
    // faces
    for (int i=0;i<6;i++) {
      total_dof_per_entity[2][i] = 0;
      current_dof_per_entity[2][i] = 0;
    }
    total_dof_per_entity[3][0] = 0;
    current_dof_per_entity[3][0] = 0;

    // tally dof on each facet.  none on vertex
    // edge dof
    for (int i=0;i<12;i++) {
      total_dof_per_entity[1][i] = order;
    }
    // face dof
    for (int i=0;i<6;i++) {
      total_dof_per_entity[2][i] = 2 * (order - 1) * order;
    }
    // internal dof
    total_dof_per_entity[3][0] = this->basisCardinality_ - 12 * order - 6 * total_dof_per_entity[2][0];

    int tagcur = 0;
    for (int k=0;k<closedBasis_.getCardinality();k++) {
      const int dimk = closedDofTags[k][0];
      const int entk = closedDofTags[k][1];
      for (int j=0;j<closedBasis_.getCardinality();j++) {
	const int dimj = closedDofTags[j][0];
	const int entj = closedDofTags[j][1];
	for (int i=0;i<openBasis_.getCardinality();i++) { 
	  const int dimi = openDofTags[i][0];
	  const int enti = openDofTags[i][1];
	  int dofdim;
	  int dofent;
	  ProductTopology::lineProduct3d(dimi,enti,dimj,entj,dimk,entk,dofdim,dofent);
	  tags[4*tagcur] = dofdim;
	  tags[4*tagcur+1] = dofent;
	  tags[4*tagcur+2] = current_dof_per_entity[dofdim][dofent];
	  tags[4*tagcur+3] = total_dof_per_entity[dofdim][dofent];
	  current_dof_per_entity[dofdim][dofent]++;
	  tagcur++;
	}
      }
    }
    for (int k=0;k<closedBasis_.getCardinality();k++) {
      const int dimk = closedDofTags[k][0];
      const int entk = closedDofTags[k][1];
      for (int j=0;j<openBasis_.getCardinality();j++) {
	const int dimj = openDofTags[j][0];
	const int entj = openDofTags[j][1];
	for (int i=0;i<closedBasis_.getCardinality();i++) { 
	  const int dimi = closedDofTags[i][0];
	  const int enti = closedDofTags[i][1];
	  int dofdim;
	  int dofent;
	  ProductTopology::lineProduct3d(dimi,enti,dimj,entj,dimk,entk,dofdim,dofent);
	  tags[4*tagcur] = dofdim;
	  tags[4*tagcur+1] = dofent;
	  tags[4*tagcur+2] = current_dof_per_entity[dofdim][dofent];
	  tags[4*tagcur+3] = total_dof_per_entity[dofdim][dofent];
	  current_dof_per_entity[dofdim][dofent]++;
	  tagcur++;
	}
      }
    }
    for (int k=0;k<openBasis_.getCardinality();k++) {
      const int dimk = openDofTags[k][0];
      const int entk = openDofTags[k][1];
      for (int j=0;j<closedBasis_.getCardinality();j++) {
	const int dimj = closedDofTags[j][0];
	const int entj = closedDofTags[j][1];
	for (int i=0;i<closedBasis_.getCardinality();i++) { 
	  const int dimi = closedDofTags[i][0];
	  const int enti = closedDofTags[i][1];
	  int dofdim;
	  int dofent;
	  ProductTopology::lineProduct3d(dimi,enti,dimj,entj,dimk,entk,dofdim,dofent);
	  tags[4*tagcur] = dofdim;
	  tags[4*tagcur+1] = dofent;
	  tags[4*tagcur+2] = current_dof_per_entity[dofdim][dofent];
	  tags[4*tagcur+3] = total_dof_per_entity[dofdim][dofent];
	  current_dof_per_entity[dofdim][dofent]++;
	  tagcur++;
	}
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
  void Basis_HCURL_HEX_In_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
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
    FieldContainer<Scalar> zPoints(dim0,1);
    
    for (int i=0;i<dim0;i++) {
      xPoints(i,0) = inputPoints(i,0);
      yPoints(i,0) = inputPoints(i,1);
      zPoints(i,0) = inputPoints(i,2);
    }
    
    switch (operatorType) {
    case OPERATOR_VALUE:
      {
	FieldContainer<Scalar> closedBasisValsXPts( closedBasis_.getCardinality() , dim0 );
	FieldContainer<Scalar> closedBasisValsYPts( closedBasis_.getCardinality() , dim0 );
	FieldContainer<Scalar> closedBasisValsZPts( closedBasis_.getCardinality() , dim0 );
	FieldContainer<Scalar> openBasisValsXPts( openBasis_.getCardinality() , dim0 );
	FieldContainer<Scalar> openBasisValsYPts( openBasis_.getCardinality() , dim0 );
	FieldContainer<Scalar> openBasisValsZPts( openBasis_.getCardinality() , dim0 );
	
	closedBasis_.getValues( closedBasisValsXPts , xPoints , OPERATOR_VALUE );
	closedBasis_.getValues( closedBasisValsYPts , yPoints , OPERATOR_VALUE );
	closedBasis_.getValues( closedBasisValsZPts , zPoints , OPERATOR_VALUE );
	openBasis_.getValues( openBasisValsXPts , xPoints , OPERATOR_VALUE );
	openBasis_.getValues( openBasisValsYPts , yPoints , OPERATOR_VALUE );
	openBasis_.getValues( openBasisValsZPts , zPoints , OPERATOR_VALUE );

	// first we get the "horizontal basis functions that are tangent to the x-varying edges
	int bfcur = 0;
	for (int k=0;k<closedBasis_.getCardinality();k++) {
	  for (int j=0;j<closedBasis_.getCardinality();j++) {
	    for (int i=0;i<openBasis_.getCardinality();i++) {
	      for (int l=0;l<dim0;l++) {
		outputValues(bfcur,l,0) = openBasisValsXPts(i,l) * closedBasisValsYPts(j,l) * closedBasisValsZPts(k,l);
		outputValues(bfcur,l,1) = 0.0;
		outputValues(bfcur,l,2) = 0.0;
	      }
	      bfcur++;
	    }
	  }
	}
	
	// second we get the basis functions in the direction of the y-varying edges
	for (int k=0;k<closedBasis_.getCardinality();k++) {
	  for (int j=0;j<openBasis_.getCardinality();j++) {
	    for (int i=0;i<closedBasis_.getCardinality();i++) {
	      for (int l=0;l<dim0;l++) {
		outputValues(bfcur,l,0) = 0.0;
		outputValues(bfcur,l,1) = closedBasisValsXPts(i,l) * openBasisValsYPts(j,l) * closedBasisValsZPts(k,l);
		outputValues(bfcur,l,2) = 0.0;
	      }
	      bfcur++;
	    }
	  }
	}	

	// third we get the basis functions in the direction of the y-varying edges
	for (int k=0;k<openBasis_.getCardinality();k++) {
	  for (int j=0;j<closedBasis_.getCardinality();j++) {
	    for (int i=0;i<closedBasis_.getCardinality();i++) {
	      for (int l=0;l<dim0;l++) {
		outputValues(bfcur,l,0) = 0.0;
		outputValues(bfcur,l,1) = 0.0;
		outputValues(bfcur,l,2) = closedBasisValsXPts(i,l) * closedBasisValsYPts(j,l) * openBasisValsZPts(k,l);
	      }
	      bfcur++;
	    }
	  }
	}
      }
      break;
    case OPERATOR_CURL:
      {
	FieldContainer<Scalar> closedBasisValsXPts( closedBasis_.getCardinality() , dim0 );
	FieldContainer<Scalar> closedBasisValsYPts( closedBasis_.getCardinality() , dim0 );
	FieldContainer<Scalar> closedBasisValsZPts( closedBasis_.getCardinality() , dim0 );
	FieldContainer<Scalar> closedBasisDerivsXPts( closedBasis_.getCardinality() , dim0 , 1 );
	FieldContainer<Scalar> closedBasisDerivsYPts( closedBasis_.getCardinality() , dim0 , 1 );
	FieldContainer<Scalar> closedBasisDerivsZPts( closedBasis_.getCardinality() , dim0 , 1 );
	FieldContainer<Scalar> openBasisValsXPts( openBasis_.getCardinality() , dim0 );
	FieldContainer<Scalar> openBasisValsYPts( openBasis_.getCardinality() , dim0 );
	FieldContainer<Scalar> openBasisValsZPts( openBasis_.getCardinality() , dim0 );
	
	closedBasis_.getValues( closedBasisValsXPts , xPoints , OPERATOR_VALUE );
	closedBasis_.getValues( closedBasisValsYPts , yPoints , OPERATOR_VALUE );
	closedBasis_.getValues( closedBasisValsZPts , zPoints , OPERATOR_VALUE );
	closedBasis_.getValues( closedBasisDerivsXPts , xPoints , OPERATOR_D1 );
	closedBasis_.getValues( closedBasisDerivsYPts , yPoints , OPERATOR_D1 );
	closedBasis_.getValues( closedBasisDerivsZPts , zPoints , OPERATOR_D1 );
	openBasis_.getValues( openBasisValsXPts , xPoints , OPERATOR_VALUE );
	openBasis_.getValues( openBasisValsYPts , yPoints , OPERATOR_VALUE );
	openBasis_.getValues( openBasisValsZPts , zPoints , OPERATOR_VALUE );
	
	int bfcur = 0;

	// first we get the basis functions that are tangent to the x-varying edges
	for (int k=0;k<closedBasis_.getCardinality();k++) {
	  for (int j=0;j<closedBasis_.getCardinality();j++) {
	    for (int i=0;i<openBasis_.getCardinality();i++) {
	      for (int l=0;l<dim0;l++) {
		outputValues(bfcur,l,0) = 0.0;
		outputValues(bfcur,l,1) = openBasisValsXPts(i,l) * closedBasisValsYPts(j,l) * closedBasisDerivsZPts(k,l,0);
		outputValues(bfcur,l,2) = -openBasisValsXPts(i,l) * closedBasisDerivsYPts(j,l,0) * closedBasisValsZPts(k,l);
	      }
	      bfcur++;
	    }
	  }
	}
	
	// second we get the basis functions in the direction of the y-varying edges
	for (int k=0;k<closedBasis_.getCardinality();k++) {
	  for (int j=0;j<openBasis_.getCardinality();j++) {
	    for (int i=0;i<closedBasis_.getCardinality();i++) {
	      for (int l=0;l<dim0;l++) {
		outputValues(bfcur,l,0) = -closedBasisValsXPts(i,l) * openBasisValsYPts(j,l) * closedBasisDerivsZPts(k,l,0);
		outputValues(bfcur,l,1) = 0.0;
		outputValues(bfcur,l,2) = closedBasisDerivsXPts(i,l,0) * openBasisValsYPts(j,l) * closedBasisValsZPts(k,l);
	      }
	      bfcur++;
	    }
	  }
	}	

	// third we get the basis functions in the direction of the y-varying edges
	for (int k=0;k<openBasis_.getCardinality();k++) {
	  for (int j=0;j<closedBasis_.getCardinality();j++) {
	    for (int i=0;i<closedBasis_.getCardinality();i++) {
	      for (int l=0;l<dim0;l++) {
		outputValues(bfcur,l,0) = closedBasisValsXPts(i,l) * closedBasisDerivsYPts(j,l,0) * openBasisValsZPts(k,l);
		outputValues(bfcur,l,1) = -closedBasisDerivsXPts(i,l,0) * closedBasisValsYPts(j,l) * openBasisValsZPts(k,l);
		outputValues(bfcur,l,2) = 0.0;
	      }
	      bfcur++;
	    }
	  }
	}		
      }
      break;
    case OPERATOR_DIV:
      TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
			  ">>> ERROR (Basis_HCURL_HEX_In_FEM): DIV is invalid operator for HCURL Basis Functions");
      break;
      
    case OPERATOR_GRAD:
      TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_GRAD), std::invalid_argument,
			  ">>> ERROR (Basis_HCURL_HEX_In_FEM): GRAD is invalid operator for HCURL Basis Functions");
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
      TEUCHOS_TEST_FOR_EXCEPTION( ( (operatorType == OPERATOR_D1)    ||
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
			  ">>> ERROR (Basis_HCURL_HEX_In_FEM): Invalid operator type");
      break;
      
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( (operatorType != OPERATOR_VALUE) &&
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
			  ">>> ERROR (Basis_HCURL_HEX_In_FEM): Invalid operator type");
    }
  }
  
  
  
  template<class Scalar, class ArrayScalar>
  void Basis_HCURL_HEX_In_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
							      const ArrayScalar &    inputPoints,
							      const ArrayScalar &    cellVertices,
							      const EOperator        operatorType) const {
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
			">>> ERROR (Basis_HCURL_HEX_In_FEM): FEM Basis calling an FVD member function");
  }

 template<class Scalar, class ArrayScalar>
 void Basis_HCURL_HEX_In_FEM<Scalar,ArrayScalar>::getDofCoords(ArrayScalar & DofCoords) const
  {
    // x-component basis functions
    int cur = 0;

    for (int k=0;k<closedPts_.dimension(0);k++)
      {
        for (int j=0;j<closedPts_.dimension(0);j++)
          {
            for (int i=0;i<openPts_.dimension(0);i++)
              {
                DofCoords(cur,0) = openPts_(i,0);
                DofCoords(cur,1) = closedPts_(j,0);
                DofCoords(cur,2) = closedPts_(k,0);
                cur++;
              }
          }
      }
   // y-component basis functions
    for (int k=0;k<closedPts_.dimension(0);k++)
      {
        for (int j=0;j<openPts_.dimension(0);j++)
          {
            for (int i=0;i<closedPts_.dimension(0);i++)
              {
                DofCoords(cur,0) = closedPts_(i,0);
                DofCoords(cur,1) = openPts_(j,0);
                DofCoords(cur,2) = closedPts_(k,0);
                cur++;
              }
          }
      }

    // z-component basis functions
    for (int k=0;k<openPts_.dimension(0);k++)
      {
        for (int j=0;j<closedPts_.dimension(0);j++)
          {
            for (int i=0;i<closedPts_.dimension(0);i++)
              {
                DofCoords(cur,0) = closedPts_(i,0);
                DofCoords(cur,1) = closedPts_(j,0);
                DofCoords(cur,2) = openPts_(k,0);
                cur++;
              }
          }
      }

    return;
  }
  
}// namespace Intrepid
