#ifndef INTREPID_HGRAD_HEX_CN_FEMDEF_HPP
#define INTREPID_HGRAD_HEX_CN_FEMDEF_HPP
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
//                    Denis Ridzal (dridzal@sandia.gov) or
//                    Robert Kirby (robert.c.kirby@ttu.edu)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HGRAD_HEX_Cn_FEMDef.hpp
    \brief  Definition file for the Intrepid::HGRAD_HEX_Cn_FEM class.
    \author Created by R. Kirby.
 */

namespace Intrepid {
  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_HEX_Cn_FEM<Scalar,ArrayScalar>::Basis_HGRAD_HEX_Cn_FEM( const int orderx , 
								      const int ordery ,
								      const int orderz ,
								      const ArrayScalar &pts_x ,
								      const ArrayScalar &pts_y ,
								      const ArrayScalar &pts_z ):
    bases_( 3 )
  {

    bases_[0] = Teuchos::rcp( new Basis_HGRAD_LINE_Cn_FEM< Scalar , ArrayScalar >( orderx , pts_x ) );
    bases_[1] = Teuchos::rcp( new Basis_HGRAD_LINE_Cn_FEM< Scalar , ArrayScalar >( ordery , pts_y ) );
    bases_[2] = Teuchos::rcp( new Basis_HGRAD_LINE_Cn_FEM< Scalar , ArrayScalar >( orderz , pts_z ) );

    this->basisCardinality_ = (orderx+1)*(ordery+1)*(orderz+1);
    if (orderx >= ordery && orderx >= orderz ) {
      this->basisDegree_ = orderx;
    }
    else if (ordery >= orderx && ordery >= orderz) {
      this->basisDegree_ = ordery;
    }
    else {
      this->basisDegree_ = orderz;
    }
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
    this -> basisType_         = BASIS_FEM_FIAT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;
  }

  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_HEX_Cn_FEM<Scalar,ArrayScalar>::Basis_HGRAD_HEX_Cn_FEM( const int order , 
								      const EPointType & pointType ):
    bases_( 3 )
  {
    bases_[0] = Teuchos::rcp( new Basis_HGRAD_LINE_Cn_FEM< Scalar , ArrayScalar >( order , pointType ) );
    // basis is same in each direction, so I only need to instantiate it once!
    bases_[1] = bases_[0];
    bases_[2] = bases_[0];
    this->basisCardinality_ = (order+1)*(order+1)*(order+1);
    this->basisDegree_ = order;
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
    this -> basisType_         = BASIS_FEM_FIAT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;
  }

  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_HEX_Cn_FEM<Scalar,ArrayScalar>::initializeTags()
  {
    // Basis-dependent initializations
    int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
    int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
    int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
    int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
  
    // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 

    std::vector<int> tags( tagSize * this->getCardinality() );

    // temporarily just put everything on the cell itself
    for (int i=0;i<this->getCardinality();i++) {
       tags[tagSize*i] = 3;
       tags[tagSize*i+1] = 0;
       tags[tagSize*i+2] = i;
       tags[tagSize*i+3] = this->getCardinality();
     }

    Basis_HGRAD_LINE_Cn_FEM<Scalar,ArrayScalar> &xBasis_ = *bases_[0];
    Basis_HGRAD_LINE_Cn_FEM<Scalar,ArrayScalar> &yBasis_ = *bases_[1];
    Basis_HGRAD_LINE_Cn_FEM<Scalar,ArrayScalar> &zBasis_ = *bases_[2];


    // now let's try to do it "right"
    // let's get the x, y, z bases and their dof
    const std::vector<std::vector<int> >& xdoftags = xBasis_.getAllDofTags();
    const std::vector<std::vector<int> >& ydoftags = yBasis_.getAllDofTags();
    const std::vector<std::vector<int> >& zdoftags = zBasis_.getAllDofTags();

    std::map<int,std::map<int,int> > total_dof_per_entity;
    std::map<int,std::map<int,int> > current_dof_per_entity;

    // vertex dof
    for (int i=0;i<8;i++) {
      total_dof_per_entity[0][i] = 0;
      current_dof_per_entity[0][i] = 0;
    }
    // edge dof (12 edges)
    for (int i=0;i<12;i++) {
      total_dof_per_entity[1][i] = 0;
      current_dof_per_entity[1][i] = 0;
    }
    // face dof (6 faces)
    for (int i=0;i<6;i++) {
      total_dof_per_entity[2][i] = 0;
      current_dof_per_entity[2][i] = 0;
    }
    // internal dof
    total_dof_per_entity[3][0] = 0;
    //total_dof_per_entity[3][1] = 0;

    // let's tally the total degrees of freedom on each entity
    for (int k=0;k<zBasis_.getCardinality();k++) {
      const int zdim = zdoftags[k][0];
      const int zent = zdoftags[k][1];
      for (int j=0;j<yBasis_.getCardinality();j++) {
	const int ydim = ydoftags[j][0];
	const int yent = ydoftags[j][1];
	for (int i=0;i<xBasis_.getCardinality();i++) {
	  const int xdim = xdoftags[i][0];
	  const int xent = xdoftags[i][1];
	  int dofdim;
	  int dofent;
	  ProductTopology::lineProduct3d( xdim , xent , ydim , yent , zdim , zent , dofdim , dofent );
	  std::cout << i << " " << j << " " << k << "\t" << dofdim << " " << dofent << std::endl;
	  total_dof_per_entity[dofdim][dofent] += 1;

	}
      }
    }

    int tagcur = 0;
    for (int k=0;k<zBasis_.getCardinality();k++) {
      const int zdim = zdoftags[k][0];
      const int zent = zdoftags[k][1];
      for (int j=0;j<yBasis_.getCardinality();j++) {
	const int ydim = ydoftags[j][0];
	const int yent = ydoftags[j][1];
	for (int i=0;i<xBasis_.getCardinality();i++) {
	  const int xdim = xdoftags[i][0];
	  const int xent = xdoftags[i][1];
	  int dofdim;
	  int dofent;
	  ProductTopology::lineProduct3d( xdim , xent , ydim , yent , zdim , zent , dofdim , dofent );
	  tags[4*tagcur] = dofdim;
	  tags[4*tagcur+1] = dofent;
	  tags[4*tagcur+2] = current_dof_per_entity[dofdim][dofent];
	  current_dof_per_entity[dofdim][dofent]++;
	  tags[4*tagcur+3] = total_dof_per_entity[dofdim][dofent];
	  tagcur++;
	}
      }
    }

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
  void Basis_HGRAD_HEX_Cn_FEM<Scalar,ArrayScalar>::getValues( ArrayScalar &outputValues ,
							       const ArrayScalar &inputPoints ,
							       const EOperator operatorType ) const 
  {
#ifdef HAVE_INTREPID_DEBUG
  Intrepid::getValues_HGRAD_Args<Scalar, ArrayScalar>(outputValues,
                                                      inputPoints,
                                                      operatorType,
                                                      this -> getBaseCellTopology(),
                                                      this -> getCardinality() );
#endif

    Basis_HGRAD_LINE_Cn_FEM<Scalar,ArrayScalar> &xBasis_ = *bases_[0];
    Basis_HGRAD_LINE_Cn_FEM<Scalar,ArrayScalar> &yBasis_ = *bases_[1];
    Basis_HGRAD_LINE_Cn_FEM<Scalar,ArrayScalar> &zBasis_ = *bases_[2];


    FieldContainer<Scalar> xInputPoints(inputPoints.dimension(0),1);
    FieldContainer<Scalar> yInputPoints(inputPoints.dimension(0),1);
    FieldContainer<Scalar> zInputPoints(inputPoints.dimension(0),1);

    for (int i=0;i<inputPoints.dimension(0);i++) {
      xInputPoints(i,0) = inputPoints(i,0);
      yInputPoints(i,0) = inputPoints(i,1);
      zInputPoints(i,0) = inputPoints(i,2);
    }

    switch (operatorType) {
    case OPERATOR_VALUE:
      {
	FieldContainer<Scalar> xBasisValues(xBasis_.getCardinality(),xInputPoints.dimension(0));
	FieldContainer<Scalar> yBasisValues(yBasis_.getCardinality(),yInputPoints.dimension(0));
	FieldContainer<Scalar> zBasisValues(zBasis_.getCardinality(),zInputPoints.dimension(0));

	xBasis_.getValues(xBasisValues,xInputPoints,OPERATOR_VALUE);
	yBasis_.getValues(yBasisValues,yInputPoints,OPERATOR_VALUE);
	zBasis_.getValues(zBasisValues,zInputPoints,OPERATOR_VALUE);

	int bfcur = 0;
	for (int k=0;k<zBasis_.getCardinality();k++) {
	  for (int j=0;j<yBasis_.getCardinality();j++) {
	    for (int i=0;i<xBasis_.getCardinality();i++) {
	      for (int l=0;l<inputPoints.dimension(0);l++) {
		outputValues(bfcur,l) = xBasisValues(i,l) * yBasisValues(j,l) * zBasisValues(k,l);
	      }
	      bfcur++;
	    }
	  }
	}
      }
      break;
    case OPERATOR_D1:
    case OPERATOR_GRAD:
      {
        FieldContainer<Scalar> xBasisValues(xBasis_.getCardinality(),xInputPoints.dimension(0));
        FieldContainer<Scalar> yBasisValues(yBasis_.getCardinality(),yInputPoints.dimension(0));
        FieldContainer<Scalar> zBasisValues(zBasis_.getCardinality(),zInputPoints.dimension(0));
        FieldContainer<Scalar> xBasisDerivs(xBasis_.getCardinality(),xInputPoints.dimension(0),1);
        FieldContainer<Scalar> yBasisDerivs(yBasis_.getCardinality(),yInputPoints.dimension(0),1);
        FieldContainer<Scalar> zBasisDerivs(zBasis_.getCardinality(),zInputPoints.dimension(0),1);

        xBasis_.getValues(xBasisValues,xInputPoints,OPERATOR_VALUE);
        yBasis_.getValues(yBasisValues,yInputPoints,OPERATOR_VALUE);
        zBasis_.getValues(zBasisValues,zInputPoints,OPERATOR_VALUE);
        xBasis_.getValues(xBasisDerivs,xInputPoints,OPERATOR_D1);
        yBasis_.getValues(yBasisDerivs,yInputPoints,OPERATOR_D1);    	
        zBasis_.getValues(zBasisDerivs,zInputPoints,OPERATOR_D1);    	

	int bfcur = 0;
	for (int k=0;k<zBasis_.getCardinality();k++) {
	  for (int j=0;j<yBasis_.getCardinality();j++) {
	    for (int i=0;i<xBasis_.getCardinality();i++) {
	      for (int l=0;l<inputPoints.dimension(0);l++) {
		outputValues(bfcur,l,0) = xBasisDerivs(i,l,0) * yBasisValues(j,l) * zBasisValues(k,l);
		outputValues(bfcur,l,1) = xBasisValues(i,l) * yBasisDerivs(j,l,0) * zBasisValues(k,l);
		outputValues(bfcur,l,2) = xBasisValues(i,l) * yBasisValues(j,l) * zBasisDerivs(k,l,0);
	      }
	      bfcur++;
	    }
	  }
	}
      }
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
        FieldContainer<Scalar> xBasisValues(xBasis_.getCardinality(),xInputPoints.dimension(0));
        FieldContainer<Scalar> yBasisValues(yBasis_.getCardinality(),yInputPoints.dimension(0));
        FieldContainer<Scalar> zBasisValues(yBasis_.getCardinality(),zInputPoints.dimension(0));

        Teuchos::Array<int> partialMult;
	
        for (int d=0;d<getDkCardinality(operatorType,3);d++) {
          getDkMultiplicities( partialMult , d , operatorType , 3 );
          if (partialMult[0] == 0) {
            xBasisValues.resize(xBasis_.getCardinality(),xInputPoints.dimension(0));
            xBasis_.getValues( xBasisValues , xInputPoints, OPERATOR_VALUE );
          }
          else {
            xBasisValues.resize(xBasis_.getCardinality(),xInputPoints.dimension(0),1);
            EOperator xop = (EOperator) ( (int) OPERATOR_D1 + partialMult[0] - 1 );
            xBasis_.getValues( xBasisValues , xInputPoints, xop );
            xBasisValues.resize(xBasis_.getCardinality(),xInputPoints.dimension(0));
          }
          if (partialMult[1] == 0) {
            yBasisValues.resize(yBasis_.getCardinality(),yInputPoints.dimension(0));
            yBasis_.getValues( yBasisValues , yInputPoints, OPERATOR_VALUE );
          }
          else {
            yBasisValues.resize(yBasis_.getCardinality(),yInputPoints.dimension(0),1);
            EOperator yop = (EOperator) ( (int) OPERATOR_D1 + partialMult[1] - 1 );
            yBasis_.getValues( yBasisValues , yInputPoints, yop );
            yBasisValues.resize(yBasis_.getCardinality(),yInputPoints.dimension(0));
          }
          if (partialMult[2] == 0) {
            zBasisValues.resize(zBasis_.getCardinality(),zInputPoints.dimension(0));
            zBasis_.getValues( zBasisValues , zInputPoints, OPERATOR_VALUE );
          }
          else {
            zBasisValues.resize(zBasis_.getCardinality(),zInputPoints.dimension(0),1);
            EOperator zop = (EOperator) ( (int) OPERATOR_D1 + partialMult[2] - 1 );
            zBasis_.getValues( zBasisValues , zInputPoints, zop );
            zBasisValues.resize(zBasis_.getCardinality(),zInputPoints.dimension(0));
          }


          int bfcur = 0;
	  for (int k=0;k<zBasis_.getCardinality();k++) {
	    for (int j=0;j<yBasis_.getCardinality();j++) {
	      for (int i=0;i<xBasis_.getCardinality();i++) {
		for (int l=0;l<inputPoints.dimension(0);l++) {
		  outputValues(bfcur,l,d) = xBasisValues(i,l) * yBasisValues(j,l) * zBasisValues(k,l);
		}
		bfcur++;
	      }
	    }
	  }
	}
      }
      break;
    default:
        TEST_FOR_EXCEPTION( true , std::invalid_argument,
                            ">>> ERROR (Basis_HGRAD_HEX_Cn_FEM): Operator type not implemented");
        break;
    }
  }

  template<class Scalar,class ArrayScalar>
void Basis_HGRAD_HEX_Cn_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                             const ArrayScalar &    inputPoints,
                                                             const ArrayScalar &    cellVertices,
                                                             const EOperator        operatorType) const {
  TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (Basis_HGRAD_HEX_Cn_FEM): FEM Basis calling an FVD member function");
}

  
}// namespace Intrepid

#endif
