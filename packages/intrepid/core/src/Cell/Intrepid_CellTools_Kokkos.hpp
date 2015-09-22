#ifndef INTREPID_CELTOOLS_KOKKOS_HPP
#define INTREPID_CELTOOLS_KOKKOS_HPP

#ifdef HAVE_INTREPID_KOKKOSCORE

#include "Intrepid_CellTools.hpp"
/*
namespace Intrepid{

	template<class Scalar>
	template<class Scalar1,class Scalar2,class Scalar3,class Layout,class MemorySpace>
	void CellTools<Scalar>::setJacobianTemp(Kokkos::View<Scalar1,Layout,MemorySpace> &   jacobian,
                            const Kokkos::View<Scalar2,Layout,MemorySpace> &           points,
                            const Kokkos::View<Scalar3,Layout,MemorySpace>  &           cellWorkset,
                            const shards::CellTopology & cellTopo,
                            const int &                  whichCell){
  
 	CellTools<Scalar>::setJacobianTempSpecKokkos<Kokkos::View<Scalar1,Layout,MemorySpace>,Kokkos::View<Scalar2,Layout,MemorySpace>, Kokkos::View<Scalar3,Layout,MemorySpace>, Rank<Kokkos::View<Scalar2,Layout,MemorySpace> >::value ,Rank<Kokkos::View<Scalar1,Layout,MemorySpace> >::value>(jacobian, points, cellWorkset, cellTopo, whichCell);
     }

template<class ViewType1, class ViewType2>
struct points_temppointsCopy2 {
  ViewType1 tempPoints;
  typename ViewType2::const_type points;
 
  points_temppointsCopy2(ViewType1 tempPoints_, ViewType2 points_):tempPoints(tempPoints_),points(points_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int pt) const {
    for(int dm = 0; dm < points.dimension(1) ; dm++){
              tempPoints(pt, dm) = points(pt, dm);
            }
            }
};
template<class ViewType1, class ViewType2,class ViewType3>
struct calculateJacobian2_4 {
  ViewType1 jacobian;
  typename ViewType2::const_type cellWorkset;
  typename ViewType3::const_type basisGrads;
  calculateJacobian2_4(ViewType1 jacobian_, ViewType2 cellWorkset_, ViewType3 basisGrads_):jacobian(jacobian_),cellWorkset(cellWorkset_),basisGrads(basisGrads_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int cellOrd, int numPoints, int spaceDim, int basisCardinality) const {
             for(int pointOrd = 0; pointOrd < numPoints; pointOrd++) {
                for(int row = 0; row < spaceDim; row++){
                  for(int col = 0; col < spaceDim; col++){
                    for(int bfOrd = 0; bfOrd < basisCardinality; bfOrd++){
                      jacobian(cellOrd, pointOrd, row, col) += cellWorkset(cellOrd, bfOrd, row)*basisGrads(bfOrd, pointOrd, col);
                    } // bfOrd
                  } // col
                } // row
              } // pointOrd
            }
};
	template<class Scalar>
	template< class ArrayJac, class ArrayPoint, class ArrayCell>
	struct CellTools<Scalar>::setJacobianTempSpecKokkos<ArrayJac,ArrayPoint, ArrayCell, 2 ,4>{
	setJacobianTempSpecKokkos(ArrayJac &   jacobian,
                            const ArrayPoint &           points,
                            const ArrayCell  &           cellWorkset,
                            const shards::CellTopology & cellTopo,
                            const int &                  whichCell){
      INTREPID_VALIDATE( validateArguments_setJacobian(jacobian, points, cellWorkset, whichCell,  cellTopo) );

    int spaceDim  = (int)cellTopo.getDimension();
    int numCells  = cellWorkset.dimension(0);
    //points can be rank-2 (P,D), or rank-3 (C,P,D)
    int numPoints = points.dimension(0);
    
    // Jacobian is computed using gradients of an appropriate H(grad) basis function: define RCP to the base class
    Teuchos::RCP< Basis< Scalar, FieldContainer<Scalar> > > HGRAD_Basis;
    
    // Choose the H(grad) basis depending on the cell topology. \todo define maps for shells and beams
    switch( cellTopo.getKey() ){
      
      // Standard Base topologies (number of cellWorkset = number of vertices)
      case shards::Line<2>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_LINE_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Triangle<3>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TRI_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Quadrilateral<4>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_QUAD_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Tetrahedron<4>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Hexahedron<8>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Wedge<6>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Pyramid<5>::key:
	    HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_PYR_C1_FEM<Scalar, FieldContainer<Scalar> >() );
	    break;
        
      // Standard Extended topologies
      case shards::Triangle<6>::key:    
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TRI_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
      case shards::Quadrilateral<9>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_QUAD_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Tetrahedron<10>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Tetrahedron<11>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_COMP12_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Hexahedron<20>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_I2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Hexahedron<27>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Wedge<15>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_I2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Wedge<18>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Pyramid<13>::key:
	HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_PYR_I2_FEM<Scalar, FieldContainer<Scalar> >() );
	break;
        
        // These extended topologies are not used for mapping purposes
      case shards::Quadrilateral<8>::key:
        TEUCHOS_TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported. ");
        break;
        
        // Base and Extended Line, Beam and Shell topologies  
      case shards::Line<3>::key:
      case shards::Beam<2>::key:
      case shards::Beam<3>::key:
      case shards::ShellLine<2>::key:
      case shards::ShellLine<3>::key:
      case shards::ShellTriangle<3>::key:
      case shards::ShellTriangle<6>::key:
      case shards::ShellQuadrilateral<4>::key:
      case shards::ShellQuadrilateral<8>::key:
      case shards::ShellQuadrilateral<9>::key:
        TEUCHOS_TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported. ");
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported.");        
    }// switch  
    
    // Temp (F,P,D) array for the values of basis functions gradients at the reference points
    int basisCardinality = HGRAD_Basis -> getCardinality();
  //  FieldContainer<Scalar> basisGrads(basisCardinality, numPoints, spaceDim);
    Kokkos::View<Scalar***> basisGrads(basisCardinality, numPoints, spaceDim);
    // Initialize jacobian

    Kokkos::deep_copy( jacobian, 0.0); 
	
          // getValues requires rank-2 (P,D) input array, but points cannot be passed directly as argument because they are a user type
          Kokkos::View<Scalar**>tempPoints(points.dimension(0), points.dimension(1));
         // FieldContainer<Scalar> tempPoints( points.dimension(0), points.dimension(1) );
          // Copy point set corresponding to this cell oridinal to the temp (P,D) array
        
         parallel_for(points.dimension(0),points_temppointsCopy2<Kokkos::View<Scalar**>,ArrayPoint>(tempPoints,points));
          HGRAD_Basis -> getValues(basisGrads, tempPoints, OPERATOR_GRAD);
          
          // The outer loops select the multi-index of the Jacobian entry: cell, point, row, col
          // If whichCell = -1, all jacobians are computed, otherwise a single cell jacobian is computed
          int cellLoop = numCells  ;
          
          parallel_for(cellLoop,calculateJacobian2_4<ArrayJac,ArrayCell, Kokkos::View<Scalar***> >(jacobian, cellWorkset, basisGrads),numPoints,spaceDim,basisCardinality);
         }
         };

template<class ViewType1, class ViewType2,class ViewType3>
struct calculateJacobian2_3 {
  ViewType1 jacobian;
  typename ViewType2::const_type cellWorkset;
  typename ViewType3::const_type basisGrads;
  calculateJacobian2_3(ViewType1 jacobian_, ViewType2 cellWorkset_, ViewType3 basisGrads_):jacobian(jacobian_),cellWorkset(cellWorkset_),basisGrads(basisGrads_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int pointOrd, int numPoints, int spaceDim, int basisCardinality,int whichCell) const {
             
                for(int row = 0; row < spaceDim; row++){
                  for(int col = 0; col < spaceDim; col++){
                    for(int bfOrd = 0; bfOrd < basisCardinality; bfOrd++){
                      jacobian(pointOrd, row, col) += cellWorkset(whichCell, bfOrd, row)*basisGrads(bfOrd, pointOrd, col);
                    } // bfOrd
                  } // col
                } // row
           
            }
};
	template<class Scalar>
	template< class ArrayJac, class ArrayPoint, class ArrayCell>
	struct CellTools<Scalar>::setJacobianTempSpecKokkos<ArrayJac,ArrayPoint, ArrayCell, 2 ,3>{
	setJacobianTempSpecKokkos(ArrayJac &   jacobian,
                            const ArrayPoint &           points,
                            const ArrayCell  &           cellWorkset,
                            const shards::CellTopology & cellTopo,
                            const int &                  whichCell){
      INTREPID_VALIDATE( validateArguments_setJacobian(jacobian, points, cellWorkset, whichCell,  cellTopo) );

    int spaceDim  = (int)cellTopo.getDimension();
    int numCells  = cellWorkset.dimension(0);
    //points can be rank-2 (P,D), or rank-3 (C,P,D)
    int numPoints = points.dimension(0);
    
    // Jacobian is computed using gradients of an appropriate H(grad) basis function: define RCP to the base class
    Teuchos::RCP< Basis< Scalar, FieldContainer<Scalar> > > HGRAD_Basis;
    
    // Choose the H(grad) basis depending on the cell topology. \todo define maps for shells and beams
    switch( cellTopo.getKey() ){
      
      // Standard Base topologies (number of cellWorkset = number of vertices)
      case shards::Line<2>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_LINE_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Triangle<3>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TRI_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Quadrilateral<4>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_QUAD_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Tetrahedron<4>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Hexahedron<8>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Wedge<6>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Pyramid<5>::key:
	    HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_PYR_C1_FEM<Scalar, FieldContainer<Scalar> >() );
	    break;
        
      // Standard Extended topologies
      case shards::Triangle<6>::key:    
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TRI_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
      case shards::Quadrilateral<9>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_QUAD_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Tetrahedron<10>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Tetrahedron<11>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_COMP12_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Hexahedron<20>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_I2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Hexahedron<27>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Wedge<15>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_I2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Wedge<18>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Pyramid<13>::key:
	HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_PYR_I2_FEM<Scalar, FieldContainer<Scalar> >() );
	break;
        
        // These extended topologies are not used for mapping purposes
      case shards::Quadrilateral<8>::key:
        TEUCHOS_TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported. ");
        break;
        
        // Base and Extended Line, Beam and Shell topologies  
      case shards::Line<3>::key:
      case shards::Beam<2>::key:
      case shards::Beam<3>::key:
      case shards::ShellLine<2>::key:
      case shards::ShellLine<3>::key:
      case shards::ShellTriangle<3>::key:
      case shards::ShellTriangle<6>::key:
      case shards::ShellQuadrilateral<4>::key:
      case shards::ShellQuadrilateral<8>::key:
      case shards::ShellQuadrilateral<9>::key:
        TEUCHOS_TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported. ");
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported.");        
    }// switch  
    
    // Temp (F,P,D) array for the values of basis functions gradients at the reference points
    int basisCardinality = HGRAD_Basis -> getCardinality();
  //  FieldContainer<Scalar> basisGrads(basisCardinality, numPoints, spaceDim);
    Kokkos::View<Scalar***> basisGrads(basisCardinality, numPoints, spaceDim);
    // Initialize jacobian

    Kokkos::deep_copy( jacobian, 0.0); 
	
          // getValues requires rank-2 (P,D) input array, but points cannot be passed directly as argument because they are a user type
          Kokkos::View<Scalar**>tempPoints(points.dimension(0), points.dimension(1));
         // FieldContainer<Scalar> tempPoints( points.dimension(0), points.dimension(1) );
          // Copy point set corresponding to this cell oridinal to the temp (P,D) array
        
         parallel_for(points.dimension(0),points_temppointsCopy2<Kokkos::View<Scalar**>,ArrayPoint>(tempPoints,points));
          HGRAD_Basis -> getValues(basisGrads, tempPoints, OPERATOR_GRAD);
          
          // The outer loops select the multi-index of the Jacobian entry: point, row, col
         parallel_for(numPoints,calculateJacobian2_3<ArrayJac,ArrayCell, Kokkos::View<Scalar***> >(jacobian, cellWorkset, basisGrads),numPoints,spaceDim,basisCardinality,whichCell);
         }
         }; 
                            
template<class ViewType1, class ViewType2,class ViewType3>
struct calculateJacobian3_4 {
  ViewType1 jacobian;
  typename ViewType2::const_type cellWorkset;
  typename ViewType3::const_type basisGrads;
  calculateJacobian3_4(ViewType1 jacobian_, ViewType2 cellWorkset_, ViewType3 basisGrads_):jacobian(jacobian_),cellWorkset(cellWorkset_),basisGrads(basisGrads_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int pointOrd, int numPoints, int spaceDim, int basisCardinality,int whichCell) const {
             
                for(int row = 0; row < spaceDim; row++){
                  for(int col = 0; col < spaceDim; col++){
                    for(int bfOrd = 0; bfOrd < basisCardinality; bfOrd++){
                      jacobian(pointOrd, row, col) += cellWorkset(whichCell, bfOrd, row)*basisGrads(bfOrd, pointOrd, col);
                    } // bfOrd
                  } // col
                } // row
           
            }
};              
    template<class Scalar>             
	template<class ArrayJac, class ArrayPoint, class ArrayCell>
	struct CellTools<Scalar>::setJacobianTempSpecKokkos<ArrayJac,ArrayPoint, ArrayCell, 3 , 4>{
	setJacobianTempSpecKokkos(ArrayJac &   jacobian,
                            const ArrayPoint &           points,
                            const ArrayCell  &           cellWorkset,
                            const shards::CellTopology & cellTopo,
                            const int &                  whichCell){
            INTREPID_VALIDATE( validateArguments_setJacobian(jacobian, points, cellWorkset, whichCell,  cellTopo) );
    
    int spaceDim  = (int)cellTopo.getDimension();
    int numCells  = cellWorkset.dimension(0);

    int numPoints = points.dimension(1);
    
    // Jacobian is computed using gradients of an appropriate H(grad) basis function: define RCP to the base class
    Teuchos::RCP< Basis< Scalar, FieldContainer<Scalar> > > HGRAD_Basis;
    
    // Choose the H(grad) basis depending on the cell topology. \todo define maps for shells and beams
    switch( cellTopo.getKey() ){
      
      // Standard Base topologies (number of cellWorkset = number of vertices)
      case shards::Line<2>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_LINE_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Triangle<3>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TRI_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Quadrilateral<4>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_QUAD_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Tetrahedron<4>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Hexahedron<8>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Wedge<6>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_C1_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Pyramid<5>::key:
	    HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_PYR_C1_FEM<Scalar, FieldContainer<Scalar> >() );
	    break;
        
      // Standard Extended topologies
      case shards::Triangle<6>::key:    
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TRI_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
      case shards::Quadrilateral<9>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_QUAD_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Tetrahedron<10>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Tetrahedron<11>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_COMP12_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Hexahedron<20>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_I2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Hexahedron<27>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Wedge<15>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_I2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;
        
      case shards::Wedge<18>::key:
        HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_C2_FEM<Scalar, FieldContainer<Scalar> >() );
        break;

      case shards::Pyramid<13>::key:
	HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_PYR_I2_FEM<Scalar, FieldContainer<Scalar> >() );
	break;
        
        // These extended topologies are not used for mapping purposes
      case shards::Quadrilateral<8>::key:
        TEUCHOS_TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported. ");
        break;
        
        // Base and Extended Line, Beam and Shell topologies  
      case shards::Line<3>::key:
      case shards::Beam<2>::key:
      case shards::Beam<3>::key:
      case shards::ShellLine<2>::key:
      case shards::ShellLine<3>::key:
      case shards::ShellTriangle<3>::key:
      case shards::ShellTriangle<6>::key:
      case shards::ShellQuadrilateral<4>::key:
      case shards::ShellQuadrilateral<8>::key:
      case shards::ShellQuadrilateral<9>::key:
        TEUCHOS_TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported. ");
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                            ">>> ERROR (Intrepid::CellTools::setJacobian): Cell topology not supported.");        
    }// switch  
    
    // Temp (F,P,D) array for the values of basis functions gradients at the reference points
    int basisCardinality = HGRAD_Basis -> getCardinality();
    FieldContainer<Scalar> basisGrads(basisCardinality, numPoints, spaceDim);
    
    // Initialize jacobian
    for(int i = 0; i < jacobian.size(); i++){
      jacobian[i] = 0.0;
    }
          
          // getValues requires rank-2 (P,D) input array, refPoints cannot be used as argument: need temp (P,D) array
          FieldContainer<Scalar> tempPoints( points.dimension(1), points.dimension(2) );
          
          for(int cellOrd = 0; cellOrd < numCells; cellOrd++) {
            
            // Copy point set corresponding to this cell oridinal to the temp (P,D) array
            for(int pt = 0; pt < points.dimension(1); pt++){
              for(int dm = 0; dm < points.dimension(2) ; dm++){
                tempPoints(pt, dm) = points(cellOrd, pt, dm);
              }//dm
            }//pt
            
            // Compute gradients of basis functions at this set of ref. points
            HGRAD_Basis -> getValues(basisGrads, tempPoints, OPERATOR_GRAD);
            
            // Compute jacobians for the point set corresponding to the current cellordinal
            for(int pointOrd = 0; pointOrd < numPoints; pointOrd++) {
              for(int row = 0; row < spaceDim; row++){
                for(int col = 0; col < spaceDim; col++){
                  
                  // The entry is computed by contracting the basis index. Number of basis functions and vertices must be the same
                  for(int bfOrd = 0; bfOrd < basisCardinality; bfOrd++){
                    jacobian(cellOrd, pointOrd, row, col) += cellWorkset(cellOrd, bfOrd, row)*basisGrads(bfOrd, pointOrd, col);
                  } // bfOrd
                } // col
              } // row
            } // pointOrd
          }//cellOrd
 	  }//end function
       };                     
}//end namespace Intrepid_Kokkos  
*/
#endif

#endif
