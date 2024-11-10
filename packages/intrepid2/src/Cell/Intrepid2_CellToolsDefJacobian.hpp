// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file   Intrepid2_CellToolsDefJacobian.hpp
    \brief  Definition file for the Jacobian functions in the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/
#ifndef __INTREPID2_CELLTOOLS_DEF_HPP__
#define __INTREPID2_CELLTOOLS_DEF_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

#include "Intrepid2_Kernels.hpp"
#include "Intrepid2_DataTools.hpp"

namespace Intrepid2 {


  //============================================================================================//
  //                                                                                            //
  //                     Jacobian, inverse Jacobian and Jacobian determinant                    //
  //                                                                                            //
  //============================================================================================//

  namespace FunctorCellTools {
    /**
     \brief Functor for calculation of Jacobian on cell workset see Intrepid2::CellTools for more
    */
    template<typename jacobianViewType, 
             typename worksetCellType, 
             typename basisGradType>
    struct F_setJacobian {
            jacobianViewType _jacobian;
      const worksetCellType  _worksetCells;
      const basisGradType    _basisGrads;
      const int              _startCell;
      const int              _endCell;

      KOKKOS_INLINE_FUNCTION
      F_setJacobian( jacobianViewType jacobian_,
                     worksetCellType  worksetCells_,
                     basisGradType    basisGrads_,
                     const int        startCell_,
                     const int        endCell_)
        : _jacobian(jacobian_), _worksetCells(worksetCells_), _basisGrads(basisGrads_),
          _startCell(startCell_), _endCell(endCell_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cell,
                      const ordinal_type point) const {
        
        const ordinal_type dim = _jacobian.extent(2); // dim2 and dim3 should match
        
        const ordinal_type gradRank = rank(_basisGrads);
        if ( gradRank == 3)
        {
          const ordinal_type cardinality = _basisGrads.extent(0);
          for (ordinal_type i=0;i<dim;++i)
            for (ordinal_type j=0;j<dim;++j) {
              _jacobian(cell, point, i, j) = 0;
              for (ordinal_type bf=0;bf<cardinality;++bf)
                _jacobian(cell, point, i, j) += _worksetCells(cell+_startCell, bf, i) * _basisGrads(bf, point, j);
            }
        }
        else
        {
          const ordinal_type cardinality = _basisGrads.extent(1);
          for (ordinal_type i=0;i<dim;++i)
          for (ordinal_type j=0;j<dim;++j) {
            _jacobian(cell, point, i, j) = 0;
            for (ordinal_type bf=0;bf<cardinality;++bf)
              _jacobian(cell, point, i, j) += _worksetCells(cell+_startCell, bf, i) * _basisGrads(cell, bf, point, j);
          }
        }
      }
    };
  }

  template<typename DeviceType>
  template<class PointScalar>
  Data<PointScalar,DeviceType> CellTools<DeviceType>::allocateJacobianDet( const Data<PointScalar,DeviceType> & jacobian )
  {
    auto extents           = jacobian.getExtents(); // C,P,D,D, which we reduce to C,P
    auto variationTypes    = jacobian.getVariationTypes();
    const bool cellVaries  = (variationTypes[0] != CONSTANT);
    const bool pointVaries = (variationTypes[1] != CONSTANT);
    
    extents[2]         = 1;
    extents[3]         = 1;
    variationTypes[2]  = CONSTANT;
    variationTypes[3]  = CONSTANT;
    
    if ( cellVaries && pointVaries )
    {
      auto data = jacobian.getUnderlyingView4();
      auto detData = getMatchingViewWithLabel(data, "Jacobian det data", data.extent_int(0), data.extent_int(1));
      return Data<PointScalar,DeviceType>(detData,2,extents,variationTypes);
    }
    else if (cellVaries || pointVaries)
    {
      auto data = jacobian.getUnderlyingView3();
      auto detData = getMatchingViewWithLabel(data, "Jacobian det data", data.extent_int(0));
      return Data<PointScalar,DeviceType>(detData,2,extents,variationTypes);
    }
    else
    {
      auto data = jacobian.getUnderlyingView1();
      auto detData = getMatchingViewWithLabel(data, "Jacobian det data", 1);
      return Data<PointScalar,DeviceType>(detData,2,extents,variationTypes);
    }
  }

  template<typename DeviceType>
  template<class PointScalar>
  Data<PointScalar,DeviceType> CellTools<DeviceType>::allocateJacobianInv( const Data<PointScalar,DeviceType> & jacobian )
  {
    auto extents        = jacobian.getExtents(); // C,P,D,D
    auto variationTypes = jacobian.getVariationTypes();
    int jacDataRank     = jacobian.getUnderlyingViewRank();
      
    if ( jacDataRank == 4 )
    {
      auto jacData = jacobian.getUnderlyingView4();
      auto invData = getMatchingViewWithLabel(jacData, "Jacobian inv data",jacData.extent(0),jacData.extent(1),jacData.extent(2),jacData.extent(3));
      return Data<PointScalar,DeviceType>(invData,4,extents,variationTypes);
    }
    else if (jacDataRank == 3)
    {
      auto jacData = jacobian.getUnderlyingView3();
      auto invData = getMatchingViewWithLabel(jacData, "Jacobian inv data",jacData.extent(0),jacData.extent(1),jacData.extent(2));
      return Data<PointScalar,DeviceType>(invData,4,extents,variationTypes);
    }
    else if (jacDataRank == 2)
    {
      auto jacData = jacobian.getUnderlyingView2();
      auto invData = getMatchingViewWithLabel(jacData, "Jacobian inv data",jacData.extent(0),jacData.extent(1));
      return Data<PointScalar,DeviceType>(invData,4,extents,variationTypes);
    }
    else if (jacDataRank == 1)
    {
      auto jacData = jacobian.getUnderlyingView1();
      auto invData = getMatchingViewWithLabel(jacData, "Jacobian inv data",jacData.extent(0));
      return Data<PointScalar,DeviceType>(invData,4,extents,variationTypes);
    }
    else
    {
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "allocateJacobianInv requires jacobian to vary in *some* dimension…");
      return Data<PointScalar,DeviceType>(); // unreachable statement; this line added to avoid compiler warning on CUDA
    }
  }

  template<typename DeviceType>
  template<class PointScalar>
  void CellTools<DeviceType>::setJacobianDet( Data<PointScalar,DeviceType> &jacobianDet, const Data<PointScalar,DeviceType> & jacobian )
  {
    auto variationTypes = jacobian.getVariationTypes();
    
    const int CELL_DIM  = 0;
    const int POINT_DIM = 1;
    const int D1_DIM    = 2;
    const bool cellVaries  = (variationTypes[CELL_DIM]  != CONSTANT);
    const bool pointVaries = (variationTypes[POINT_DIM] != CONSTANT);
    
    auto det2 = KOKKOS_LAMBDA (const PointScalar &a, const PointScalar &b, const PointScalar &c, const PointScalar &d) -> PointScalar
    {
      return a * d - b * c;
    };
    
    auto det3 = KOKKOS_LAMBDA (const PointScalar &a, const PointScalar &b, const PointScalar &c,
                               const PointScalar &d, const PointScalar &e, const PointScalar &f,
                               const PointScalar &g, const PointScalar &h, const PointScalar &i) -> PointScalar
    {
      return a * det2(e,f,h,i) - b * det2(d,f,g,i) + c * det2(d,e,g,h);
    };
    
    auto det4 = KOKKOS_LAMBDA (const PointScalar &a, const PointScalar &b, const PointScalar &c, const PointScalar &d,
                               const PointScalar &e, const PointScalar &f, const PointScalar &g, const PointScalar &h,
                               const PointScalar &i, const PointScalar &j, const PointScalar &k, const PointScalar &l,
                               const PointScalar &m, const PointScalar &n, const PointScalar &o, const PointScalar &p) -> PointScalar
    {
      return a * det3(f,g,h,j,k,l,n,o,p) - b * det3(e,g,h,i,k,l,m,o,p) + c * det3(e,f,h,i,j,l,m,n,p) - d * det3(e,f,g,i,j,k,m,n,o);
    };
    
    if (variationTypes[D1_DIM] == BLOCK_PLUS_DIAGONAL)
    {
      if (cellVaries && pointVaries)
      {
        auto data    = jacobian.getUnderlyingView3();
        auto detData = jacobianDet.getUnderlyingView2();
        
        Kokkos::parallel_for(
        Kokkos::MDRangePolicy<ExecSpaceType,Kokkos::Rank<2>>({0,0},{data.extent_int(0),data.extent_int(1)}),
        KOKKOS_LAMBDA (int cellOrdinal, int pointOrdinal) {
          const int blockWidth   = jacobian.blockPlusDiagonalLastNonDiagonal() + 1;
          const int numDiagonals = data.extent_int(2) - blockWidth * blockWidth;
          const int spaceDim     = blockWidth + numDiagonals;
          
          PointScalar det = 1.0;
          switch (blockWidth)
          {
            case 0:
              det = 1.0;
              break;
            case 1:
            {
              det = data(cellOrdinal,pointOrdinal,0);
              break;
            }
            case 2:
            {
              const auto & a = data(cellOrdinal,pointOrdinal,0);
              const auto & b = data(cellOrdinal,pointOrdinal,1);
              const auto & c = data(cellOrdinal,pointOrdinal,2);
              const auto & d = data(cellOrdinal,pointOrdinal,3);
              det = det2(a,b,c,d);
              
              break;
            }
            case 3:
            {
              const auto & a = data(cellOrdinal,pointOrdinal,0);
              const auto & b = data(cellOrdinal,pointOrdinal,1);
              const auto & c = data(cellOrdinal,pointOrdinal,2);
              const auto & d = data(cellOrdinal,pointOrdinal,3);
              const auto & e = data(cellOrdinal,pointOrdinal,4);
              const auto & f = data(cellOrdinal,pointOrdinal,5);
              const auto & g = data(cellOrdinal,pointOrdinal,6);
              const auto & h = data(cellOrdinal,pointOrdinal,7);
              const auto & i = data(cellOrdinal,pointOrdinal,8);
              det = det3(a,b,c,d,e,f,g,h,i);
              
              break;
            }
            case 4:
            {
              const auto & a = data(cellOrdinal,pointOrdinal,0);
              const auto & b = data(cellOrdinal,pointOrdinal,1);
              const auto & c = data(cellOrdinal,pointOrdinal,2);
              const auto & d = data(cellOrdinal,pointOrdinal,3);
              const auto & e = data(cellOrdinal,pointOrdinal,4);
              const auto & f = data(cellOrdinal,pointOrdinal,5);
              const auto & g = data(cellOrdinal,pointOrdinal,6);
              const auto & h = data(cellOrdinal,pointOrdinal,7);
              const auto & i = data(cellOrdinal,pointOrdinal,8);
              const auto & j = data(cellOrdinal,pointOrdinal,9);
              const auto & k = data(cellOrdinal,pointOrdinal,10);
              const auto & l = data(cellOrdinal,pointOrdinal,11);
              const auto & m = data(cellOrdinal,pointOrdinal,12);
              const auto & n = data(cellOrdinal,pointOrdinal,13);
              const auto & o = data(cellOrdinal,pointOrdinal,14);
              const auto & p = data(cellOrdinal,pointOrdinal,15);
              det = det4(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p);
              
              break;
            }
            default:
              INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "blocks with block width > 4 not supported for BLOCK_PLUS_DIAGONAL");
          }
          // incorporate the remaining, diagonal entries
          const int offset = blockWidth * blockWidth;
          
          for (int d=blockWidth; d<spaceDim; d++)
          {
            const int index = d-blockWidth+offset;
            det *= data(cellOrdinal,pointOrdinal,index);
          }
          detData(cellOrdinal,pointOrdinal) = det;
        });
      }
      else if (cellVaries || pointVaries) // exactly one of cell,point vary -- whichever it is will be in rank 0 of data, invData
      {
        auto data    = jacobian.getUnderlyingView2();
        auto detData = jacobianDet.getUnderlyingView1();
        
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,data.extent_int(0)),
        KOKKOS_LAMBDA (const int &cellPointOrdinal) {
          const int blockWidth   = jacobian.blockPlusDiagonalLastNonDiagonal() + 1;
          const int numDiagonals = data.extent_int(2) - blockWidth * blockWidth;
          const int spaceDim     = blockWidth + numDiagonals;
          
          PointScalar det = 1.0;
          switch (blockWidth)
          {
            case 0:
              det = 1.0;
              break;
            case 1:
            {
              det = data(cellPointOrdinal,0);
              break;
            }
            case 2:
            {
              const auto & a = data(cellPointOrdinal,0);
              const auto & b = data(cellPointOrdinal,1);
              const auto & c = data(cellPointOrdinal,2);
              const auto & d = data(cellPointOrdinal,3);
              det = det2(a,b,c,d);
              
              break;
            }
            case 3:
            {
              const auto & a = data(cellPointOrdinal,0);
              const auto & b = data(cellPointOrdinal,1);
              const auto & c = data(cellPointOrdinal,2);
              const auto & d = data(cellPointOrdinal,3);
              const auto & e = data(cellPointOrdinal,4);
              const auto & f = data(cellPointOrdinal,5);
              const auto & g = data(cellPointOrdinal,6);
              const auto & h = data(cellPointOrdinal,7);
              const auto & i = data(cellPointOrdinal,8);
              det = det3(a,b,c,d,e,f,g,h,i);
              
              break;
            }
            default:
              INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "blocks with block width > 3 not supported for BLOCK_PLUS_DIAGONAL");
          }
          // incorporate the remaining, diagonal entries
          const int offset = blockWidth * blockWidth;
          
          for (int d=blockWidth; d<spaceDim; d++)
          {
            const int index = d-blockWidth+offset;
            det *= data(cellPointOrdinal,index);
          }
          detData(cellPointOrdinal) = det;
        });
      }
      else // neither cell nor point varies
      {
        auto data    = jacobian.getUnderlyingView1();
        auto detData = jacobianDet.getUnderlyingView1();
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,1),
        KOKKOS_LAMBDA (const int &dummyArg) {
          const int blockWidth   = jacobian.blockPlusDiagonalLastNonDiagonal() + 1;
          const int numDiagonals = jacobian.extent_int(2) - blockWidth * blockWidth;
          const int spaceDim     = blockWidth + numDiagonals;
                  
          PointScalar det = 1.0;
          switch (blockWidth)
          {
            case 0:
              det = 1.0;
              break;
            case 1:
            {
              det = data(0);
              break;
            }
            case 2:
            {
              const auto & a = data(0);
              const auto & b = data(1);
              const auto & c = data(2);
              const auto & d = data(3);
              det = det2(a,b,c,d);
              
              break;
            }
            case 3:
            {
              const auto & a = data(0);
              const auto & b = data(1);
              const auto & c = data(2);
              const auto & d = data(3);
              const auto & e = data(4);
              const auto & f = data(5);
              const auto & g = data(6);
              const auto & h = data(7);
              const auto & i = data(8);
              det = det3(a,b,c,d,e,f,g,h,i);
              
              break;
            }
            default:
              INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "blocks with block width > 3 not supported for BLOCK_PLUS_DIAGONAL");
          }
          // incorporate the remaining, diagonal entries
          const int offset = blockWidth * blockWidth;
          
          for (int d=blockWidth; d<spaceDim; d++)
          {
            const int index = d-blockWidth+offset;
            det *= data(index);
          }
          detData(0) = det;
        });
      }
    }
    else if ( jacobian.getUnderlyingViewRank() == 4 )
    {
      auto data    = jacobian.getUnderlyingView4();
      auto detData = jacobianDet.getUnderlyingView2();
      Intrepid2::RealSpaceTools<ExecSpaceType>::det(detData, data);
    }
    else if ( jacobian.getUnderlyingViewRank() == 3 )
    {
      auto data    = jacobian.getUnderlyingView3();
      auto detData = jacobianDet.getUnderlyingView1();
      Intrepid2::RealSpaceTools<ExecSpaceType>::det(detData, data);
    }
    else if ( jacobian.getUnderlyingViewRank() == 2 )
    {
      auto data    = jacobian.getUnderlyingView2();
      auto detData = jacobianDet.getUnderlyingView1();
      Kokkos::parallel_for("fill jacobian det", Kokkos::RangePolicy<ExecSpaceType>(0,1), KOKKOS_LAMBDA(const int &i)
      {
        detData(0) = Intrepid2::Kernels::Serial::determinant(data);
      });
    }
    else
    {
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "jacobian's underlying view must have rank 2,3, or 4");
    }
  }

  template<typename DeviceType>
  template<class PointScalar>
  void CellTools<DeviceType>::setJacobianDetInv( Data<PointScalar,DeviceType> &jacobianDetInv, const Data<PointScalar,DeviceType> & jacobian )
  {
    setJacobianDet(jacobianDetInv, jacobian);
    
    auto unitData = jacobianDetInv.allocateConstantData(1.0);
    jacobianDetInv.storeInPlaceQuotient(unitData, jacobianDetInv);
  }

  template<typename DeviceType>
  template<class PointScalar>
  void CellTools<DeviceType>::setJacobianInv( Data<PointScalar,DeviceType> &jacobianInv, const Data<PointScalar,DeviceType> & jacobian )
  {
    auto variationTypes  = jacobian.getVariationTypes();
    
    const int CELL_DIM  = 0;
    const int POINT_DIM = 1;
    const int D1_DIM    = 2;
    
    auto det2 = KOKKOS_LAMBDA (const PointScalar &a, const PointScalar &b, const PointScalar &c, const PointScalar &d) -> PointScalar
    {
      return a * d - b * c;
    };
    
    auto det3 = KOKKOS_LAMBDA (const PointScalar &a, const PointScalar &b, const PointScalar &c,
                               const PointScalar &d, const PointScalar &e, const PointScalar &f,
                               const PointScalar &g, const PointScalar &h, const PointScalar &i) -> PointScalar
    {
      return a * det2(e,f,h,i) - b * det2(d,f,g,i) + c * det2(d,e,g,h);
    };
    
    if (variationTypes[D1_DIM] == BLOCK_PLUS_DIAGONAL)
    {
      const bool cellVaries  = variationTypes[CELL_DIM]  != CONSTANT;
      const bool pointVaries = variationTypes[POINT_DIM] != CONSTANT;
      if (cellVaries && pointVaries)
      {
        auto data    = jacobian.getUnderlyingView3();
        auto invData = jacobianInv.getUnderlyingView3();
        
        Kokkos::parallel_for(
        Kokkos::MDRangePolicy<ExecSpaceType,Kokkos::Rank<2>>({0,0},{data.extent_int(0),data.extent_int(1)}),
        KOKKOS_LAMBDA (int cellOrdinal, int pointOrdinal) {
          const int blockWidth   = jacobian.blockPlusDiagonalLastNonDiagonal() + 1;
          const int numDiagonals = data.extent_int(2) - blockWidth * blockWidth;
          const int spaceDim     = blockWidth + numDiagonals;
          
          switch (blockWidth)
          {
            case 0:
              break;
            case 1:
            {
              invData(cellOrdinal,pointOrdinal,0) = 1.0 / data(cellOrdinal,pointOrdinal,0);
              break;
            }
            case 2:
            {
              const auto & a = data(cellOrdinal,pointOrdinal,0);
              const auto & b = data(cellOrdinal,pointOrdinal,1);
              const auto & c = data(cellOrdinal,pointOrdinal,2);
              const auto & d = data(cellOrdinal,pointOrdinal,3);
              const PointScalar det = det2(a,b,c,d);
              
              invData(cellOrdinal,pointOrdinal,0) =   d/det;
              invData(cellOrdinal,pointOrdinal,1) = - b/det;
              invData(cellOrdinal,pointOrdinal,2) = - c/det;
              invData(cellOrdinal,pointOrdinal,3) =   a/det;
              break;
            }
            case 3:
            {
              const auto & a = data(cellOrdinal,pointOrdinal,0);
              const auto & b = data(cellOrdinal,pointOrdinal,1);
              const auto & c = data(cellOrdinal,pointOrdinal,2);
              const auto & d = data(cellOrdinal,pointOrdinal,3);
              const auto & e = data(cellOrdinal,pointOrdinal,4);
              const auto & f = data(cellOrdinal,pointOrdinal,5);
              const auto & g = data(cellOrdinal,pointOrdinal,6);
              const auto & h = data(cellOrdinal,pointOrdinal,7);
              const auto & i = data(cellOrdinal,pointOrdinal,8);
              const PointScalar det = det3(a,b,c,d,e,f,g,h,i);
              
              {
                const auto val0 =   e*i - h*f;
                const auto val1 = - d*i + g*f;
                const auto val2 =   d*h - g*e;
                
                invData(cellOrdinal,pointOrdinal,0) = val0/det;
                invData(cellOrdinal,pointOrdinal,1) = val1/det;
                invData(cellOrdinal,pointOrdinal,2) = val2/det;
              }
              {
                const auto val0 =   h*c - b*i;
                const auto val1 =   a*i - g*c;
                const auto val2 = - a*h + g*b;
                
                invData(cellOrdinal,pointOrdinal,3) = val0/det;
                invData(cellOrdinal,pointOrdinal,4) = val1/det;
                invData(cellOrdinal,pointOrdinal,5) = val2/det;
              }
              {
                const auto val0 =   b*f - e*c;
                const auto val1 = - a*f + d*c;
                const auto val2 =   a*e - d*b;

                invData(cellOrdinal,pointOrdinal,6) = val0/det;
                invData(cellOrdinal,pointOrdinal,7) = val1/det;
                invData(cellOrdinal,pointOrdinal,8) = val2/det;
              }
              break;
            }
            default:
              INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "blocks with block width > 3 not supported for BLOCK_PLUS_DIAGONAL");
          }
          // handle the remaining, diagonal entries
          const int offset = blockWidth * blockWidth;
          
          for (int d=blockWidth; d<spaceDim; d++)
          {
            const int index = d-blockWidth+offset;
            invData(cellOrdinal,pointOrdinal,index) = 1.0 / data(cellOrdinal,pointOrdinal,index);
          }
        });
      }
      else if (cellVaries || pointVaries) // exactly one of cell,point vary -- whichever it is will be in rank 0 of data, invData
      {
        auto data    = jacobian.getUnderlyingView2();
        auto invData = jacobianInv.getUnderlyingView2();
        
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,data.extent_int(0)),
        KOKKOS_LAMBDA (const int &cellPointOrdinal) {
          const int blockWidth   = jacobian.blockPlusDiagonalLastNonDiagonal() + 1;
          const int numDiagonals = data.extent_int(1) - blockWidth * blockWidth;
          const int spaceDim     = blockWidth + numDiagonals;
          
          switch (blockWidth)
          {
            case 0:
              break;
            case 1:
            {
              invData(cellPointOrdinal,0) = 1.0 / data(cellPointOrdinal,0);
              break;
            }
            case 2:
            {
              const auto & a = data(cellPointOrdinal,0);
              const auto & b = data(cellPointOrdinal,1);
              const auto & c = data(cellPointOrdinal,2);
              const auto & d = data(cellPointOrdinal,3);
              const PointScalar det = det2(a,b,c,d);
              
              invData(cellPointOrdinal,0) =   d/det;
              invData(cellPointOrdinal,1) = - b/det;
              invData(cellPointOrdinal,2) = - c/det;
              invData(cellPointOrdinal,3) =   a/det;
              break;
            }
            case 3:
            {
              const auto & a = data(cellPointOrdinal,0);
              const auto & b = data(cellPointOrdinal,1);
              const auto & c = data(cellPointOrdinal,2);
              const auto & d = data(cellPointOrdinal,3);
              const auto & e = data(cellPointOrdinal,4);
              const auto & f = data(cellPointOrdinal,5);
              const auto & g = data(cellPointOrdinal,6);
              const auto & h = data(cellPointOrdinal,7);
              const auto & i = data(cellPointOrdinal,8);
              const PointScalar det = det3(a,b,c,d,e,f,g,h,i);
              
              {
                const auto val0 =   e*i - h*f;
                const auto val1 = - d*i + g*f;
                const auto val2 =   d*h - g*e;
                
                invData(cellPointOrdinal,0) = val0/det;
                invData(cellPointOrdinal,1) = val1/det;
                invData(cellPointOrdinal,2) = val2/det;
              }
              {
                const auto val0 =   h*c - b*i;
                const auto val1 =   a*i - g*c;
                const auto val2 = - a*h + g*b;
                
                invData(cellPointOrdinal,3) = val0/det;
                invData(cellPointOrdinal,4) = val1/det;
                invData(cellPointOrdinal,5) = val2/det;
              }
              {
                const auto val0 =   b*f - e*c;
                const auto val1 = - a*f + d*c;
                const auto val2 =   a*e - d*b;

                invData(cellPointOrdinal,6) = val0/det;
                invData(cellPointOrdinal,7) = val1/det;
                invData(cellPointOrdinal,8) = val2/det;
              }
              break;
            }
            default:
              INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "blocks with block width > 3 not supported for BLOCK_PLUS_DIAGONAL in setJacobianInv()");
          }
          // handle the remaining, diagonal entries
          const int offset = blockWidth * blockWidth;
          
          for (int d=blockWidth; d<spaceDim; d++)
          {
            const int index = d-blockWidth+offset;
            invData(cellPointOrdinal,index) = 1.0 / data(cellPointOrdinal,index);
          }
        });
      }
      else // neither cell nor point varies
      {
        auto data    = jacobian.getUnderlyingView1();
        auto invData = jacobianInv.getUnderlyingView1();
        
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,1),
        KOKKOS_LAMBDA (const int &dummyArg) {
          const int blockWidth   = jacobian.blockPlusDiagonalLastNonDiagonal() + 1;
          const int numDiagonals = data.extent_int(0) - blockWidth * blockWidth;
          const int spaceDim     = blockWidth + numDiagonals;
          
          switch (blockWidth)
          {
            case 0:
              break;
            case 1:
            {
              invData(0) = 1.0 / data(0);
              break;
            }
            case 2:
            {
              const auto & a = data(0);
              const auto & b = data(1);
              const auto & c = data(2);
              const auto & d = data(3);
              const PointScalar det = det2(a,b,c,d);
              
              invData(0) =   d/det;
              invData(1) = - b/det;
              invData(2) = - c/det;
              invData(3) =   a/det;
              break;
            }
            case 3:
            {
              const auto & a = data(0);
              const auto & b = data(1);
              const auto & c = data(2);
              const auto & d = data(3);
              const auto & e = data(4);
              const auto & f = data(5);
              const auto & g = data(6);
              const auto & h = data(7);
              const auto & i = data(8);
              const PointScalar det = det3(a,b,c,d,e,f,g,h,i);
              
              {
                const auto val0 =   e*i - h*f;
                const auto val1 = - d*i + g*f;
                const auto val2 =   d*h - g*e;
                
                invData(0) = val0/det;
                invData(1) = val1/det;
                invData(2) = val2/det;
              }
              {
                const auto val0 =   h*c - b*i;
                const auto val1 =   a*i - g*c;
                const auto val2 = - a*h + g*b;
                
                invData(3) = val0/det;
                invData(4) = val1/det;
                invData(5) = val2/det;
              }
              {
                const auto val0 =   b*f - e*c;
                const auto val1 = - a*f + d*c;
                const auto val2 =   a*e - d*b;

                invData(6) = val0/det;
                invData(7) = val1/det;
                invData(8) = val2/det;
              }
              break;
            }
            default:
              INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "blocks with block width > 3 not supported for BLOCK_PLUS_DIAGONAL in setJacobianInv()");
          }
          // handle the remaining, diagonal entries
          const int offset = blockWidth * blockWidth;
          
          for (int d=blockWidth; d<spaceDim; d++)
          {
            const int index = d-blockWidth+offset;
            invData(index) = 1.0 / data(index);
          }
        });
      }
    }
    else if ( jacobian.getUnderlyingViewRank() == 4 )
    {
      auto data    = jacobian.getUnderlyingView4();
      auto invData = jacobianInv.getUnderlyingView4();
      
      Intrepid2::RealSpaceTools<ExecSpaceType>::inverse(invData, data);
    }
    else if ( jacobian.getUnderlyingViewRank() == 3 )
    {
      auto data    = jacobian.getUnderlyingView3();
      auto invData = jacobianInv.getUnderlyingView3();
      
      Intrepid2::RealSpaceTools<ExecSpaceType>::inverse(invData, data);
    }
    else if ( jacobian.getUnderlyingViewRank() == 2 ) // Cell, point constant; D1, D2 vary normally
    {
      auto data    = jacobian.getUnderlyingView2();
      auto invData = jacobianInv.getUnderlyingView2();
      
      Kokkos::parallel_for("fill jacobian inverse", Kokkos::RangePolicy<ExecSpaceType>(0,1), KOKKOS_LAMBDA(const int &i)
      {
        Intrepid2::Kernels::Serial::inverse(invData,data);
      });
    }
    else
    {
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "jacobian's underlying view must have rank 2,3, or 4");
    }
  }

  template<typename DeviceType>
  template<typename JacobianViewType,
           typename BasisGradientsType,
           typename WorksetType>
  void
  CellTools<DeviceType>::
  setJacobian(       JacobianViewType jacobian,
               const WorksetType      worksetCell,
               const BasisGradientsType gradients, const int startCell, const int endCell)
  {
    constexpr bool is_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(jacobian)::memory_space>::accessible;
    static_assert(is_accessible, "CellTools<DeviceType>::setJacobian(..): output view's memory space is not compatible with DeviceType");

    using FunctorType      = FunctorCellTools::F_setJacobian<JacobianViewType,WorksetType,BasisGradientsType> ;
    
    // resolve the -1 default argument for endCell into the true end cell index
    int endCellResolved = (endCell == -1) ? worksetCell.extent_int(0) : endCell;
    
    using range_policy_type = Kokkos::MDRangePolicy
      < ExecSpaceType, Kokkos::Rank<2>, Kokkos::IndexType<ordinal_type> >;
    range_policy_type policy( { 0, 0 },
                              { jacobian.extent(0), jacobian.extent(1) } );
    Kokkos::parallel_for( policy, FunctorType(jacobian, worksetCell, gradients, startCell, endCellResolved) );
  }

  template<typename DeviceType>
  template<typename JacobianViewType,
           typename PointViewType,
           typename WorksetType,
           typename HGradBasisType>
  void
  CellTools<DeviceType>::
  setJacobian(       JacobianViewType             jacobian,
               const PointViewType                points,
               const WorksetType                  worksetCell,
               const Teuchos::RCP<HGradBasisType> basis,
               const int startCell, const int endCell) {
    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(jacobian)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(points)::memory_space>::accessible;
    static_assert(are_accessible, "CellTools<DeviceType>::setJacobian(..): input/output views' memory spaces are not compatible with DeviceType");

#ifdef HAVE_INTREPID2_DEBUG    
    CellTools_setJacobianArgs(jacobian, points, worksetCell, basis->getBaseCellTopology(), startCell, endCell);
    //static_assert(std::is_same( pointValueType, decltype(basis->getDummyOutputValue()) ));
#endif
    const auto cellTopo = basis->getBaseCellTopology();
    const ordinal_type spaceDim = cellTopo.getDimension();
    const ordinal_type numCells = jacobian.extent(0);
    
    //points can be rank-2 (P,D), or rank-3 (C,P,D)
    const ordinal_type pointRank = points.rank();
    const ordinal_type numPoints = (pointRank == 2 ? points.extent(0) : points.extent(1));
    const ordinal_type basisCardinality = basis->getCardinality();
    
    // the following does not work for RCP; its * operator returns reference not the object
    //typedef typename decltype(*basis)::output_value_type gradValueType;
    //typedef Kokkos::DynRankView<decltype(basis->getDummyOutputValue()),DeviceType> gradViewType;

    auto vcprop = Kokkos::common_view_alloc_prop(points);
    using GradViewType = Kokkos::DynRankView<typename decltype(vcprop)::value_type,DeviceType>;

    GradViewType grads;

    switch (pointRank) {
    case 2: {
      // For most FEMs
      grads = GradViewType(Kokkos::view_alloc("CellTools::setJacobian::grads", vcprop),basisCardinality, numPoints, spaceDim);
      basis->getValues(grads, 
                       points, 
                       OPERATOR_GRAD);
      break;
    }
    case 3: { 
      // For CVFEM
      grads = GradViewType(Kokkos::view_alloc("CellTools::setJacobian::grads", vcprop), numCells, basisCardinality, numPoints, spaceDim);
      for (ordinal_type cell=0;cell<numCells;++cell) 
        basis->getValues(Kokkos::subview( grads,  cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() ),  
                         Kokkos::subview( points, cell, Kokkos::ALL(), Kokkos::ALL() ),  
                         OPERATOR_GRAD);
      break;
    }
    }
    
    setJacobian(jacobian, worksetCell, grads, startCell, endCell);
  }

  template<typename DeviceType>
  template<typename JacobianInvViewType,                                   
           typename JacobianViewType>                                      
  void                                                                                               
  CellTools<DeviceType>::
  setJacobianInv(       JacobianInvViewType jacobianInv,     
                  const JacobianViewType    jacobian ) {
#ifdef HAVE_INTREPID2_DEBUG
    CellTools_setJacobianInvArgs(jacobianInv, jacobian);
#endif
    RealSpaceTools<DeviceType>::inverse(jacobianInv, jacobian);
  }
  
  template<typename DeviceType>
  template<typename JacobianDetViewType,                                   
           typename JacobianViewType>                                      
  void                                                                                               
  CellTools<DeviceType>::
  setJacobianDet(       JacobianDetViewType jacobianDet,    
                  const JacobianViewType    jacobian ) {
#ifdef HAVE_INTREPID2_DEBUG
    CellTools_setJacobianDetArgs(jacobianDet, jacobian);
#endif
    RealSpaceTools<DeviceType>::det(jacobianDet, jacobian);
  }

  template<typename DeviceType>
  template<typename Scalar>
  void
  CellTools<DeviceType>::
  setJacobianDividedByDet( Data<Scalar,DeviceType> & jacobianDividedByDet,
                          const Data<Scalar,DeviceType> & jacobian,
                          const Data<Scalar,DeviceType> & jacobianDetInv)
  {
    DataTools::multiplyByCPWeights(jacobianDividedByDet,jacobian,jacobianDetInv);
  }
}

#endif
