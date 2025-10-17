// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Tests against structured integration facilities - "synthetic" test cases (i.e., no geometry specified).
    \author Nathan V. Roberts
*/

#include "Kokkos_Core.hpp"

#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_UnitTestRepository.hpp>

#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_IntegrationTools.hpp>
#include <Intrepid2_TestUtils.hpp>

#include "Intrepid2_Data.hpp"
#include "Intrepid2_TensorData.hpp"
#include "Intrepid2_TensorPoints.hpp"
#include "Intrepid2_TransformedBasisValues.hpp"
#include "Intrepid2_VectorData.hpp"

#include "Intrepid2_ScalarView.hpp"

namespace
{
 using namespace Intrepid2;

template<typename DeviceType,class Scalar>
void IT_integrate(Data<Scalar,DeviceType> integrals, const TransformedBasisValues<Scalar,DeviceType> & basisValuesLeft,
                  const TensorData<Scalar,DeviceType> & cellMeasures,
                  const TransformedBasisValues<Scalar,DeviceType> & basisValuesRight, const bool sumInto = false,
                  double* approximateFlops = NULL)
{
  using ExecutionSpace = typename DeviceType::execution_space;

  // ordinal filter is used for Serendipity basis; we don't yet support Serendipity for sum factorization.
  // (when we do, the strategy will likely be to do sum factorized assembly and then filter the result).
  const bool  leftHasOrdinalFilter =  basisValuesLeft.basisValues().ordinalFilter().extent_int(0) > 0;
  const bool rightHasOrdinalFilter = basisValuesRight.basisValues().ordinalFilter().extent_int(0) > 0;
  TEUCHOS_TEST_FOR_EXCEPTION(leftHasOrdinalFilter || rightHasOrdinalFilter, std::invalid_argument, "Ordinal filters for BasisValues are not yet supported by IntegrationTools");
  
  if (approximateFlops != NULL)
  {
    *approximateFlops = 0;
  }
    
  // integral data may have shape (C,F1,F2) or (if the variation type is CONSTANT in the cell dimension) shape (F1,F2)
  const int integralViewRank = integrals.getUnderlyingViewRank();
  
  if (!sumInto)
  {
    integrals.clear();
  }
  
  const int cellDataExtent = integrals.getDataExtent(0);
  const int numFieldsLeft  = basisValuesLeft.numFields();
  const int numFieldsRight = basisValuesRight.numFields();
  const int spaceDim       = basisValuesLeft.spaceDim();
  
  INTREPID2_TEST_FOR_EXCEPTION(basisValuesLeft.spaceDim() != basisValuesRight.spaceDim(), std::invalid_argument, "vectorDataLeft and vectorDataRight must agree on the space dimension");
  
  const int leftFamilyCount  =  basisValuesLeft.basisValues().numFamilies();
  const int rightFamilyCount = basisValuesRight.basisValues().numFamilies();
  
  // we require that the number of tensor components in the vectors are the same for each vector entry
  // this is not strictly necessary, but it makes implementation easier, and we don't at present anticipate other use cases
  int numTensorComponentsLeft = -1;
  const bool leftIsVectorValued = basisValuesLeft.vectorData().isValid();
  
  if (leftIsVectorValued)
  {
    const auto &refVectorLeft   = basisValuesLeft.vectorData();
    int numFamiliesLeft         = refVectorLeft.numFamilies();
    int numVectorComponentsLeft = refVectorLeft.numComponents();
    Kokkos::Array<int,7> maxFieldsForComponentLeft  {0,0,0,0,0,0,0};
    for (int familyOrdinal=0; familyOrdinal<numFamiliesLeft; familyOrdinal++)
    {
      for (int vectorComponent=0; vectorComponent<numVectorComponentsLeft; vectorComponent++)
      {
        const TensorData<Scalar,DeviceType> &tensorData = refVectorLeft.getComponent(familyOrdinal,vectorComponent);
        if (tensorData.numTensorComponents() > 0)
        {
          if (numTensorComponentsLeft == -1)
          {
            numTensorComponentsLeft = tensorData.numTensorComponents();
          }
          INTREPID2_TEST_FOR_EXCEPTION(numVectorComponentsLeft != tensorData.numTensorComponents(), std::invalid_argument, "Each valid entry in vectorDataLeft must have the same number of tensor components as every other");
          for (int r=0; r<numTensorComponentsLeft; r++)
          {
            maxFieldsForComponentLeft[r] = std::max(tensorData.getTensorComponent(r).extent_int(0), maxFieldsForComponentLeft[r]);
          }
        }
      }
    }
  }
  else
  {
    numTensorComponentsLeft = basisValuesLeft.basisValues().tensorData(0).numTensorComponents(); // family ordinal 0
    for (int familyOrdinal = 0; familyOrdinal < leftFamilyCount; familyOrdinal++)
    {
      INTREPID2_TEST_FOR_EXCEPTION(basisValuesLeft.basisValues().tensorData(familyOrdinal).numTensorComponents() != numTensorComponentsLeft, std::invalid_argument, "All families must match in the number of tensor components");
    }
  }
  int numTensorComponentsRight = -1;
  const bool rightIsVectorValued = basisValuesRight.vectorData().isValid();
  
  if (rightIsVectorValued)
  {
    const auto &refVectorRight   = basisValuesRight.vectorData();
    int numFamiliesRight         = refVectorRight.numFamilies();
    int numVectorComponentsRight = refVectorRight.numComponents();
    Kokkos::Array<int,7> maxFieldsForComponentRight {0,0,0,0,0,0,0};
    for (int familyOrdinal=0; familyOrdinal<numFamiliesRight; familyOrdinal++)
    {
      for (int vectorComponent=0; vectorComponent<numVectorComponentsRight; vectorComponent++)
      {
        const auto &tensorData = refVectorRight.getComponent(familyOrdinal,vectorComponent);
        if (tensorData.numTensorComponents() > 0)
        {
          if (numTensorComponentsRight == -1)
          {
            numTensorComponentsRight = tensorData.numTensorComponents();
          }
          INTREPID2_TEST_FOR_EXCEPTION(numVectorComponentsRight != tensorData.numTensorComponents(), std::invalid_argument, "Each valid entry in vectorDataRight must have the same number of tensor components as every other");
          for (int r=0; r<numTensorComponentsRight; r++)
          {
            maxFieldsForComponentRight[r] = std::max(tensorData.getTensorComponent(r).extent_int(0), maxFieldsForComponentRight[r]);
          }
        }
      }
    }
    INTREPID2_TEST_FOR_EXCEPTION(numTensorComponentsRight != numTensorComponentsLeft, std::invalid_argument, "Right families must match left in the number of tensor components");
  }
  else
  {
    // check that right tensor component count agrees with left
    for (int familyOrdinal=0; familyOrdinal< rightFamilyCount; familyOrdinal++)
    {
      INTREPID2_TEST_FOR_EXCEPTION(basisValuesRight.basisValues().tensorData(familyOrdinal).numTensorComponents() != numTensorComponentsLeft, std::invalid_argument, "Right families must match left in the number of tensor components");
    }
  }
  const int numPointTensorComponents = cellMeasures.numTensorComponents() - 1;
    
  if ((numPointTensorComponents == numTensorComponentsLeft) && basisValuesLeft.axisAligned() && basisValuesRight.axisAligned())
  {
    // cellMeasures is a non-trivial tensor product, and the pullbacks are all diagonals.
    
    // in this case, the integrals in each tensorial direction are entirely separable
    // allocate some temporary storage for the integrals in each tensorial direction
    
    // if cellMeasures is a nontrivial tensor product, that means that all cells have the same shape, up to scaling.

    Kokkos::Array<int,Parameters::MaxTensorComponents> pointDimensions;
    for (int r=0; r<numPointTensorComponents; r++)
    {
      // first tensorial component of cell measures is the cell dimension; after that we have (P1,P2,…)
      pointDimensions[r] = cellMeasures.getTensorComponent(r+1).extent_int(0);
    }

    // only one of these will be a valid container:
    Kokkos::View<Scalar**,  DeviceType> integralView2;
    Kokkos::View<Scalar***, DeviceType> integralView3;
    if (integralViewRank == 3)
    {
      integralView3 = integrals.getUnderlyingView3();
    }
    else
    {
      integralView2 = integrals.getUnderlyingView2();
    }
    for (int leftFamilyOrdinal=0; leftFamilyOrdinal<leftFamilyCount; leftFamilyOrdinal++)
    {
      int a_offset = 0; // left vector component offset
      int leftFieldOffset = basisValuesLeft.basisValues().familyFieldOrdinalOffset(leftFamilyOrdinal);
      
      const int leftVectorComponentCount = leftIsVectorValued ? basisValuesLeft.vectorData().numComponents() : 1;
      for (int leftVectorComponentOrdinal = 0; leftVectorComponentOrdinal < leftVectorComponentCount; leftVectorComponentOrdinal++)
      {
        TensorData<Scalar,DeviceType> leftComponent = leftIsVectorValued ? basisValuesLeft.vectorData().getComponent(leftFamilyOrdinal, leftVectorComponentOrdinal)
                                                                         : basisValuesLeft.basisValues().tensorData(leftFamilyOrdinal);
        if (!leftComponent.isValid())
        {
          a_offset++; // empty components are understood to take up one dimension
          continue;
        }
        const int leftDimSpan = leftComponent.extent_int(2);
          
        const int leftComponentFieldCount = leftComponent.extent_int(0);
        
        for (int rightFamilyOrdinal=0; rightFamilyOrdinal<rightFamilyCount; rightFamilyOrdinal++)
        {
          int b_offset = 0; // right vector component offset
          int rightFieldOffset = basisValuesRight.vectorData().familyFieldOrdinalOffset(rightFamilyOrdinal);

          const int rightVectorComponentCount = rightIsVectorValued ? basisValuesRight.vectorData().numComponents() : 1;
          for (int rightVectorComponentOrdinal = 0; rightVectorComponentOrdinal < rightVectorComponentCount; rightVectorComponentOrdinal++)
          {
            TensorData<Scalar,DeviceType> rightComponent = rightIsVectorValued ? basisValuesRight.vectorData().getComponent(rightFamilyOrdinal, rightVectorComponentOrdinal)
                                                                               : basisValuesRight.basisValues().tensorData(rightFamilyOrdinal);
            if (!rightComponent.isValid())
            {
              b_offset++; // empty components are understood to take up one dimension
              continue;
            }
            const int rightDimSpan = rightComponent.extent_int(2);
              
            const int rightComponentFieldCount = rightComponent.extent_int(0);
                    
            // we only accumulate for a == b (since this is the axis-aligned case).  Do the a, b spans intersect?
            if ((a_offset + leftDimSpan <= b_offset) || (b_offset + rightDimSpan <= a_offset)) // no intersection
            {
              b_offset += rightDimSpan;
              
              continue;
            }
            
            // if the a, b spans intersect, by assumption they should align exactly
            INTREPID2_TEST_FOR_EXCEPTION(( a_offset != b_offset), std::logic_error, "left and right dimension offsets should match.");
            INTREPID2_TEST_FOR_EXCEPTION(( leftDimSpan != rightDimSpan), std::invalid_argument, "left and right components must span the same number of dimensions.");
            
            const int d_start = a_offset;
            const int d_end   = d_start + leftDimSpan;
            
            using ComponentIntegralsArray = Kokkos::Array< Kokkos::Array<ScalarView<Scalar,DeviceType>, Parameters::MaxTensorComponents>, Parameters::MaxTensorComponents>;
            ComponentIntegralsArray componentIntegrals;
            for (int r=0; r<numPointTensorComponents; r++)
            {
              /*
               Four vector data cases to consider:
               1. Both vector data containers are filled with axial components - first component in 3D has form (f,0,0), second (0,f,0), third (0,0,f).
               2. Both vector data containers have arbitrary components - in 3D: (f1,f2,f3) where f1 is given by the first component, f2 by the second, f3 by the third.
               3. First container is axial, second arbitrary.
               4. First is arbitrary, second axial.
               
               But note that in all four cases, the structure of the integral is the same: you have three vector component integrals that get summed.  The actual difference between
               the cases does not show up in the reference-space integrals here, but in the accumulation in physical space below, where the tensor field numbering comes into play.
               
               The choice between axial and arbitrary affects the way the fields are numbered; the arbitrary components' indices refer to the same vector function, so they correspond,
               while the axial components refer to distinct scalar functions, so their numbering in the data container is cumulative.
              */
                
              Data<Scalar,DeviceType>  quadratureWeights = cellMeasures.getTensorComponent(r+1);
              const int numPoints = pointDimensions[r];
                
                // It may be worth considering the possibility that some of these components point to the same data -- if so, we could possibly get better data locality by making the corresponding componentIntegral entries point to the same location as well.  (And we could avoid some computations here.)
        
              Data<Scalar,DeviceType>  leftTensorComponent =  leftComponent.getTensorComponent(r);
              Data<Scalar,DeviceType> rightTensorComponent = rightComponent.getTensorComponent(r);
              
              const int  leftTensorComponentDimSpan =  leftTensorComponent.extent_int(2);
              const int  leftTensorComponentFields  =  leftTensorComponent.extent_int(0);
              const int rightTensorComponentDimSpan = rightTensorComponent.extent_int(2);
              const int rightTensorComponentFields  = rightTensorComponent.extent_int(0);
              
              INTREPID2_TEST_FOR_EXCEPTION(( leftTensorComponentDimSpan != rightTensorComponentDimSpan), std::invalid_argument, "left and right components must span the same number of dimensions.");
            
              for (int d=d_start; d<d_end; d++)
              {
                ScalarView<Scalar,DeviceType> componentIntegralView;
                
                const bool allocateFadStorage = !(std::is_standard_layout<Scalar>::value && std::is_trivial<Scalar>::value);
                if (allocateFadStorage)
                {
                  auto fad_size_output = dimension_scalar(integrals.getUnderlyingView());
                  componentIntegralView = ScalarView<Scalar,DeviceType>("componentIntegrals for tensor component " + std::to_string(r) + ", in dimension " + std::to_string(d), leftTensorComponentFields, rightTensorComponentFields, fad_size_output);
                }
                else
                {
                  componentIntegralView = ScalarView<Scalar,DeviceType>("componentIntegrals for tensor component " + std::to_string(r) + ", in dimension " + std::to_string(d), leftTensorComponentFields, rightTensorComponentFields);
                }
            
                auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{leftTensorComponentFields,rightTensorComponentFields});
                
                Impl::F_RefSpaceIntegral<Scalar, DeviceType> refSpaceIntegralFunctor(componentIntegralView, leftTensorComponent, rightTensorComponent, quadratureWeights,
                                                                                     leftTensorComponentDimSpan);
                Kokkos::parallel_for("compute componentIntegrals", policy, refSpaceIntegralFunctor);
                
                componentIntegrals[r][d] = componentIntegralView;
                
                if (approximateFlops != NULL)
                {
                  *approximateFlops += leftTensorComponentFields*rightTensorComponentFields*numPoints*(3); // two multiplies, one add in innermost loop
                }
              } // d
            } // r
            
            ExecutionSpace().fence();
            Impl::F_ComputeIntegral<Scalar,DeviceType> computeIntegralFunctor(integralView2, integralView3, leftComponent, rightComponent, cellMeasures, basisValuesLeft, basisValuesRight, componentIntegrals, d_start, d_end, numPointTensorComponents, leftFieldOffset, rightFieldOffset, integralViewRank);
            Kokkos::Array<int,3> upperBounds {cellDataExtent,leftComponentFieldCount,rightComponentFieldCount}; // separately declared in effort to get around Intel 17.0.1 compiler weirdness.
            Kokkos::Array<int,3> lowerBounds {0,0,0};
            auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>(lowerBounds, upperBounds);
            // TODO: note that for best performance, especially with Fad types, we should replace this parallel for with a Functor and use hierarchical parallelism
            Kokkos::parallel_for("compute field integrals", policy, computeIntegralFunctor);
            b_offset += rightDimSpan;
          } // rightVectorComponentOrdinal
        } // rightFamilyOrdinal
        a_offset += leftDimSpan;
      } // leftVectorComponentOrdinal
    } // leftFamilyOrdinal
    
    if (approximateFlops != NULL)
    {
      // TODO: check the accuracy of this
      *approximateFlops += (2 + spaceDim * (3 + numPointTensorComponents)) * cellDataExtent * numFieldsLeft * numFieldsRight;
    }
  }
  else // general case (not axis-aligned + affine tensor-product structure)
  {
    // prepare composed transformation matrices
    const Data<Scalar,DeviceType> & leftTransform  = basisValuesLeft.transform();
    const Data<Scalar,DeviceType> & rightTransform = basisValuesRight.transform();
    const bool transposeLeft  = true;
    const bool transposeRight = false;
//    auto timer = Teuchos::TimeMonitor::getNewTimer("mat-mat");
//    timer->start();
    // transforms can be matrices -- (C,P,D,D): rank 4 -- or scalar weights -- (C,P): rank 2 -- or vector weights -- (C,P,D): rank 3
    Data<Scalar,DeviceType> composedTransform;
    // invalid/empty transforms are used when the identity is intended.
    const int leftRank  = leftTransform.rank();
    const int rightRank = rightTransform.rank();
    
    if (leftTransform.isValid() && rightTransform.isValid())
    {
      const bool bothRank4 = (leftRank == 4) && (rightRank == 4);
      const bool bothRank3 = (leftRank == 3) && (rightRank == 3);
      const bool bothRank2 = (leftRank == 2) && (rightRank == 2);
      const bool ranks32   = ((leftRank == 3) && (rightRank == 2)) || ((leftRank == 2) && (rightRank == 3));
      const bool ranks42   = ((leftRank == 4) && (rightRank == 2)) || ((leftRank == 2) && (rightRank == 4));
      
      if (bothRank4) // (C,P,D,D)
      {
        composedTransform = Data<Scalar,DeviceType>::allocateMatMatResult(transposeLeft, leftTransform, transposeRight, rightTransform);
        composedTransform.storeMatMat(transposeLeft, leftTransform, transposeRight, rightTransform);
        
        // if the composedTransform matrices are full, the following is a good estimate.  If they have some diagonal portions, this will overcount.
        if (approximateFlops != NULL)
        {
          *approximateFlops += composedTransform.getUnderlyingViewSize() * (spaceDim - 1) * 2;
        }
      }
      else if (bothRank3) // (C,P,D)
      {
        // re-cast leftTransform as a rank 4 (C,P,1,D) object -- a 1 x D matrix at each (C,P).
        const int newRank   = 4;
        auto extents        = leftTransform.getExtents();
        auto variationTypes = leftTransform.getVariationTypes();
        extents[3]               = extents[2];
        extents[2]               = 1;
        variationTypes[3]        = variationTypes[2];
        variationTypes[2]        = CONSTANT;
        auto leftTransformMatrix = leftTransform.shallowCopy(newRank, extents, variationTypes);
        
        // re-cast rightTransform as a rank 4 (C,P,1,D) object -- a 1 x D matrix at each (C,P)
        extents                  = rightTransform.getExtents();
        variationTypes           = rightTransform.getVariationTypes();
        extents[3]               = extents[2];
        extents[2]               = 1;
        variationTypes[3]        = variationTypes[2];
        variationTypes[2]        = CONSTANT;
        auto rightTransformMatrix = rightTransform.shallowCopy(newRank, extents, variationTypes);
        
        composedTransform = Data<Scalar,DeviceType>::allocateMatMatResult(transposeLeft, leftTransformMatrix, transposeRight, rightTransformMatrix); // false: don't transpose
        composedTransform.storeMatMat(transposeLeft, leftTransformMatrix, transposeRight, rightTransformMatrix);
                
        if (approximateFlops != NULL)
        {
          *approximateFlops += composedTransform.getUnderlyingViewSize(); // one multiply per entry
        }
      }
      else if (bothRank2)
      {
        composedTransform = leftTransform.allocateInPlaceCombinationResult(leftTransform, rightTransform);
        composedTransform.storeInPlaceProduct(leftTransform, rightTransform);
        
        // re-cast composedTranform as a rank 4 (C,P,1,1) object -- a 1 x 1 matrix at each (C,P).
        const int newRank   = 4;
        auto extents        = composedTransform.getExtents();
        auto variationTypes = composedTransform.getVariationTypes();
        composedTransform = composedTransform.shallowCopy(newRank, extents, variationTypes);
        if (approximateFlops != NULL)
        {
          *approximateFlops += composedTransform.getUnderlyingViewSize(); // one multiply per entry
        }
      }
      else if (ranks32) // rank 2 / rank 3 combination.
      {
        const auto & rank3Transform = (leftRank == 3) ? leftTransform : rightTransform;
        const auto & rank2Transform = (leftRank == 2) ? leftTransform : rightTransform;
        
        composedTransform = DataTools::multiplyByCPWeights(rank3Transform, rank2Transform);
        
        // re-cast composedTransform as a rank 4 object:
        // logically, the original rank-3 transform can be understood as a 1xD matrix.  The composed transform is leftTransform^T * rightTransform, so:
        // - if left  has the rank-3 transform, composedTransform should be a (C,P,D,1) object -- a D x 1 matrix at each (C,P).
        // - if right has the rank-3 transform, composedTransform should be a (C,P,1,D) object -- a 1 x D matrix at each (C,P).
        const int newRank   = 4;
        auto extents        = composedTransform.getExtents();
        auto variationTypes = composedTransform.getVariationTypes();
        if (leftRank == 3)
        {
          // extents[3] and variationTypes[3] will already be 1 and CONSTANT, respectively
          // extents[3]               = 1;
          // variationTypes[3]        = CONSTANT;
        }
        else
        {
          extents[3]               = extents[2];
          extents[2]               = 1;
          variationTypes[3]        = variationTypes[2];
          variationTypes[2]        = CONSTANT;
        }
        composedTransform = composedTransform.shallowCopy(newRank, extents, variationTypes);
      }
      else if (ranks42) // rank 4 / rank 2 combination.
      {
        if (leftRank == 4)
        {
          // want to transpose left matrix, and multiply by the values from rightTransform
          // start with the multiplication:
          auto composedTransformTransposed = DataTools::multiplyByCPWeights(leftTransform, rightTransform);
          composedTransform = DataTools::transposeMatrix(composedTransformTransposed);
        }
        else // (leftRank == 2)
        {
          composedTransform = DataTools::multiplyByCPWeights(rightTransform, leftTransform);
        }
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported transform combination");
      }
    }
    else if (leftTransform.isValid())
    {
      // rightTransform is the identity
      switch (leftRank)
      {
        case 4: composedTransform = DataTools::transposeMatrix(leftTransform); break;
        case 3:
        {
          // - if left  has the rank-3 transform, composedTransform should be a (C,P,D,1) object -- a D x 1 matrix at each (C,P).
          const int newRank   = 4;
          auto extents        = leftTransform.getExtents();
          auto variationTypes = leftTransform.getVariationTypes();
          
          composedTransform = leftTransform.shallowCopy(newRank, extents, variationTypes);
        }
          break;
        case 2: composedTransform = leftTransform; break;
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported transform combination");
      }
    }
    else if (rightTransform.isValid())
    {
      // leftTransform is the identity
      composedTransform = rightTransform;
      switch (rightRank)
      {
        case 4: composedTransform = rightTransform; break;
        case 3:
        {
          // - if right has the rank-3 transform, composedTransform should be a (C,P,1,D) object -- a 1 x D matrix at each (C,P).
          const int newRank   = 4;
          auto extents        = rightTransform.getExtents();
          auto variationTypes = rightTransform.getVariationTypes();
          extents[3]          = extents[2];
          variationTypes[3]   = variationTypes[2];
          extents[2]          = 1;
          variationTypes[2]   = CONSTANT;
          
          composedTransform = rightTransform.shallowCopy(newRank, extents, variationTypes);
        }
          break;
        case 2: composedTransform = rightTransform; break;
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported transform combination");
      }
    }
    else
    {
      // both left and right transforms are identity
      Kokkos::Array<ordinal_type,4> extents {basisValuesLeft.numCells(),basisValuesLeft.numPoints(),spaceDim,spaceDim};
      Kokkos::Array<DataVariationType,4> variationTypes {CONSTANT,CONSTANT,BLOCK_PLUS_DIAGONAL,BLOCK_PLUS_DIAGONAL};
      
      Kokkos::View<Scalar*,DeviceType> identityUnderlyingView("Intrepid2::FST::integrate() - identity view",spaceDim);
      Kokkos::deep_copy(identityUnderlyingView, 1.0);
      composedTransform = Data<Scalar,DeviceType>(identityUnderlyingView,extents,variationTypes);
    }
    
//    timer->stop();
//    std::cout << "Completed mat-mat in " << timer->totalElapsedTime() << " seconds.\n";
//    timer->reset();
    
    const int leftComponentCount  = leftIsVectorValued ? basisValuesLeft. vectorData().numComponents() : 1;
    const int rightComponentCount = rightIsVectorValued ? basisValuesRight.vectorData().numComponents() : 1;
    
    int leftFieldOrdinalOffset = 0; // keeps track of the number of fields in prior families
    for (int leftFamilyOrdinal=0; leftFamilyOrdinal<leftFamilyCount; leftFamilyOrdinal++)
    {
      // "a" keeps track of the spatial dimension over which we are integrating in the left vector
      // components are allowed to span several dimensions; we keep track of the offset for the component in a_offset
      int a_offset = 0;
      bool haveLaunchedContributionToCurrentFamilyLeft = false; // helps to track whether we need a Kokkos::fence before launching a kernel.
      for (int leftComponentOrdinal=0; leftComponentOrdinal<leftComponentCount; leftComponentOrdinal++)
      {
        TensorData<Scalar,DeviceType> leftComponent = leftIsVectorValued ? basisValuesLeft.vectorData().getComponent(leftFamilyOrdinal, leftComponentOrdinal)
                                                                         : basisValuesLeft.basisValues().tensorData(leftFamilyOrdinal);
        if (!leftComponent.isValid())
        {
           // represents zero
          a_offset += basisValuesLeft.vectorData().numDimsForComponent(leftComponentOrdinal);
          continue;
        }
           
        int rightFieldOrdinalOffset = 0; // keeps track of the number of fields in prior families
        for (int rightFamilyOrdinal=0; rightFamilyOrdinal<rightFamilyCount; rightFamilyOrdinal++)
        {
          // "b" keeps track of the spatial dimension over which we are integrating in the right vector
          // components are allowed to span several dimensions; we keep track of the offset for the component in b_offset
          bool haveLaunchedContributionToCurrentFamilyRight = false; // helps to track whether we need a Kokkos::fence before launching a kernel.
          int b_offset = 0;
          for (int rightComponentOrdinal=0; rightComponentOrdinal<rightComponentCount; rightComponentOrdinal++)
          {
            TensorData<Scalar,DeviceType> rightComponent = rightIsVectorValued ? basisValuesRight.vectorData().getComponent(rightFamilyOrdinal, rightComponentOrdinal)
                                                                               : basisValuesRight.basisValues().tensorData(rightFamilyOrdinal);
            if (!rightComponent.isValid())
            {
               // represents zero
              b_offset += basisValuesRight.vectorData().numDimsForComponent(rightComponentOrdinal);
              continue;
            }
            
            INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(leftComponent.numTensorComponents() != rightComponent.numTensorComponents(), std::invalid_argument, "left TensorData and right TensorData have different number of tensor components.  This is not supported.");
            
            const int vectorSize = getVectorSizeForHierarchicalParallelism<Scalar>();
            Kokkos::TeamPolicy<ExecutionSpace> policy = Kokkos::TeamPolicy<ExecutionSpace>(cellDataExtent,Kokkos::AUTO(),vectorSize);
            
            // TODO: expose the options for forceNonSpecialized and usePointCacheForRank3Tensor through an IntegrationAlgorithm enumeration.
            // AUTOMATIC: let Intrepid2 choose an algorithm based on the inputs (and, perhaps, the execution space)
            // STANDARD: don't use sum factorization or axis alignment -- just do the simple contraction, a (p+1)^9 algorithm in 3D
            // SUM_FACTORIZATION                              // (p+1)^7 algorithm in 3D
            // SUM_FACTORIZATION_AXIS_ALIGNED                 // (p+1)^6 algorithm in 3D
            // SUM_FACTORIZATION_FORCE_GENERIC_IMPLEMENTATION // mainly intended for testing purposes (specialized implementations perform better when they are provided)
            // SUM_FACTORIZATION_WITH_POINT_CACHE             // novel (p+1)^7 (in 3D) algorithm in Intrepid2; unclear as yet when and whether this may be a superior approach
            bool forceNonSpecialized = false; // We might expose this in the integrate() arguments in the future.  We *should* default to false in the future.
            bool usePointCacheForRank3Tensor = true; // EXPERIMENTAL; has better performance under CUDA, but slightly worse performance than standard on serial CPU
            
            // in one branch or another below, we will launch a parallel kernel that contributes to (leftFamily, rightFamily) field ordinal pairs.
            // if we have already launched something that contributes to that part of the integral container, we need a Kokkos fence() to ensure that these do not interfere with each other.
            if (haveLaunchedContributionToCurrentFamilyLeft && haveLaunchedContributionToCurrentFamilyRight)
            {
              ExecutionSpace().fence();
            }
            haveLaunchedContributionToCurrentFamilyLeft  = true;
            haveLaunchedContributionToCurrentFamilyRight = true;
            
            if (integralViewRank == 2)
            {
              if (usePointCacheForRank3Tensor && (leftComponent.numTensorComponents() == 3))
              {
                auto functor = Impl::F_IntegratePointValueCache<Scalar, DeviceType, 2>(integrals, leftComponent, composedTransform, rightComponent, cellMeasures, a_offset, b_offset, leftFieldOrdinalOffset, rightFieldOrdinalOffset);
                
                const int recommendedTeamSize = policy.team_size_recommended(functor,Kokkos::ParallelForTag());
                const int teamSize            = functor.teamSize(recommendedTeamSize);
                
                policy = Kokkos::TeamPolicy<DeviceType>(cellDataExtent,teamSize,vectorSize);
                
                Kokkos::parallel_for("F_IntegratePointValueCache rank 2", policy, functor);
                
                if (approximateFlops != NULL)
                {
                  *approximateFlops += functor.approximateFlopCountPerCell() * integrals.getDataExtent(0);
                }
              }
              else
              {
                auto functor = Impl::F_Integrate<Scalar, DeviceType, 2>(integrals, leftComponent, composedTransform, rightComponent, cellMeasures, a_offset, b_offset, leftFieldOrdinalOffset, rightFieldOrdinalOffset, forceNonSpecialized);
                
                const int recommendedTeamSize = policy.team_size_recommended(functor,Kokkos::ParallelForTag());
                const int teamSize            = functor.teamSize(recommendedTeamSize);
                
                policy = Kokkos::TeamPolicy<ExecutionSpace>(cellDataExtent,teamSize,vectorSize);
                
                Kokkos::parallel_for("F_Integrate rank 2", policy, functor);
                
                if (approximateFlops != NULL)
                {
                  *approximateFlops += functor.approximateFlopCountPerCell() * integrals.getDataExtent(0);
                }
              }
            }
            else if (integralViewRank == 3)
            {
              if (usePointCacheForRank3Tensor && (leftComponent.numTensorComponents() == 3))
              {
                auto functor = Impl::F_IntegratePointValueCache<Scalar, DeviceType, 3>(integrals, leftComponent, composedTransform, rightComponent, cellMeasures, a_offset, b_offset, leftFieldOrdinalOffset, rightFieldOrdinalOffset);
                
                const int recommendedTeamSize = policy.team_size_recommended(functor,Kokkos::ParallelForTag());
                const int teamSize            = functor.teamSize(recommendedTeamSize);
                
                policy = Kokkos::TeamPolicy<ExecutionSpace>(cellDataExtent,teamSize,vectorSize);
                
                Kokkos::parallel_for("F_IntegratePointValueCache rank 3", policy, functor);
                
                if (approximateFlops != NULL)
                {
                  *approximateFlops += functor.approximateFlopCountPerCell() * integrals.getDataExtent(0);
                }
              }
              else
              {
                auto functor = Impl::F_Integrate<Scalar, DeviceType, 3>(integrals, leftComponent, composedTransform, rightComponent, cellMeasures, a_offset, b_offset, leftFieldOrdinalOffset, rightFieldOrdinalOffset, forceNonSpecialized);
                
                const int recommendedTeamSize = policy.team_size_recommended(functor,Kokkos::ParallelForTag());
                const int teamSize            = functor.teamSize(recommendedTeamSize);
                
                policy = Kokkos::TeamPolicy<DeviceType>(cellDataExtent,teamSize,vectorSize);
                
                Kokkos::parallel_for("F_Integrate rank 3", policy, functor);
                
                if (approximateFlops != NULL)
                {
                  *approximateFlops += functor.approximateFlopCountPerCell() * integrals.getDataExtent(0);
                }
              }
            }
            b_offset += rightIsVectorValued ? basisValuesRight.vectorData().numDimsForComponent(rightComponentOrdinal) : 1;
          }
          rightFieldOrdinalOffset += rightIsVectorValued ? basisValuesRight.vectorData().numFieldsInFamily(rightFamilyOrdinal) : basisValuesRight.basisValues().numFieldsInFamily(rightFamilyOrdinal);
        }
        a_offset += leftIsVectorValued ? basisValuesLeft.vectorData().numDimsForComponent(leftComponentOrdinal) : 1;
      }
      leftFieldOrdinalOffset += leftIsVectorValued ? basisValuesLeft.vectorData().numFieldsInFamily(leftFamilyOrdinal) : basisValuesLeft.basisValues().numFieldsInFamily(leftFamilyOrdinal);
    }
  }
//  if (approximateFlops != NULL)
//  {
//    std::cout << "Approximate flop count (new): " << *approximateFlops << std::endl;
//  }
  ExecutionSpace().fence(); // make sure we've finished writing to integrals container before we return
}

TEUCHOS_UNIT_TEST( PR14546, Distill14546SegFault )
{
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  using D  = Data<DataScalar,DeviceType>;
  using TD = TensorData<DataScalar,DeviceType>;
  using VD = VectorData<DataScalar,DeviceType>;
  using TBV = TransformedBasisValues<DataScalar,DeviceType>;
  
  const int spaceDim = 1;
    
  auto oneElementView = getFixedRankView<DataScalar>("oneElementView", 1);
  Kokkos::deep_copy(oneElementView, 1.0);
  
  Kokkos::Array<int,2> extents {1,1};
  Kokkos::Array<DataVariationType,2> variationTypes {GENERAL,GENERAL};
  D fieldComponentData(oneElementView,extents,variationTypes);

  TD  tensorData(std::vector<D>{fieldComponentData,fieldComponentData});
  
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim, spaceDim);
  Kokkos::deep_copy(identityMatrixView, 1.0);
  
  Kokkos::Array<int,4> transformExtents {1, 1, 1, 1};
  Kokkos::Array<DataVariationType,4> transformationVariationType {GENERAL, GENERAL, GENERAL, GENERAL};
  
  D explicitIdentityMatrix(identityMatrixView, transformExtents, transformationVariationType);
  
  const int numFamilies = 2;
  TD nullTD;
  Kokkos::Array<TD, spaceDim > firstFamilyLeft  {tensorData};
  Kokkos::Array< Kokkos::Array<TD, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft, firstFamilyLeft};
  
  VD vectorDataLeft(vectorComponentsLeft);
  
  Kokkos::Array<TD, spaceDim > firstFamilyRight  {tensorData};
  Kokkos::Array< Kokkos::Array<TD, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight, firstFamilyRight};
  
  // imitate VD construction vectorDataRight(vectorComponentsRight)
  Kokkos::Array< Kokkos::Array<TD, spaceDim>, numFamilies> vectorComponents = vectorComponentsRight;
  using VectorArray = Kokkos::Array< TD, Parameters::MaxVectorComponents >; // for axis-aligned case, these correspond entry-wise to the axis with which the vector values align
  using FamilyVectorArray = Kokkos::Array< VectorArray, Parameters::MaxTensorComponents>;

  FamilyVectorArray vectorComponents_; // outer: family ordinal; inner: component/spatial dimension ordinal
  vectorComponents_[0][0] = vectorComponents[0][0];
  vectorComponents_[1][0] = vectorComponents[1][0];
  
  TBV  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
  TBV transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataLeft);
  
  D constantCellMeasuresData(1.0, Kokkos::Array<int,2>{1,1});
  TD constantCellMeasures(constantCellMeasuresData);
  
  // these assignments imitate a function call with arguments (tbvLeft, cellMeasures, tbvRight)
  const TBV      tbvLeft = transformedUnitVectorDataLeft;
  const TD  cellMeasures = constantCellMeasures;
  const TBV     tbvRight = transformedUnitVectorDataLeft;
  
  using IT = Intrepid2::IntegrationTools<DeviceType>;
  
//  D integralsBaseline  = IT::allocateIntegralData(tbvLeft, cellMeasures, tbvRight);
  // imitate call to allocateIntegralData() for integralsBaseline
  TBV tbvLeft_aid     = tbvLeft;
  TD cellMeasures_aid = cellMeasures;
  TBV tbvRight_aid    = tbvRight;
  
  D integralsIntegrate = IT::allocateIntegralData(tbvLeft, cellMeasures, tbvRight);
  
  // these assignments imitate a function call (to integrate_baseline) with arguments (integralsBaseline, tbvLeft, cellMeasures, tbvRight)
//  D integrals_bl = integralsBaseline;
  const TBV tbvLeft_bl  = tbvLeft;
//  const TD cellMeasures_bl = cellMeasures;
  const TBV tbvRight_bl = tbvRight;
  
  IT_integrate<DeviceType,DataScalar>(integralsIntegrate, tbvLeft, cellMeasures, tbvRight);
}

TEUCHOS_UNIT_TEST( PR14546, AllocationIssue )
{
  // not sure if this is the same issue -- I think so, but in case it isn't, I'm creating a separate test
  using DataScalar  = double;
  using DeviceType = DefaultTestDeviceType;
  using D  = Data<DataScalar,DeviceType>;
  using TD = TensorData<DataScalar,DeviceType>;
  using VD = VectorData<DataScalar,DeviceType>;
  using TBV = TransformedBasisValues<DataScalar,DeviceType>;
  
  const int spaceDim = 1;
  
  auto oneElementView = getFixedRankView<DataScalar>("oneElementView", 1);
  Kokkos::deep_copy(oneElementView, 1.0);
  
  Kokkos::Array<int,2> extents {1,1};
  Kokkos::Array<DataVariationType,2> variationTypes {GENERAL,GENERAL};
  D fieldComponentData(oneElementView,extents,variationTypes);

  TD tensorData(std::vector<D>{fieldComponentData});
  
  auto identityMatrixView = getFixedRankView<DataScalar>("identity matrix", spaceDim, spaceDim);
  Kokkos::deep_copy(identityMatrixView, 1.0);
  
  Kokkos::Array<int,4> transformExtents {1, 1, 1, 1};
  Kokkos::Array<DataVariationType,4> transformationVariationType {GENERAL, GENERAL, GENERAL, GENERAL};
  
  D explicitIdentityMatrix(identityMatrixView, transformExtents, transformationVariationType);
  {
    auto data = getMatchingViewWithLabel(identityMatrixView, "Data mat-mat result", 1, 1, 1, 1);
    std::cout << "data.size(): " << data.size() << std::endl;
    std::cout << "Got to line " << __LINE__ << std::endl;
  }
  
  const int numFamilies = 1;
  TD nullTD;
  Kokkos::Array<TD, spaceDim > firstFamilyLeft  {tensorData};
//  Kokkos::Array<TD, spaceDim > secondFamilyLeft {tensorData};
  Kokkos::Array< Kokkos::Array<TD, spaceDim>, numFamilies> vectorComponentsLeft {firstFamilyLeft};//, secondFamilyLeft};
  
  VD vectorDataLeft(vectorComponentsLeft);
  
  Kokkos::Array<TD, spaceDim > firstFamilyRight  {tensorData};
  Kokkos::Array<TD, spaceDim > secondFamilyRight {tensorData};
  Kokkos::Array< Kokkos::Array<TD, spaceDim>, numFamilies> vectorComponentsRight {firstFamilyRight}; //, secondFamilyRight};
  
  VD vectorDataRight(vectorComponentsRight);
  
  TBV  transformedUnitVectorDataLeft(explicitIdentityMatrix,vectorDataLeft);
  TBV transformedUnitVectorDataRight(explicitIdentityMatrix,vectorDataLeft);
  
  D constantCellMeasuresData(1.0, Kokkos::Array<int,2>{1,1});
  TD constantCellMeasures(constantCellMeasuresData);
  
  // these assignments imitate a function call with arguments (tbvLeft, cellMeasures, tbvRight)
  const TBV      tbvLeft = transformedUnitVectorDataLeft;
  const TD  cellMeasures = constantCellMeasures;
  const TBV     tbvRight = transformedUnitVectorDataLeft;
  
  using IT = Intrepid2::IntegrationTools<DeviceType>;
  
  auto integralsBaseline  = IT::allocateIntegralData(tbvLeft, cellMeasures, tbvRight);
  auto integralsIntegrate = IT::allocateIntegralData(tbvLeft, cellMeasures, tbvRight);
  
  {
    auto data = getMatchingViewWithLabel(identityMatrixView, "Data mat-mat result", 1, 1, 1, 1);
    std::cout << "data.size(): " << data.size() << std::endl;
    std::cout << "Got to line " << __LINE__ << std::endl;
  }
  
  // these assignments imitate a function call (to integrate_baseline) with arguments (integralsBaseline, tbvLeft, cellMeasures, tbvRight)
  D integrals_bl = integralsBaseline;
  const TBV tbvLeft_bl  = tbvLeft;
  const TD cellMeasures_bl = cellMeasures;
  const TBV tbvRight_bl = tbvRight;
    
  {
    auto data = getMatchingViewWithLabel(identityMatrixView, "Data mat-mat result", 1, 1, 1, 1);
    std::cout << "data.size(): " << data.size() << std::endl;
    std::cout << "Got to line " << __LINE__ << std::endl;
  }
  
  IT_integrate<DeviceType,DataScalar>(integralsIntegrate, tbvLeft, cellMeasures, tbvRight);
}

} // anonymous namespace
