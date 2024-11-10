// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_TensorBasis.hpp
 \brief  Implementation of bases that are tensor products of two or three component bases.
 \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_TensorBasis_h
#define Intrepid2_TensorBasis_h

#include <Kokkos_DynRankView.hpp>

#include <Teuchos_RCP.hpp>

#include <Intrepid2_config.h>

#include <map>
#include <set>
#include <vector>

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_DeviceAssert.hpp"
#include "Intrepid2_TensorTopologyMap.hpp"
#include "Intrepid2_TensorViewIterator.hpp"
#include "Intrepid2_Utils.hpp" // defines FAD_VECTOR_SIZE, VECTOR_SIZE

#include "Intrepid2_CellTopology.hpp"

namespace Intrepid2
{
  template<ordinal_type spaceDim>
  KOKKOS_INLINE_FUNCTION
  ordinal_type getDkEnumeration(const Kokkos::Array<int,spaceDim> &entries);

  template<ordinal_type spaceDim>
  KOKKOS_INLINE_FUNCTION
  void getDkEnumerationInverse(Kokkos::Array<int,spaceDim> &entries, const ordinal_type dkEnum, const ordinal_type operatorOrder);
  
  template<>
  KOKKOS_INLINE_FUNCTION
  void getDkEnumerationInverse<1>(Kokkos::Array<int,1> &entries, const ordinal_type dkEnum, const ordinal_type operatorOrder)
  {
    entries[0] = operatorOrder;
  }
  
  template<>
  KOKKOS_INLINE_FUNCTION
  void getDkEnumerationInverse<2>(Kokkos::Array<int,2> &entries, const ordinal_type dkEnum, const ordinal_type operatorOrder)
  {
    entries[0] = operatorOrder - dkEnum;
    entries[1] = dkEnum;
  }
  
  template<>
  KOKKOS_INLINE_FUNCTION
  void getDkEnumerationInverse<3>(Kokkos::Array<int,3> &entries, const ordinal_type dkEnum, const ordinal_type operatorOrder)
  {
    // formula is zMult + (yMult+zMult)*(yMult+zMult+1)/2; where xMult+yMult+zMult = operatorOrder
    // it seems complicated to work out a formula that will invert this.  For the present we just take a brute force approach,
    // using getDkEnumeration() to check each possibility
    for (ordinal_type yMult=0; yMult<=operatorOrder; yMult++)
    {
      for (ordinal_type zMult=0; zMult<=operatorOrder-yMult; zMult++)
      {
        const ordinal_type xMult = operatorOrder-(zMult+yMult);
        if (dkEnum == getDkEnumeration<3>(xMult,yMult,zMult))
        {
          entries[0] = xMult;
          entries[1] = yMult;
          entries[2] = zMult;
          return;
        }
      }
    }
  }
  
  template<ordinal_type spaceDim>
  KOKKOS_INLINE_FUNCTION
  void getDkEnumerationInverse(Kokkos::Array<int,spaceDim> &entries, const ordinal_type dkEnum, const ordinal_type operatorOrder)
  {
    // for operator order k, the recursive formula defining getDkEnumeration is:
    // getDkEnumeration(k0,k1,…,k_{n-1}) = getDkCardinality(k - k0) + getDkEnumeration(k1,…,k_{n-1})
    // The entries are in reverse lexicographic order.  We search for k0, by incrementing k0 until getDkEnumeration(k0,0,…,0) <= dkEnum
    // Then we recursively call getDkEnumerationInverse<spaceDim-1>({k1,…,k_{n-1}}, dkEnum - getDkEnumeration(k0,0,…,0) - 1)
    
    for (int k0=0; k0<=operatorOrder; k0++)
    {
      entries[0] = k0;
      for (int d=1; d<spaceDim-1; d++)
      {
        entries[d] = 0;
      }
      // sum of entries must be equal to operatorOrder
      if (spaceDim > 1) entries[spaceDim-1] = operatorOrder - k0;
      else if (k0 != operatorOrder) continue; // if spaceDim == 1, then the only way the sum of the entries is operatorOrder is if k0 == operatorOrder
      const ordinal_type dkEnumFor_k0 = getDkEnumeration<spaceDim>(entries);
      
      if      (dkEnumFor_k0 > dkEnum) continue; // next k0
      else if (dkEnumFor_k0 == dkEnum) return;  // entries has (k0,0,…,0), and this has dkEnum as its enumeration value
      else
      {
        // (k0,0,…,0) is prior to the dkEnum entry, which means that the dkEnum entry starts with k0-1.
        entries[0] = k0 - 1;
        
        // We determine the rest of the entries through a recursive call to getDkEnumerationInverse<spaceDim - 1>().

        // ensure that we don't try to allocate an empty array…
        constexpr ordinal_type sizeForSubArray = (spaceDim > 2) ? spaceDim - 1 : 1;
        Kokkos::Array<int,sizeForSubArray> subEntries = {};
        
        // the -1 in sub-entry enumeration value accounts for the fact that the entry is the one *after* (k0,0,…,0)
        getDkEnumerationInverse<spaceDim-1>(subEntries, dkEnum - dkEnumFor_k0 - 1, operatorOrder - entries[0]);
        
        for (int i=1; i<spaceDim; i++)
        {
          entries[i] = subEntries[i-1];
        }
        return;
      }
    }
    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "entries corresponding to dkEnum not found");
  }
  
  template<>
  KOKKOS_INLINE_FUNCTION
  ordinal_type getDkEnumeration<1>(const Kokkos::Array<int,1> &entries)
  {
    return getDkEnumeration<1>(entries[0]);
  }
  
  template<ordinal_type spaceDim>
  KOKKOS_INLINE_FUNCTION
  ordinal_type getDkEnumeration(const Kokkos::Array<int,spaceDim> &entries)
  {
    ordinal_type k_minus_k0 = 0; // sum of all the entries but the first
    
    // recursive formula in general is: getDkEnumeration(k0,k1,…,k_{n-1}) = getDkCardinality(k - k0) + getDkEnumeration(k1,…,k_{n-1})
    // ensure that we don't try to allocate an empty array…
    constexpr ordinal_type sizeForSubArray = (spaceDim > 2) ? spaceDim - 1 : 1;
    Kokkos::Array<int,sizeForSubArray> remainingEntries;
    for (int i=1; i<spaceDim; i++)
    {
      k_minus_k0 += entries[i];
      remainingEntries[i-1] = entries[i];
    }
    
    if (k_minus_k0 == 0)
    {
      return 0;
    }
    else
    {
      EOperator opFor_k_minus_k0_minus_1 = (k_minus_k0 > 1) ? EOperator(OPERATOR_D1 + k_minus_k0 - 2) : EOperator(OPERATOR_VALUE);
      const ordinal_type dkCardinality = getDkCardinality(opFor_k_minus_k0_minus_1, spaceDim);
      const ordinal_type dkEnum = dkCardinality + getDkEnumeration<sizeForSubArray>(remainingEntries);
      return dkEnum;
    }
  }
  
  template<ordinal_type spaceDim1, ordinal_type spaceDim2>
  KOKKOS_INLINE_FUNCTION
  ordinal_type getDkTensorIndex(const ordinal_type dkEnum1, const ordinal_type operatorOrder1,
                                const ordinal_type dkEnum2, const ordinal_type operatorOrder2)
  {
    Kokkos::Array<int,spaceDim1> entries1 = {};
    getDkEnumerationInverse<spaceDim1>(entries1, dkEnum1, operatorOrder1);
    
    Kokkos::Array<int,spaceDim2> entries2 = {};
    getDkEnumerationInverse<spaceDim2>(entries2, dkEnum2, operatorOrder2);
    
    const int spaceDim = spaceDim1 + spaceDim2;
    Kokkos::Array<int,spaceDim> entries;
    
    for (ordinal_type d=0; d<spaceDim1; d++)
    {
      entries[d] = entries1[d];
    }
    
    for (ordinal_type d=0; d<spaceDim2; d++)
    {
      entries[d+spaceDim1] = entries2[d];
    }
    
    return getDkEnumeration<spaceDim>(entries);
  }

template<typename BasisBase>
class Basis_TensorBasis;

/** \struct  Intrepid2::OperatorTensorDecomposition
    \brief   For a multi-component tensor basis, specifies the operators to be applied to the components to produce the composite operator on the tensor basis.
*/
struct OperatorTensorDecomposition
{
  // if we want to make this usable on device, we could switch to Kokkos::Array instead of std::vector.  But this is not our immediate use case.
  std::vector< std::vector<EOperator> > ops; // outer index: vector entry ordinal; inner index: basis component ordinal. (scalar-valued operators have a single entry in outer vector)
  std::vector<double> weights; // weights for each vector entry
  ordinal_type numBasisComponents_;
  bool rotateXYNinetyDegrees_ = false; // if true, indicates that something that otherwise would have values (f_x, f_y, …) should be mapped to (-f_y, f_x, …).  This is used for H(curl) wedges (specifically, for OPERATOR_CURL).
  
  OperatorTensorDecomposition(const std::vector<EOperator> &opsBasis1, const std::vector<EOperator> &opsBasis2, const std::vector<double> vectorComponentWeights)
  :
  weights(vectorComponentWeights),
  numBasisComponents_(2)
  {
    const ordinal_type size = opsBasis1.size();
    const ordinal_type opsBasis2Size = opsBasis2.size();
    const ordinal_type weightsSize = weights.size();
    INTREPID2_TEST_FOR_EXCEPTION(size != opsBasis2Size, std::invalid_argument, "opsBasis1.size() != opsBasis2.size()");
    INTREPID2_TEST_FOR_EXCEPTION(size != weightsSize,   std::invalid_argument, "opsBasis1.size() != weights.size()");
    
    for (ordinal_type i=0; i<size; i++)
    {
      ops.push_back(std::vector<EOperator>{opsBasis1[i],opsBasis2[i]});
    }
  }
  
  OperatorTensorDecomposition(const std::vector< std::vector<EOperator> > &vectorEntryOps, const std::vector<double> &vectorComponentWeights)
  :
  ops(vectorEntryOps),
  weights(vectorComponentWeights)
  {
    const ordinal_type numVectorComponents = ops.size();
    const ordinal_type weightsSize = weights.size();
    INTREPID2_TEST_FOR_EXCEPTION(numVectorComponents != weightsSize,   std::invalid_argument, "opsBasis1.size() != weights.size()");
    
    INTREPID2_TEST_FOR_EXCEPTION(numVectorComponents == 0,   std::invalid_argument, "must have at least one entry!");
    
    ordinal_type numBases = 0;
    for (ordinal_type i=0; i<numVectorComponents; i++)
    {
      if (numBases == 0)
      {
        numBases = ops[i].size();
      }
      else if (ops[i].size() != 0)
      {
        const ordinal_type opsiSize = ops[i].size();
        INTREPID2_TEST_FOR_EXCEPTION(numBases != opsiSize, std::invalid_argument, "must have one operator for each basis in each nontrivial entry in vectorEntryOps");
      }
    }
    INTREPID2_TEST_FOR_EXCEPTION(numBases == 0, std::invalid_argument, "at least one vectorEntryOps entry must be non-trivial");
    numBasisComponents_ = numBases;
  }
  
  OperatorTensorDecomposition(const std::vector<EOperator> &basisOps, const double weight = 1.0)
  :
  ops({basisOps}),
  weights({weight}),
  numBasisComponents_(basisOps.size())
  {}
  
  OperatorTensorDecomposition(const EOperator &opBasis1, const EOperator &opBasis2, double weight = 1.0)
  :
  ops({ std::vector<EOperator>{opBasis1, opBasis2} }),
  weights({weight}),
  numBasisComponents_(2)
  {}
  
  OperatorTensorDecomposition(const EOperator &opBasis1, const EOperator &opBasis2, const EOperator &opBasis3, double weight = 1.0)
  :
  ops({ std::vector<EOperator>{opBasis1, opBasis2, opBasis3} }),
  weights({weight}),
  numBasisComponents_(3)
  {}
  
  ordinal_type numVectorComponents() const
  {
    return ops.size(); // will match weights.size()
  }
  
  ordinal_type numBasisComponents() const
  {
    return numBasisComponents_;
  }
  
  double weight(const ordinal_type &vectorComponentOrdinal) const
  {
    return weights[vectorComponentOrdinal];
  }
  
  bool identicallyZeroComponent(const ordinal_type &vectorComponentOrdinal) const
  {
    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(vectorComponentOrdinal < 0,                      std::invalid_argument, "vectorComponentOrdinal is out of bounds");
    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(vectorComponentOrdinal >= numVectorComponents(), std::invalid_argument, "vectorComponentOrdinal is out of bounds");
    return ops[vectorComponentOrdinal].size() == 0;
  }
  
  EOperator op(const ordinal_type &vectorComponentOrdinal, const ordinal_type &basisOrdinal) const
  {
    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(vectorComponentOrdinal < 0,                      std::invalid_argument, "vectorComponentOrdinal is out of bounds");
    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(vectorComponentOrdinal >= numVectorComponents(), std::invalid_argument, "vectorComponentOrdinal is out of bounds");
    if (identicallyZeroComponent(vectorComponentOrdinal))
    {
      return OPERATOR_MAX; // by convention: zero in this component
    }
    else
    {
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(basisOrdinal < 0,                    std::invalid_argument, "basisOrdinal is out of bounds");
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(basisOrdinal >= numBasisComponents_, std::invalid_argument, "basisOrdinal is out of bounds");
      return ops[vectorComponentOrdinal][basisOrdinal];
    }
  }
  
  //! takes as argument bases that are components in this decomposition, and decomposes them further if they are tensor bases.  Returns a fully expanded decomposition.
  template<typename DeviceType, typename OutputValueType, class PointValueType>
  OperatorTensorDecomposition expandedDecomposition(std::vector< Teuchos::RCP<Basis<DeviceType,OutputValueType,PointValueType> > > &bases)
  {
    const ordinal_type basesSize = bases.size();
    INTREPID2_TEST_FOR_EXCEPTION(basesSize != numBasisComponents_, std::invalid_argument, "The number of bases provided must match the number of basis components in this decomposition");
    
    ordinal_type numExpandedBasisComponents = 0;
    using BasisBase   = Basis<DeviceType,OutputValueType,PointValueType>;
    using TensorBasis = Basis_TensorBasis<BasisBase>;
    std::vector<TensorBasis*> basesAsTensorBasis(numBasisComponents_);
    for (ordinal_type basisComponentOrdinal=0; basisComponentOrdinal<numBasisComponents_; basisComponentOrdinal++)
    {
      TensorBasis* basisAsTensorBasis = dynamic_cast<TensorBasis*>(bases[basisComponentOrdinal].get());
      basesAsTensorBasis[basisComponentOrdinal] = basisAsTensorBasis;
      if (basisAsTensorBasis)
      {
        numExpandedBasisComponents += basisAsTensorBasis->getTensorBasisComponents().size();
      }
      else
      {
        numExpandedBasisComponents += 1;
      }
    }
    
    std::vector< std::vector<EOperator> > expandedOps; // outer index: vector entry ordinal; inner index: basis component ordinal.
    std::vector<double> expandedWeights;
    const ordinal_type opsSize = ops.size();
    for (ordinal_type simpleVectorEntryOrdinal=0; simpleVectorEntryOrdinal<opsSize; simpleVectorEntryOrdinal++)
    {
      if (identicallyZeroComponent(simpleVectorEntryOrdinal))
      {
        expandedOps.push_back(std::vector<EOperator>{});
        expandedWeights.push_back(0.0);
        continue;
      }
      
      std::vector< std::vector<EOperator> > expandedBasisOpsForSimpleVectorEntry(1); // start out with one outer entry; expands if a component is vector-valued
      
      // this lambda appends an op to each of the vector components
      auto addExpandedOp = [&expandedBasisOpsForSimpleVectorEntry](const EOperator &op)
      {
        const ordinal_type size = expandedBasisOpsForSimpleVectorEntry.size();
        for (ordinal_type i=0; i<size; i++)
        {
          expandedBasisOpsForSimpleVectorEntry[i].push_back(op);
        }
      };
      
      // this lambda takes a scalar-valued (single outer entry) expandedBasisOps and expands it
      // according to the number of vector entries coming from the vector-valued component basis
      auto vectorizeExpandedOps = [&expandedBasisOpsForSimpleVectorEntry](const int &numSubVectors)
      {
        // we require that this only gets called once per simpleVectorEntryOrdinal -- i.e., only one basis component gets to be vector-valued.
        INTREPID2_TEST_FOR_EXCEPTION(expandedBasisOpsForSimpleVectorEntry.size() != 1, std::invalid_argument, "multiple basis components may not be vector-valued!");
        for (ordinal_type i=1; i<numSubVectors; i++)
        {
          expandedBasisOpsForSimpleVectorEntry.push_back(expandedBasisOpsForSimpleVectorEntry[0]);
        }
      };
      
      std::vector<EOperator> subVectorOps;     // only used if one of the components is vector-valued
      std::vector<double> subVectorWeights {weights[simpleVectorEntryOrdinal]};
      for (ordinal_type basisComponentOrdinal=0; basisComponentOrdinal<numBasisComponents_; basisComponentOrdinal++)
      {
        const auto &op = ops[simpleVectorEntryOrdinal][basisComponentOrdinal];
        
        if (! basesAsTensorBasis[basisComponentOrdinal])
        {
          addExpandedOp(op);
        }
        else
        {
          OperatorTensorDecomposition basisOpDecomposition = basesAsTensorBasis[basisComponentOrdinal]->getOperatorDecomposition(op);
          if (basisOpDecomposition.numVectorComponents() > 1)
          {
            // We don't currently support a use case where we have multiple component bases that are vector-valued:
            INTREPID2_TEST_FOR_EXCEPTION(subVectorWeights.size() > 1, std::invalid_argument, "Unhandled case: multiple component bases are vector-valued");
            // We do support a single vector-valued case, though; this splits the current simpleVectorEntryOrdinal into an appropriate number of components that come in order in the expanded vector
            ordinal_type numSubVectors = basisOpDecomposition.numVectorComponents();
            vectorizeExpandedOps(numSubVectors);
            
            double weightSoFar = subVectorWeights[0];
            for (ordinal_type subVectorEntryOrdinal=1; subVectorEntryOrdinal<numSubVectors; subVectorEntryOrdinal++)
            {
              subVectorWeights.push_back(weightSoFar * basisOpDecomposition.weight(subVectorEntryOrdinal));
            }
            subVectorWeights[0] *= basisOpDecomposition.weight(0);
            for (ordinal_type subVectorEntryOrdinal=0; subVectorEntryOrdinal<numSubVectors; subVectorEntryOrdinal++)
            {
              for (ordinal_type subComponentBasis=0; subComponentBasis<basisOpDecomposition.numBasisComponents(); subComponentBasis++)
              {
                const auto &basisOp = basisOpDecomposition.op(subVectorEntryOrdinal, subComponentBasis);
                expandedBasisOpsForSimpleVectorEntry[subVectorEntryOrdinal].push_back(basisOp);
              }
            }
          }
          else
          {
            double componentWeight = basisOpDecomposition.weight(0);
            const ordinal_type size = subVectorWeights.size();
            for (ordinal_type i=0; i<size; i++)
            {
              subVectorWeights[i] *= componentWeight;
            }
            ordinal_type subVectorEntryOrdinal = 0;
            const ordinal_type numBasisComponents = basisOpDecomposition.numBasisComponents();
            for (ordinal_type subComponentBasis=0; subComponentBasis<numBasisComponents; subComponentBasis++)
            {
              const auto &basisOp = basisOpDecomposition.op(subVectorEntryOrdinal, basisComponentOrdinal);
              addExpandedOp( basisOp );
            }
          }
        }
      }
      
      // sanity check on the new expandedOps entries:
      for (ordinal_type i=0; i<static_cast<ordinal_type>(expandedBasisOpsForSimpleVectorEntry.size()); i++)
      {
        const ordinal_type size = expandedBasisOpsForSimpleVectorEntry[i].size();
        INTREPID2_TEST_FOR_EXCEPTION(size != numExpandedBasisComponents, std::logic_error, "each vector in expandedBasisOpsForSimpleVectorEntry should have as many entries as there are expanded basis components");
      }
      
      expandedOps.insert(expandedOps.end(), expandedBasisOpsForSimpleVectorEntry.begin(), expandedBasisOpsForSimpleVectorEntry.end());
      expandedWeights.insert(expandedWeights.end(), subVectorWeights.begin(), subVectorWeights.end());
    }
    // check that vector lengths agree:
    INTREPID2_TEST_FOR_EXCEPTION(expandedOps.size() != expandedWeights.size(), std::logic_error, "expandedWeights and expandedOps do not agree on the number of vector components");
    
    OperatorTensorDecomposition result(expandedOps, expandedWeights);
    result.setRotateXYNinetyDegrees(rotateXYNinetyDegrees_);
    return result;
  }
  
  //! If true, this flag indicates that a vector component that spans the first two dimensions should be rotated by 90 degrees clockwise, mapping (x,y) to (-y,x).  If there is no such vector component, the flag should be ignored.  As of this writing, this is used only by the "derived" H(curl) basis on the wedge.
  bool rotateXYNinetyDegrees() const
  {
    return rotateXYNinetyDegrees_;
  }
  
  void setRotateXYNinetyDegrees(const bool &value)
  {
    rotateXYNinetyDegrees_ = value;
  }
};

  /** \class  Intrepid2::TensorViewFunctor
      \brief  Functor for computing values for the TensorBasis class.
   
   This functor is not intended for use outside of \ref Intrepid2::Basis_TensorBasis.
  */
  template<class ExecutionSpace, class OutputScalar, class OutputFieldType>
  class TensorViewFunctor
  {
    using ScratchSpace       = typename ExecutionSpace::scratch_memory_space;
    using OutputScratchView  = Kokkos::View<OutputScalar*,ScratchSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    
    using TeamPolicy = Kokkos::TeamPolicy<ExecutionSpace>;
    using TeamMember = typename TeamPolicy::member_type;
    
    using TensorViewIteratorType = ::Intrepid2::TensorViewIterator<OutputFieldType, OutputFieldType, OutputFieldType, OutputScalar>;
    using RankCombinationType = typename TensorViewIteratorType::RankCombinationType;
    
    OutputFieldType  output_; // F,P[,D…]
    OutputFieldType  input1_; // F1,P[,D…] or F1,P1[,D…]
    OutputFieldType  input2_; // F2,P[,D…] or F2,P2[,D…]
    
    int numFields_, numPoints_;
    int numFields1_, numPoints1_;
    int numFields2_, numPoints2_;
    
    bool tensorPoints_; // if true, input1 and input2 refer to values at decomposed points, and P = P1 * P2.  If false, then the two inputs refer to points in the full-dimensional space, and their point lengths are the same as that of the final output.
    
    using RankCombinationViewType = typename TensorViewIteratorType::RankCombinationViewType;
    RankCombinationViewType rank_combinations_;// indicates the policy by which the input views will be combined in output view
    
    double weight_;
    
  public:
    
    TensorViewFunctor(OutputFieldType output, OutputFieldType inputValues1, OutputFieldType inputValues2,
                      bool tensorPoints, double weight)
    : output_(output), input1_(inputValues1), input2_(inputValues2), tensorPoints_(tensorPoints), weight_(weight)
    {
      numFields_ = output.extent_int(0);
      numPoints_ = output.extent_int(1);
      
      numFields1_ = inputValues1.extent_int(0);
      numPoints1_ = inputValues1.extent_int(1);
      
      numFields2_ = inputValues2.extent_int(0);
      numPoints2_ = inputValues2.extent_int(1);
      
      if (!tensorPoints_)
      {
        // then the point counts should all match
        INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != numPoints1_, std::invalid_argument, "incompatible point counts");
        INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != numPoints2_, std::invalid_argument, "incompatible point counts");
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != numPoints1_ * numPoints2_, std::invalid_argument, "incompatible point counts");
      }
      
      INTREPID2_TEST_FOR_EXCEPTION(numFields_ != numFields1_ * numFields2_, std::invalid_argument, "incompatible field sizes");
      
      const ordinal_type max_rank = std::max(inputValues1.rank(),inputValues2.rank());
      // at present, no supported case will result in an output rank greater than both input ranks
      
      const ordinal_type outputRank = output.rank();
      INTREPID2_TEST_FOR_EXCEPTION(outputRank > max_rank, std::invalid_argument, "Unsupported view combination.");
      rank_combinations_ = RankCombinationViewType("Rank_combinations_", max_rank);
      auto rank_combinations_host = Kokkos::create_mirror_view(rank_combinations_);
      
      rank_combinations_host[0] = TensorViewIteratorType::TENSOR_PRODUCT; // field combination is always tensor product
      rank_combinations_host[1] = tensorPoints ? TensorViewIteratorType::TENSOR_PRODUCT : TensorViewIteratorType::DIMENSION_MATCH; // tensorPoints controls interpretation of the point dimension
      for (ordinal_type d=2; d<max_rank; d++)
      {
        // d >= 2 have the interpretation of spatial dimensions (gradients, etc.)
        // we let the extents of the containers determine what we're doing here
        if ((inputValues1.extent_int(d) == inputValues2.extent_int(d)) && (output.extent_int(d) == 1))
        {
          rank_combinations_host[d] = TensorViewIteratorType::TENSOR_CONTRACTION;
        }
        else if (((inputValues1.extent_int(d) == output.extent_int(d)) && (inputValues2.extent_int(d) == 1))
                 || ((inputValues2.extent_int(d) == output.extent_int(d)) && (inputValues1.extent_int(d) == 1))
                 )
        {
          // this looks like multiplication of a vector by a scalar, resulting in a vector
          // this can be understood as a tensor product
          rank_combinations_host[d] = TensorViewIteratorType::TENSOR_PRODUCT;
        }
        else if ((inputValues1.extent_int(d) == inputValues2.extent_int(d)) && (output.extent_int(d) == inputValues1.extent_int(d) * inputValues2.extent_int(d)))
        {
          // this is actually a generalization of the above case: a tensor product, something like a vector outer product
          rank_combinations_host[d] = TensorViewIteratorType::TENSOR_PRODUCT;
        }
        else if ((inputValues1.extent_int(d) == inputValues2.extent_int(d)) && (output.extent_int(d) == inputValues1.extent_int(d)))
        {
          // it's a bit weird (I'm not aware of the use case, in the present context), but we can handle this case by adopting DIMENSION_MATCH here
          // this is something like MATLAB's .* and .+ operators, which operate entry-wise
          rank_combinations_host[d] = TensorViewIteratorType::DIMENSION_MATCH;
        }
        else
        {
          std::cout << "inputValues1.extent_int(" << d << ") = " << inputValues1.extent_int(d) << std::endl;
          std::cout << "inputValues2.extent_int(" << d << ") = " << inputValues2.extent_int(d) << std::endl;
          std::cout << "output.extent_int("       << d << ") = " << output.extent_int(d) << std::endl;
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "unable to find an interpretation for this combination of views");
        }
      }
      Kokkos::deep_copy(rank_combinations_,rank_combinations_host);
    }
    
    KOKKOS_INLINE_FUNCTION
    void operator()( const TeamMember & teamMember ) const
    {
      auto fieldOrdinal1 = teamMember.league_rank();
      
      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
        TensorViewIteratorType it(output_,input1_,input2_,rank_combinations_);
        const int FIELD_ORDINAL_DIMENSION = 0;
        it.setLocation({fieldOrdinal1,0,0,0,0,0,0},{fieldOrdinal2,0,0,0,0,0,0});
        int next_increment_rank = FIELD_ORDINAL_DIMENSION; // used to initialize prev_increment_rank at the start of the do/while loop.  Notionally, we last incremented in the field ordinal rank to get to the {fieldOrdinal1,0,0,0,0,0,0},{fieldOrdinal2,0,0,0,0,0,0} location.
        OutputScalar accumulator = 0;
        
        do
        {
          accumulator += weight_ * it.getView1Entry() * it.getView2Entry();
          next_increment_rank = it.nextIncrementRank();
          
          if ((next_increment_rank < 0) || (rank_combinations_[next_increment_rank] != TensorViewIteratorType::TENSOR_CONTRACTION))
          {
            // then we've finished the accumulation and should set the value
            it.set(accumulator);
            // reset the accumulator:
            accumulator = 0;
          }
        } while (it.increment() > FIELD_ORDINAL_DIMENSION);
      });
    }
  };
  
  /** \class  Intrepid2::Basis_TensorBasis
      \brief  Basis defined as the tensor product of two component bases.
   
   The cell topology for the tensor basis is the tensor product of the cell topologies on which the component bases are defined;
   \see Intrepid2::TensorTopologyMap.
   
   The basis is ordered such that the Basis1 field ordinals are the fastest-moving index; the formula for the composite field ordinal is:
     compositeFieldOrdinal = componentFieldOrdinal2 * basis1Cardinality + componentFieldOrdinal1
   This is done so that we can consider, e.g. Basis1 as the "x" dimension and Basis2 as the "y" dimension, and have the basis ordered in
   the same way that existing Intrepid2 bases on the quadrilateral are ordered, namely, one moves along the x dimension first, moving
   across the quadrilateral dofs "row-wise".
   
  */
  template<typename BasisBaseClass = void>
  class Basis_TensorBasis
  :
  public BasisBaseClass
  {
  public:
    using BasisBase = BasisBaseClass;
    using BasisPtr  = Teuchos::RCP<BasisBase>;
  
  protected:
    BasisPtr basis1_;
    BasisPtr basis2_;
    
    std::vector<BasisPtr> tensorComponents_;
    
    std::string name_; // name of the basis
    
    int numTensorialExtrusions_; // relative to cell topo returned by getBaseCellTopology().
  public:
    using DeviceType = typename BasisBase::DeviceType;
    using ExecutionSpace  = typename BasisBase::ExecutionSpace;
    using OutputValueType = typename BasisBase::OutputValueType;
    using PointValueType  = typename BasisBase::PointValueType;
    
    using OrdinalTypeArray1DHost = typename BasisBase::OrdinalTypeArray1DHost;
    using OrdinalTypeArray2DHost = typename BasisBase::OrdinalTypeArray2DHost;
    using OutputViewType         = typename BasisBase::OutputViewType;
    using PointViewType          = typename BasisBase::PointViewType;
    using TensorBasis            = Basis_TensorBasis<BasisBaseClass>;
  public:
    /** \brief  Constructor.
        \param [in] basis1 - the first component basis
        \param [in] basis2 - the second component basis
        \param [in] functionSpace - the function space to which the composite basis belongs (use FUNCTION_SPACE_MAX for unknown/unspecified function space)
        \param [in] useShardsCellTopologyAndTags - if true, attempt to assign a shards CellTopology corresponding to the tensor topology (shards Quad and Hex do not have tensor structure; this will map dofs appropriately) -- supported for 2D and 3D hypercubes
     */
    Basis_TensorBasis(BasisPtr basis1, BasisPtr basis2, EFunctionSpace functionSpace = FUNCTION_SPACE_MAX,
                      const bool useShardsCellTopologyAndTags = false)
    :
    basis1_(basis1),basis2_(basis2)
    {
      this->functionSpace_ = functionSpace;
      
      Basis_TensorBasis* basis1AsTensor = dynamic_cast<Basis_TensorBasis*>(basis1_.get());
      if (basis1AsTensor)
      {
        auto basis1Components = basis1AsTensor->getTensorBasisComponents();
        tensorComponents_.insert(tensorComponents_.end(), basis1Components.begin(), basis1Components.end());
      }
      else
      {
        tensorComponents_.push_back(basis1_);
      }
      
      Basis_TensorBasis* basis2AsTensor = dynamic_cast<Basis_TensorBasis*>(basis2_.get());
      if (basis2AsTensor)
      {
        auto basis2Components = basis2AsTensor->getTensorBasisComponents();
        tensorComponents_.insert(tensorComponents_.end(), basis2Components.begin(), basis2Components.end());
      }
      else
      {
        tensorComponents_.push_back(basis2_);
      }
      
      this->basisCardinality_  = basis1->getCardinality() * basis2->getCardinality();
      this->basisDegree_       = std::max(basis1->getDegree(), basis2->getDegree());
      
      {
        std::ostringstream basisName;
        basisName << basis1->getName() << " x " << basis2->getName();
        name_ = basisName.str();
      }
      
      // set cell topology
      this->basisCellTopologyKey_ = tensorComponents_[0]->getBaseCellTopology().getKey();
      this->numTensorialExtrusions_ = tensorComponents_.size() - 1;
      
      this->basisType_         = basis1_->getBasisType();
      this->basisCoordinates_  = COORDINATES_CARTESIAN;
      
      ordinal_type spaceDim1 = basis1_->getDomainDimension();
      ordinal_type spaceDim2 = basis2_->getDomainDimension();
      
      INTREPID2_TEST_FOR_EXCEPTION(spaceDim2 != 1, std::invalid_argument, "TensorBasis only supports 1D bases in basis2_ position");
      
      if (this->getBasisType() == BASIS_FEM_HIERARCHICAL)
      {
        // fill in degree lookup:
        int degreeSize = basis1_->getPolynomialDegreeLength() + basis2_->getPolynomialDegreeLength();
        this->fieldOrdinalPolynomialDegree_   = OrdinalTypeArray2DHost("TensorBasis - field ordinal polynomial degree", this->basisCardinality_, degreeSize);
        this->fieldOrdinalH1PolynomialDegree_ = OrdinalTypeArray2DHost("TensorBasis - field ordinal polynomial H^1 degree", this->basisCardinality_, degreeSize);
        
        const ordinal_type basis1Cardinality = basis1_->getCardinality();
        const ordinal_type basis2Cardinality = basis2_->getCardinality();
        
        int degreeLengthField1 = basis1_->getPolynomialDegreeLength();
        int degreeLengthField2 = basis2_->getPolynomialDegreeLength();
        
        for (ordinal_type fieldOrdinal1 = 0; fieldOrdinal1 < basis1Cardinality; fieldOrdinal1++)
        {
          OrdinalTypeArray1DHost degreesField1   = basis1_->getPolynomialDegreeOfField(fieldOrdinal1);
          OrdinalTypeArray1DHost h1DegreesField1 = basis1_->getH1PolynomialDegreeOfField(fieldOrdinal1);
          for (ordinal_type fieldOrdinal2 = 0; fieldOrdinal2 < basis2Cardinality; fieldOrdinal2++)
          {
            OrdinalTypeArray1DHost degreesField2   = basis2_->getPolynomialDegreeOfField(fieldOrdinal2);
            OrdinalTypeArray1DHost h1DegreesField2 = basis2_->getH1PolynomialDegreeOfField(fieldOrdinal2);
            const ordinal_type tensorFieldOrdinal = fieldOrdinal2 * basis1Cardinality + fieldOrdinal1;
            
            for (int d3=0; d3<degreeLengthField1; d3++)
            {
              this->fieldOrdinalPolynomialDegree_  (tensorFieldOrdinal,d3) =   degreesField1(d3);
              this->fieldOrdinalH1PolynomialDegree_(tensorFieldOrdinal,d3) = h1DegreesField1(d3);
            }
            for (int d3=0; d3<degreeLengthField2; d3++)
            {
              this->fieldOrdinalPolynomialDegree_  (tensorFieldOrdinal,d3+degreeLengthField1) =   degreesField2(d3);
              this->fieldOrdinalH1PolynomialDegree_(tensorFieldOrdinal,d3+degreeLengthField1) = h1DegreesField2(d3);
            }
          }
        }
      }
      
      if (useShardsCellTopologyAndTags)
      {
        setShardsTopologyAndTags();
      }
      else
      {
        // we build tags recursively, making reference to basis1_ and basis2_'s tags to produce the tensor product tags.
  //      // initialize tags
        const auto & cardinality = this->basisCardinality_;
  
        // Basis-dependent initializations
        const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
        const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
        const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
        const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
        const ordinal_type posDfCnt = 3;        // position in the tag, counting from 0, of DoF count for the subcell
  
        OrdinalTypeArray1DHost tagView("tag view", cardinality*tagSize);
  
        // we assume that basis2_ is defined on a line, and that basis1_ is defined on a domain that is once-extruded in by that line.
        auto cellTopo = CellTopology::cellTopology(tensorComponents_[0]->getBaseCellTopology(), numTensorialExtrusions_);
        auto basis1Topo = cellTopo->getTensorialComponent();
        
        const ordinal_type spaceDim = spaceDim1 + spaceDim2;
        const ordinal_type sideDim   = spaceDim - 1;
        
        const OrdinalTypeArray2DHost ordinalToTag1 = basis1_->getAllDofTags();
        const OrdinalTypeArray2DHost ordinalToTag2 = basis2_->getAllDofTags();
                
        for (int fieldOrdinal1=0; fieldOrdinal1<basis1_->getCardinality(); fieldOrdinal1++)
        {
          ordinal_type subcellDim1   = ordinalToTag1(fieldOrdinal1,posScDim);
          ordinal_type subcellOrd1   = ordinalToTag1(fieldOrdinal1,posScOrd);
          ordinal_type subcellDfCnt1 = ordinalToTag1(fieldOrdinal1,posDfCnt);
          for (int fieldOrdinal2=0; fieldOrdinal2<basis2_->getCardinality(); fieldOrdinal2++)
          {
            ordinal_type subcellDim2   = ordinalToTag2(fieldOrdinal2,posScDim);
            ordinal_type subcellOrd2   = ordinalToTag2(fieldOrdinal2,posScOrd);
            ordinal_type subcellDfCnt2 = ordinalToTag2(fieldOrdinal2,posDfCnt);
            
            ordinal_type subcellDim = subcellDim1 + subcellDim2;
            ordinal_type subcellOrd;
            if (subcellDim2 == 0)
            {
              // vertex node in extrusion; the subcell is not extruded but belongs to one of the two "copies"
              // of the basis1 topology
              ordinal_type sideOrdinal = cellTopo->getTensorialComponentSideOrdinal(subcellOrd2); // subcellOrd2 is a "side" of the line topology
              subcellOrd = CellTopology::getSubcellOrdinalMap(cellTopo, sideDim, sideOrdinal,
                                                              subcellDim1, subcellOrd1);
            }
            else
            {
              // line subcell in time; the subcell *is* extruded in final dimension
              subcellOrd = cellTopo->getExtrudedSubcellOrdinal(subcellDim1, subcellOrd1);
              if (subcellOrd == -1)
              {
                std::cout << "ERROR: -1 subcell ordinal.\n";
                subcellOrd = cellTopo->getExtrudedSubcellOrdinal(subcellDim1, subcellOrd1);
              }
            }
            ordinal_type tensorFieldOrdinal = fieldOrdinal2 * basis1_->getCardinality() + fieldOrdinal1;
      //        cout << "(" << fieldOrdinal1 << "," << fieldOrdinal2 << ") --> " << i << endl;
            ordinal_type dofOffsetOrdinal1 = ordinalToTag1(fieldOrdinal1,posDfOrd);
            ordinal_type dofOffsetOrdinal2 = ordinalToTag2(fieldOrdinal2,posDfOrd);
            ordinal_type dofsForSubcell1   = ordinalToTag1(fieldOrdinal1,posDfCnt);
            ordinal_type dofOffsetOrdinal  = dofOffsetOrdinal2 * dofsForSubcell1 + dofOffsetOrdinal1;
            tagView(tagSize*tensorFieldOrdinal + posScDim) = subcellDim; // subcellDim
            tagView(tagSize*tensorFieldOrdinal + posScOrd) = subcellOrd; // subcell ordinal
            tagView(tagSize*tensorFieldOrdinal + posDfOrd) = dofOffsetOrdinal;  // ordinal of the specified DoF relative to the subcell
            tagView(tagSize*tensorFieldOrdinal + posDfCnt) = subcellDfCnt1 * subcellDfCnt2; // total number of DoFs associated with the subcell
          }
        }
        
        //        // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
        //        // tags are constructed on host
        this->setOrdinalTagData(this->tagToOrdinal_,
                                this->ordinalToTag_,
                                tagView,
                                this->basisCardinality_,
                                tagSize,
                                posScDim,
                                posScOrd,
                                posDfOrd);
      }
    }
    
    void setShardsTopologyAndTags()
    {
      shards::CellTopology cellTopo1 = basis1_->getBaseCellTopology();
      shards::CellTopology cellTopo2 = basis2_->getBaseCellTopology();
      
      auto cellKey1 = basis1_->getBaseCellTopology().getKey();
      auto cellKey2 = basis2_->getBaseCellTopology().getKey();
      
      const int numTensorialExtrusions = basis1_->getNumTensorialExtrusions() + basis2_->getNumTensorialExtrusions();
      if ((cellKey1 == shards::Line<2>::key) && (cellKey2 == shards::Line<2>::key) && (numTensorialExtrusions == 0))
      {
        this->basisCellTopologyKey_ = shards::Quadrilateral<4>::key;
      }
      else if (   ((cellKey1 == shards::Quadrilateral<4>::key) && (cellKey2 == shards::Line<2>::key))
               || ((cellKey2 == shards::Quadrilateral<4>::key) && (cellKey1 == shards::Line<2>::key))
               || ((cellKey1 == shards::Line<2>::key) && (cellKey2 == shards::Line<2>::key) && (numTensorialExtrusions == 1))
              )
      {
        this->basisCellTopologyKey_ = shards::Hexahedron<8>::key;
      }
      else if ((cellKey1 == shards::Triangle<3>::key) && (cellKey2 == shards::Line<2>::key))
      {
        this->basisCellTopologyKey_ = shards::Wedge<6>::key;
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Cell topology combination not yet supported");
      }
      
      // numTensorialExtrusions_ is relative to the baseCellTopology; what we've just done is found a cell topology of the same spatial dimension as the extruded topology, so now numTensorialExtrusions_ should be 0.
      numTensorialExtrusions_ = 0;
      
      // initialize tags
      {
        const auto & cardinality = this->basisCardinality_;
        
        // Basis-dependent initializations
        const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
        const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
        const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
        const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
        
        OrdinalTypeArray1DHost tagView("tag view", cardinality*tagSize);
        
        shards::CellTopology cellTopo = this->getBaseCellTopology();
        
        ordinal_type tensorSpaceDim  = cellTopo.getDimension();
        ordinal_type spaceDim1       = cellTopo1.getDimension();
        ordinal_type spaceDim2       = cellTopo2.getDimension();
        
        TensorTopologyMap topoMap(cellTopo1, cellTopo2);
        
        for (ordinal_type d=0; d<=tensorSpaceDim; d++) // d: tensorial dimension
        {
          ordinal_type d2_max = std::min(spaceDim2, d);
          for (ordinal_type d2=0; d2<=d2_max; d2++)
          {
            ordinal_type d1 = d-d2;
            if (d1 > spaceDim1) continue;
            
            ordinal_type subcellCount2 = cellTopo2.getSubcellCount(d2);
            ordinal_type subcellCount1 = cellTopo1.getSubcellCount(d1);
            for (ordinal_type subcellOrdinal2=0; subcellOrdinal2<subcellCount2; subcellOrdinal2++)
            {
              ordinal_type subcellDofCount2 = basis2_->getDofCount(d2, subcellOrdinal2);
              for (ordinal_type subcellOrdinal1=0; subcellOrdinal1<subcellCount1; subcellOrdinal1++)
              {
                ordinal_type subcellDofCount1 = basis1_->getDofCount(d1, subcellOrdinal1);
                ordinal_type tensorLocalDofCount = subcellDofCount1 * subcellDofCount2;
                for (ordinal_type localDofID2 = 0; localDofID2<subcellDofCount2; localDofID2++)
                {
                  ordinal_type fieldOrdinal2 = basis2_->getDofOrdinal(d2, subcellOrdinal2, localDofID2);
                  OrdinalTypeArray1DHost degreesField2;
                  if (this->basisType_ == BASIS_FEM_HIERARCHICAL) degreesField2 = basis2_->getPolynomialDegreeOfField(fieldOrdinal2);
                  for (ordinal_type localDofID1 = 0; localDofID1<subcellDofCount1; localDofID1++)
                  {
                    ordinal_type fieldOrdinal1 = basis1_->getDofOrdinal(d1, subcellOrdinal1, localDofID1);
                    ordinal_type tensorLocalDofID = localDofID2 * subcellDofCount1 + localDofID1;
                    ordinal_type tensorFieldOrdinal = fieldOrdinal2 * basis1_->getCardinality() + fieldOrdinal1;
                    tagView(tensorFieldOrdinal*tagSize+0) = d; // subcell dimension
                    tagView(tensorFieldOrdinal*tagSize+1) = topoMap.getCompositeSubcellOrdinal(d1, subcellOrdinal1, d2, subcellOrdinal2);
                    tagView(tensorFieldOrdinal*tagSize+2) = tensorLocalDofID;
                    tagView(tensorFieldOrdinal*tagSize+3) = tensorLocalDofCount;
                  } // localDofID1
                } // localDofID2
              } // subcellOrdinal1
            }   // subcellOrdinal2
          }
        }
        
        //        // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
        //        // tags are constructed on host
        this->setOrdinalTagData(this->tagToOrdinal_,
                                this->ordinalToTag_,
                                tagView,
                                this->basisCardinality_,
                                tagSize,
                                posScDim,
                                posScOrd,
                                posDfOrd);
      }
    }
    
    virtual int getNumTensorialExtrusions() const override
    {
      return numTensorialExtrusions_;
    }
    
    /** \brief  Given "Dk" enumeration indices for the component bases, returns a Dk enumeration index for the composite basis.
        \param [in] dkEnum1         - Dk enumeration index for first component basis
        \param [in] operatorOrder1  - operator order for the first component basis
        \param [in] dkEnum2         - Dk enumeration index for second component basis
        \param [in] operatorOrder2  - operator order for the second component basis
     
        \return Dk enumeration index for the composite basis, corresponding to operator order operatorOrder1 + operatorOrder2.
     */
    ordinal_type getTensorDkEnumeration(ordinal_type dkEnum1, ordinal_type operatorOrder1,
                                        ordinal_type dkEnum2, ordinal_type operatorOrder2) const
    {
      ordinal_type spaceDim1 = basis1_->getDomainDimension();
      ordinal_type spaceDim2 = basis2_->getDomainDimension();
      
      // We support total spaceDim <= 7.
      switch (spaceDim1)
      {
        case 0:
        {
          INTREPID2_TEST_FOR_EXCEPTION(operatorOrder1 > 0, std::invalid_argument, "For spaceDim1 = 0, operatorOrder1 must be 0.");
          return dkEnum2;
        }
        case 1:
          switch (spaceDim2)
        {
          case 1: return getDkTensorIndex<1, 1>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 2: return getDkTensorIndex<1, 2>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 3: return getDkTensorIndex<1, 3>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 4: return getDkTensorIndex<1, 4>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 5: return getDkTensorIndex<1, 5>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 6: return getDkTensorIndex<1, 6>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported dimension combination");
        }
        case 2:
          switch (spaceDim2)
        {
          case 1: return getDkTensorIndex<2, 1>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 2: return getDkTensorIndex<2, 2>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 3: return getDkTensorIndex<2, 3>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 4: return getDkTensorIndex<2, 4>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 5: return getDkTensorIndex<2, 5>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported dimension combination");
        }
        case 3:
          switch (spaceDim2)
        {
          case 1: return getDkTensorIndex<3, 1>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 2: return getDkTensorIndex<3, 2>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 3: return getDkTensorIndex<3, 3>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 4: return getDkTensorIndex<3, 4>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported dimension combination");
        }
        case 4:
          switch (spaceDim2)
        {
          case 1: return getDkTensorIndex<4, 1>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 2: return getDkTensorIndex<4, 2>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 3: return getDkTensorIndex<4, 3>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported dimension combination");
        }
        case 5:
          switch (spaceDim2)
        {
          case 1: return getDkTensorIndex<5, 1>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 2: return getDkTensorIndex<5, 2>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported dimension combination");
        }
        case 6:
          switch (spaceDim2)
        {
          case 1: return getDkTensorIndex<6, 1>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported dimension combination");
        }
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported dimension combination");
      }
    }
    
    /** \brief Returns a simple decomposition of the specified operator: what operator(s) should be applied to basis1, and what operator(s) to basis2.  A one-element OperatorTensorDecomposition corresponds to a single TensorData entry; a multiple-element OperatorTensorDecomposition corresponds to a VectorData object with axialComponents = false.
     
     Subclasses must override this method.
    */
    virtual OperatorTensorDecomposition getSimpleOperatorDecomposition(const EOperator &operatorType) const
    {
      const int spaceDim  = this->getDomainDimension();
      
      const EOperator VALUE = Intrepid2::OPERATOR_VALUE;
      
      std::vector< std::vector<EOperator> > opsVALUE{{VALUE, VALUE}};
      
      std::vector< std::vector<EOperator> > ops(spaceDim);
      
      switch (operatorType)
      {
        case VALUE:
          ops = opsVALUE;
          break;
        case OPERATOR_DIV:
        case OPERATOR_CURL:
          // DIV and CURL are multi-family bases; subclasses are required to override
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported operator type - TensorBasis subclass should override");
          break;
        case OPERATOR_GRAD:
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
        case OPERATOR_Dn:
        {
          auto opOrder = getOperatorOrder(operatorType); // number of derivatives that we take in total
          const int dkCardinality = getDkCardinality(operatorType, 2); // 2 because we have two tensor component bases, basis1_ and basis2_
          
          ops = std::vector< std::vector<EOperator> >(dkCardinality);
          
          // the Dk enumeration happens in lexicographic order (reading from left to right: x, y, z, etc.)
          // this governs the nesting order of the dkEnum1, dkEnum2 for loops below: dkEnum2 should increment fastest.
          for (int derivativeCountComp2=0; derivativeCountComp2<=opOrder; derivativeCountComp2++)
          {
            int derivativeCountComp1=opOrder-derivativeCountComp2;
            EOperator op1 = (derivativeCountComp1 == 0) ? OPERATOR_VALUE : EOperator(OPERATOR_D1 + (derivativeCountComp1 - 1));
            EOperator op2 = (derivativeCountComp2 == 0) ? OPERATOR_VALUE : EOperator(OPERATOR_D1 + (derivativeCountComp2 - 1));
            
            int dkCardinality1 = getDkCardinality(op1, 1); // use dim = 1 because this is a "simple" decomposition -- full decomposition will expand within the dimensions of basis1_
            int dkCardinality2 = getDkCardinality(op2, 1); // use dim = 1 because this is a "simple" decomposition -- full decomposition will expand within the dimensions of basis2_
            
            for (int dkEnum1=0; dkEnum1<dkCardinality1; dkEnum1++)
            {
              for (int dkEnum2=0; dkEnum2<dkCardinality2; dkEnum2++)
              {
                ordinal_type dkTensorIndex = getDkTensorIndex<1, 1>(dkEnum1, derivativeCountComp1, dkEnum2, derivativeCountComp2);
                ops[dkTensorIndex] = std::vector<EOperator>{op1, op2};
              }
            }
          }
        }
          break;
      }
      
      std::vector<double> weights(ops.size(), 1.0);
      return OperatorTensorDecomposition(ops, weights);
    }
    
    /** \brief Returns a full decomposition of the specified operator.  (Full meaning that all TensorBasis components are expanded into their non-TensorBasis components.)
      */
    virtual OperatorTensorDecomposition getOperatorDecomposition(const EOperator operatorType) const
    {
      if (((operatorType >= OPERATOR_D1) && (operatorType <= OPERATOR_D10)) || (operatorType == OPERATOR_GRAD))
      {
        // ordering of the operators is reverse-lexicographic, reading left to right (highest-dimension is fastest-moving).
        // first entry will be (operatorType, VALUE, …, VALUE)
        // next will be (operatorType - 1, OP_D1, VALUE, …, VALUE)
        // then         (operatorType - 1, VALUE, OP_D1, …, VALUE)
        
        ordinal_type numBasisComponents = tensorComponents_.size();
        
        auto opOrder = getOperatorOrder(operatorType); // number of derivatives that we take in total
        const int dkCardinality = getDkCardinality(operatorType, numBasisComponents);
        
        std::vector< std::vector<EOperator> > ops(dkCardinality);
        
        std::vector<EOperator> prevEntry(numBasisComponents, OPERATOR_VALUE);
        prevEntry[0] = operatorType;
        
        ops[0] = prevEntry;
        
        for (ordinal_type dkOrdinal=1; dkOrdinal<dkCardinality; dkOrdinal++)
        {
          std::vector<EOperator> entry = prevEntry;
          
          // decrement to follow reverse lexicographic ordering:
          /*
           How to tell when it is time to decrement the nth entry:
           1. Let a be the sum of the opOrders for entries 0 through n-1.
           2. Let b be the sum of the nth entry and the final entry.
           3. If opOrder == a + b, then the nth entry should be decremented.
           */
          ordinal_type cumulativeOpOrder = 0;
          ordinal_type finalOpOrder = getOperatorOrder(entry[numBasisComponents-1]);
          for (ordinal_type compOrdinal=0; compOrdinal<numBasisComponents; compOrdinal++)
          {
            const ordinal_type thisOpOrder = getOperatorOrder(entry[compOrdinal]);
            cumulativeOpOrder += thisOpOrder;
            if (cumulativeOpOrder + finalOpOrder == opOrder)
            {
              // decrement this
              EOperator decrementedOp;
              if (thisOpOrder == 1)
              {
                decrementedOp = OPERATOR_VALUE;
              }
              else
              {
                decrementedOp = static_cast<EOperator>(OPERATOR_D1 + ((thisOpOrder - 1) - 1));
              }
              entry[compOrdinal]   = decrementedOp;
              const ordinal_type remainingOpOrder = opOrder - cumulativeOpOrder + 1;
              entry[compOrdinal+1] = static_cast<EOperator>(OPERATOR_D1 + (remainingOpOrder - 1));
              for (ordinal_type i=compOrdinal+2; i<numBasisComponents; i++)
              {
                entry[i] = OPERATOR_VALUE;
              }
              break;
            }
          }
          ops[dkOrdinal] = entry;
          prevEntry = entry;
        }
        std::vector<double> weights(dkCardinality, 1.0);
        
        return OperatorTensorDecomposition(ops, weights);
      }
      else
      {
        OperatorTensorDecomposition opSimpleDecomposition = this->getSimpleOperatorDecomposition(operatorType);
        std::vector<BasisPtr> componentBases {basis1_, basis2_};
        return opSimpleDecomposition.expandedDecomposition(componentBases);
      }
    }
    
    /** \brief Allocate BasisValues container suitable for passing to the getValues() variant that takes a TensorPoints container as argument.
     
        The basic exact-sequence operators are supported (VALUE, GRAD, DIV, CURL), as are the Dn operators (OPERATOR_D1 through OPERATOR_D10).
     */
    virtual BasisValues<OutputValueType,DeviceType> allocateBasisValues( TensorPoints<PointValueType,DeviceType> points, const EOperator operatorType = OPERATOR_VALUE) const override
    {
      const bool operatorIsDk = (operatorType >= OPERATOR_D1) && (operatorType <= OPERATOR_D10);
      const bool operatorSupported = (operatorType == OPERATOR_VALUE) || (operatorType == OPERATOR_GRAD) || (operatorType == OPERATOR_CURL) || (operatorType == OPERATOR_DIV) || operatorIsDk;
      INTREPID2_TEST_FOR_EXCEPTION(!operatorSupported, std::invalid_argument, "operator is not supported by allocateBasisValues");
      
      // check that points's spatial dimension matches the basis
      const int spaceDim = this->getDomainDimension();
      INTREPID2_TEST_FOR_EXCEPTION(spaceDim != points.extent_int(1), std::invalid_argument, "points must be shape (P,D), with D equal to the dimension of the basis domain");
      
      // check that points has enough tensor components
      ordinal_type numBasisComponents = tensorComponents_.size();
      if (numBasisComponents > points.numTensorComponents())
      {
        // Then we require points to have a trivial tensor structure.  (Subclasses could be more sophisticated.)
        // (More sophisticated approaches are possible here, too, but likely the most common use case in which there is not a one-to-one correspondence
        //  between basis components and point components will involve trivial tensor structure in the points...)
        INTREPID2_TEST_FOR_EXCEPTION(points.numTensorComponents() != 1, std::invalid_argument, "If points does not have the same number of tensor components as the basis, then it should have trivial tensor structure.");
        const ordinal_type numPoints = points.extent_int(0);
        auto outputView = this->allocateOutputView(numPoints, operatorType);
        
        Data<OutputValueType,DeviceType> outputData(outputView);
        TensorData<OutputValueType,DeviceType> outputTensorData(outputData);
        
        return BasisValues<OutputValueType,DeviceType>(outputTensorData);
      }
      INTREPID2_TEST_FOR_EXCEPTION(numBasisComponents > points.numTensorComponents(), std::invalid_argument, "points must have at least as many tensorial components as basis.");
      
      OperatorTensorDecomposition opDecomposition = getOperatorDecomposition(operatorType);
            
      ordinal_type numVectorComponents = opDecomposition.numVectorComponents();
      const bool useVectorData = numVectorComponents > 1;
      
      std::vector<ordinal_type> componentPointCounts(numBasisComponents);
      ordinal_type pointComponentNumber = 0;
      for (ordinal_type r=0; r<numBasisComponents; r++)
      {
        const ordinal_type compSpaceDim = tensorComponents_[r]->getDomainDimension();
        ordinal_type dimsSoFar = 0;
        ordinal_type numPointsForBasisComponent = 1;
        while (dimsSoFar < compSpaceDim)
        {
          INTREPID2_TEST_FOR_EXCEPTION(pointComponentNumber >= points.numTensorComponents(), std::invalid_argument, "Error in processing points container; perhaps it is mis-sized?");
          const int numComponentPoints = points.componentPointCount(pointComponentNumber);
          const int numComponentDims = points.getTensorComponent(pointComponentNumber).extent_int(1);
          numPointsForBasisComponent *= numComponentPoints;
          dimsSoFar += numComponentDims;
          INTREPID2_TEST_FOR_EXCEPTION(dimsSoFar > points.numTensorComponents(), std::invalid_argument, "Error in processing points container; perhaps it is mis-sized?");
          pointComponentNumber++;
        }
        componentPointCounts[r] = numPointsForBasisComponent;
      }
      
      if (useVectorData)
      {
        const int numFamilies = 1;
        std::vector< std::vector<TensorData<OutputValueType,DeviceType> > > vectorComponents(numFamilies, std::vector<TensorData<OutputValueType,DeviceType> >(numVectorComponents));
        
        const int familyOrdinal = 0;
        for (ordinal_type vectorComponentOrdinal=0; vectorComponentOrdinal<numVectorComponents; vectorComponentOrdinal++)
        {
          if (!opDecomposition.identicallyZeroComponent(vectorComponentOrdinal))
          {
            std::vector< Data<OutputValueType,DeviceType> > componentData;
            for (ordinal_type r=0; r<numBasisComponents; r++)
            {
              const int numComponentPoints = componentPointCounts[r];
              const EOperator op = opDecomposition.op(vectorComponentOrdinal, r);
              auto componentView = tensorComponents_[r]->allocateOutputView(numComponentPoints, op);
              componentData.push_back(Data<OutputValueType,DeviceType>(componentView));
            }
            vectorComponents[familyOrdinal][vectorComponentOrdinal] = TensorData<OutputValueType,DeviceType>(componentData);
          }
        }
        VectorData<OutputValueType,DeviceType> vectorData(vectorComponents);
        return BasisValues<OutputValueType,DeviceType>(vectorData);
      }
      else
      {
        // TensorData: single tensor product
        std::vector< Data<OutputValueType,DeviceType> > componentData;
        
        const ordinal_type vectorComponentOrdinal = 0;
        for (ordinal_type r=0; r<numBasisComponents; r++)
        {
          const int numComponentPoints = componentPointCounts[r];
          const EOperator op = opDecomposition.op(vectorComponentOrdinal, r);
          auto componentView = tensorComponents_[r]->allocateOutputView(numComponentPoints, op);
          
          const int rank = 2; // (F,P) -- TensorData-only BasisValues are always scalar-valued.  Use VectorData for vector-valued.
          // (we need to be explicit about the rank argument because GRAD, even in 1D, elevates to rank 3), so e.g. DIV of HDIV uses a componentView that is rank 3;
          //  we want Data to insulate us from that fact)
          const Kokkos::Array<int,7> extents {componentView.extent_int(0), componentView.extent_int(1), 1,1,1,1,1};
          Kokkos::Array<DataVariationType,7> variationType {GENERAL, GENERAL, CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT };
          componentData.push_back(Data<OutputValueType,DeviceType>(componentView, rank, extents, variationType));
        }
        
        TensorData<OutputValueType,DeviceType> tensorData(componentData);
        
        std::vector< TensorData<OutputValueType,DeviceType> > tensorDataEntries {tensorData};
        return BasisValues<OutputValueType,DeviceType>(tensorDataEntries);
      }
    }
    
    // since the getValues() below only overrides the FEM variant, we specify that
    // we use the base class's getValues(), which implements the FVD variant by throwing an exception.
    // (It's an error to use the FVD variant on this basis.)
    using BasisBase::getValues;
    
    /** \brief  Method to extract component points from composite points.
        \param [in]  inputPoints                  - points defined on the composite cell topology
        \param [in]  attemptTensorDecomposition   - if true, attempt to find a tensor decomposition.
        \param [out] inputPoints1                 - points defined on the first component cell topology
        \param [out] inputPoints2                 - points defined on the second component cell topology
        \param [out] tensorDecompositionSucceeded - if true, the attempt to find a tensor decomposition succeeded.
     
     At present, attemptTensorDecomposition is ignored, and tensorDecompositionSucceeded will always return false.
     However, we intend to support the tensor decomposition in the future, which will allow substantial optimizations
     in computation of tensor bases.
     */
    void getComponentPoints(const PointViewType inputPoints, const bool attemptTensorDecomposition,
                            PointViewType & inputPoints1, PointViewType & inputPoints2, bool &tensorDecompositionSucceeded) const
    {
      INTREPID2_TEST_FOR_EXCEPTION(attemptTensorDecomposition, std::invalid_argument, "tensor decomposition not yet supported");
      
      // for inputPoints that are actually tensor-product of component quadrature points (say),
      // having just the one input (which will have a lot of redundant point data) is suboptimal
      // The general case can have unique x/y/z coordinates at every point, though, so we have to support that
      // when this interface is used.  But we may try detecting that the data is tensor-product and compressing
      // from there...  Ultimately, we should also add a getValues() variant that takes multiple input point containers,
      // one for each tensorial dimension.
      
      // this initial implementation is intended to simplify development of 2D and 3D bases, while also opening
      // the possibility of higher-dimensional bases.  It is not necessarily optimized for speed/memory.  There
      // are things we can do in this regard, which may become important for matrix-free computations wherein
      // basis values don't get stored but are computed dynamically.
      
      int spaceDim1 = basis1_->getDomainDimension();
      int spaceDim2 = basis2_->getDomainDimension();
      
      int totalSpaceDim   = inputPoints.extent_int(1);
      
      TEUCHOS_ASSERT(spaceDim1 + spaceDim2 == totalSpaceDim);
      
      // first pass: just take subviews to get input points -- this will result in redundant computations when points are themselves tensor product (i.e., inputPoints itself contains redundant data)
      
      inputPoints1 = Kokkos::subview(inputPoints,Kokkos::ALL(),std::make_pair(0,spaceDim1));
      inputPoints2 = Kokkos::subview(inputPoints,Kokkos::ALL(),std::make_pair(spaceDim1,totalSpaceDim));
      
      //      std::cout << "inputPoints : " << inputPoints.extent(0) << " x " << inputPoints.extent(1) << std::endl;
      //      std::cout << "inputPoints1 : " << inputPoints1.extent(0) << " x " << inputPoints1.extent(1) << std::endl;
      //      std::cout << "inputPoints2 : " << inputPoints2.extent(0) << " x " << inputPoints2.extent(1) << std::endl;
      
      tensorDecompositionSucceeded = false;
    }
    
    /** \brief  Fills in spatial locations (coordinates) of degrees of freedom (nodes) on the reference cell
        \param [out] dofCoords - the container into which to place the degrees of freedom.
     
     dofCoords should have shape (F,D), where the field dimension matches the cardinality of the basis, and D is the
     spatial dimension of the topology on which the basis is defined.
     
     Note that getDofCoords() is not supported by all bases; in particular, hierarchical bases do not generally support this.
     */
    virtual void getDofCoords( typename BasisBase::ScalarViewType dofCoords ) const override
    {
      int spaceDim1 = basis1_->getBaseCellTopology().getDimension();
      int spaceDim2 = basis2_->getBaseCellTopology().getDimension();
      
      using ValueType    = typename BasisBase::ScalarViewType::value_type;
      using ResultLayout = typename DeduceLayout< typename BasisBase::ScalarViewType >::result_layout;
      using ViewType     = Kokkos::DynRankView<ValueType, ResultLayout, DeviceType >;
      
      const ordinal_type basisCardinality1 = basis1_->getCardinality();
      const ordinal_type basisCardinality2 = basis2_->getCardinality();

      ViewType dofCoords1("dofCoords1",basisCardinality1,spaceDim1);
      ViewType dofCoords2("dofCoords2",basisCardinality2,spaceDim2);
      
      basis1_->getDofCoords(dofCoords1);
      basis2_->getDofCoords(dofCoords2);
      
      Kokkos::RangePolicy<ExecutionSpace> policy(0, basisCardinality2);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const int fieldOrdinal2)
       {
         for (int fieldOrdinal1=0; fieldOrdinal1<basisCardinality1; fieldOrdinal1++)
         {
           const ordinal_type fieldOrdinal = fieldOrdinal1 + fieldOrdinal2 * basisCardinality1;
           for (int d1=0; d1<spaceDim1; d1++)
           {
             dofCoords(fieldOrdinal,d1) = dofCoords1(fieldOrdinal1,d1);
           }
           for (int d2=0; d2<spaceDim2; d2++)
           {
             dofCoords(fieldOrdinal,spaceDim1+d2) = dofCoords2(fieldOrdinal2,d2);
           }
         }
       });
    }
    

    /** \brief  Fills in coefficients of degrees of freedom on the reference cell
        \param [out] dofCoeffs - the container into which to place the degrees of freedom.

     dofCoeffs has the rank 1 or 2 and the first dimension equals the cardinality of the basis.
     dofCoeffs has rank 2 if either basis1 or basis2 is a vector basis, in which case the second dimension equals the number of vector components
     We do not support the case where both basis1 and basis2 are vector bases

     Note that getDofCoeffs() is not supported by all bases; in particular, hierarchical bases do not generally support this.
     */
    virtual void getDofCoeffs( typename BasisBase::ScalarViewType dofCoeffs ) const override
    {
      using ValueType    = typename BasisBase::ScalarViewType::value_type;
      using ResultLayout = typename DeduceLayout< typename BasisBase::ScalarViewType >::result_layout;
      using ViewType     = Kokkos::DynRankView<ValueType, ResultLayout, DeviceType >;

      const ordinal_type basisCardinality1 = basis1_->getCardinality();
      const ordinal_type basisCardinality2 = basis2_->getCardinality();

      bool isVectorBasis1 = getFieldRank(basis1_->getFunctionSpace()) == 1;
      bool isVectorBasis2 = getFieldRank(basis2_->getFunctionSpace()) == 1;

      INTREPID2_TEST_FOR_EXCEPTION(isVectorBasis1 && isVectorBasis2, std::invalid_argument, "the case in which basis1 and basis2 are vector bases is not supported");

      int basisDim1 = isVectorBasis1 ? basis1_->getBaseCellTopology().getDimension() : 1;
      int basisDim2 = isVectorBasis2 ? basis2_->getBaseCellTopology().getDimension() : 1;

      auto dofCoeffs1 = isVectorBasis1 ? ViewType("dofCoeffs1",basis1_->getCardinality(), basisDim1) : ViewType("dofCoeffs1",basis1_->getCardinality());
      auto dofCoeffs2 = isVectorBasis2 ? ViewType("dofCoeffs2",basis2_->getCardinality(), basisDim2) : ViewType("dofCoeffs2",basis2_->getCardinality());
      
      basis1_->getDofCoeffs(dofCoeffs1);
      basis2_->getDofCoeffs(dofCoeffs2);

      Kokkos::RangePolicy<ExecutionSpace> policy(0, basisCardinality2);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const int fieldOrdinal2)
       {
         for (int fieldOrdinal1=0; fieldOrdinal1<basisCardinality1; fieldOrdinal1++)
         {
           const ordinal_type fieldOrdinal = fieldOrdinal1 + fieldOrdinal2 * basisCardinality1;
           for (int d1 = 0; d1 <basisDim1; d1++) {
             for (int d2 = 0; d2 <basisDim2; d2++) {
               dofCoeffs.access(fieldOrdinal,d1+d2)  = dofCoeffs1.access(fieldOrdinal1,d1);
               dofCoeffs.access(fieldOrdinal,d1+d2) *= dofCoeffs2.access(fieldOrdinal2,d2);
             }
           }
         }
       });
    }

    /** \brief  Returns basis name
     
     \return the name of the basis
     */
    virtual
    const char*
    getName() const override {
      return name_.c_str();
    }
    
    std::vector<BasisPtr> getTensorBasisComponents() const
    {
      return tensorComponents_;
    }
    
    /** \brief  Evaluation of a FEM basis on a <strong>reference cell</strong>, using point and output value containers that allow preservation of tensor-product structure.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values.  Should be allocated using Basis::allocateBasisValues().
        \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points.  This should be allocated using Cubature::allocateCubaturePoints() and filled using Cubature::getCubature().
        \param  operatorType      [in]  - the operator acting on the basis function
     
        This is the preferred getValues() method for TensorBasis and DirectSumBasis and their subclasses.  It allows a reduced memory footprint and optimized integration, etc.
    */
    virtual
    void
    getValues(       BasisValues<OutputValueType,DeviceType> outputValues,
               const TensorPoints<PointValueType,DeviceType>  inputPoints,
               const EOperator operatorType = OPERATOR_VALUE ) const override
    {
      const ordinal_type numTensorComponents = tensorComponents_.size();
      if (inputPoints.numTensorComponents() < numTensorComponents)
      {
        // then we require that both inputPoints and outputValues trivial tensor structure
        INTREPID2_TEST_FOR_EXCEPTION( inputPoints.numTensorComponents() != 1, std::invalid_argument, "If inputPoints differs from the tensor basis in component count, then inputPoints must have trivial tensor product structure" );
        INTREPID2_TEST_FOR_EXCEPTION( outputValues.numFamilies() != 1, std::invalid_argument, "If inputPoints differs from the tensor basis in component count, outputValues must have a single family with trivial tensor product structure" );
        INTREPID2_TEST_FOR_EXCEPTION( outputValues.tensorData().numTensorComponents() != 1, std::invalid_argument, "If inputPoints differs from the tensor basis in component count, outputValues must have a single family with trivial tensor product structure" );
        
        OutputViewType outputView = outputValues.tensorData().getTensorComponent(0).getUnderlyingView();
        PointViewType   pointView = inputPoints.getTensorComponent(0);
        this->getValues(outputView, pointView, operatorType);
        return;
      }
      
      OperatorTensorDecomposition operatorDecomposition = getOperatorDecomposition(operatorType);
      
      const ordinal_type numVectorComponents = operatorDecomposition.numVectorComponents();
      const bool               useVectorData = numVectorComponents > 1;
      const ordinal_type  numBasisComponents = operatorDecomposition.numBasisComponents();
      
      for (ordinal_type vectorComponentOrdinal=0; vectorComponentOrdinal<numVectorComponents; vectorComponentOrdinal++)
      {
        const double weight = operatorDecomposition.weight(vectorComponentOrdinal);
        ordinal_type pointComponentOrdinal = 0;
        for (ordinal_type basisOrdinal=0; basisOrdinal<numBasisComponents; basisOrdinal++, pointComponentOrdinal++)
        {
          const EOperator op = operatorDecomposition.op(vectorComponentOrdinal, basisOrdinal);
          // by convention, op == OPERATOR_MAX signals a zero component; skip
          if (op != OPERATOR_MAX)
          {
            const int vectorFamily = 0; // TensorBasis always has just a single family; multiple families arise in DirectSumBasis
            auto tensorData = useVectorData ? outputValues.vectorData().getComponent(vectorFamily,vectorComponentOrdinal) : outputValues.tensorData();
            INTREPID2_TEST_FOR_EXCEPTION( ! tensorData.getTensorComponent(basisOrdinal).isValid(), std::invalid_argument, "Invalid output component encountered");
                        
            const Data<OutputValueType,DeviceType> & outputData = tensorData.getTensorComponent(basisOrdinal);
            
            auto basisValueView = outputData.getUnderlyingView();
            PointViewType  pointView = inputPoints.getTensorComponent(pointComponentOrdinal);
            const ordinal_type basisDomainDimension = tensorComponents_[basisOrdinal]->getDomainDimension();
            if (pointView.extent_int(1) == basisDomainDimension)
            {
              tensorComponents_[basisOrdinal]->getValues(basisValueView, pointView, op);
            }
            else
            {
              // we need to wrap the basisValueView in a BasisValues container, and to wrap the point components in a TensorPoints container.
              
              // combine point components to build up to basisDomainDimension
              ordinal_type dimsSoFar = 0;
              std::vector< ScalarView<PointValueType,DeviceType> > basisPointComponents;
              while (dimsSoFar < basisDomainDimension)
              {
                INTREPID2_TEST_FOR_EXCEPTION(pointComponentOrdinal >= inputPoints.numTensorComponents(), std::invalid_argument, "Error in processing points container; perhaps it is mis-sized?");
                const auto & pointComponent = inputPoints.getTensorComponent(pointComponentOrdinal);
                const ordinal_type numComponentDims   = pointComponent.extent_int(1);
                dimsSoFar += numComponentDims;
                INTREPID2_TEST_FOR_EXCEPTION(dimsSoFar > inputPoints.numTensorComponents(), std::invalid_argument, "Error in processing points container; perhaps it is mis-sized?");
                basisPointComponents.push_back(pointComponent);
                if (dimsSoFar < basisDomainDimension)
                {
                  // we will pass through this loop again, so we should increment the point component ordinal
                  pointComponentOrdinal++;
                }
              }
              
              TensorPoints<PointValueType, DeviceType> basisPoints(basisPointComponents);
              
              bool useVectorData2 = (basisValueView.rank() == 3);

              BasisValues<OutputValueType,DeviceType> basisValues;
              if (useVectorData2)
              {
                VectorData<OutputValueType,DeviceType> vectorData(outputData);
                basisValues = BasisValues<OutputValueType,DeviceType>(vectorData);
              }
              else
              {
                TensorData<OutputValueType,DeviceType> tensorData2(outputData);
                basisValues = BasisValues<OutputValueType,DeviceType>(tensorData2);
              }
              
              tensorComponents_[basisOrdinal]->getValues(basisValues, basisPoints, op);
            }
            
            // op.rotateXYNinetyDegrees() is set to true for one of the H(curl) wedge families
            // (due to the fact that Intrepid2::EOperator does not allow us to extract individual vector components
            //  via, e.g., OPERATOR_X, OPERATOR_Y, etc., we don't have a way of expressing the decomposition
            //  just in terms of EOperator and component-wise scalar weights; we could also do this via component-wise
            //  matrix weights, but this would involve a more intrusive change to the implementation).
            const bool spansXY = (vectorComponentOrdinal == 0) && (basisValueView.extent_int(2) == 2);
            if (spansXY && operatorDecomposition.rotateXYNinetyDegrees())
            {
              // map from (f_x,f_y) --> (-f_y,f_x)
              auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{basisValueView.extent_int(0),basisValueView.extent_int(1)});
              Kokkos::parallel_for("rotateXYNinetyDegrees", policy,
              KOKKOS_LAMBDA (const int &fieldOrdinal, const int &pointOrdinal) {
                const auto  f_x = basisValueView(fieldOrdinal,pointOrdinal,0); // copy
                const auto &f_y = basisValueView(fieldOrdinal,pointOrdinal,1); // reference
                basisValueView(fieldOrdinal,pointOrdinal,0) = -f_y;
                basisValueView(fieldOrdinal,pointOrdinal,1) =  f_x;
              });
            }
            
            // if weight is non-trivial (not 1.0), then we need to multiply one of the component views by weight.
            // we do that for the first basisOrdinal's values
            if ((weight != 1.0) && (basisOrdinal == 0))
            {
              if (basisValueView.rank() == 2)
              {
                auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{basisValueView.extent_int(0),basisValueView.extent_int(1)});
                Kokkos::parallel_for("multiply basisValueView by weight", policy,
                KOKKOS_LAMBDA (const int &fieldOrdinal, const int &pointOrdinal) {
                  basisValueView(fieldOrdinal,pointOrdinal) *= weight;
                });
              }
              else if (basisValueView.rank() == 3)
              {
                auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{basisValueView.extent_int(0),basisValueView.extent_int(1),basisValueView.extent_int(2)});
                Kokkos::parallel_for("multiply basisValueView by weight", policy,
                KOKKOS_LAMBDA (const int &fieldOrdinal, const int &pointOrdinal, const int &d) {
                  basisValueView(fieldOrdinal,pointOrdinal,d) *= weight;
                });
              }
              else
              {
                INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported rank for basisValueView");
              }
            }
          }
        }
      }
    }
    
    /** \brief  Evaluation of a FEM basis on a <strong>reference cell</strong>.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values
        \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points
        \param  operatorType      [in]  - the operator acting on the basis functions

        \remark For rank and dimension specifications of the output array see Section
        \ref basis_md_array_sec.  Dimensions of <var>ArrayScalar</var> arguments are checked
        at runtime if HAVE_INTREPID2_DEBUG is defined.

        \remark A FEM basis spans a COMPLETE or INCOMPLETE polynomial space on the reference cell
        which is a smooth function space. Thus, all operator types that are meaningful for the
        approximated function space are admissible. When the order of the operator exceeds the
        degree of the basis, the output array is filled with the appropriate number of zeros.
    */
    void getValues( OutputViewType outputValues, const PointViewType  inputPoints,
                   const EOperator operatorType = OPERATOR_VALUE ) const override
    {
      bool tensorPoints;  // true would mean that we take the tensor product of inputPoints1 and inputPoints2 (and that this would be equivalent to inputPoints as given -- i.e., inputPoints1 and inputPoints2 would be a tensor decomposition of inputPoints)
      bool attemptTensorDecomposition = false; // support for this not yet implemented
      PointViewType inputPoints1, inputPoints2;
      getComponentPoints(inputPoints, attemptTensorDecomposition, inputPoints1, inputPoints2, tensorPoints);
      
      const auto functionSpace = this->getFunctionSpace();
      
      if ((functionSpace == FUNCTION_SPACE_HVOL) || (functionSpace == FUNCTION_SPACE_HGRAD))
      {
        // then we can handle VALUE, GRAD, and Op_Dn without reference to subclass
        switch (operatorType)
        {
          case OPERATOR_VALUE:
          case OPERATOR_GRAD:
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
          {
            auto opOrder = getOperatorOrder(operatorType); // number of derivatives that we take in total
            // the Dk enumeration happens in lexicographic order (reading from left to right: x, y, z, etc.)
            // this governs the nesting order of the dkEnum1, dkEnum2 for loops below: dkEnum2 should increment fastest.
            for (int derivativeCountComp2=0; derivativeCountComp2<=opOrder; derivativeCountComp2++)
            {
              int derivativeCountComp1=opOrder-derivativeCountComp2;
              EOperator op1 = (derivativeCountComp1 == 0) ? OPERATOR_VALUE : EOperator(OPERATOR_D1 + (derivativeCountComp1 - 1));
              EOperator op2 = (derivativeCountComp2 == 0) ? OPERATOR_VALUE : EOperator(OPERATOR_D1 + (derivativeCountComp2 - 1));
              
              int spaceDim1 = inputPoints1.extent_int(1);
              int spaceDim2 = inputPoints2.extent_int(1);
              
              int dkCardinality1 = (op1 != OPERATOR_VALUE) ? getDkCardinality(op1, spaceDim1) : 1;
              int dkCardinality2 = (op2 != OPERATOR_VALUE) ? getDkCardinality(op2, spaceDim2) : 1;
              
              int basisCardinality1 = basis1_->getCardinality();
              int basisCardinality2 = basis2_->getCardinality();
              
              int totalPointCount = tensorPoints ? inputPoints1.extent_int(0) * inputPoints2.extent_int(0) : inputPoints1.extent_int(0);
              
              int pointCount1, pointCount2;
              if (tensorPoints)
              {
                pointCount1 = inputPoints1.extent_int(0);
                pointCount2 = inputPoints2.extent_int(0);
              }
              else
              {
                pointCount1 = totalPointCount;
                pointCount2 = totalPointCount;
              }
              
              OutputViewType outputValues1, outputValues2;
              if (op1 == OPERATOR_VALUE)
                outputValues1 = getMatchingViewWithLabel(outputValues, "output values - basis 1",basisCardinality1,pointCount1);
              else
                outputValues1 = getMatchingViewWithLabel(outputValues, "output values - basis 1",basisCardinality1,pointCount1,dkCardinality1);
              
              if (op2 == OPERATOR_VALUE)
                outputValues2 = getMatchingViewWithLabel(outputValues, "output values - basis 2",basisCardinality2,pointCount2);
              else
                outputValues2 = getMatchingViewWithLabel(outputValues, "output values - basis 2",basisCardinality2,pointCount2,dkCardinality2);
                
              basis1_->getValues(outputValues1,inputPoints1,op1);
              basis2_->getValues(outputValues2,inputPoints2,op2);
              
              const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputValueType>();
              const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointValueType>();
              const int vectorSize = std::max(outputVectorSize,pointVectorSize);
              
              auto policy = Kokkos::TeamPolicy<ExecutionSpace>(basisCardinality1,Kokkos::AUTO(),vectorSize);
              
              double weight = 1.0;
              using FunctorType = TensorViewFunctor<ExecutionSpace, OutputValueType, OutputViewType>;
              
              for (int dkEnum1=0; dkEnum1<dkCardinality1; dkEnum1++)
              {
                auto outputValues1_dkEnum1 = (op1 != OPERATOR_VALUE) ? Kokkos::subview(outputValues1,Kokkos::ALL(),Kokkos::ALL(),dkEnum1)
                : Kokkos::subview(outputValues1,Kokkos::ALL(),Kokkos::ALL());
                for (int dkEnum2=0; dkEnum2<dkCardinality2; dkEnum2++)
                {
                  auto outputValues2_dkEnum2 = (op2 != OPERATOR_VALUE) ? Kokkos::subview(outputValues2,Kokkos::ALL(),Kokkos::ALL(),dkEnum2)
                  : Kokkos::subview(outputValues2,Kokkos::ALL(),Kokkos::ALL());
                  
                  ordinal_type dkTensorIndex = getTensorDkEnumeration(dkEnum1, derivativeCountComp1, dkEnum2, derivativeCountComp2);
                  auto outputValues_dkTensor = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),dkTensorIndex);
                  // Note that there may be performance optimizations available here:
                  // - could eliminate interior for loop in favor of having a vector-valued outputValues1_dk
                  // - could add support to TensorViewFunctor (and probably TensorViewIterator) for this kind of tensor Dk type of traversal
                  //   (this would allow us to eliminate both for loops here)
                  // At the moment, we defer such optimizations on the idea that this may not ever become a performance bottleneck.
                  FunctorType functor(outputValues_dkTensor, outputValues1_dkEnum1, outputValues2_dkEnum2, tensorPoints, weight);
                  Kokkos::parallel_for("TensorViewFunctor", policy , functor);
                }
              }
            }
          }
            break;
          default: // non-OPERATOR_Dn case must be handled by subclass.
            this->getValues(outputValues, operatorType, inputPoints1, inputPoints2, tensorPoints);
        }
      }
      else
      {
        // not HVOL or HGRAD; subclass must handle
        this->getValues(outputValues, operatorType, inputPoints1, inputPoints2, tensorPoints);
      }
    }
    
    /** \brief  Evaluation of a tensor FEM basis on a <strong>reference cell</strong>; subclasses should override this.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values
        \param  operatorType      [in]  - the operator acting on the basis functions
        \param  inputPoints1      [in]  - rank-2 array (P1,D1) with the evaluation points for basis1
        \param  inputPoints2      [in]  - rank-2 array (P2,D2) with the evaluation points for basis2
        \param  tensorPoints      [in]  - whether the points should be interpreted as tensor components of the evaluation points, or in a one-to-one correspondence

       Subclasses should override this method; this gives them an opportunity to specify how operatorType should be decomposed into operators on the component bases.
     
       If tensorPoints is true, then the points dimension of outputValues should be (P1*P2).
       If tensorPoints is false, then P1 should equal P2, and these should match the points dimension of outputValues.
     
     There are three variants of getValues:
     1. The three-argument version defined by Intrepid2::Basis.  TensorBasis provides an implementation of this, which calls the five-argument version (this one).
     2. The five-argument version (this method), which provides separate point sets for the component bases, and must be specified by subclasses.  Typical implementations call the seven-argument version.
     3. The seven-argument version (below), implemented by TensorBasis, which provides separate point sets and operators for the component bases, as well as an optional weight.
     
     The intent is that subclasses implement this five-argument version; in that implementation, they need to do little else than call the seven-argument version below.
     
     Note that the three-argument implementation handles the OPERATOR_Dn operators directly; that is, subclasses can omit any consideration of OPERATOR_Dn operators in their implementation of the five-argument version.
    */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType  inputPoints1, const PointViewType  inputPoints2,
                           bool tensorPoints) const
    {
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "one-operator, two-inputPoints getValues should be overridden by TensorBasis subclasses");
    }
    
    /** \brief  Evaluation of a tensor FEM basis on a <strong>reference cell</strong>.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values
        \param  inputPoints1      [in]  - rank-2 array (P1,D1) with the evaluation points for basis1
        \param  operatorType1     [in]  - the operator acting on basis1
        \param  inputPoints2      [in]  - rank-2 array (P2,D2) with the evaluation points for basis2
        \param  operatorType2     [in]  - the operator acting on basis2
        \param  tensorPoints      [in]  - whether the points should be interpreted as tensor components of the evaluation points, or in a one-to-one correspondence
        \param  weight            [in]  - optional weight (typically 1.0 or -1.0)
     
       If tensorPoints is true, then the points dimension of outputValues should be (P1*P2).
       If tensorPoints is false, then P1 should equal P2, and these should match the points dimension of outputValues.
     
     There are three variants of getValues:
     1. The three-argument version defined by Intrepid2::Basis.  TensorBasis provides an implementation of this, which calls the five-argument version (this one).
     2. The five-argument version (above), which provides separate point sets for the component bases, and must be specified by subclasses.  Typical implementations call the seven-argument version.
     3. The seven-argument version (this method), implemented by TensorBasis, which provides separate point sets and operators for the component bases, as well as an optional weight.
     
     Subclasses should override the five-argument version above; in their implementation, they need to do little else than call this seven-argument version.
    */
    void getValues( OutputViewType outputValues,
                   const PointViewType  inputPoints1, const EOperator operatorType1,
                   const PointViewType  inputPoints2, const EOperator operatorType2,
                   bool tensorPoints, double weight=1.0) const
    {
      int basisCardinality1 = basis1_->getCardinality();
      int basisCardinality2 = basis2_->getCardinality();
      
      int totalPointCount = tensorPoints ? inputPoints1.extent_int(0) * inputPoints2.extent_int(0) : inputPoints1.extent_int(0);
      
      int pointCount1, pointCount2;
      if (tensorPoints)
      {
        pointCount1 = inputPoints1.extent_int(0);
        pointCount2 = inputPoints2.extent_int(0);
      }
      else
      {
        pointCount1 = totalPointCount;
        pointCount2 = totalPointCount;
      }
      
      const ordinal_type spaceDim1 = inputPoints1.extent_int(1);
      const ordinal_type spaceDim2 = inputPoints2.extent_int(1);
      
      INTREPID2_TEST_FOR_EXCEPTION(!tensorPoints && (totalPointCount != inputPoints2.extent_int(0)),
                                   std::invalid_argument, "If tensorPoints is false, the point counts must match!");
            
      const ordinal_type opRank1 = getOperatorRank(basis1_->getFunctionSpace(), operatorType1, spaceDim1);
      const ordinal_type opRank2 = getOperatorRank(basis2_->getFunctionSpace(), operatorType2, spaceDim2);
      
      const ordinal_type outputRank1 = opRank1 + getFieldRank(basis1_->getFunctionSpace());
      const ordinal_type outputRank2 = opRank2 + getFieldRank(basis2_->getFunctionSpace());
      
      OutputViewType outputValues1, outputValues2;
      if (outputRank1 == 0)
      {
        outputValues1 = getMatchingViewWithLabel(outputValues,"output values - basis 1",basisCardinality1,pointCount1);
      }
      else if (outputRank1 == 1)
      {
        outputValues1 = getMatchingViewWithLabel(outputValues,"output values - basis 1",basisCardinality1,pointCount1,spaceDim1);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported opRank1");
      }
      
      if (outputRank2 == 0)
      {
        outputValues2 = getMatchingViewWithLabel(outputValues,"output values - basis 2",basisCardinality2,pointCount2);
      }
      else if (outputRank2 == 1)
      {
        outputValues2 = getMatchingViewWithLabel(outputValues,"output values - basis 2",basisCardinality2,pointCount2,spaceDim2);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported opRank2");
      }
      
      basis1_->getValues(outputValues1,inputPoints1,operatorType1);
      basis2_->getValues(outputValues2,inputPoints2,operatorType2);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputValueType>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointValueType>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);
      
      auto policy = Kokkos::TeamPolicy<ExecutionSpace>(basisCardinality1,Kokkos::AUTO(),vectorSize);
      
      using FunctorType = TensorViewFunctor<ExecutionSpace, OutputValueType, OutputViewType>;
      
      FunctorType functor(outputValues, outputValues1, outputValues2, tensorPoints, weight);
      Kokkos::parallel_for("TensorViewFunctor", policy , functor);
    }
    
    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual HostBasisPtr<OutputValueType, PointValueType>
    getHostBasis() const override {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "TensorBasis subclasses must override getHostBasis");
    }
  }; // Basis_TensorBasis
  
  /** \struct Intrepid2::TensorBasis3_Functor
      \brief  Functor for computing values for the TensorBasis3 class.
   
   This functor is not intended for use outside of \ref Intrepid2::Basis_TensorBasis3.
   
   We may replace usage of this functor with TensorViewFunctor in the future.  This would likely allow more TensorBasis3 use cases.
  */
  template<class ExecutionSpace, class OutputScalar, class OutputFieldType>
  struct TensorBasis3_Functor
  {
    using ScratchSpace       = typename ExecutionSpace::scratch_memory_space;
    using OutputScratchView  = Kokkos::View<OutputScalar*,ScratchSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    
    using TeamPolicy = Kokkos::TeamPolicy<ExecutionSpace>;
    using TeamMember = typename TeamPolicy::member_type;
    
    OutputFieldType  output_; // F,P
    OutputFieldType  input1_; // F1,P[,D] or F1,P1[,D]
    OutputFieldType  input2_; // F2,P[,D] or F2,P2[,D]
    OutputFieldType  input3_; // F2,P[,D] or F2,P2[,D]
    
    int numFields_, numPoints_;
    int numFields1_, numPoints1_;
    int numFields2_, numPoints2_;
    int numFields3_, numPoints3_;
    
    bool tensorPoints_; // if true, input1, input2, input3 refer to values at decomposed points, and P = P1 * P2 * P3.  If false, then the three inputs refer to points in the full-dimensional space, and their point lengths are the same as that of the final output.
    
    double weight_;
    
    TensorBasis3_Functor(OutputFieldType output, OutputFieldType inputValues1, OutputFieldType inputValues2, OutputFieldType inputValues3,
                         bool tensorPoints, double weight)
    : output_(output), input1_(inputValues1), input2_(inputValues2), input3_(inputValues3), tensorPoints_(tensorPoints), weight_(weight)
    {
      numFields_ = output.extent_int(0);
      numPoints_ = output.extent_int(1);
      
      numFields1_ = inputValues1.extent_int(0);
      numPoints1_ = inputValues1.extent_int(1);
      
      numFields2_ = inputValues2.extent_int(0);
      numPoints2_ = inputValues2.extent_int(1);
      
      numFields3_ = inputValues3.extent_int(0);
      numPoints3_ = inputValues3.extent_int(1);
      /*
       We don't yet support tensor-valued bases here (only vector and scalar).  The main design question is how the layouts
       of the input containers relates to the layout of the output container.  The work we've done in TensorViewIterator basically
       shows the choices that can be made.  It does appear that in most cases (at least (most of?) those supported by TensorViewIterator),
       we can infer from the dimensions of input/output containers what choice should be made in each dimension.
       */
      INTREPID2_TEST_FOR_EXCEPTION(inputValues1.rank() > 3, std::invalid_argument, "ranks greater than 3 not yet supported");
      INTREPID2_TEST_FOR_EXCEPTION(inputValues2.rank() > 3, std::invalid_argument, "ranks greater than 3 not yet supported");
      INTREPID2_TEST_FOR_EXCEPTION(inputValues3.rank() > 3, std::invalid_argument, "ranks greater than 3 not yet supported");
      INTREPID2_TEST_FOR_EXCEPTION((inputValues1.rank() == 3) && (inputValues2.rank() == 3), std::invalid_argument, "two vector-valued input ranks not yet supported");
      INTREPID2_TEST_FOR_EXCEPTION((inputValues1.rank() == 3) && (inputValues3.rank() == 3), std::invalid_argument, "two vector-valued input ranks not yet supported");
      INTREPID2_TEST_FOR_EXCEPTION((inputValues2.rank() == 3) && (inputValues3.rank() == 3), std::invalid_argument, "two vector-valued input ranks not yet supported");
      
      if (!tensorPoints_)
      {
        // then the point counts should all match
        INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != numPoints1_, std::invalid_argument, "incompatible point counts");
        INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != numPoints2_, std::invalid_argument, "incompatible point counts");
        INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != numPoints3_, std::invalid_argument, "incompatible point counts");
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != numPoints1_ * numPoints2_ * numPoints3_, std::invalid_argument, "incompatible point counts");
      }
      
      INTREPID2_TEST_FOR_EXCEPTION(numFields_ != numFields1_ * numFields2_ * numFields3_, std::invalid_argument, "incompatible field sizes");
    }
    
    KOKKOS_INLINE_FUNCTION
    void operator()( const TeamMember & teamMember ) const
    {
      auto fieldOrdinal1 = teamMember.league_rank();
      
      if (!tensorPoints_)
      {
        if ((input1_.rank() == 2) && (input2_.rank() == 2) && (input3_.rank() == 2))
        {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal=0; pointOrdinal<numPoints_; pointOrdinal++)
              {
                output_(fieldOrdinal,pointOrdinal) = weight_ * input1_(fieldOrdinal1,pointOrdinal) * input2_(fieldOrdinal2,pointOrdinal) * input3_(fieldOrdinal3,pointOrdinal);
              }
            }
          });
        }
        else if (input1_.rank() == 3)
        {
          int spaceDim = input1_.extent_int(2);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal=0; pointOrdinal<numPoints_; pointOrdinal++)
              {
                for (int d=0; d<spaceDim; d++)
                {
                  output_(fieldOrdinal,pointOrdinal,d) = weight_ * input1_(fieldOrdinal1,pointOrdinal,d) * input2_(fieldOrdinal2,pointOrdinal) * input3_(fieldOrdinal3,pointOrdinal);
                }
              }
            }
          });
        }
        else if (input2_.rank() == 3)
        {
          int spaceDim = input2_.extent_int(2);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal=0; pointOrdinal<numPoints_; pointOrdinal++)
              {
                for (int d=0; d<spaceDim; d++)
                {
                  output_(fieldOrdinal,pointOrdinal,d) = weight_ * input1_(fieldOrdinal1,pointOrdinal) * input2_(fieldOrdinal2,pointOrdinal,d) * input3_(fieldOrdinal3,pointOrdinal);
                }
              }
            }
          });
        }
        else if (input3_.rank() == 3)
        {
          int spaceDim = input3_.extent_int(2);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal=0; pointOrdinal<numPoints_; pointOrdinal++)
              {
                for (int d=0; d<spaceDim; d++)
                {
                  output_(fieldOrdinal,pointOrdinal,d) = weight_ * input1_(fieldOrdinal1,pointOrdinal) * input2_(fieldOrdinal2,pointOrdinal) * input3_(fieldOrdinal3,pointOrdinal,d);
                }
              }
            }
          });
        }
        else
        {
          // unsupported rank combination -- enforced in constructor
        }
      }
      else
      {
        if ((input1_.rank() == 2) && (input2_.rank() == 2) && (input3_.rank() == 2) )
        {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal3=0; pointOrdinal3<numPoints3_; pointOrdinal3++)
              {
                for (int pointOrdinal2=0; pointOrdinal2<numPoints2_; pointOrdinal2++)
                {
                  for (int pointOrdinal1=0; pointOrdinal1<numPoints1_; pointOrdinal1++)
                  {
                    int pointOrdinal = (pointOrdinal3 * numPoints2_ + pointOrdinal2) * numPoints1_ + pointOrdinal1;
                    output_(fieldOrdinal,pointOrdinal) = weight_ * input1_(fieldOrdinal1,pointOrdinal1) * input2_(fieldOrdinal2,pointOrdinal2) * input3_(fieldOrdinal3,pointOrdinal3);
                  }
                }
              }
            }
          });
        }
        else if (input1_.rank() == 3) // based on constructor requirements, this means the others are rank 2
        {
          int spaceDim = input1_.extent_int(2);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal3=0; pointOrdinal3<numPoints3_; pointOrdinal3++)
              {
                for (int pointOrdinal2=0; pointOrdinal2<numPoints2_; pointOrdinal2++)
                {
                  for (int pointOrdinal1=0; pointOrdinal1<numPoints1_; pointOrdinal1++)
                  {
                    int pointOrdinal = (pointOrdinal3 * numPoints2_ + pointOrdinal2) * numPoints1_ + pointOrdinal1;
                    for (int d=0; d<spaceDim; d++)
                    {
                      output_(fieldOrdinal,pointOrdinal,d) = weight_ * input1_(fieldOrdinal1,pointOrdinal1,d) * input2_(fieldOrdinal2,pointOrdinal2) * input3_(fieldOrdinal3,pointOrdinal3);
                    }
                  }
                }
              }
            }
          });
        }
        else if (input2_.rank() == 3) // based on constructor requirements, this means the others are rank 2
        {
          int spaceDim = input2_.extent_int(2);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal3=0; pointOrdinal3<numPoints3_; pointOrdinal3++)
              {
                for (int pointOrdinal2=0; pointOrdinal2<numPoints2_; pointOrdinal2++)
                {
                  for (int pointOrdinal1=0; pointOrdinal1<numPoints1_; pointOrdinal1++)
                  {
                    int pointOrdinal = (pointOrdinal3 * numPoints2_ + pointOrdinal2) * numPoints1_ + pointOrdinal1;
                    for (int d=0; d<spaceDim; d++)
                    {
                      output_(fieldOrdinal,pointOrdinal,d) = weight_ * input1_(fieldOrdinal1,pointOrdinal1) * input2_(fieldOrdinal2,pointOrdinal2,d) * input3_(fieldOrdinal3,pointOrdinal3);
                    }
                  }
                }
              }
            }
          });
        }
        else if (input3_.rank() == 3) // based on constructor requirements, this means the others are rank 2
        {
          int spaceDim = input3_.extent_int(2);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal3=0; pointOrdinal3<numPoints3_; pointOrdinal3++)
              {
                for (int pointOrdinal2=0; pointOrdinal2<numPoints2_; pointOrdinal2++)
                {
                  for (int pointOrdinal1=0; pointOrdinal1<numPoints1_; pointOrdinal1++)
                  {
                    int pointOrdinal = (pointOrdinal3 * numPoints2_ + pointOrdinal2) * numPoints1_ + pointOrdinal1;
                    for (int d=0; d<spaceDim; d++)
                    {
                      output_(fieldOrdinal,pointOrdinal,d) = weight_ * input1_(fieldOrdinal1,pointOrdinal1) * input2_(fieldOrdinal2,pointOrdinal2) * input3_(fieldOrdinal3,pointOrdinal3,d);
                    }
                  }
                }
              }
            }
          });
        }
        else
        {
          // unsupported rank combination -- enforced in constructor
        }
      }
    }
  }; // TensorBasis3_Functor
  
  
  template<typename BasisBaseClass = void>
  class Basis_TensorBasis3
  : public Basis_TensorBasis<BasisBaseClass>
  {
    using BasisBase   = BasisBaseClass;
    using TensorBasis = Basis_TensorBasis<BasisBase>;
  public:
    using typename BasisBase::OutputViewType;
    using typename BasisBase::PointViewType;
    using typename BasisBase::ScalarViewType;

    using typename BasisBase::OutputValueType;
    using typename BasisBase::PointValueType;

    using typename BasisBase::ExecutionSpace;

    using BasisPtr  = Teuchos::RCP<BasisBase>;
  protected:
    BasisPtr basis1_;
    BasisPtr basis2_;
    BasisPtr basis3_;
  public:
    Basis_TensorBasis3(BasisPtr basis1, BasisPtr basis2, BasisPtr basis3, const bool useShardsCellTopologyAndTags = false)
    :
    TensorBasis(Teuchos::rcp( new TensorBasis(basis1,basis2,FUNCTION_SPACE_MAX,useShardsCellTopologyAndTags)),
                basis3,
                FUNCTION_SPACE_MAX,useShardsCellTopologyAndTags),
    basis1_(basis1),
    basis2_(basis2),
    basis3_(basis3)
    {}
    
    /** \brief Returns a full decomposition of the specified operator.  (Full meaning that all TensorBasis components are expanded into their non-TensorBasis components.)
     
     The outer container has a length that corresponds to the number of vector components in output; the inner container has one entry per basis component.
     
     A one-element outer container corresponds to a single TensorData entry; a multiple-element outer container corresponds to a VectorData object with axialComponents = false.
    */
    virtual OperatorTensorDecomposition getOperatorDecomposition(const EOperator operatorType) const override
    {
      OperatorTensorDecomposition opSimpleDecomposition = this->getSimpleOperatorDecomposition(operatorType);
      std::vector<BasisPtr> componentBases {basis1_, basis2_, basis3_};
      return opSimpleDecomposition.expandedDecomposition(componentBases);
    }
    
    using TensorBasis::getValues;
    
    /** \brief  Evaluation of a tensor FEM basis on a <strong>reference cell</strong>.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values
        \param  operatorType      [in]  - the operator acting on the basis functions
        \param  inputPoints12     [in]  - rank-2 array (P12,D12) with the evaluation points for basis12
        \param  inputPoints3      [in]  - rank-2 array (P3,D3) with the evaluation points for basis3
        \param  tensorPoints      [in]  - whether the points should be interpreted as tensor components of the evaluation points, or in a one-to-one correspondence
     
       If tensorPoints is true, then the points dimension of outputValues should be (P12*P3).
       If tensorPoints is false, then P12 should equal P3, and these should match the points dimension of outputValues.
     
     There are four variants of getValues:
     1. The three-argument version defined by Intrepid2::Basis.  TensorBasis provides an implementation of this, which calls the five-argument version (this one).
     2. The five-argument version (this method), which provides partially separated point sets for the component bases.  TensorBasis3 provides an implementation of this, which calls the six-argument version.
     3. The six-argument version, which fully separates the point sets for the component bases.  Subclasses should implement this; they essentially specify the decomposition of the operator.
     4. The nine-argument version (below), implemented by TensorBasis3, which provides separate point sets and operators for the component bases, as well as an optional weight.
     
     The intent is that subclasses implement the six-argument version; in that implementation, they need to do little else than call the nine-argument version below.
     
     Note that the three-argument implementation handles the OPERATOR_Dn operators directly; that is, subclasses can omit any consideration of OPERATOR_Dn operators in their implementation of the six-argument version.
    */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType inputPoints12, const PointViewType  inputPoints3,
                           bool tensorPoints) const override
    {
      // TODO: rework this to use superclass's getComponentPoints.
      
      int spaceDim1 = basis1_->getDomainDimension();
      int spaceDim2 = basis2_->getDomainDimension();
      
      int totalSpaceDim12 = inputPoints12.extent_int(1);
      
      TEUCHOS_ASSERT(spaceDim1 + spaceDim2 == totalSpaceDim12);
      
      if (!tensorPoints)
      {
        auto inputPoints1 = Kokkos::subview(inputPoints12,Kokkos::ALL(),std::make_pair(0,spaceDim1));
        auto inputPoints2 = Kokkos::subview(inputPoints12,Kokkos::ALL(),std::make_pair(spaceDim1,totalSpaceDim12));
        
        this->getValues(outputValues, operatorType, inputPoints1, inputPoints2, inputPoints3, tensorPoints);
      }
      else
      {
        // superclass doesn't (yet) have a clever way to detect tensor points in a single container
        // we'd need something along those lines here to detect them in inputPoints12.
        // if we do add such a mechanism to superclass, it should be simple enough to call that from here
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "This method does not yet handle tensorPoints=true");
      }
    }
    
    /** \brief  Evaluation of a tensor FEM basis on a <strong>reference cell</strong>; subclasses should override this.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values
        \param  operatorType      [in]  - the operator acting on the basis functions
        \param  inputPoints1      [in]  - rank-2 array (P1,D1) with the evaluation points for basis1
        \param  inputPoints1      [in]  - rank-2 array (P2,D2) with the evaluation points for basis2
        \param  inputPoints3      [in]  - rank-2 array (P3,D3) with the evaluation points for basis3
        \param  tensorPoints      [in]  - whether the points should be interpreted as tensor components of the evaluation points, or in a one-to-one correspondence

       Subclasses should override this method; this gives them an opportunity to specify how operatorType should be decomposed into operators on the component bases.
     
       If tensorPoints is true, then the points dimension of outputValues should be (P1*P2*P3).
       If tensorPoints is false, then P1 should equal P2 and P2 should equal P3, and these should match the points dimension of outputValues.
     
     There are four variants of getValues:
     1. The three-argument version defined by Intrepid2::Basis.  TensorBasis provides an implementation of this, which calls the five-argument version (this one).
     2. The five-argument version (above), which provides partially separated point sets for the component bases.  TensorBasis3 provides an implementation of this, which calls the six-argument version.
     3. The six-argument version (this method), which fully separates the point sets for the component bases.  Subclasses should implement this; they essentially specify the decomposition of the operator.
     4. The nine-argument version (below), implemented by TensorBasis3, which provides separate point sets and operators for the component bases, as well as an optional weight.
     
     The intent is that subclasses implement this six-argument version; in that implementation, they need to do little else than call the nine-argument version below.
     
     Note that the three-argument implementation handles the OPERATOR_Dn operators directly; that is, subclasses can omit any consideration of OPERATOR_Dn operators in their implementation of the five-argument version.
    */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType  inputPoints1, const PointViewType  inputPoints2, const PointViewType inputPoints3,
                           bool tensorPoints) const = 0;
    
    /** \brief  Evaluation of a tensor FEM basis on a <strong>reference cell</strong>; subclasses should override this.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values
        \param  inputPoints1      [in]  - rank-2 array (P1,D1) with the evaluation points for basis1
        \param  operatorType1     [in]  - the operator acting on basis1
        \param  inputPoints2      [in]  - rank-2 array (P2,D2) with the evaluation points for basis2
        \param  operatorType2     [in]  - the operator acting on basis2
        \param  inputPoints3      [in]  - rank-2 array (P3,D3) with the evaluation points for basis3
        \param  operatorType3     [in]  - the operator acting on basis3
        \param  tensorPoints      [in]  - whether the points should be interpreted as tensor components of the evaluation points, or in a one-to-one correspondence
     
       If tensorPoints is true, then the points dimension of outputValues should be (P1*P2*P3).
       If tensorPoints is false, then P1 should equal P2 and P2 should equal P3, and these should match the points dimension of outputValues.
     
     There are four variants of getValues:
     1. The three-argument version defined by Intrepid2::Basis.  TensorBasis provides an implementation of this, which calls the five-argument version (this one).
     2. The five-argument version (above), which provides partially separated point sets for the component bases.  TensorBasis3 provides an implementation of this, which calls the six-argument version.
     3. The six-argument version (this method), which fully separates the point sets for the component bases.  Subclasses should implement this; they essentially specify the decomposition of the operator.
     4. The nine-argument version (below), implemented by TensorBasis3, which provides separate point sets and operators for the component bases, as well as an optional weight.
     
     The intent is that subclasses implement this six-argument version; in that implementation, they need to do little else than call the nine-argument version below.
     
     Note that the three-argument implementation handles the OPERATOR_Dn operators directly; that is, subclasses can omit any consideration of OPERATOR_Dn operators in their implementation of the five-argument version.
    */
    void getValues( OutputViewType outputValues,
                   const PointViewType  inputPoints1, const EOperator operatorType1,
                   const PointViewType  inputPoints2, const EOperator operatorType2,
                   const PointViewType  inputPoints3, const EOperator operatorType3,
                   bool tensorPoints, double weight=1.0) const
    {
      int basisCardinality1 = basis1_->getCardinality();
      int basisCardinality2 = basis2_->getCardinality();
      int basisCardinality3 = basis3_->getCardinality();
      
      int spaceDim1 = inputPoints1.extent_int(1);
      int spaceDim2 = inputPoints2.extent_int(1);
      int spaceDim3 = inputPoints3.extent_int(1);
      
      int totalPointCount;
      int pointCount1, pointCount2, pointCount3;
      if (tensorPoints)
      {
        pointCount1 = inputPoints1.extent_int(0);
        pointCount2 = inputPoints2.extent_int(0);
        pointCount3 = inputPoints3.extent_int(0);
        totalPointCount = pointCount1 * pointCount2 * pointCount3;
      }
      else
      {
        totalPointCount = inputPoints1.extent_int(0);
        pointCount1 = totalPointCount;
        pointCount2 = totalPointCount;
        pointCount3 = totalPointCount;
        
        INTREPID2_TEST_FOR_EXCEPTION((totalPointCount != inputPoints2.extent_int(0)) || (totalPointCount != inputPoints3.extent_int(0)),
                                     std::invalid_argument, "If tensorPoints is false, the point counts must match!");
      }
      
      // structure of this implementation:
      /*
       - allocate output1, output2, output3 containers
       - either:
       1. split off the tensor functor call into its own method in TensorBasis, and
       - call it once with output1, output2, placing these in another newly allocated output12, then
       - call it again with output12, output3
       OR
       2. create a 3-argument tensor functor and call it with output1,output2,output3
       
       At the moment, the 3-argument functor seems like a better approach.  It's likely more code, but somewhat
       more efficient and easier to understand/debug.  And the code is fairly straightforward to produce.
       */
      
      // copied from the 2-argument TensorBasis implementation:
      
      OutputViewType outputValues1, outputValues2, outputValues3;

      //Note: the gradient of HGRAD basis on a line has an output vector of rank 3, the last dimension being of size 1.
      //      in particular this holds even when computing the divergence of an HDIV basis, which is scalar and has rank 2.
      if ((spaceDim1 == 1) && (operatorType1 == OPERATOR_VALUE))
      {
        // use a rank 2 container for basis1
        outputValues1 = getMatchingViewWithLabel(outputValues,"output values - basis 1",basisCardinality1,pointCount1);
      }
      else
      {
        outputValues1 = getMatchingViewWithLabel(outputValues,"output values - basis 1",basisCardinality1,pointCount1,spaceDim1);
      }
      if ((spaceDim2 == 1) && (operatorType2 == OPERATOR_VALUE))
      {
        // use a rank 2 container for basis2
        outputValues2 = getMatchingViewWithLabel(outputValues,"output values - basis 2",basisCardinality2,pointCount2);
      }
      else
      {
        outputValues2 = getMatchingViewWithLabel(outputValues,"output values - basis 2",basisCardinality2,pointCount2,spaceDim2);
      }
      if ((spaceDim3 == 1) && (operatorType3 == OPERATOR_VALUE))
      {
        // use a rank 2 container for basis2
        outputValues3 = getMatchingViewWithLabel(outputValues,"output values - basis 3",basisCardinality3,pointCount3);
      }
      else
      {
        outputValues3 = getMatchingViewWithLabel(outputValues,"output values - basis 3",basisCardinality3,pointCount3,spaceDim3);
      }

      basis1_->getValues(outputValues1,inputPoints1,operatorType1);
      basis2_->getValues(outputValues2,inputPoints2,operatorType2);
      basis3_->getValues(outputValues3,inputPoints3,operatorType3);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputValueType>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointValueType>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);

      auto policy = Kokkos::TeamPolicy<ExecutionSpace>(basisCardinality1,Kokkos::AUTO(),vectorSize);
      
      using FunctorType = TensorBasis3_Functor<ExecutionSpace, OutputValueType, OutputViewType>;
      FunctorType functor(outputValues, outputValues1, outputValues2, outputValues3, tensorPoints, weight);
      Kokkos::parallel_for("TensorBasis3_Functor", policy , functor);
    }

    /** \brief  Fills in coefficients of degrees of freedom on the reference cell
      \param [out] dofCoeffs - the container into which to place the degrees of freedom.

     dofCoeffs is a rank 1 with dimension equal to the cardinality of the basis.

     Note that getDofCoeffs() is not supported by all bases; in particular, hierarchical bases do not generally support this.
     */
    virtual void getDofCoeffs( typename BasisBase::ScalarViewType dofCoeffs ) const override
    {
      using ValueType    = typename BasisBase::ScalarViewType::value_type;
      using ResultLayout = typename DeduceLayout< typename BasisBase::ScalarViewType >::result_layout;
      using ViewType     = Kokkos::DynRankView<ValueType, ResultLayout, typename TensorBasis::DeviceType >;

      const ordinal_type basisCardinality1 = basis1_->getCardinality();
      const ordinal_type basisCardinality2 = basis2_->getCardinality();
      const ordinal_type basisCardinality3 = basis3_->getCardinality();

      auto dofCoeffs1 = ViewType("dofCoeffs1",basisCardinality1);
      auto dofCoeffs2 = ViewType("dofCoeffs2",basisCardinality2);
      auto dofCoeffs3 = ViewType("dofCoeffs3",basisCardinality3);

      basis1_->getDofCoeffs(dofCoeffs1);
      basis2_->getDofCoeffs(dofCoeffs2);
      basis3_->getDofCoeffs(dofCoeffs3);

      Kokkos::RangePolicy<ExecutionSpace> policy(0, basisCardinality3);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const int fieldOrdinal3)
       {
         for (int fieldOrdinal2=0; fieldOrdinal2<basisCardinality2; fieldOrdinal2++)
          for (int fieldOrdinal1=0; fieldOrdinal1<basisCardinality1; fieldOrdinal1++)
          {
            const ordinal_type fieldOrdinal = fieldOrdinal1 + fieldOrdinal2 * basisCardinality1 + fieldOrdinal3 * (basisCardinality1*basisCardinality2);
            dofCoeffs(fieldOrdinal)  = dofCoeffs1(fieldOrdinal1);
            dofCoeffs(fieldOrdinal) *= dofCoeffs2(fieldOrdinal2) * dofCoeffs3(fieldOrdinal3);
          }
       });
    }
    
    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual HostBasisPtr<OutputValueType, PointValueType>
    getHostBasis() const override {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "TensorBasis3 subclasses must override getHostBasis");
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_TensorBasis_h */
