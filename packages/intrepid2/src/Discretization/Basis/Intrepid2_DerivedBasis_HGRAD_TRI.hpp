// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_DerivedBasis_HGRAD_TRI.hpp
    \brief  Implementation of H(grad) basis on the triangle that is templated on H(grad) on the line.
    \author Created by N.V. Roberts.
 
 H(grad) on the triangle, defined following the construction of Fuentes et al.
 
 This class requires the template argument to be a hierarchical basis on the line; the construction here
 will not be nodal even if the line basis provided is nodal.  (It would be nodal on the edges, but even
 the edge functions would not be nodal on the interior of the triangle.)  Therefore, an exception is thrown
 at construction if a non-hierarchical basis is provided.
 
 In contrast to the derived bases defined on the quadrilateral and hexahedron, the interior (face) functions
 defined here are defined independently of the line basis provided as a template argument.  Following Fuentes
 et al., the interior functions are defined as Jacobi polynomials.
 
 Eventually, we hope to develop a derived basis where the interior functions are defined in terms of the
 functions defined on the line.  Ideally, we could define a generic "simplicial extension" that could be
 used to define bases on the triangle, tetrahedron, and pyramid, as well as higher-dimensional simplices.
 Ideally, this would also preserve hierarchical/nodal status.  As far as we know, however, such an extension
 has not yet been constructed in the literature; for now, it remains an aspiration for us.
 
 */

#ifndef Intrepid2_DerivedBasis_HGRAD_TRI_h
#define Intrepid2_DerivedBasis_HGRAD_TRI_h

namespace Intrepid2
{
  // TODO: place SimplexTopologyMap in its own file
  class SimplexTopologyMap
  {
    shards::CellTopology baseTopo_;
  public:
    SimplexTopologyMap(const shards::CellTopology &baseTopo)
    :
    baseTopo_(baseTopo)
    {}
    
    /** \brief  Map from component subcell ordinals to the corresponding composite subcell ordinal.
     
     \param  subcell1Dim     [in] - spatial dimension of the subcell in cellTopo1
     \param  subcell1Ordinal [in] - ordinal of the subcell in cellTopo1
     \param  subcell2Dim     [in] - spatial dimension of the line subcell (0 or 1)
     \param  subcell2Ordinal [in] - ordinal of the subcell in cellTopo2 (0 or 1)
     
     \return the subcell ordinal of the corresponding subcell in the composite cell topology.
     
     The dimension of the composite subcell is subcell1Dim + subcell2Dim.
     */
    unsigned getCompositeSubcellOrdinal(unsigned subcell1Dim, unsigned subcell1Ordinal, unsigned subcell2Dim, unsigned subcell2Ordinal)
    {
      INTREPID2_TEST_FOR_EXCEPTION(subcell1Dim > baseTopo_.getDimension(), std::invalid_argument, "Invalid base subcell dimension");
      INTREPID2_TEST_FOR_EXCEPTION(subcell2Dim > 1, std::invalid_argument, "Invalid line subcell dimension");
      const unsigned lineSubcellCount = (subcell2Dim == 0) ? 2 : 1;
      INTREPID2_TEST_FOR_EXCEPTION(subcell1Ordinal >= baseTopo_.getSubcellCount(subcell1Dim), std::invalid_argument, "Invalid base subcell ordinal");
      INTREPID2_TEST_FOR_EXCEPTION(subcell2Ordinal >= lineSubcellCount, std::invalid_argument, "Invalid line subcell ordinal");
      
      // the subcell numbering in dimension d is as follows:
      // - base subcells (0 to subcellCount(d)-1)
      // - extended subcells (subcellCount(d) to subcellCount(d) + subcellCount(d-1))
      
      if ((subcell2Dim == 0) && (subcell2Ordinal == 0))
      {
        // base subcell
        return subcell1Ordinal;
      }
      else if (subcell2Ordinal == 1)
      {
        // we do this test to ensure that even the apex satisfies the rule that the dimension of the composite subcell is subcell1Dim + subcell2Dim
        // (geometrically, one can argue that *all* the subcells in the base are mapped to the apex when selecting vertex 1 in the line)
        INTREPID2_TEST_FOR_EXCEPTION(subcell1Dim != 0, std::invalid_argument, "For line vertex 1, subcell in base must be a vertex");
        // this is the apex; numbered after the vertices in the base topo
        return baseTopo_.getVertexCount();
      }
      else
      {
        // subcell1 is being extended in the new dimension
        // we number these extended subcells following the same-dimensional subcells that exist in the base
        const int d = subcell1Dim + subcell2Dim;
        return baseTopo_.getSubcellCount(d) + subcell1Ordinal;
      }
    }
  };
  
  // TODO: place SimplicialExtensionBasis in its own file.
  template<class BaseTopoBasis, class LineBasis>
  class SimplicialExtensionBasis
  :
  public Basis<typename LineBasis::ExecutionSpace,typename LineBasis::OutputValueType,typename LineBasis::PointValueType>
  {
  public:
    using ExecutionSpace  = typename LineBasis::ExecutionSpace;
    using OutputValueType = typename LineBasis::OutputValueType;
    using PointValueType  = typename LineBasis::PointValueType;
    
    using OutputViewType = typename LineBasis::outputViewType;
    using PointViewType  = typename LineBasis::pointViewType ;
    using ScalarViewType = typename LineBasis::scalarViewType;
    
    using OrdinalTypeArray1DHost = typename LineBasis::OrdinalTypeArray1DHost;
    using OrdinalTypeArray2DHost = typename LineBasis::OrdinalTypeArray2DHost;
  protected:
    LineBasis     lineBasis_;
    BaseTopoBasis baseBasis_;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order
     */
    SimplicialExtensionBasis(int polyOrder)
    :
    lineBasis_(LineBasis(polyOrder)),
    baseBasis_(BaseTopoBasis(polyOrder))
    {
      this->functionSpace_ = FUNCTION_SPACE_HGRAD;
      
      auto lineCardinality     = lineBasis_.getCardinality();
      this->basisCardinality_  = lineCardinality * (lineCardinality + 1) / 2; // TODO: set this correctly in general (this is the formula for triangles)
      this->basisDegree_       = lineBasis_.getDegree();
      
      // this basis class only supports straight-edged cell topologies -- this is shards::Triangle<3>, same as shards::Triangle<>
      this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Triangle <> >());
      
      // assert that the input basis is defined on a (straight) line
      auto lineTopo    = lineBasis_.getBaseCellTopology();
      auto lineTopoKey = lineTopo.getKey();
      INTREPID2_TEST_FOR_EXCEPTION(lineTopoKey != shards::Line<2>::key, std::invalid_argument, "Input basis must be defined on a (straight) line");
      
      auto baseTopo    = baseBasis_.getBaseCellTopology();
      
      // because the basis we construct here will not be nodal even if the input line basis is nodal,
      // we require that the input basis be hierarchical.
      TEUCHOS_TEST_FOR_EXCEPTION(lineBasis_.getBasisType() != BASIS_FEM_HIERARCHICAL, std::invalid_argument, "non-hierarchical input bases are not (yet) supported");
      
      this->basisType_         = lineBasis_.getBasisType();
      this->basisCoordinates_  = COORDINATES_CARTESIAN;
      
      // initialize tags
      {
        const auto & cardinality = this->basisCardinality_;
        
        // Basis-dependent initializations
        const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
        const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
        const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
        const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
        
        OrdinalTypeArray1DHost tagView("tag view", cardinality*tagSize);
        
        shards::CellTopology cellTopo = this->basisCellTopology_;
        
        int spaceDim  = cellTopo.getDimension();

        if (this->getBasisType() == BASIS_FEM_HIERARCHICAL)
        {
          const int degreeSize = lineBasis_.getPolynomialDegreeLength();
          this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("Simplicial extension basis - field ordinal polynomial degree", this->basisCardinality_, degreeSize);
        }
        
        // for a generic simplicial extension, we would need to define a topology map, something like our tensor topology map, which
        // maps from the the subcell enumeration on the input topology combined with a subcell on the (extending) line topology
        // to the subcell enumeration on the simplicial extension topology.  The combination of the first point of the line topology
        // and any subcell on the input topology would be that subcell on the base; the combination of the second point of the line topology
        // with any subcell will be the topmost point.  The combination of the interior of the line topology and any subcell of the base
        // will be the simplicial extension of that subcell.
        
        SimplexTopologyMap topoMap(baseTopo);
        
        // we order the basis members according to the dimension of the subcell to which they belong:
        // vertex degrees of freedom come first, then edges, faces, 3D interiors
        int fieldOrdinal = 0;
        for (int d=0; d<=spaceDim; d++)
        {
          // d: dimension of the subcell
          // d2: dimension of the subcell in the extended dimension
          const int d2_max = std::min(1,d);
          for (int d2=0; d2<=d2_max; d2++)
          {
            const int d1 = d-d2; // d1: dimension of subcell in the base dimension
            const unsigned subcellCount1 = baseTopo.getSubcellCount(d1);
            for (unsigned subcellOrdinal1=0; subcellOrdinal1<subcellCount1; subcellOrdinal1++)
            {
              const unsigned subcellOrdinal2 = 0; // for d2=0, only use the vertex that touches the base; for d2=1, only one interior
              // (we handle the vertex 1 case for the extending line topo below -- that's the apex)
              ordinal_type subcellDofCount = baseBasis_.getDofCount(d1, subcellOrdinal1) * lineBasis_.getDofCount(d2, subcellOrdinal2);
              
              const int compositeSubcellOrdinal = topoMap.getCompositeSubcellOrdinal(d1, subcellOrdinal1, d2, subcellOrdinal2);
              for (ordinal_type localDofID = 0; localDofID < subcellDofCount; localDofID++)
              {
                tagView(fieldOrdinal*tagSize+0) = d; // subcell dimension
                tagView(fieldOrdinal*tagSize+1) = compositeSubcellOrdinal;
                tagView(fieldOrdinal*tagSize+2) = localDofID;
                tagView(fieldOrdinal*tagSize+3) = subcellDofCount;
              }
              fieldOrdinal += subcellDofCount;
            }
          }
          if (d==0)
          {
            // handle the apex here
            const int apexOrdinal = topoMap.getCompositeSubcellOrdinal(0, 0, 0, 1);
            ordinal_type apexDofCount = lineBasis_.getDofCount(0,1);
            for (ordinal_type localDofID = 0; localDofID < apexDofCount; localDofID++)
            {
              tagView(fieldOrdinal*tagSize+0) = d; // subcell dimension
              tagView(fieldOrdinal*tagSize+1) = apexOrdinal;
              tagView(fieldOrdinal*tagSize+2) = localDofID;
              tagView(fieldOrdinal*tagSize+3) = apexDofCount;
            }
            fieldOrdinal += apexDofCount;
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
    virtual void getValues( OutputViewType outputValues, const PointViewType  inputPoints,
                            const EOperator operatorType = OPERATOR_VALUE ) const override
    {
      // TODO: implement this (call the 4-argument version below)
    }
    
    /** \brief  multi-component getValues() method (required/called by SimplicialExtensionBasis)
        \param [out] outputValues - the view into which to place the output values
        \param [in] operatorType - the operator on the basis
        \param [in] inputPoints1 - input points in the base topology dimension(s)
        \param [in] inputPoints2 - input points in the extension dimension
     
     inputPoints1 and inputPoints2 should correspond entrywise to the evaluation points.
     This generally leads to redundant computations; in the future, we may offer a way to
     express a simplicial point set (analogous to the "tensorPoints" option in TensorBasis).
     */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType inputPoints1, const PointViewType inputPoints2) const = 0;
    
    /** \brief  Evaluation of a simplicial-extension FEM basis on a <strong>reference cell</strong>.
     
     Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
     points in the <strong>reference cell</strong> for which the basis is defined.
     
     \param  outputValues      [out] - variable rank array with the basis values
     \param  inputPoints1      [in]  - rank-2 array (P1,D1) with the evaluation points for basis1
     \param  operatorType1     [in]  - the operator acting on basis1
     \param  inputPoints2      [in]  - rank-2 array (P2,D2) with the evaluation points for basis2
     \param  operatorType2     [in]  - the operator acting on basis2
     \param  weight            [in]  - optional weight (typically 1.0 or -1.0)
     
     P1 should equal P2, and these should match the points dimension of outputValues.  In the future, we may offer a way to
     express a simplicial point set (analogous to the "tensorPoints" option in TensorBasis), at which point P1 and P2 could
     differ, and the point dimension of outputValues would be some function of P1 and P2.
     
     There are three variants of getValues:
     1. The three-argument version defined by Intrepid2::Basis.  SimplicialExtensionBasis provides an implementation of this, which calls the five-argument version (this one).
     2. The four-argument version (above), which provides separate point sets for the component bases, and must be specified by subclasses.  Typical implementations call the six-argument version.
     3. The six-argument version (this method), implemented by SimplicialExtensionBasis, which provides separate point sets and operators for the component bases, as well as an optional weight.
     
     Subclasses should override the four-argument version above; in their implementation, they need to do little else than determine the operator splitting and call this six-argument version.
     */
    void getValues( OutputViewType outputValues,
                   const PointViewType  inputPoints1, const EOperator operatorType1,
                   const PointViewType  inputPoints2, const EOperator operatorType2,
                   double weight=1.0) const
    {
      int baseBasisCardinality = baseBasis_.getCardinality();
      int lineBasisCardinality = lineBasis_.getCardinality();
      
      const int totalPointCount = inputPoints1.extent_int(0);
      
      // if/when we offer
      const int pointCount1 = totalPointCount;
      const int pointCount2 = totalPointCount;
      
      int spaceDim1 = inputPoints1.extent_int(1);
      int spaceDim2 = inputPoints2.extent_int(1);
      
      INTREPID2_TEST_FOR_EXCEPTION(totalPointCount != inputPoints2.extent_int(0),
                                   std::invalid_argument, "If tensorPoints is false, the point counts must match!");
      
      int opRankBase = getOperatorRank(baseBasis_.getFunctionSpace(), operatorType1, spaceDim1);
      int opRankLine = getOperatorRank(lineBasis_.getFunctionSpace(), operatorType2, spaceDim2);
      
      OutputViewType outputValues1, outputValues2;
      if (opRankBase == 0)
      {
        outputValues1 = getMatchingViewWithLabel(outputValues,"output values - base basis",baseBasisCardinality,pointCount1);
      }
      else if (opRankBase == 1)
      {
        outputValues1 = getMatchingViewWithLabel(outputValues,"output values - base basis",baseBasisCardinality,pointCount1,spaceDim1);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported opRankBase");
      }
      
      if (opRankLine == 0)
      {
        outputValues2 = getMatchingViewWithLabel(outputValues,"output values - line basis",lineBasisCardinality,pointCount2);
      }
      else if (opRankLine == 1)
      {
        outputValues2 = getMatchingViewWithLabel(outputValues,"output values - line basis",lineBasisCardinality,pointCount2,spaceDim2);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported opRank2");
      }
      
      baseBasis_.getValues(outputValues1,inputPoints1,operatorType1);
      lineBasis_.getValues(outputValues2,inputPoints2,operatorType2);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputValueType>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointValueType>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);
      
      auto policy = Kokkos::TeamPolicy<ExecutionSpace>(baseBasisCardinality,Kokkos::AUTO(),vectorSize);
      
      using FunctorType = TensorViewFunctor<ExecutionSpace, OutputValueType, OutputViewType>;
      
      const bool tensorPoints = false;
      FunctorType functor(outputValues, outputValues1, outputValues2, tensorPoints, weight);
      Kokkos::parallel_for( policy , functor, "TensorViewFunctor");
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HGRAD_TRI_h */
