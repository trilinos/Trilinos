// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_IntrepidManager_hpp
#define percept_IntrepidManager_hpp

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

#include <math.h>

#include <percept/Percept.hpp>

#include "Intrepid2_DefaultCubatureFactory.hpp"

#include <percept/function/StringFunction.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/MDArray.hpp>
#include <percept/PerceptMesh.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include "Teuchos_RCP.hpp"


// Shards includes
#include "Shards_CellTopology.hpp"
#include <percept/function/internal/ComputeBases.hpp>

#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_ArrayView.hpp"
// #include "Shards_Array.hpp"
#include "Teuchos_RCP.hpp"
// #include "Teuchos_BLAS.hpp"
#include "Teuchos_oblackholestream.hpp"
//#include "Teuchos_Assert.hpp"

#ifdef SHARDS_ARRAY_BOUNDS_CHECKING
xxx error 
#endif


using namespace Intrepid2;

#define IM_TAG( ADT ) ADT ## _TAG

#define IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( ADT )  \
  class ADT  : public shards::ArrayDimTag {         \
  public:                                           \
    const char * name() const ;                     \
    static const ADT & tag();                       \
    int num;                                        \
    ADT(int n);                                     \
    ~ADT();                                         \
    ADT( const ADT & adt);                          \
    ADT & operator = ( const ADT & adt);            \
  private:                                          \
    ADT();                                          \
  }

/** \brief  Macro for implementing the body of a simple ArrayDimTag
 *  \param ADT  name of the tag.
 */
#define IM_SHARDS_ARRAY_DIM_TAG_IMPLEMENTATION( ADT )                   \
  const char * ADT::name() const                                        \
  { static const char n[] = # ADT ; return n; }                         \
  const ADT & ADT::tag() { static const ADT self ; return self ; }      \
  ADT::ADT(int n) { num =n;}                                            \
  ADT::  ~ADT() {}                                                      \
  ADT::ADT( const ADT & adt) { num=adt.num;}                            \
  ADT & ADT::operator = ( const ADT & adt) {num=adt.num; return *this;} \
  ADT::ADT() {}                                                         


  namespace percept
  {

    using shards::CellTopology;

    IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Elements_Tag );
    IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Cub_Points_Tag );
    IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( NodesPerElem_Tag );
    IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Spatial_Dim_Tag );
    IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( DOFs_Tag );
    IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( BasisFields_Tag ); 

    void tni(void); 

    /**
     * |-------------------------------------------------------------------------------------------------|
     * |   Index type              | Dimension |  Description                                            |
     * |---------------------------|-----------|---------------------------------------------------------|
     * |   point                   |    [P]    |  number of points stored in an MD array                 |
     * |   vertex                  |    [V]    |  number of nodes stored in an MD aray                   |
     * |   field                   |    [F]    |  number of fields stored in an MD array                 |
     * |   basis field             |    [B]    |  number of basis fields stored in an MD array           |
     * |   cell                    |    [C]    |  number of cells stored in an MD array                  |
     * |   field coordinate        |    [D]    |  space dimension                                        |
     * |   derivative ordinal      |    [K]    |  cardinality of the set of kth derivatives              |
     * |                           |           |                                                         |
     * |   dof                     |   [DOF]   |  number of DOFs   stored in an MD array                 |
     * |-------------------------------------------------------------------------------------------------|
     *
     *  Note: Intrepid2 really doesn't have a concept of "DOF" at a node.  It's either a single variable,
     *    or a vector- or tensor-valued variable.  So, no DOF-related arrays as used herein can be used
     *    with Intrepid2 - you must call Intrepd2 one DOF at a time.
     * 
     *  MPerego: this is correct, however using subviews it should be possible to select a particular DOF 
     *           avoiding deep copying views.
     *
     * DynRankView<double> cub_points(numCubPoints, spaceDim);
     * DynRankView<double> cub_weights(numCubPoints);
     * 
     * DynRankView<double> cell_nodes(numCells, numNodes, spaceDim);
     * 
     * DynRankView<double> jacobian(numCells, numCubPoints, spaceDim, spaceDim);
     * DynRankView<double> jacobian_inv(numCells, numCubPoints, spaceDim, spaceDim);
     * DynRankView<double> jacobian_det(numCells, numCubPoints);
     * DynRankView<double> weighted_measure(numCells, numCubPoints);
     * 
     * DynRankView<double> grad_at_cub_points(numFields, numCubPoints, spaceDim);
     * DynRankView<double> transformed_grad_at_cub_points(numCells, numFields, numCubPoints, spaceDim);
     * DynRankView<double> weighted_transformed_grad_at_cub_points(numCells, numFields, numCubPoints, spaceDim);
     * DynRankView<double> stiffness_matrices(numCells, numFields, numFields);
     * 
     * 
     */

    class IntrepidManager
    {
    private:

#if (defined(__PGI) && defined(USE_PGI_7_1_COMPILER_BUG_WORKAROUND))
      // workaround for PGI compiler bug

      using BasisType = Intrepid2::Basis<Kokkos::HostSpace, double, double >;
      using BasisTypeRCP = Intrepid2::BasisPtr<Kokkos::HostSpace, double, double >;

      static void bootstrap();
#endif

    public:
      typedef IntrepidManager IM;
  
#define NUM(AClass) im.m_ ## AClass  . num

      /// ([P],[D])
      class CubaturePoints : public MDArray
      {
        const IntrepidManager& m_im;
      public:
        typedef MDArray BaseType;

        CubaturePoints(const IM& im) ;

        void copyTo(MDArray& mda)
        {
          Kokkos::resize(mda, m_im.m_Cub_Points_Tag.num, m_im.m_Spatial_Dim_Tag.num);
          Kokkos::deep_copy(mda, *this);
        }

        using BaseType::operator();
      };

      /// ([P])
      class CubatureWeights : public MDArray
      {

      public:
        typedef MDArray BaseType;

        CubatureWeights(const IM& im);
        using BaseType::operator();

      };

      /// ([C], [V], [D])
      class CellWorkSet : public MDArray
      {
      public:
        typedef MDArray BaseType;

        CellWorkSet(const IM& im);
        using BaseType::operator();
      };

      /// ([C], [P], [D])
      class PhysicalCoords : public MDArray
      {

        const IntrepidManager& m_im;
      public:
        typedef MDArray BaseType;

        PhysicalCoords(const IM& im);

        void operator()(CellWorkSet& c, CubaturePoints& xi);
        
        void copyTo(MDArray& mda)
        {
          Kokkos::resize(mda, m_im.m_Elements_Tag.num, m_im.m_Cub_Points_Tag.num, m_im.m_Spatial_Dim_Tag.num);
          Kokkos::deep_copy(mda, *this);
        }
        
	using BaseType::operator();

      };

      /// ([C], [P], [D], [D])
      class Jacobian : public MDArray
      {
        const IntrepidManager& m_im;
      public:
        typedef MDArray BaseType;

        Jacobian(const IM& im);
        void operator()(CubaturePoints& xi, CellWorkSet& c, CellTopology& topo);

        void copyTo(MDArray& mda)
        {
          Kokkos::resize(mda, m_im.m_Elements_Tag.num, m_im.m_Cub_Points_Tag.num, m_im.m_Spatial_Dim_Tag.num, m_im.m_Spatial_Dim_Tag.num);
          Kokkos::deep_copy(mda, *this);
        }

        using BaseType::operator();
      };

      /// ([C], [P], [D])
      class FaceNormal : public MDArray
      {
      public:
        typedef MDArray BaseType;

        FaceNormal(const IM& im);
        void operator()(Jacobian& jac, int i_face, CellTopology& topo);

        using BaseType::operator();
      };

      /// ([C], [P], [D], [D])
      class JacobianInverse : public MDArray
      {

      public:
        typedef MDArray BaseType;

        JacobianInverse(const IM& im);
        void operator()(Jacobian& jac);
        using BaseType::operator();
      };

      /// ([C], [P])
      class JacobianDet : public MDArray
      {

      public:
        typedef MDArray BaseType;

        JacobianDet(const IM& im);

        void operator()(Jacobian& jac);
        using BaseType::operator();
      };


      /// weights multiplied by Jacobian det at cubature points
      /// ([C], [P])
      class WeightedMeasure : public MDArray
      {

      public:
        typedef MDArray BaseType;

        WeightedMeasure(const IM& im);

        void operator()(CubatureWeights& w, JacobianDet& dJ);
        using BaseType::operator();

      };

      /// ([C], [P], [DOF])
      class IntegrandValuesDOF : public MDArray
      {

      public:
        typedef MDArray BaseType;

        IntegrandValuesDOF(const IM& im);

        void copyFrom(MDArray& mda);

        using BaseType::operator();
      };

      /// ([C], [P])
      class IntegrandValues : public MDArray
      {

      public:
        typedef MDArray BaseType;

        IntegrandValues(const IM& im);

        void copyFrom(MDArray& mda);
        void copyFrom(IM& im, MDArray& mda, int iDof);

        using BaseType::operator();
      };

      /// ([C], [DOF])
      class IntegralDOF : public MDArray
      {

      public:
        typedef MDArray BaseType;

        IntegralDOF(const IM& im);

        /// wXdOmega: ([C], [P])
        /// iv:       ([C], [P], [DOF])
        /// this:     ([C], [DOF])
        void operator()(IntegrandValuesDOF& iv, WeightedMeasure& wXdOmega);

        using BaseType::operator();
      };

      /// ([C])
      class Integral : public MDArray
      {

      public:
        typedef MDArray BaseType;

        Integral(const IM& im);

        /// wXdOmega: ([C], [P])
        /// iv:       ([C], [P])
        /// this:     ([C])
        void operator()(IntegrandValues& iv, WeightedMeasure& wXdOmega);
        using BaseType::operator();
      };


      /// ([C],[P],[DOF]): evaluated field values at each integration point in each cell: 
      class FieldValues : public MDArray
      {
      public:
        typedef MDArray BaseType;
        FieldValues(const IM& im);

        void operator()(const stk::mesh::BulkData & bulk, const stk::mesh::Entity element, MDArray& transformed_basis_values, stk::mesh::FieldBase* field);
        void operator()(const stk::mesh::BulkData & bulk, const stk::mesh::Entity element, MDArray& transformed_basis_values, stk::mesh::FieldBase* field, MDArray& output_field_values);
      };

      /// these are the "transformed_basis_values" at each point in each cell in the work set
      /// ([C],[B],[P]), or ([C],[B],[P],[D]) for GRAD
      /// here we assume that [B] is equivalent to [V]

      class Bases : public MDArray
      {
        ComputeBases m_cb;
      public:
        typedef MDArray BaseType;
        Bases(const IM& im);

        using BaseType::operator();

        void operator()(const stk::mesh::BulkData& bulk, const stk::mesh::Entity element, const MDArray& parametric_coordinates);
        void operator()(const stk::mesh::BulkData& bulk, const stk::mesh::Bucket& bucket, const MDArray& parametric_coordinates);
        
      };


      //--------------------------------------------------------------------------------------------------------------------------------------------------
      //--------------------------------------------------------------------------------------------------------------------------------------------------
      //--------------------------------------------------------------------------------------------------------------------------------------------------
      IntrepidManager(Elements_Tag el, Cub_Points_Tag ct, NodesPerElem_Tag nc, Spatial_Dim_Tag st, DOFs_Tag dt);

      IntrepidManager(Elements_Tag el, CellTopology& cellTopo, unsigned cubDegree = 2);

      void setupCubature(CellTopology& cellTopo, unsigned cubDegree=2);

      static void isInElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates, unsigned& found_it, const stk::mesh::Entity element,
                              const stk::mesh::BulkData& bulkData);
      static void more_template_instantiations();

      //void print() {}
      Elements_Tag m_Elements_Tag;
      Cub_Points_Tag m_Cub_Points_Tag;
      NodesPerElem_Tag m_NodesPerElem_Tag;
      Spatial_Dim_Tag m_Spatial_Dim_Tag;
      DOFs_Tag m_DOFs_Tag;

      // BasisFields_Tag m_BasisFields_Tag;

      CellTopology *m_topo;
      Teuchos::RCP<Cubature<Kokkos::HostSpace, double, double > > m_cub;

    };

  }

//MPerego: Is Rank used at all?

template<class entity>
struct Rank;

template<>
struct Rank<percept::IntrepidManager::CubaturePoints>{
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Cub_Points_Tag );
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Spatial_Dim_Tag );
static const int value=shards::ArrayVector<double, shards::NaturalOrder, Cub_Points_Tag, Spatial_Dim_Tag>::Rank;
//static const int value=2;
//percept::IntrepidManager::CubaturePoints::Rank;
};
template<>
struct Rank<percept::IntrepidManager::WeightedMeasure>{
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Elements_Tag );
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Cub_Points_Tag );
static const int value=shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag>::Rank;		
//static const int value=2;
//percept::IntrepidManager::CubaturePoints::Rank;
};
template<>
struct Rank<percept::IntrepidManager::Jacobian>{
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Elements_Tag );
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Cub_Points_Tag );
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Spatial_Dim_Tag );
static const int value=shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag, Spatial_Dim_Tag, Spatial_Dim_Tag>::Rank;	
//static const int value=3;
//percept::IntrepidManager::CubaturePoints::Rank;
};
template<>
struct Rank<percept::IntrepidManager::CellWorkSet>{
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Elements_Tag );
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( NodesPerElem_Tag );
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Spatial_Dim_Tag );
static const int value=shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, NodesPerElem_Tag, Spatial_Dim_Tag>::Rank;
//static const int value=3;
//percept::IntrepidManager::CubaturePoints::Rank;
};
template<>
struct Rank<percept::IntrepidManager::JacobianDet>{
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Elements_Tag );
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Cub_Points_Tag );
static const int value=shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag>::Rank;	
//static const int value=2;
//percept::IntrepidManager::CubaturePoints::Rank;
};
template<>
struct Rank<percept::IntrepidManager::CubatureWeights>{

IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Cub_Points_Tag );

static const int value=shards::ArrayVector<double, shards::NaturalOrder, Cub_Points_Tag>::Rank;	
//static const int value=2;
//percept::IntrepidManager::CubaturePoints::Rank;
};
template<>
struct Rank<percept::IntrepidManager::JacobianInverse>{
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Cub_Points_Tag );
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Elements_Tag );
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Spatial_Dim_Tag );
static const int value=shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag, Spatial_Dim_Tag, Spatial_Dim_Tag>::Rank;
};
template<>
struct Rank<percept::IntrepidManager::PhysicalCoords>{
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Cub_Points_Tag );
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Elements_Tag );
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Spatial_Dim_Tag );
static const int value=shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag, Spatial_Dim_Tag>::Rank;
};
template<>
struct Rank<percept::IntrepidManager::Integral>{
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Cub_Points_Tag );
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Elements_Tag );
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Spatial_Dim_Tag );
static const int value=shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag>::Rank;
};
template<>
struct Rank<percept::IntrepidManager::IntegrandValues>{
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Cub_Points_Tag );
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Elements_Tag );
IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( Spatial_Dim_Tag );
static const int value=shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag>::Rank;
};
#endif
