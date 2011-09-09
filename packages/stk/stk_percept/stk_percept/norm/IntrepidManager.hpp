/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef stk_percept_IntrepidManager_hpp
#define stk_percept_IntrepidManager_hpp

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

#include <math.h>

#include <stk_percept/Percept.hpp>

#include "Intrepid_DefaultCubatureFactory.hpp"

#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/MDArray.hpp>
#include <stk_percept/PerceptMesh.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include "Teuchos_RCP.hpp"


// Shards includes
#include "Shards_CellTopology.hpp"
#include <stk_percept/function/internal/ComputeBases.hpp>

#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_ArrayView.hpp"
// #include "Shards_Array.hpp"
#include "Teuchos_RCP.hpp"
// #include "Teuchos_BLAS.hpp"
#include "Teuchos_oblackholestream.hpp"
//#include "Teuchos_TestForException.hpp"


using namespace std;
using namespace Intrepid;

#define IM_TAG( ADT ) ADT ## _TAG

#define IM_SHARDS_ARRAY_DIM_TAG_DECLARATION( ADT )  \
  class ADT  : public shards::ArrayDimTag {         \
  public:                                           \
    const char * name() const ;                     \
    static const ADT & tag();                       \
    static int num;                                 \
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
  int ADT::num = 0;                                                     \
  const char * ADT::name() const                                        \
  { static const char n[] = # ADT ; return n; }                         \
  const ADT & ADT::tag() { static const ADT self ; return self ; }      \
  ADT::ADT(int n) { num =n;}                                            \
  ADT::  ~ADT() {}                                                      \
  ADT::ADT( const ADT & adt) { num=adt.num;}                            \
  ADT & ADT::operator = ( const ADT & adt) {num=adt.num; return *this;} \
  ADT::ADT() {}                                                         


namespace stk
{
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
     *  Note: Intrepid really doesn't have a concept of "DOF" at a node.  It's either a single variable,
     *    or a vector- or tensor-valued variable.  So, no DOF-related arrays as used herein can be used
     *    with Intrepid - you must call Intrepd one DOF at a time.
     *
     * FieldContainer<double> cub_points(numCubPoints, spaceDim);
     * FieldContainer<double> cub_weights(numCubPoints);
     * 
     * FieldContainer<double> cell_nodes(numCells, numNodes, spaceDim);
     * 
     * FieldContainer<double> jacobian(numCells, numCubPoints, spaceDim, spaceDim);
     * FieldContainer<double> jacobian_inv(numCells, numCubPoints, spaceDim, spaceDim);
     * FieldContainer<double> jacobian_det(numCells, numCubPoints);
     * FieldContainer<double> weighted_measure(numCells, numCubPoints);
     * 
     * FieldContainer<double> grad_at_cub_points(numFields, numCubPoints, spaceDim);
     * FieldContainer<double> transformed_grad_at_cub_points(numCells, numFields, numCubPoints, spaceDim);
     * FieldContainer<double> weighted_transformed_grad_at_cub_points(numCells, numFields, numCubPoints, spaceDim);
     * FieldContainer<double> stiffness_matrices(numCells, numFields, numFields);
     * 
     * 
     */

    class IntrepidManager
    {
    private:

#if (defined(__PGI) && defined(USE_PGI_7_1_COMPILER_BUG_WORKAROUND))
      // workaround for PGI compiler bug

      typedef Intrepid::Basis<double, MDArray > BasisType;
      typedef Teuchos::RCP<BasisType>           BasisTypeRCP;

      static void bootstrap();
#endif

    public:
      typedef IntrepidManager IM;
  
#define NUM(AClass) im.m_ ## AClass  . num

      //hexJacobian(numCells, numCubPoints, spaceDim, spaceDim);


      class temp
      {
        double m_dummy;
      public:

        double&       operator()(int i1, int i2, int i3)        { tni(); return m_dummy; }
        const double& operator()(int i1, int i2, int i3) const  { tni(); return m_dummy; }

        double&       operator()(int i1, int i2, int i3, int i4)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2, int i3, int i4) const { tni(); return m_dummy;}

        double&       operator()(int i1, int i2, int i3, int i4, int i5)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2, int i3, int i4, int i5) const { tni(); return m_dummy;}

      };

      /// ([P],[D])
      class CubaturePoints : public shards::ArrayVector<double, shards::NaturalOrder, Cub_Points_Tag, Spatial_Dim_Tag>
      {
        IntrepidManager& m_im;
      public:
        typedef shards::ArrayVector<double, shards::NaturalOrder, Cub_Points_Tag, Spatial_Dim_Tag> BaseType;
        typedef shards::Array<double, shards::NaturalOrder, Cub_Points_Tag, Spatial_Dim_Tag> BaseBaseType;

        CubaturePoints(IM& im) ;

        void copyTo(MDArray& mda)
        {
          mda.resize(m_im.m_Cub_Points_Tag.num, m_im.m_Spatial_Dim_Tag.num);
          mda.setValues(contiguous_data(), size());
        }
#if 1
        //using BaseBaseType::operator();

        double m_dummy;
        double&       operator()(int i1, int i2)       ; //{ tni(); return m_dummy; }
        const double& operator()(int i1, int i2) const ; //{ tni(); return m_dummy; }

        double&       operator()(int i1, int i2, int i3); //       { tni(); return m_dummy; }
        const double& operator()(int i1, int i2, int i3) const; // { tni(); return m_dummy; }
#endif
      };

      /// ([P])
      class CubatureWeights : public shards::ArrayVector<double, shards::NaturalOrder, Cub_Points_Tag>
      {

      public:
        typedef shards::ArrayVector<double, shards::NaturalOrder, Cub_Points_Tag> BaseType;

        CubatureWeights(IM& im);

#if 1
        double m_dummy;

        double&       operator()(int i1, int i2)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2) const { tni(); return m_dummy;}

        double&       operator()(int i1, int i2, int i3)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2, int i3) const { tni(); return m_dummy;}

        double&       operator()(int i1, int i2, int i3, int i4)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2, int i3, int i4) const { tni(); return m_dummy;}

        using BaseType::operator();
#endif
      };

      /// ([C], [V], [D])
      class CellWorkSet : public shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, NodesPerElem_Tag, Spatial_Dim_Tag>
      {
      public:
        typedef shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, NodesPerElem_Tag, Spatial_Dim_Tag> BaseType;

        CellWorkSet(IM& im) ;
        using BaseType::operator();
        //void operator()(BulkData& bulkData, Bucket& bucket);
      };

      /// ([C], [P], [D])
      class PhysicalCoords : public shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag, Spatial_Dim_Tag>
      {

        IntrepidManager& m_im;
      public:
        typedef shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag, Spatial_Dim_Tag> BaseType;
        typedef shards::Array<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag, Spatial_Dim_Tag> BaseBaseType;

        PhysicalCoords(IM& im);

        void operator()(CellWorkSet& c, CubaturePoints& xi);
        
        void copyTo(MDArray& mda)
        {
          mda.resize(m_im.m_Elements_Tag.num, m_im.m_Cub_Points_Tag.num, m_im.m_Spatial_Dim_Tag.num);
          mda.setValues(contiguous_data(), size());
        }
#if 1
        //using BaseBaseType::operator();

        double m_dummy;
        double&       operator()(int i1, int i2)       { tni(); return m_dummy; }
        const double& operator()(int i1, int i2) const { tni(); return m_dummy; }

        double&       operator()(int i1, int i2, int i3)      ;// { tni(); return m_dummy; }
        const double& operator()(int i1, int i2, int i3) const ;//{ tni(); return m_dummy; }

#endif

      };

      /// ([C], [P], [D], [D])
      class Jacobian : public shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag, Spatial_Dim_Tag, Spatial_Dim_Tag>
      {

      public:
        typedef shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag, Spatial_Dim_Tag, Spatial_Dim_Tag> BaseType;

        Jacobian(IM& im);
        void operator()(CubaturePoints& xi, CellWorkSet& c, CellTopology& topo);

        double m_dummy;
        double&       operator()(int i1, int i2, int i3);
        const double& operator()(int i1, int i2, int i3) const;
        using BaseType::operator();
      };

      /// ([C], [P], [D], [D])
      class JacobianInverse : public shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag, Spatial_Dim_Tag, Spatial_Dim_Tag>
      {

      public:
        typedef shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag, Spatial_Dim_Tag, Spatial_Dim_Tag> BaseType;

        JacobianInverse(IM& im);
        void operator()(Jacobian& jac);
#if 0
        double m_dummy;
        double&       operator()(int i1, int i2, int i3) { tni(); return m_dummy;}
        const double& operator()(int i1, int i2, int i3) const { tni(); return m_dummy;}
        using BaseType::operator();
#endif
      };

      /// ([C], [P])
      class JacobianDet : public shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag>
      {

      public:
        typedef shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag> BaseType;

        //Jacobian(shards::shards::Array<Cub_Points_Tag>& xi, shards::Array<NodesPerElem_Tag>& c, Topology& topo) 
        JacobianDet(IM& im);

        void operator()(Jacobian& jac);
        using BaseType::operator();

        double m_dummy;
        double&       operator()(int i1, int i2, int i3);
        const double& operator()(int i1, int i2, int i3) const;

      };


      /// weights multiplied by Jacobian det at cubature points
      /// ([C], [P])
      class WeightedMeasure : public shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag>
      {

      public:
        typedef shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag> BaseType;

        WeightedMeasure(IM& im);

        void operator()(CubatureWeights& w, JacobianDet& dJ);
        using BaseType::operator();

#if 1
        double m_dummy;

        double&       operator()(int i1, int i2, int i3)        { tni(); return m_dummy; }
        const double& operator()(int i1, int i2, int i3) const  { tni(); return m_dummy; }

        double&       operator()(int i1, int i2, int i3, int i4)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2, int i3, int i4) const { tni(); return m_dummy;}

        double&       operator()(int i1, int i2, int i3, int i4, int i5)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2, int i3, int i4, int i5) const { tni(); return m_dummy;}
#endif
      };

      /// ([C], [P], [DOF])
      class IntegrandValuesDOF : public shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag, DOFs_Tag>
      {

      public:
        typedef shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag, DOFs_Tag> BaseType;

        IntegrandValuesDOF(IM& im);

        void copyFrom(MDArray& mda);

        using BaseType::operator();
        double m_dummy;

        double&       operator()(int i1, int i2)       { tni(); return m_dummy; }
        const double& operator()(int i1, int i2) const { tni(); return m_dummy; }

        //double&       operator()(int i1, int i2, int i3)       { tni(); return m_dummy; }
        //const double& operator()(int i1, int i2, int i3) const { tni(); return m_dummy; }

        double&       operator()(int i1, int i2, int i3, int i4)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2, int i3, int i4) const { tni(); return m_dummy;}

        double&       operator()(int i1, int i2, int i3, int i4, int i5)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2, int i3, int i4, int i5) const { tni(); return m_dummy;}

      };

      /// ([C], [P])
      class IntegrandValues : public shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag>
      {

      public:
        typedef shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag> BaseType;

        IntegrandValues(IM& im);

        void copyFrom(MDArray& mda);
        void copyFrom(IntrepidManager& im, MDArray& mda, int iDof);

        using BaseType::operator();
        double m_dummy;

        //double&       operator()(int i1, int i2)       { tni(); return m_dummy; }
        //const double& operator()(int i1, int i2) const { tni(); return m_dummy; }

        double&       operator()(int i1, int i2, int i3)       { tni(); return m_dummy; }
        const double& operator()(int i1, int i2, int i3) const { tni(); return m_dummy; }

        double&       operator()(int i1, int i2, int i3, int i4)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2, int i3, int i4) const { tni(); return m_dummy;}

        double&       operator()(int i1, int i2, int i3, int i4, int i5)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2, int i3, int i4, int i5) const { tni(); return m_dummy;}

      };

      /// ([C], [DOF])
      class IntegralDOF : public shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, DOFs_Tag>
      {

      public:
        typedef shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, DOFs_Tag > BaseType;

        IntegralDOF(IM& im);

        /// wXdOmega: ([C], [P])
        /// iv:       ([C], [P], [DOF])
        /// this:     ([C], [DOF])
        void operator()(IntegrandValuesDOF& iv, WeightedMeasure& wXdOmega, int comp_type);

        using BaseType::operator();

        double m_dummy;
        double&       operator()(int i1)       { tni(); return m_dummy;}
        const double& operator()(int i1) const { tni(); return m_dummy;}

//         double&       operator()(int i1, int i2)       { tni(); return m_dummy;}
//         const double& operator()(int i1, int i2) const { tni(); return m_dummy;}

        double&       operator()(int i1, int i2, int i3)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2, int i3) const { tni(); return m_dummy;}

        double&       operator()(int i1, int i2, int i3, int i4, int i5)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2, int i3, int i4, int i5) const { tni(); return m_dummy;}

      };

      /// ([C])
      class Integral : public shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag>
      {

      public:
        typedef shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag > BaseType;

        Integral(IM& im);

        /// wXdOmega: ([C], [P])
        /// iv:       ([C], [P])
        /// this:     ([C])
        void operator()(IntegrandValues& iv, WeightedMeasure& wXdOmega, int comp_type);
        using BaseType::operator();

        double m_dummy;
        //double&       operator()(int i1)       { tni(); return m_dummy;}
        //const double& operator()(int i1) const { tni(); return m_dummy;}

        double&       operator()(int i1, int i2)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2) const { tni(); return m_dummy;}

        double&       operator()(int i1, int i2, int i3)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2, int i3) const { tni(); return m_dummy;}

        double&       operator()(int i1, int i2, int i3, int i4, int i5)       { tni(); return m_dummy;}
        const double& operator()(int i1, int i2, int i3, int i4, int i5) const { tni(); return m_dummy;}

      };

      // FIXME - change to shards array
      //class FieldValues : public shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag>

      /// ([C],[P],[DOF]): evaluated field values at each integration point in each cell: 
      class FieldValues : public MDArray
      {
      public:
        typedef MDArray BaseType;
        FieldValues(IM& im);
        void operator()(const stk::mesh::Entity& element, MDArray& transformed_basis_values, mesh::FieldBase* field);
        void operator()(const stk::mesh::Entity& element, MDArray& transformed_basis_values, mesh::FieldBase* field, MDArray& output_field_values);
      };

      // FIXME - change to shards array
      /// these are the "transformed_basis_values" at each point in each cell in the work set
      /// ([C],[B],[P]), or ([C],[B],[P],[D]) for GRAD
      /// here we assume that [B] is equivalent to [V]

      class Bases : public MDArray
      {
        ComputeBases m_cb;
      public:
        typedef MDArray BaseType;
        Bases(IM& im);

        using BaseType::operator();

        void operator()(const stk::mesh::Entity& element, const MDArray& parametric_coordinates);
        void operator()(const stk::mesh::Bucket& bucket, const MDArray& parametric_coordinates);
        
      };


      //--------------------------------------------------------------------------------------------------------------------------------------------------
      //--------------------------------------------------------------------------------------------------------------------------------------------------
      //--------------------------------------------------------------------------------------------------------------------------------------------------
      IntrepidManager(Elements_Tag el, Cub_Points_Tag ct, NodesPerElem_Tag nc, Spatial_Dim_Tag st, DOFs_Tag dt);

      IntrepidManager(Elements_Tag el, CellTopology& cellTopo, unsigned cubDegree = 2);

      void setupCubature(CellTopology& cellTopo, unsigned cubDegree=2);

      static void isInElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates, unsigned& found_it, const mesh::Entity& element,
                              const mesh::BulkData& bulkData);
      static void more_template_instantiations();

      //void print() {}
      Elements_Tag m_Elements_Tag;
      Cub_Points_Tag m_Cub_Points_Tag;
      NodesPerElem_Tag m_NodesPerElem_Tag;
      Spatial_Dim_Tag m_Spatial_Dim_Tag;
      DOFs_Tag m_DOFs_Tag;

      // BasisFields_Tag m_BasisFields_Tag;

      CellTopology *m_topo;
      Teuchos::RCP<Cubature<double, CubaturePoints, CubatureWeights > > m_cub;

    };


    //std::ostream& operator<<(std::ostream& os, const shards::Array<double, shards::NaturalOrder>& container) ;
    

  }
}


#endif
