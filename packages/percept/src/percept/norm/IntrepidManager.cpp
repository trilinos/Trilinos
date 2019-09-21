// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <stdexcept>

// Intrepid includes

#define PGI_INSTANTIATION_FILE
#include <percept/Percept.hpp>


#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"

#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid_Basis.hpp"
#include <percept/Intrepid_HGRAD_WEDGE_C2_Serendipity_FEM.hpp>
#include <percept/Intrepid_HGRAD_QUAD_C2_Serendipity_FEM.hpp>
#include <percept/Intrepid_HGRAD_HEX_C2_Serendipity_FEM.hpp>


#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

#include "IntrepidManager.hpp"

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <Shards_CellTopology.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <percept/Util.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/function/Function.hpp>
#include <percept/function/ElementOp.hpp>

#include <percept/function/internal/ComputeBases.hpp>
#include <percept/function/internal/ComputeFieldValues.hpp>

  namespace percept
  {

    IM_SHARDS_ARRAY_DIM_TAG_IMPLEMENTATION( Elements_Tag )
    IM_SHARDS_ARRAY_DIM_TAG_IMPLEMENTATION( Cub_Points_Tag )
    IM_SHARDS_ARRAY_DIM_TAG_IMPLEMENTATION( NodesPerElem_Tag )
    IM_SHARDS_ARRAY_DIM_TAG_IMPLEMENTATION( Spatial_Dim_Tag )
    IM_SHARDS_ARRAY_DIM_TAG_IMPLEMENTATION( DOFs_Tag )   // [F]  FIXME
    IM_SHARDS_ARRAY_DIM_TAG_IMPLEMENTATION( BasisFields_Tag )   // [B]  FIXME

#if (defined(__PGI) && defined(USE_PGI_7_1_COMPILER_BUG_WORKAROUND))
      // workaround for PGI compiler bug
    void IntrepidManager::bootstrap()
    {
      static BasisTypeRCP a1 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_HEX_C1_FEM<double, MDArray >() );

      // FIXME
      static BasisTypeRCP a2 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_LINE_C1_FEM<double, MDArray >() );
      static BasisTypeRCP a3 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TRI_C1_FEM<double, MDArray >() );
      static BasisTypeRCP a4 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TRI_C2_FEM<double, MDArray >() );

      static BasisTypeRCP a5 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, MDArray >() );
      static BasisTypeRCP a6 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM<double, MDArray >() );
      static BasisTypeRCP a7 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_QUAD_C2_FEM<double, MDArray >() );

      static BasisTypeRCP a8 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_HEX_C1_FEM<double, MDArray >() );
      static BasisTypeRCP a9 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_HEX_C2_Serendipity_FEM<double, MDArray >() );
      static BasisTypeRCP a10 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_HEX_C2_FEM<double, MDArray >() );

      static BasisTypeRCP a11 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TET_C1_FEM<double, MDArray >() );
      static BasisTypeRCP a12 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TET_C2_FEM<double, MDArray >() );

      static BasisTypeRCP a13 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_WEDGE_C1_FEM<double, MDArray >() );

      // Intrepid doesn't support wedge 15
      static BasisTypeRCP a14 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_WEDGE_C2_Serendipity_FEM<double, MDArray >() );


      // Shells
      static BasisTypeRCP a15 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TRI_C1_FEM<double, MDArray >() );
      static BasisTypeRCP a16 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TRI_C2_FEM<double, MDArray >() );
      static BasisTypeRCP a17 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, MDArray >() );
      static BasisTypeRCP a18 = Teuchos::rcp ( new Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM<double, MDArray >() );
    }
#endif

    void tni(void)
    {
      throw std::runtime_error("not implemented");
    }

    void IntrepidManager::setupCubature( CellTopology& cell_topo, unsigned cubDegree)
    {
      DefaultCubatureFactory<double, CubaturePoints, CubatureWeights> cubFactory;         // create cubature factory
      m_cub = cubFactory.create(cell_topo, cubDegree);                                    // create default cubature
      unsigned numCubPoints = m_cub->getNumPoints();                                      // retrieve number of cubature points
      m_Cub_Points_Tag = Cub_Points_Tag(numCubPoints);
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([P],[D])
    IntrepidManager::CubaturePoints::CubaturePoints(IM& im) : BaseType( NUM(Cub_Points_Tag), NUM(Spatial_Dim_Tag) ), m_im(im)
    {
    }

#if 1
    //using BaseBaseType::operator();

    double& IntrepidManager::CubaturePoints::operator()(int i1, int i2, int i3)
    {
      throw std::runtime_error("CubaturePoints:: operator()(int i1, int i2, int i3) not implemented");
      return m_dummy;
    }
    const double& IntrepidManager::CubaturePoints::operator()(int i1, int i2, int i3) const {
      throw std::runtime_error("CubaturePoints:: operator()(int i1, int i2, int i3) not implemented");
      return m_dummy;
    }

    double& IntrepidManager::CubaturePoints::operator()(int i1, int i2) { return BaseBaseType::operator()(i1, i2); }
    const double& IntrepidManager::CubaturePoints::operator()(int i1, int i2) const { return BaseBaseType::operator()(i1, i2); }
#endif
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([P])
    IntrepidManager::CubatureWeights::   CubatureWeights(IM& im) : BaseType( NUM(Cub_Points_Tag))
    {
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([C], [V], [D])
    IntrepidManager::CellWorkSet::
    CellWorkSet(IM& im) : BaseType(NUM(Elements_Tag), NUM(NodesPerElem_Tag), NUM(Spatial_Dim_Tag))
    {

    }

#if 0
    void
    IntrepidManager::CellWorkSet::
    operator()(BulkData& bulkData, Bucket& bucket)
    {
      // FIXME
    }
#endif
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([C], [P], [D])
    IntrepidManager::PhysicalCoords::
    PhysicalCoords(IM& im) : BaseType( NUM(Elements_Tag), NUM(Cub_Points_Tag), NUM(Spatial_Dim_Tag) ), m_im(im)
    {
    }

    void
    IntrepidManager::PhysicalCoords::
    operator()(CellWorkSet& c, CubaturePoints& xi)
    {
      //FieldContainer<double> images(1, cellDim );
      //Intrepid::CellTools<double>::mapToPhysicalFrame(images, preImages, triNodes, triangle_3, whichCell);
      CellTools<double>::mapToPhysicalFrame(*this, xi, c, *m_im.m_topo, -1);
    }
#if 1
    //using BaseBaseType::operator();

    //double m_dummy;
#if 0
    double&
    IntrepidManager::PhysicalCoords::
    operator()(int i1, int i2, int i3)
    {
      throw std::runtime_error("PhysicalCoords:: operator()(int i1, int i2, int i3) not implemented");
      return m_dummy;
    }
    const double&
    IntrepidManager::PhysicalCoords::

    operator()(int i1, int i2, int i3) const {
      throw std::runtime_error("PhysicalCoords:: operator()(int i1, int i2, int i3) not implemented");
      return m_dummy;
    }
#endif

    double&     IntrepidManager::PhysicalCoords::
    operator()(int i1, int i2, int i3) { return BaseBaseType::operator()(i1, i2, i3); }
    const double&     IntrepidManager::PhysicalCoords::
    operator()(int i1, int i2, int i3) const { return BaseBaseType::operator()(i1, i2, i3); }
#endif

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
      /// ([C], [P], [D], [D])
    IntrepidManager::Jacobian::
    Jacobian(IM& im) : BaseType(NUM(Elements_Tag), NUM(Cub_Points_Tag), NUM(Spatial_Dim_Tag), NUM(Spatial_Dim_Tag)), m_im(im)
    {

    }
    void
    IntrepidManager::Jacobian::

    operator()(CubaturePoints& xi, CellWorkSet& c, CellTopology& topo)
    {
      CellTools<double>::setJacobian(*this, xi, c, topo);           // compute cell Jacobians
    }
    //double m_dummy;
    double&     IntrepidManager::Jacobian::
    operator()(int i1, int i2, int i3) { return m_dummy;}
    const double&     IntrepidManager::Jacobian::
    operator()(int i1, int i2, int i3) const { return m_dummy;}


    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
      /// ([C], [P], [D])
    IntrepidManager::FaceNormal::
    FaceNormal(IM& im) : BaseType(NUM(Elements_Tag), NUM(Cub_Points_Tag), NUM(Spatial_Dim_Tag)), m_im(im)
    {

    }
    void
    IntrepidManager::FaceNormal::
    operator()(Jacobian& J, int i_face, CellTopology& topo)
    {
      // FIXME - don't instantiate this, since Intrepid isn't completely free of its own MDArray, FieldContainer
#if 0
      MDArray J_mda;
      J.copyTo(J_mda);
      MDArray fn_mda(m_im.m_Elements_Tag.num, m_im.m_Cub_Points_Tag.num, m_im.m_Spatial_Dim_Tag.num);
      CellTools<double>::getPhysicalFaceNormals(fn_mda, J_mda, i_face, cell_topo);
      this->copyFrom(fn_mda);
#endif
      //CellTools<double>::getPhysicalFaceNormals(*this, jac, i_face, topo);
      throw std::runtime_error(" don't instantiate this, since Intrepid isn't completely free of its own MDArray, FieldContainer");
    }
    //double m_dummy;
    double&     IntrepidManager::FaceNormal::
    operator()(int i1, int i2) { return m_dummy;}
    const double&     IntrepidManager::FaceNormal::
    operator()(int i1, int i2) const { return m_dummy;}



    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([C], [P], [D], [D])
    IntrepidManager::JacobianInverse::
    JacobianInverse(IM& im) : BaseType(NUM(Elements_Tag), NUM(Cub_Points_Tag), NUM(Spatial_Dim_Tag), NUM(Spatial_Dim_Tag))
    {

    }
    void    IntrepidManager::JacobianInverse::
    operator()(Jacobian& jac)
    {
      CellTools<double>::setJacobianInv(*this, jac);
    }


    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
      /// ([C], [P])
    IntrepidManager::JacobianDet::
    JacobianDet(IM& im) : BaseType(NUM(Elements_Tag), NUM(Cub_Points_Tag))
    {
    }

    void
    IntrepidManager::JacobianDet::
    operator()(Jacobian& jac)
    {
      CellTools<double>::setJacobianDet(*this, jac);
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
      /// ([C], [P])
    IntrepidManager::WeightedMeasure::
    WeightedMeasure(IM& im) : BaseType(NUM(Elements_Tag), NUM(Cub_Points_Tag))
    {
    }

    void
    IntrepidManager::WeightedMeasure::
    operator()(CubatureWeights& w, JacobianDet& dJ)
    {
      FunctionSpaceTools::computeCellMeasure<double>(*this, dJ, w);
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
      /// ([C], [P], [DOF])
    IntrepidManager::IntegrandValuesDOF::
    IntegrandValuesDOF(IM& im) : BaseType(NUM(Elements_Tag), NUM(Cub_Points_Tag), NUM(DOFs_Tag))
    {
    }
    void
    IntrepidManager::IntegrandValuesDOF::
    copyFrom(MDArray& mda)
    {
      copy(&mda[0], &mda[0] + mda.size(), this->contiguous_data());
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
      /// ([C], [P])
    IntrepidManager::IntegrandValues::
    IntegrandValues(IM& im) : BaseType(NUM(Elements_Tag), NUM(Cub_Points_Tag))
    {
    }
    void
    IntrepidManager::IntegrandValues::
    copyFrom(MDArray& mda)
    {
      copy(&mda[0], &mda[0] + mda.size(), this->contiguous_data());
    }
    void
    IntrepidManager::IntegrandValues::
    copyFrom(IntrepidManager& im, MDArray& mda, int iDof)
    {
      for (int iCell = 0; iCell < im.m_Elements_Tag.num; iCell++)
        {
          for (int iPoint = 0; iPoint < im.m_Cub_Points_Tag.num; iPoint++)
            {
              (*this)(iCell, iPoint) = mda(iCell, iPoint, iDof);
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([C], [DOF])
    IntrepidManager::IntegralDOF::
    IntegralDOF(IM& im) : BaseType(NUM(Elements_Tag), NUM(DOFs_Tag))
    {
    }

    /// wXdOmega: ([C], [P])
    /// iv:       ([C], [P], [DOF])
    /// this:     ([C], [DOF])
    void
    IntrepidManager::IntegralDOF::
    operator()(IntegrandValuesDOF& iv, WeightedMeasure& wXdOmega, int comp_type)
    {
      VERIFY_OP(iv.rank(),       == , 3,                     "IntrepidManager::Integral::operator() bad iv rank");
      VERIFY_OP(wXdOmega.rank(), == , 2,                     "IntrepidManager::Integral::operator() bad wXdOmega rank");
      VERIFY_OP((*this).rank(),  == , 2,                     "IntrepidManager::Integral::operator() bad (*this) rank");
      VERIFY_OP(iv.dimension(0), == , this->dimension(0),    "IntrepidManager::Integral::operator() bad");
      VERIFY_OP(iv.dimension(1), == , wXdOmega.dimension(1), "IntrepidManager::Integral::operator() bad");
      VERIFY_OP(iv.dimension(0), == , wXdOmega.dimension(0), "IntrepidManager::Integral::operator() bad");
      VERIFY_OP(iv.dimension(2), == , this->dimension(1),    "IntrepidManager::Integral::operator() bad");

#if 0
      for (int iDof = 0; iDof < iv.dimension(2); iDof++)
        {
          Integral Is(...);
          Is(iv, wXdOmega, COMP_BLAS);
          FunctionSpaceTools::integrate<double>(*this, iv, wXdOmega,  COMP_BLAS);
        }
#endif
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([C])
    IntrepidManager::Integral::
    Integral(IM& im) : BaseType(NUM(Elements_Tag))
    {
    }

    /// wXdOmega: ([C], [P])
    /// iv:       ([C], [P])
    /// this:     ([C])
    void
    IntrepidManager::Integral::
    operator()(IntegrandValues& iv, WeightedMeasure& wXdOmega, int comp_type)
    {
      VERIFY_OP(iv.rank(),       == , 2,                     "IntrepidManager::Integral::operator() bad iv rank");
      VERIFY_OP(wXdOmega.rank(), == , 2,                     "IntrepidManager::Integral::operator() bad wXdOmega rank");
      VERIFY_OP((*this).rank(),  == , 1,                     "IntrepidManager::Integral::operator() bad (*this) rank");
      VERIFY_OP(iv.dimension(0), == , this->dimension(0),    "IntrepidManager::Integral::operator() bad");
      VERIFY_OP(iv.dimension(1), == , wXdOmega.dimension(1), "IntrepidManager::Integral::operator() bad");
      VERIFY_OP(iv.dimension(0), == , wXdOmega.dimension(0), "IntrepidManager::Integral::operator() bad");

      FunctionSpaceTools::integrate<double>(*this, iv, wXdOmega,  COMP_BLAS);
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([C],[B],[P]), or ([C],[B],[P],[D]) for GRAD
    static ComputeBases s_compute_bases;

    IntrepidManager::Bases::
    Bases(IM& im) : BaseType(NUM(Elements_Tag), NUM(NodesPerElem_Tag), NUM(Cub_Points_Tag))
                  , m_cb(ComputeBases()) 
         //, m_cb(new ComputeBases()) 
    {
      //m_cb =&s_compute_bases;
    }

    void
    IntrepidManager::Bases::
    operator()(const stk::mesh::BulkData& bulk, const stk::mesh::Entity element, const MDArray& parametric_coordinates)
    {
      m_cb.getBases(bulk, element, parametric_coordinates, *this);
    }

    void
    IntrepidManager::Bases::
    operator()(const stk::mesh::BulkData& bulk, const stk::mesh::Bucket& bucket, const MDArray& parametric_coordinates)
    {
      m_cb.getBases(bucket, parametric_coordinates, *this);
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([C],[P],[DOF]): evaluated field values at each integration point in each cell:
    IntrepidManager::FieldValues::
    FieldValues(IM& im) : BaseType(NUM(Elements_Tag), NUM(Cub_Points_Tag), NUM(DOFs_Tag)) {}

    void IntrepidManager::FieldValues::operator()(const stk::mesh::BulkData& bulk, const stk::mesh::Entity element, MDArray& transformed_basis_values, stk::mesh::FieldBase* field)
    {
      ComputeFieldValues cfv;
      cfv.get_fieldValues(bulk, element, transformed_basis_values, field, *this);
    }
    void IntrepidManager::FieldValues::operator()(const stk::mesh::BulkData& bulk, const stk::mesh::Entity element, MDArray& transformed_basis_values, stk::mesh::FieldBase* field, MDArray& output_field_values)
    {
      ComputeFieldValues cfv;
      cfv.get_fieldValues(bulk, element, transformed_basis_values, field, output_field_values);
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------

    IntrepidManager::IntrepidManager(Elements_Tag el, Cub_Points_Tag ct, NodesPerElem_Tag nc, Spatial_Dim_Tag st,
                                     DOFs_Tag dt) :
      m_Elements_Tag(el.num), m_Cub_Points_Tag(ct.num), m_NodesPerElem_Tag(nc.num), m_Spatial_Dim_Tag(st.num),  m_DOFs_Tag(dt.num),
      m_topo(0)
    {
    }

    IntrepidManager::IntrepidManager(Elements_Tag el, CellTopology& cellTopo, unsigned cubDegree)
      : m_Elements_Tag(el.num), m_Cub_Points_Tag(1), m_NodesPerElem_Tag(1), m_Spatial_Dim_Tag(1), m_DOFs_Tag(1),
        m_topo(&cellTopo)
    {
      setupCubature(cellTopo, cubDegree);
      m_NodesPerElem_Tag.num = cellTopo.getNodeCount();
      m_Spatial_Dim_Tag.num = cellTopo.getDimension();
    }





#if 0
    void IntrepidManager::getCubature(CubaturePoints& cp, CubatureWeights& cw)
    {
      m_cub->getCubature(cp, cw);
    }
#endif

    void test()
    {
      typedef IntrepidManager IM;
      IntrepidManager im(Elements_Tag(1), Cub_Points_Tag(4), NodesPerElem_Tag(4), Spatial_Dim_Tag(3), DOFs_Tag(1));

      CellTopology* topo=0;

      IM::Jacobian          J  (im);
      IM::JacobianInverse  Ji  (im);
      IM::JacobianDet      dJ  (im);

      IM::CubaturePoints   xi  (im);
      IM::CellWorkSet       c  (im);
      IM::CubatureWeights   w  (im);

      IM::PhysicalCoords   pc  (im);

      im.m_cub->getCubature(xi, w);
      pc(c, xi);

#if 0
      IM::WeightedMeasure wXdOmega(im);

      wXdOmega = w * dJ;
#endif

#if 0
      IM::GradField gradP("pressure");
      IM::TransformedGradField tgradP(gradP);
      IM::WeightedTransformedGrad wtgradP(tgradP);

      IM::GradField gradN_i("basis");
      IM::TransformedGradField tgradN_i(gradN_i);
      IM::WeightedTransformedGrad wtgradN_i(tgradN_i);
      IM::StiffnessMatrix K;

      IM::WeightedMeasure wXdOmega;
      IM::HGRADPhysicalSpace tgradN_i;
      IM::MultipliedByMeasure wXdOmegaTGradN_i;
      IM::IntegratedOp K( COMP_CPP );
#endif

      //IM::JacobianInv ji(im);
      //std::cout << "j.size()= " << j.size() << std::endl;
      //im.print();


#if 0
      xi.set();
      w.set();
      gradN_i(xi);
#endif
      //Step 6: Apply cell tools
      J(xi, c, *topo);
      Ji(J);
      dJ(J);

#if 0
      //Step 7: Apply function space tools

      wXdOmega = w * dJ;

      tgradN_i = Ji * gradN_i;

      wXdOmegaTGradN_i = wXdOmega * tgradN_i;

      K =  tgradN_i * wXdOmegaTGradN_i;
#endif

    }

    void IntrepidManager::isInElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates, unsigned& found_it, const stk::mesh::Entity element,
                                      const stk::mesh::BulkData& bulkData)
    {
      found_it = 0;

      // FIXME consider caching the coords_field in FieldFunction
      const stk::mesh::MetaData& metaData = stk::mesh::MetaData::get(bulkData);
      CoordinatesFieldType *coords_field = metaData.get_field<CoordinatesFieldType >(stk::topology::NODE_RANK, "coordinates");

      const stk::mesh::Bucket & bucket = bulkData.bucket(element);
      const CellTopologyData * const bucket_cell_topo_data = stk::mesh::get_cell_topology(bucket.topology()).getCellTopologyData();
      if (!bucket_cell_topo_data)
        {
          stk::mesh::EntityRank bucket_rank = bucket.entity_rank();
          std::cout << "bucket_rank = " << bucket_rank << std::endl;
          throw std::runtime_error("IntrepidManager::bogus topology");
        }

      unsigned numCells = 1; // FIXME

      shards::CellTopology topo(bucket_cell_topo_data);
      unsigned numNodes = topo.getNodeCount();
      unsigned cellDim  = topo.getDimension();
      MDArray cellWorkset(numCells, numNodes, cellDim);

      /// FIXME -- fill cellWorkset
      const MyPairIterRelation elem_nodes(bulkData, element, stk::topology::NODE_RANK );


      for (unsigned iCell = 0; iCell < numCells; iCell++)
        {
          for (unsigned iNode = 0; iNode < numNodes; iNode++)
            {
              stk::mesh::Entity node = elem_nodes[iNode].entity();
              double * node_coord_data = stk::mesh::field_data( *coords_field , node);
              for (unsigned iDim=0; iDim < cellDim; iDim++)
                {
                  cellWorkset(iCell, iNode, iDim) = node_coord_data[iDim];
                }
            }
        }

      // FIXME for multiple points
      if (input_phy_points.rank() == 1)
        {
          VERIFY_1("IsInElement::isInElement bad rank of input_phy_points");
        }
      VERIFY_OP(input_phy_points.dimension(0), == , 1, "IsInElement::isInElement bad input_phy_points 1st dim");
      if (input_phy_points.dimension(1) < (int)cellDim)
        {
          std::cout << "IsInElement::isInElement bad input_phy_points 2nd dim";
        }
      VERIFY_OP(input_phy_points.dimension(1), >= , (int)cellDim, "IsInElement::isInElement bad input_phy_points 2nd dim");

      if (found_parametric_coordinates.rank() == 1)
        {
          VERIFY_1("IsInElement::isInElement bad rank of found_parametric_coordinates");
        }
      VERIFY_OP(found_parametric_coordinates.dimension(0), == , 1, "IsInElement::isInElement bad found_parametric_coordinates 1st dim");
      VERIFY_OP(found_parametric_coordinates.dimension(1), == , (int)cellDim,
                "IsInElement::isInElement bad found_parametric_coordinates 2nd dim");

      unsigned cellOrd = 0;  // FIXME
      Intrepid::CellTools<double>::mapToReferenceFrame(found_parametric_coordinates, input_phy_points, cellWorkset, topo, cellOrd);
      MDArrayUInt inclusion_results(1);  // FIXME
      double threshold = 1.e-4; // (INTREPID_THRESHOLD default = 10*double_eps ~ 20e-16)
      Intrepid::CellTools<double>::checkPointwiseInclusion(inclusion_results, found_parametric_coordinates, topo, threshold);
      found_it = inclusion_results(0);
      if (found_it)
        {
          // for testing only
          if (0)
            {
              FieldContainer<double> images(1, cellDim );
              //Intrepid::CellTools<double>::mapToPhysicalFrame(images, preImages, triNodes, triangle_3, whichCell);
              Intrepid::CellTools<double>::mapToPhysicalFrame(images, found_parametric_coordinates, cellWorkset, topo, cellOrd);
            }
        }
    }

    void IntrepidManager::more_template_instantiations()
    {
      throw std::runtime_error("NEVER invoke this method - for REDSTORM template instantiation only");
      CellTopology *cell_topo_ptr = 0;
      CellTopology& cell_topo = *cell_topo_ptr;
      unsigned numCells = 1;
      unsigned numNodes = 1;
      unsigned spaceDim = 1;

      // Rank-3 array with dimensions (C,N,D) for the node coordinates of 3 traingle cells
      FieldContainer<double> cellNodes(numCells, numNodes, spaceDim);

      DefaultCubatureFactory<double> cubFactory;                                              // create cubature factory
      unsigned cubDegree = 2;                                                                      // set cubature degree, e.g. 2
      Teuchos::RCP<Cubature<double> > myCub = cubFactory.create(cell_topo, cubDegree);         // create default cubature

      unsigned numCubPoints = myCub->getNumPoints();                                               // retrieve number of cubature points

      FieldContainer<double> cub_points(numCubPoints, spaceDim);
      FieldContainer<double> cub_weights(numCubPoints);

      // Rank-4 array (C,P,D,D) for the Jacobian and its inverse and Rank-2 array (C,P) for its determinant
      FieldContainer<double> jacobian(numCells, numCubPoints, spaceDim, spaceDim);
      FieldContainer<double> jacobian_inv(numCells, numCubPoints, spaceDim, spaceDim);
      FieldContainer<double> jacobian_det(numCells, numCubPoints);

      myCub->getCubature(cub_points, cub_weights);                                          // retrieve cubature points and weights

      // Methods to compute cell Jacobians, their inverses and their determinants

      CellTools<double>::setJacobian(jacobian, cub_points, cellNodes, cell_topo);           // compute cell Jacobians
      CellTools<double>::setJacobianInv(jacobian_inv, jacobian);                            // compute inverses of cell Jacobians
      CellTools<double>::setJacobianDet(jacobian_det, jacobian);                            // compute determinants of cell Jacobians

      FieldContainer<double> weightedMeasure(numCells, numCubPoints);
      FieldContainer<double> onesLeft(numCells,  numCubPoints);
      FieldContainer<double> volume(numCells);

      // compute weighted measure
      FunctionSpaceTools::computeCellMeasure<double>(weightedMeasure, jacobian_det, cub_weights);

      // integrate to get volume
      FunctionSpaceTools::integrate<double>(volume, onesLeft, weightedMeasure,  COMP_BLAS);

    }



  }
