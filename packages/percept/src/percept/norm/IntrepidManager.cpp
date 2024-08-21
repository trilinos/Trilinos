// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <stdexcept>

// Intrepid2 includes

#define PGI_INSTANTIATION_FILE
#include <percept/Percept.hpp>


#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include <percept/element/intrepid/BasisTable.hpp>
#include "Intrepid2_RealSpaceTools.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_Utils.hpp"

#include "IntrepidManager.hpp"

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

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
      static BasisTypeRCP a1 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_HEX_C1_FEM<Kokkos::HostSpace, double, double >() );

      // FIXME
      static BasisTypeRCP a2 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_LINE_C1_FEM<Kokkos::HostSpace, double, double >() );
      static BasisTypeRCP a3 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_TRI_C1_FEM<Kokkos::HostSpace, double, double >() );
      static BasisTypeRCP a4 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_TRI_C2_FEM<Kokkos::HostSpace, double, double >() );

      static BasisTypeRCP a5 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<Kokkos::HostSpace, double, double >() );
      static BasisTypeRCP a6 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_QUAD_I2_FEM<Kokkos::HostSpace, double, double >() );
      static BasisTypeRCP a7 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_QUAD_C2_FEM<Kokkos::HostSpace, double, double >() );

      static BasisTypeRCP a8 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_HEX_C1_FEM<Kokkos::HostSpace, double, double >() );
      static BasisTypeRCP a9 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_HEX_I2_FEM<Kokkos::HostSpace, double, double >() );
      static BasisTypeRCP a10 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_HEX_C2_FEM<Kokkos::HostSpace, double, double >() );

      static BasisTypeRCP a11 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_TET_C1_FEM<Kokkos::HostSpace, double, double >() );
      static BasisTypeRCP a12 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_TET_C2_FEM<Kokkos::HostSpace, double, double >() );

      static BasisTypeRCP a13 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_WEDGE_C1_FEM<Kokkos::HostSpace, double, double >() );
      static BasisTypeRCP a14 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_WEDGE_I2_FEM<Kokkos::HostSpace, double, double >() );


      // Shells
      static BasisTypeRCP a15 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_TRI_C1_FEM<Kokkos::HostSpace, double, double >() );
      static BasisTypeRCP a16 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_TRI_C2_FEM<Kokkos::HostSpace, double, double >() );
      static BasisTypeRCP a17 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<Kokkos::HostSpace, double, double >() );
      static BasisTypeRCP a18 = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_QUAD_I2_FEM<Kokkos::HostSpace, double, double >() );
    }
#endif

    void tni(void)
    {
      throw std::runtime_error("not implemented");
    }

    void IntrepidManager::setupCubature( CellTopology& cell_topo, unsigned cubDegree)
    {
      DefaultCubatureFactory cubFactory;         // create cubature factory
      m_cub = cubFactory.create<Kokkos::HostSpace, double, double >(cell_topo, cubDegree);                                    // create default cubature
      unsigned numCubPoints = m_cub->getNumPoints();                                      // retrieve number of cubature points
      m_Cub_Points_Tag = Cub_Points_Tag(numCubPoints);
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([P],[D])
    IntrepidManager::CubaturePoints::CubaturePoints(const IM& im) : BaseType("IM::CubaturePoints", NUM(Cub_Points_Tag), NUM(Spatial_Dim_Tag) ), m_im(im)
    {
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([P])
    IntrepidManager::CubatureWeights::CubatureWeights(const IM& im) : BaseType("IM::CubatureWeights", NUM(Cub_Points_Tag))
    {
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([C], [V], [D])
    IntrepidManager::CellWorkSet::
    CellWorkSet(const IM& im) : BaseType("IM::CellWorkSet", NUM(Elements_Tag), NUM(NodesPerElem_Tag), NUM(Spatial_Dim_Tag))
    {

    }


    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([C], [P], [D])
    IntrepidManager::PhysicalCoords::
    PhysicalCoords(const IM& im) : BaseType("IM::PhysicalCoords", NUM(Elements_Tag), NUM(Cub_Points_Tag), NUM(Spatial_Dim_Tag) ), m_im(im)
    {
    }

    void
    IntrepidManager::PhysicalCoords::
    operator()(CellWorkSet& c, CubaturePoints& xi)
    {
      //Intrepid2::CellTools<Kokkos::HostSpace>::mapToPhysicalFrame(images, preImages, triNodes, triangle_3);
      auto basis = BasisTable::getInstance()->getBasis(*m_im.m_topo);
      Intrepid2::CellTools<Kokkos::HostSpace>::mapToPhysicalFrame(*this, xi, c, basis);
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
      /// ([C], [P], [D], [D])
    IntrepidManager::Jacobian::
    Jacobian(const IM& im) : BaseType("IM::Jacobian", NUM(Elements_Tag), NUM(Cub_Points_Tag), NUM(Spatial_Dim_Tag), NUM(Spatial_Dim_Tag)), m_im(im)
    {

    }
    void
    IntrepidManager::Jacobian::

    operator()(CubaturePoints& xi, CellWorkSet& c, CellTopology& topo)
    {
      CellTools<Kokkos::HostSpace>::setJacobian(*this, xi, c, topo);           // compute cell Jacobians
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
      /// ([C], [P], [D])
    IntrepidManager::FaceNormal::
    FaceNormal(const IM& im) : BaseType("IM::FaceNormal", NUM(Elements_Tag), NUM(Cub_Points_Tag), NUM(Spatial_Dim_Tag))
    {

    }
    void
    IntrepidManager::FaceNormal::
    operator()(Jacobian& jac, int i_face, CellTopology& topo)
    {
      Intrepid2::CellTools<Kokkos::HostSpace>::getPhysicalFaceNormals(*this, jac, i_face, topo);
    }


    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([C], [P], [D], [D])
    IntrepidManager::JacobianInverse::
    JacobianInverse(const IM& im) : BaseType("IM::JacobianInverse", NUM(Elements_Tag), NUM(Cub_Points_Tag), NUM(Spatial_Dim_Tag), NUM(Spatial_Dim_Tag))
    {

    }
    void    IntrepidManager::JacobianInverse::
    operator()(Jacobian& jac)
    {
      CellTools<Kokkos::HostSpace>::setJacobianInv(*this, jac);
    }


    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
      /// ([C], [P])
    IntrepidManager::JacobianDet::
    JacobianDet(const IM& im) : BaseType("IM::JacobianDet", NUM(Elements_Tag), NUM(Cub_Points_Tag))
    {
    }

    void
    IntrepidManager::JacobianDet::
    operator()(Jacobian& jac)
    {
      Intrepid2::CellTools<Kokkos::HostSpace>::setJacobianDet(*this, jac);
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
      /// ([C], [P])
    IntrepidManager::WeightedMeasure::
    WeightedMeasure(const IM& im) : BaseType("IM::WeightedMeasure", NUM(Elements_Tag), NUM(Cub_Points_Tag))
    {
    }

    void
    IntrepidManager::WeightedMeasure::
    operator()(CubatureWeights& w, JacobianDet& dJ)
    {
      FunctionSpaceTools<Kokkos::HostSpace>::computeCellMeasure(*this, dJ, w);
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
      /// ([C], [P], [DOF])
    IntrepidManager::IntegrandValuesDOF::
    IntegrandValuesDOF(const IM& im) : BaseType("IM::IntegrandValuesDOF", NUM(Elements_Tag), NUM(Cub_Points_Tag), NUM(DOFs_Tag))
    {
    }
    void
    IntrepidManager::IntegrandValuesDOF::
    copyFrom(MDArray& mda)
    {
      Kokkos::deep_copy(*this, mda);
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
      /// ([C], [P])
    IntrepidManager::IntegrandValues::
    IntegrandValues(const IM& im) : BaseType("IM::IntegrandValues", NUM(Elements_Tag), NUM(Cub_Points_Tag))
    {
    }
    void
    IntrepidManager::IntegrandValues::
    copyFrom(MDArray& mda)
    {
     Kokkos::deep_copy(*this, mda);
    }
    void
    IntrepidManager::IntegrandValues::
    copyFrom(IntrepidManager& im, MDArray& mda, int iDof)
    {
      Kokkos::deep_copy(*this, Kokkos::subview(mda, std::make_pair(0, im.m_Elements_Tag.num),  std::make_pair(0, im.m_Cub_Points_Tag.num), iDof));
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([C], [DOF])
    IntrepidManager::IntegralDOF::
    IntegralDOF(const IM& im) : BaseType("IM::IntegralDOF", NUM(Elements_Tag), NUM(DOFs_Tag))
    {
    }

    /// wXdOmega: ([C], [P])
    /// iv:       ([C], [P], [DOF])
    /// this:     ([C], [DOF])
    void
    IntrepidManager::IntegralDOF::
    operator()(IntegrandValuesDOF& iv, WeightedMeasure& wXdOmega)
    {
      VERIFY_OP(iv.rank(),       == , 3,                     "IntrepidManager::Integral::operator() bad iv rank");
      VERIFY_OP(wXdOmega.rank(), == , 2,                     "IntrepidManager::Integral::operator() bad wXdOmega rank");
      VERIFY_OP((*this).rank(),  == , 2,                     "IntrepidManager::Integral::operator() bad (*this) rank");
      VERIFY_OP(iv.extent_int(0), == , this->extent_int(0),    "IntrepidManager::Integral::operator() bad");
      VERIFY_OP(iv.extent_int(1), == , wXdOmega.extent_int(1), "IntrepidManager::Integral::operator() bad");
      VERIFY_OP(iv.extent_int(0), == , wXdOmega.extent_int(0), "IntrepidManager::Integral::operator() bad");
      VERIFY_OP(iv.extent_int(2), == , this->extent_int(1),    "IntrepidManager::Integral::operator() bad");
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([C])
    IntrepidManager::Integral::
    Integral(const IM& im) : BaseType("IM::Integral", NUM(Elements_Tag))
    {
    }

    /// wXdOmega: ([C], [P])
    /// iv:       ([C], [P])
    /// this:     ([C])
    void
    IntrepidManager::Integral::
    operator()(IntegrandValues& iv, WeightedMeasure& wXdOmega)
    {
      VERIFY_OP(iv.rank(),       == , 2,                     "IntrepidManager::Integral::operator() bad iv rank");
      VERIFY_OP(wXdOmega.rank(), == , 2,                     "IntrepidManager::Integral::operator() bad wXdOmega rank");
      VERIFY_OP((*this).rank(),  == , 1,                     "IntrepidManager::Integral::operator() bad (*this) rank");
      VERIFY_OP(iv.extent_int(0), == , this->extent_int(0),    "IntrepidManager::Integral::operator() bad");
      VERIFY_OP(iv.extent_int(1), == , wXdOmega.extent_int(1), "IntrepidManager::Integral::operator() bad");
      VERIFY_OP(iv.extent_int(0), == , wXdOmega.extent_int(0), "IntrepidManager::Integral::operator() bad");

      FunctionSpaceTools<Kokkos::HostSpace>::integrate(*this, iv, wXdOmega);
    }

    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------
    /// ([C],[B],[P]), or ([C],[B],[P],[D]) for GRAD
    static ComputeBases s_compute_bases;

    IntrepidManager::Bases::
    Bases(const IM& im) : BaseType("IM::Bases", NUM(Elements_Tag), NUM(NodesPerElem_Tag), NUM(Cub_Points_Tag))
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
    FieldValues(const IM& im) : BaseType("IM::FieldValues", NUM(Elements_Tag), NUM(Cub_Points_Tag), NUM(DOFs_Tag)) {}

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


      //Step 6: Apply cell tools
      J(xi, c, *topo);
      Ji(J);
      dJ(J);

    }

    void IntrepidManager::isInElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates, unsigned& found_it, const stk::mesh::Entity element,
                                      const stk::mesh::BulkData& bulkData)
    {
      found_it = 0;

      // FIXME consider caching the coords_field in FieldFunction
      const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
      stk::mesh::FieldBase *coords_field = metaData.get_field(stk::topology::NODE_RANK, "coordinates");

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
      MDArray cellWorkset("cellWorkset", numCells, numNodes, cellDim);

      /// FIXME -- fill cellWorkset
      const MyPairIterRelation elem_nodes(bulkData, element, stk::topology::NODE_RANK );


      for (unsigned iCell = 0; iCell < numCells; iCell++)
        {
          for (unsigned iNode = 0; iNode < numNodes; iNode++)
            {
              stk::mesh::Entity node = elem_nodes[iNode].entity();
              double * node_coord_data = static_cast<double*>(stk::mesh::field_data( *coords_field , node));
              for (unsigned iDim=0; iDim < cellDim; iDim++)
                {
                  cellWorkset(iCell, iNode, iDim) = node_coord_data[iDim];
                }
            }
        }

      Intrepid2::CellTools<Kokkos::HostSpace>::mapToReferenceFrame(found_parametric_coordinates, input_phy_points, cellWorkset, topo);
      MDArrayUInt inclusion_results("inclusion_results", 1, 1);  // FIXME
      double threshold = 1.e-4; // (INTREPID_THRESHOLD default = 10*double_eps ~ 20e-16)
      Intrepid2::CellTools<Kokkos::HostSpace>::checkPointwiseInclusion(inclusion_results, found_parametric_coordinates, topo, threshold);
      found_it = inclusion_results(0,0);
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
      MDArray cellNodes("cellNodes", numCells, numNodes, spaceDim);

      DefaultCubatureFactory cubFactory;                                              // create cubature factory
      unsigned cubDegree = 2;                                                                      // set cubature degree, e.g. 2
      auto myCub = cubFactory.create<Kokkos::HostSpace, double, double>(cell_topo, cubDegree);         // create default cubature

      unsigned numCubPoints = myCub->getNumPoints();                                               // retrieve number of cubature points

      MDArray cub_points("cub_points", numCubPoints, spaceDim);
      MDArray cub_weights("cub_weights", numCubPoints);

      // Rank-4 array (C,P,D,D) for the Jacobian and its inverse and Rank-2 array (C,P) for its determinant
      MDArray jacobian("jacobian", numCells, numCubPoints, spaceDim, spaceDim);
      MDArray jacobian_inv("jacobian_inv", numCells, numCubPoints, spaceDim, spaceDim);
      MDArray jacobian_det("jacobian_det", numCells, numCubPoints);

      myCub->getCubature(cub_points, cub_weights);                                          // retrieve cubature points and weights

      // Methods to compute cell Jacobians, their inverses and their determinants

      CellTools<Kokkos::HostSpace>::setJacobian(jacobian, cub_points, cellNodes, cell_topo);           // compute cell Jacobians
      CellTools<Kokkos::HostSpace>::setJacobianInv(jacobian_inv, jacobian);                            // compute inverses of cell Jacobians
      CellTools<Kokkos::HostSpace>::setJacobianDet(jacobian_det, jacobian);                            // compute determinants of cell Jacobians

      MDArray weightedMeasure("weightedMeasure", numCells, numCubPoints);
      MDArray onesLeft("onesLeft", numCells,  numCubPoints);
      MDArray volume("volume", numCells);

      // compute weighted measure
      FunctionSpaceTools<Kokkos::HostSpace>::computeCellMeasure(weightedMeasure, jacobian_det, cub_weights);

      // integrate to get volume
      FunctionSpaceTools<Kokkos::HostSpace>::integrate(volume, onesLeft, weightedMeasure);

    }



  }
