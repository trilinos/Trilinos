// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_IntegratedOp_hpp
#define percept_IntegratedOp_hpp

#include <cmath>
#include <math.h>
#include <map>
#include <vector>

#include <typeinfo>

#include <percept/Util.hpp>

#include <percept/function/MDArray.hpp>

#include <percept/function/FunctionOperator.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/Function.hpp>
#include <percept/function/CompositeFunction.hpp>
#include <percept/function/internal/HasValue.hpp>
#include <percept/function/StringFunction.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/function/ElementOp.hpp>
#include <percept/function/BucketOp.hpp>
#include <percept/function/MultipleFieldFunction.hpp>

#include <percept/PerceptMesh.hpp>
#include <percept/ExceptionWatch.hpp>

#include <percept/norm/IntrepidManager.hpp>

  namespace percept
  {
    class IntegratedOp : public ElementOp, public BucketOp, public HasConstValue<std::vector<double> >
    {
    public:

      enum AccumulationType {
        ACCUMULATE_SUM,
        ACCUMULATE_MAX
      };

      IntegratedOp(Function& integrand,  TurboOption turboOpt=TURBO_NONE, stk::mesh::FieldBase *field=0) :
        m_nDOFs(1), m_accumulation_buffer(), m_count_elems(0), m_is_field(false), m_integrand(integrand), m_turboOpt(turboOpt),
        m_cubDegree(2), m_accumulation_type(ACCUMULATE_SUM)
      {
        if (typeid(integrand) == typeid(FieldFunction))
          {
            m_is_field = true;
          }
        if (field)
          {
            const stk::mesh::FieldBase::Restriction & r = stk::mesh::find_restriction(*field, stk::topology::NODE_RANK, stk::mesh::MetaData::get(*field).universal_part());
            unsigned stride = r.num_scalars_per_entity() ;
            m_nDOFs = stride;
          }
        else
          {
            unsigned sz = integrand.getCodomainDimensions().size();
            if (sz >= 1)
              {
                m_nDOFs = integrand.getCodomainDimensions()[ sz -1 ];
              }
          }

        m_accumulation_buffer.resize(m_nDOFs);
        init();
      }

      void setAccumulationType(AccumulationType type) { m_accumulation_type = type; }
      AccumulationType getAccumulationType() { return m_accumulation_type; }

      void setCubDegree(unsigned cubDegree) { m_cubDegree= cubDegree; }
      unsigned getCubDegree() { return m_cubDegree; }

      void init()
      {
        m_count_elems=0;
        m_accumulation_buffer.assign(m_nDOFs, 0.0);
      }

      std::vector<double>& getValue(void) { return m_accumulation_buffer; }
      unsigned getElementCount() { return m_count_elems; }

      /// innermost operation of an bucket-based loop; return value of true forces the enclosing loop to terminate and this class'
      ///   derived classes can return info back to the loop invoker
      virtual bool operator()(const stk::mesh::Bucket& bucket,  stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData);

      /// innermost operation of an element-based loop; return value of true forces the enclosing loop to terminate and this class'
      ///   derived classes can return info back to the loop invoker
      virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData);
      void init_elementOp() { init(); }
      void fini_elementOp() {}

    private:

      const stk::mesh::Bucket& mybucket(const stk::mesh::BulkData& bulkData, const stk::mesh::Bucket& bucket_or_element) const { return bucket_or_element; }
      const stk::mesh::Bucket& mybucket(const stk::mesh::BulkData& bulkData, const stk::mesh::Entity& bucket_or_element) const { return bulkData.bucket(bucket_or_element); }
      stk::mesh::Bucket& mybucket( stk::mesh::BulkData& bulkData,  stk::mesh::Bucket& bucket_or_element)  { return bucket_or_element; }
      stk::mesh::Bucket& mybucket( stk::mesh::BulkData& bulkData,  stk::mesh::Entity& bucket_or_element)  { return bulkData.bucket(bucket_or_element); }


      template<class BucketOrEntity>
      bool helper(const BucketOrEntity& bucket_or_element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
      {
        EXCEPTWATCH;

        const CellTopologyData * const cell_topo_data = stk::mesh::get_cell_topology(mybucket(bulkData, bucket_or_element).topology()).getCellTopologyData();
        CellTopology cell_topo(cell_topo_data);

        if (cell_topo_data->key == shards::Particle::key || 
            cell_topo_data->key == shards::Beam<2>::key) 
          return false;
        
        int cell_dimension = cell_topo.getDimension();
        int meta_dimension = bulkData.mesh_meta_data().spatial_dimension();

        if (cell_dimension == meta_dimension - 1)
          {
            return helperSubDim(bucket_or_element, field, bulkData);
          }

        VERIFY_OP_ON(cell_dimension, ==, meta_dimension, "Dimensions don't match");

        const stk::mesh::MetaData& meta = bulkData.mesh_meta_data();
        CoordinatesFieldType& coord_field = *meta.get_field<double>(stk::topology::NODE_RANK, "coordinates");

        // FIXME for fields not on a Node
        unsigned nDOF = m_nDOFs;

        unsigned nCells = PerceptMesh::size1(bucket_or_element);
        m_count_elems += nCells;

        typedef IntrepidManager IM;
        unsigned cubDegree = m_cubDegree;
        IM im(Elements_Tag(nCells), cell_topo, cubDegree);
        if (0)
          {
            unsigned numCubPoints = im.m_cub->getNumPoints();
            std::cout << "numCubPoints= " << numCubPoints << std::endl;
          }


        // FIXME
        im.m_DOFs_Tag.num = m_nDOFs;
        // FIXME

        IM::Jacobian              J  (im);
        IM::JacobianDet          dJ  (im);
        IM::CubaturePoints       xi  (im);
        IM::CellWorkSet          cn  (im);
        IM::CubatureWeights      wt  (im);
        IM::PhysicalCoords       pc  (im);
        IM::IntegrandValues      iv  (im);
        IM::IntegrandValuesDOF  ivD  (im);
        IM::Integral             Io  (im);
        IM::Bases                Nb  (im);

        IM::WeightedMeasure wXdJ  (im);
        IM::FieldValues       fv  (im);

        im.m_cub->getCubature(xi, wt);

        unsigned spaceDim = im.m_Spatial_Dim_Tag.num;

        PerceptMesh::fillCellNodes(bulkData, bucket_or_element,  &coord_field, cn, spaceDim);

        // get jacobian
        J(xi, cn, cell_topo);
        dJ(J);
        wXdJ(wt, dJ);


        // get physical coordinates at integration points
        pc(cn, xi);

        // get bases
#if 1
        // FIXME
        MDArray xi_mda;
        xi.copyTo(xi_mda);
        Nb(bulkData, bucket_or_element, xi_mda);
#else
        Nb(bulkData, bucket_or_element, xi);
#endif

        // apply integrand (right now we have MDArray hard-coded... FIXME - templatize on its type)
        // it should look like this (one instead of multiple lines):
#if 0
        m_integrand(pc, v);
#else
        MDArray pc_mda;
        pc.copyTo(pc_mda);
        std::vector<int>  ivDims(ivD.rank());
        for (size_t i=0; i<ivDims.size(); ++i)
          ivDims[i] = ivD.extent_int(i);


        /// NOTE: m_integrand requires the ranks of in/out MDArrays to be such that out_rank >= in_rank
        /// Thus, we use IntegrandValuesDOF with [DOF] = 1, and then copy the result to IntegrandValues
        /// which does not have the additional rightmost DOF index (Intrepid2 doesn't have the concept of
        /// DOF's, it works on scalars only for the integration routines, or at least that's how I understand
        /// it currently.

        // create an array that percept::Function will like to hold the results
        ivDims[ivDims.size()-1] = m_nDOFs;

        std::vector<size_t> dimensions(8, KOKKOS_INVALID_INDEX);
        for(size_t i=0; i< ivDims.size(); ++i) 
          dimensions[i]=ivDims[i];

        MDArray iv_mda ( "iv_mda", dimensions[0], dimensions[1], dimensions[2], dimensions[3], dimensions[4], dimensions[5], dimensions[6], dimensions[7]);

        if (m_turboOpt == TURBO_ELEMENT || m_turboOpt == TURBO_BUCKET)
          {
            m_integrand(pc_mda, iv_mda, bucket_or_element, xi_mda);
          }
        else
          {
            m_integrand(pc_mda, iv_mda);
          }

        // now, copy from the results to an array that Intrepid2::integrate will like

#endif

        for (unsigned iDof = 0; iDof < nDOF; iDof++)
          {
            iv.copyFrom(im, iv_mda, iDof);

            // get the integral
            if (m_accumulation_type == ACCUMULATE_SUM)
              {
                Io(iv, wXdJ);
              }

            //optional design:
            //
            //  Io(integrand(pc_mda, v), wXdJ(w, dJ(J(xi, c, cell_topo)), COMP_BLAS);

            for (unsigned iCell = 0; iCell < nCells; iCell++)
              {
                //                 if (Util::getFlag(0))
                //                   {
                //                     std::cout << "tmp Io(iCell)= " << Io(iCell) << std::endl;
                //                     Util::pause(true, "Io(iCell)");
                //                   }
                if (m_accumulation_type == ACCUMULATE_SUM)
                  {
                    m_accumulation_buffer[iDof] += Io(iCell);
                  }
                else if (m_accumulation_type == ACCUMULATE_MAX)
                  {
                    double valIo = 0.0;
                    for (int ivpts = 0; ivpts < iv.extent_int(1); ivpts++)
                      {
                        valIo = std::max(valIo, iv((int)iCell, ivpts));
                      }
                    //std::cout << "m_accumulation_buffer[iDof] = " << m_accumulation_buffer[iDof] << " valIO= " << valIo  << std::endl;
                    m_accumulation_buffer[iDof] = std::max(m_accumulation_buffer[iDof], valIo);
                  }
              }
          }

	
	// store field if one is given
	if (field) {
	  double * field_vals = static_cast<double *>(field_data(*field, bucket_or_element));

	  for (unsigned iDof = 0; iDof < nDOF; iDof++) {
            for (unsigned iCell = 0; iCell < nCells; iCell++) {
	      field_vals[iCell*nDOF+iDof] = Io(iCell);
	    }
	  }
	}

        return false;
      }

      /** for sub-dim elements embedded in higher dim (e.g. quad face of hex elements), we need to:
       *   1. get parent element through side relations
       *   2. get integration points and weights on quad ref element
       *   3. get integration points and weights on hex ref element
       *   4. get face normals (from Intrepid2 they are non-normalized normals, i.e. they are surface area vectors)
       *   5. get norm of face normals
       *   6. get products with cubature weights
       */

      bool helperSubDim(const stk::mesh::Bucket& bucket, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
      {
        const unsigned num_elements_in_bucket = bucket.size();

        for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
          {
            stk::mesh::Entity element = bucket[iElement];
            helperSubDim(element, field, bulkData);
          }
        return false;
      }

      bool helperSubDim(const stk::mesh::Entity child_element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
      {
        EXCEPTWATCH;

        const CellTopologyData * const child_cell_topo_data = stk::mesh::get_cell_topology(mybucket(bulkData, child_element).topology()).getCellTopologyData();
        CellTopology child_cell_topo(child_cell_topo_data);
        int child_cell_dimension = child_cell_topo.getDimension();
        int meta_dimension = bulkData.mesh_meta_data().spatial_dimension();

        // for now, only allow face (or edge)
        VERIFY_OP_ON(child_cell_dimension, ==, meta_dimension - 1, "Dimensions don't match");

        const stk::mesh::MetaData& meta = bulkData.mesh_meta_data();
        CoordinatesFieldType& coord_field = *meta.get_field<double>(stk::topology::NODE_RANK, "coordinates");

        // FIXME for fields not on a Node
        unsigned nDOF = m_nDOFs;

        unsigned nCells = PerceptMesh::size1(child_element);
        m_count_elems += nCells;

        typedef IntrepidManager IM;
        unsigned cubDegree = m_cubDegree;
        const MyPairIterRelation parent_elements(bulkData, child_element, static_cast<stk::mesh::EntityRank>(bulkData.entity_rank(child_element) + 1) );
        VERIFY_OP_ON(parent_elements.size(), ==, 1, "cant find parent");
        const stk::mesh::Entity element = parent_elements[0].entity();
        unsigned i_face = parent_elements[0].relation_ordinal();

        const CellTopologyData * const cell_topo_data = stk::mesh::get_cell_topology(mybucket(bulkData, element).topology()).getCellTopologyData();
        CellTopology cell_topo(cell_topo_data);
        int cell_dimension = cell_topo.getDimension();
        VERIFY_OP_ON(cell_dimension, ==, meta_dimension , "Dimensions don't match");

        IM im(Elements_Tag(nCells), cell_topo, cubDegree);
        IM imChild(Elements_Tag(nCells), child_cell_topo, cubDegree);
        unsigned numCubPoints_child = imChild.m_cub->getNumPoints();
        im.m_Cub_Points_Tag = Cub_Points_Tag(numCubPoints_child);

        if (0)
          {
            std::cout << "numCubPoints_child= " << numCubPoints_child
                      << " parent rank= " << bulkData.entity_rank(element)
                      << " parent topo= " << cell_topo.getName()
                      << std::endl;
          }

        // FIXME
        im.m_DOFs_Tag.num = m_nDOFs;
        // FIXME

        // _c suffix is for the child (face) element
        IM::Jacobian              J  (im);
        IM::FaceNormal           fn  (im);
        //IM::JacobianDet          dJ  (im);
        IM::CubaturePoints       xi  (im);
        IM::CubaturePoints       xi_c  (imChild);
        IM::CellWorkSet          cn  (im);
        IM::CubatureWeights      wt  (im);
        IM::CubatureWeights      wt_c  (imChild);
        IM::PhysicalCoords       pc  (im);
        IM::IntegrandValues      iv  (im);
        IM::IntegrandValuesDOF  ivD  (im);
        IM::Integral             Io  (im);
        IM::Bases                Nb  (im);

        IM::WeightedMeasure wXfn  (im);
        IM::FieldValues       fv  (im);

        imChild.m_cub->getCubature(xi_c, wt_c);

        unsigned spaceDim = im.m_Spatial_Dim_Tag.num;

        PerceptMesh::fillCellNodes(bulkData, element,  static_cast<stk::mesh::FieldBase *>(&coord_field), cn, spaceDim);

        // get parent cell integration points
        // Map Gauss points on quad to reference face: paramGaussPoints -> refGaussPoints
        CellTools<Kokkos::HostSpace>::mapToReferenceSubcell(xi,
                                                 xi_c,
                                                 2, i_face, cell_topo);  // FIXME magic

        // get jacobian
        J(xi, cn, cell_topo);
        //dJ(J);

        //shards::ArrayVector<double, shards::NaturalOrder, Elements_Tag, Cub_Points_Tag > fn_Norm;

        // FIXME
        //fn(J, i_face, cell_topo);
        MDArray J_mda;
        J.copyTo(J_mda);
        MDArray fn_mda("fn_mda", im.m_Elements_Tag.num, numCubPoints_child, spaceDim);
        CellTools<Kokkos::HostSpace>::getPhysicalFaceNormals(fn_mda, J_mda, i_face, cell_topo);

        /// get norm of fn
        for (int icell = 0; icell < im.m_Elements_Tag.num; icell++)
          {
            for (int ipt = 0; ipt < (int)numCubPoints_child; ipt++)
              {
                double sum = 0.0;
                for (int i = 0; i < (int)spaceDim; i++)
                  {
                    sum += square(fn_mda(icell, ipt, i));
                  }
                wXfn(icell, ipt) = std::sqrt(sum) * wt_c(ipt);
              }
          }

        // get physical coordinates at integration points
        pc(cn, xi);

        // get bases
#if 1
        // FIXME
        MDArray xi_mda;
        xi.copyTo(xi_mda);
        Nb(bulkData, element, xi_mda);
#else
        Nb(bulkData, element, xi);
#endif

        // apply integrand (right now we have MDArray hard-coded... FIXME - templatize on its type)
        // it should look like this (one instead of multiple lines):
#if 0
        m_integrand(pc, v);
#else
        MDArray pc_mda;
        pc.copyTo(pc_mda);
        std::vector<int>  ivDims(ivD.rank());
        for (size_t i=0; i<ivDims.size(); ++i)
          ivDims[i] = ivD.extent_int(i);


        /// NOTE: m_integrand requires the ranks of in/out MDArrays to be such that out_rank >= in_rank
        /// Thus, we use IntegrandValuesDOF with [DOF] = 1, and then copy the result to IntegrandValues
        /// which does not have the additional rightmost DOF index (Intrepid doesn't have the concept of
        /// DOF's, it works on scalars only for the integration routines, or at least that's how I understand
        /// it currently.

        // create an array that percept::Function will like to hold the results

        ivDims[ivDims.size()-1] = m_nDOFs;

        std::vector<size_t> dimensions(8, KOKKOS_INVALID_INDEX);
        for(size_t i=0; i< ivDims.size(); ++i) 
          dimensions[i]=ivDims[i];

        MDArray iv_mda ( "iv_mda", dimensions[0], dimensions[1], dimensions[2], dimensions[3], dimensions[4], dimensions[5], dimensions[6], dimensions[7]);

        if (m_turboOpt == TURBO_ELEMENT || m_turboOpt == TURBO_BUCKET)
          {
            m_integrand(pc_mda, iv_mda, element, xi_mda);
          }
        else
          {
            m_integrand(pc_mda, iv_mda);
          }

        // now, copy from the results to an array that Intrepid2::integrate will like

#endif

        for (unsigned iDof = 0; iDof < nDOF; iDof++)
          {
            iv.copyFrom(im, iv_mda, iDof);

            // get the integral
            if (m_accumulation_type == ACCUMULATE_SUM)
              {
                Io(iv, wXfn);
              }

            //optional design:
            //
            //  Io(integrand(pc_mda, v), wXdJ(w, dJ(J(xi, c, cell_topo)), COMP_BLAS);

            for (unsigned iCell = 0; iCell < nCells; iCell++)
              {
                //                 if (Util::getFlag(0))
                //                   {
                //                     std::cout << "tmp Io(iCell)= " << Io(iCell) << std::endl;
                //                     Util::pause(true, "Io(iCell)");
                //                   }
                if (m_accumulation_type == ACCUMULATE_SUM)
                  {
                    m_accumulation_buffer[iDof] += Io(iCell);
                  }
                else if (m_accumulation_type == ACCUMULATE_MAX)
                  {
                    double valIo = 0.0;
                    for (int ivpts = 0; ivpts < iv.extent_int(1); ivpts++)
                      {
                        valIo = std::max(valIo, iv((int)iCell, ivpts));
                      }
                    //std::cout << "m_accumulation_buffer[iDof] = " << m_accumulation_buffer[iDof] << " valIO= " << valIo  << std::endl;
                    m_accumulation_buffer[iDof] = std::max(m_accumulation_buffer[iDof], valIo);
                  }
              }
          }
        return false;
      }




    private:
      unsigned m_nDOFs;
      std::vector<double> m_accumulation_buffer;
      unsigned m_count_elems;
      bool m_is_field;
      Function& m_integrand;
      TurboOption m_turboOpt;
      unsigned m_cubDegree;
      AccumulationType m_accumulation_type;
    };

    //template<>


    inline
    bool IntegratedOp::operator()(const stk::mesh::Bucket& bucket,  stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
    {
      EXCEPTWATCH;
      helper(bucket, field, bulkData);
      return false;
    }

    /// innermost operation of an element-based loop; return value of true forces the enclosing loop to terminate and this class'
    ///   derived classes can return info back to the loop invoker
    inline
    bool IntegratedOp::operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
    {
      EXCEPTWATCH;
      helper(element, field, bulkData);
      return false;
    }


  }

#endif
