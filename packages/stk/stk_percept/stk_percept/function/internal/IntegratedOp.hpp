#ifndef stk_percept_IntegratedOp_hpp
#define stk_percept_IntegratedOp_hpp

#include <cmath>
#include <math.h>
#include <map>
#include <vector>

#include <typeinfo>

#include <stk_percept/Util.hpp>

#include <stk_percept/function/MDArray.hpp>

#include <stk_percept/function/FunctionOperator.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/Function.hpp>
#include <stk_percept/function/CompositeFunction.hpp>
#include <stk_percept/function/internal/HasValue.hpp>
#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/function/ElementOp.hpp>
#include <stk_percept/function/BucketOp.hpp>
#include <stk_percept/function/MultipleFieldFunction.hpp>

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/ExceptionWatch.hpp>

#include <stk_percept/norm/IntrepidManager.hpp>

namespace stk
{
  namespace percept
  {
    class IntegratedOp : public ElementOp, public BucketOp, public HasConstValue<std::vector<double> >
    {
    public:

      IntegratedOp(Function& integrand,  TurboOption turboOpt=TURBO_NONE, mesh::FieldBase *field=0) :
        m_nDOFs(1), m_accumulation_buffer(), m_count_elems(0), m_is_field(false), m_integrand(integrand), m_turboOpt(turboOpt)
      {
        if (typeid(integrand) == typeid(FieldFunction))
          {
            m_is_field = true;
          }
        if (field)
          {
            const stk::mesh::FieldBase::Restriction & r = field->restriction(stk::mesh::fem::FEMMetaData::NODE_RANK, mesh::fem::FEMMetaData::get(*field).universal_part());
            unsigned stride = r.dimension() ;
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

        init();
        m_accumulation_buffer.resize(m_nDOFs);
      }

      void init()
      {
        m_count_elems=0;
        m_accumulation_buffer.assign(m_nDOFs, 0.0);
      }

      std::vector<double>& getValue(void) { return m_accumulation_buffer; }
      unsigned getElementCount() { return m_count_elems; }

      /// innermost operation of an bucket-based loop; return value of true forces the enclosing loop to terminate and this class'
      ///   derived classes can return info back to the loop invoker
      virtual bool operator()(const stk::mesh::Bucket& bucket,  stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData);

      /// innermost operation of an element-based loop; return value of true forces the enclosing loop to terminate and this class'
      ///   derived classes can return info back to the loop invoker
      virtual bool operator()(const stk::mesh::Entity& element, stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData);
      void init_elementOp() { init(); }
      void fini_elementOp() {}

    private:


      template<class BucketOrEntity>
      bool helper(const BucketOrEntity& bucket_or_element, stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData)
      {
        EXCEPTWATCH;

        const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(bucket_or_element);
        CellTopology cell_topo(cell_topo_data);

        VectorFieldType& coord_field = *(mesh::fem::FEMMetaData::get(bulkData)).get_field<VectorFieldType>("coordinates");

        // FIXME for fields not on a Node
        unsigned nDOF = m_nDOFs;

        unsigned nCells = PerceptMesh::size1(bucket_or_element);
        m_count_elems += nCells;

        typedef IntrepidManager IM;
        unsigned cubDegree = 2;
        IM im(Elements_Tag(nCells), cell_topo, cubDegree);
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

        PerceptMesh::fillCellNodes(bucket_or_element,  &coord_field, cn, spaceDim);

        // get jacobian
        J(xi, cn, cell_topo);
        dJ(J);
        wXdJ(wt, dJ);

        if (0)
          {
	    using namespace shards;

            std::cout << "dJ= \n" << dJ << std::endl;
            std::cout << "wXdJ= \n" << wXdJ << std::endl;
            std::cout << "xi= \n" << xi << std::endl;
            std::cout << "wt= \n" << wt << std::endl;
            std::cout << "cn= \n" << cn << std::endl;
            Util::setDoPause(true);
            Util::pause(true);
          }

        // get physical coordinates at integration points
        pc(cn, xi);

        // get bases
#if 1
        // FIXME
        MDArray xi_mda;
        xi.copyTo(xi_mda);
        Nb(bucket_or_element, xi_mda);
#else
        Nb(bucket_or_element, xi);
#endif

        // apply integrand (right now we have MDArray hard-coded... FIXME - templatize on its type)
        // it should look like this (one instead of multiple lines):
#if 0
        m_integrand(pc, v);
#else
        MDArray pc_mda;
        pc.copyTo(pc_mda);
        std::vector<int>  ivDims;
        ivD.dimensions( ivDims);


        /// NOTE: m_integrand requires the ranks of in/out MDArrays to be such that out_rank >= in_rank
        /// Thus, we use IntegrandValuesDOF with [DOF] = 1, and then copy the result to IntegrandValues
        /// which does not have the additional rightmost DOF index (Intrepid doesn't have the concept of
        /// DOF's, it works on scalars only for the integration routines, or at least that's how I understand
        /// it currently.

        // create an array that stk::percept::Function will like to hold the results


        ivDims[ivDims.size()-1] = m_nDOFs;

        MDArray iv_mda ( Teuchos::Array<int>(ivDims.begin(), ivDims.end()));

        if (m_turboOpt == TURBO_ELEMENT || m_turboOpt == TURBO_BUCKET)
          {
            m_integrand(pc_mda, iv_mda, bucket_or_element, xi_mda);
          }
        else
          {
            m_integrand(pc_mda, iv_mda);
          }

        // now, copy from the results to an array that Intrepid::integrate will like

#endif

        for (unsigned iDof = 0; iDof < nDOF; iDof++)
          {
            iv.copyFrom(im, iv_mda, iDof);

            // get the integral
            Io(iv, wXdJ, COMP_BLAS);

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
                m_accumulation_buffer[iDof] += Io(iCell);
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
    };

    //template<>


    bool IntegratedOp::operator()(const stk::mesh::Bucket& bucket,  stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData)
    {
      EXCEPTWATCH;
      helper(bucket, field, bulkData);
      return false;
    }

    /// innermost operation of an element-based loop; return value of true forces the enclosing loop to terminate and this class'
    ///   derived classes can return info back to the loop invoker
    bool IntegratedOp::operator()(const stk::mesh::Entity& element, stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData)
    {
      EXCEPTWATCH;
      helper(element, field, bulkData);
      return false;
    }


  }
}
#endif
