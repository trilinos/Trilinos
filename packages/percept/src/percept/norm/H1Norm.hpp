// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_H1Norm_hpp
#define percept_H1Norm_hpp

#include <percept/norm/Norm.hpp>

  namespace percept
  {


    class H1_NormOp  : public Function, public HasFinalOp<std::vector<double> >
    {
      Function& m_integrand;
      Teuchos::RCP<Function > m_grad_integrand;
      MDArray m_grad_codomain;
      int m_spatialDim;

    public:

      /// integrand tells what fields Intrepid2 should compute, etc.
      H1_NormOp(Function& integrand, int spatialDim=3) : Function("H1_NormOp",integrand.getCodomainDimensions(), integrand.getCodomainDimensions()) , m_integrand(integrand), m_spatialDim(spatialDim) {
        m_grad_integrand = m_integrand.gradient(spatialDim);
        m_grad_codomain = m_grad_integrand->getNewCodomain();
      }

      void operator()(MDArray& pc_mda, MDArray& iv_mda, double time_value_optional=0.0)
      {
        //VERIFY_OP(pc_mda.size(), ==, iv_mda.size(), "H1_NormOp::operator() bad sizes");

        m_integrand(pc_mda, iv_mda, time_value_optional);
        std::vector<int> dims;
        iv_mda.dimensions(dims);
        dims.back() = m_grad_codomain.dimension(m_grad_codomain.rank()-1);
        VERIFY_OP_ON(dims.back(), ==, m_spatialDim, "H1Norm::operator()");
        MDArray out_grad(Teuchos::Array<int>(dims.begin(), dims.end()));
        (*m_grad_integrand)(pc_mda, out_grad, time_value_optional);
        VERIFY_OP_ON(out_grad.size(), ==, iv_mda.size()*m_spatialDim, "H1Norm::operator() 2");
        for (int i = 0; i < out_grad.size()/ m_spatialDim; i++)
          {
            double sum=0.0;
            for (int j = 0; j < m_spatialDim; j++)
              sum += SQR(out_grad[i * m_spatialDim + j]);

            sum += SQR(iv_mda(i));

            iv_mda(i) = sum;
          }
      }

      virtual void operator()(MDArray& pc_mda, MDArray& iv_mda, const stk::mesh::Entity element, const MDArray& parametric_coords, double time_value_optional=0.0)
      {
        //(*this)(domain, codomain, time_value_optional);

        m_integrand(pc_mda, iv_mda, element, parametric_coords, time_value_optional);
        std::vector<int> dims;
        iv_mda.dimensions(dims);
        dims.back() = m_grad_codomain.dimension(m_grad_codomain.rank()-1);
        VERIFY_OP_ON(dims.back(), ==, m_spatialDim, "H1Norm::operator()");
        MDArray out_grad(Teuchos::Array<int>(dims.begin(), dims.end()));
        (*m_grad_integrand)(pc_mda, out_grad, element, parametric_coords, time_value_optional);
        VERIFY_OP_ON(out_grad.size(), ==, iv_mda.size()*m_spatialDim, "H1Norm::operator() 2");
        for (int i = 0; i < out_grad.size()/ m_spatialDim; i++)
          {
            double sum=0.0;
            for (int j = 0; j < m_spatialDim; j++)
              sum += SQR(out_grad[i * m_spatialDim + j]);

            sum += SQR(iv_mda(i));

            iv_mda(i) = sum;
          }

      }

      virtual void operator()(MDArray& pc_mda, MDArray& iv_mda, const stk::mesh::Bucket& bucket, const MDArray& parametric_coords, double time_value_optional=0.0)
      {
        //(*this)(domain, codomain, time_value_optional);

        m_integrand(pc_mda, iv_mda, bucket, parametric_coords, time_value_optional);
        std::vector<int> dims;
        iv_mda.dimensions(dims);
        dims.back() = m_grad_codomain.dimension(m_grad_codomain.rank()-1);
        VERIFY_OP_ON(dims.back(), ==, m_spatialDim, "H1Norm::operator()");
        MDArray out_grad(Teuchos::Array<int>(dims.begin(), dims.end()));
        (*m_grad_integrand)(pc_mda, out_grad, bucket, parametric_coords, time_value_optional);
        VERIFY_OP_ON(out_grad.size(), ==, iv_mda.size()*m_spatialDim, "H1Norm::operator() 2");
        for (int i = 0; i < out_grad.size()/ m_spatialDim; i++)
          {
            double sum=0.0;
            for (int j = 0; j < m_spatialDim; j++)
              sum += SQR(out_grad[i * m_spatialDim + j]);

            sum += SQR(iv_mda(i));

            iv_mda(i) = sum;
          }


      }

      void finalOp(const std::vector<double>& vin, std::vector<double>& vout)
      {
        for (unsigned i = 0; i < vin.size(); i++)
          vout[i] = std::sqrt(vin[i]);
      }
    };

    /// compute the H1 norm or semi-norm
    class H1Norm : public  Norm<2>
    {
    public:
      typedef Norm<2> Base;

      H1Norm(stk::mesh::BulkData& bulkData, std::string partName, TurboOption turboOpt=TURBO_NONE, bool is_surface_norm=false) :
        Base(bulkData, partName, turboOpt, is_surface_norm) {}

      H1Norm(stk::mesh::BulkData& bulkData, MDArrayString& partNames, TurboOption turboOpt=TURBO_NONE, bool is_surface_norm=false) :
        Base(bulkData, partNames, turboOpt, is_surface_norm) {}

      H1Norm(stk::mesh::BulkData& bulkData, stk::mesh::Part *part = 0, TurboOption turboOpt=TURBO_NONE, bool is_surface_norm=false) :
        Base(bulkData, part, turboOpt, is_surface_norm) {}

#ifndef SWIG
      H1Norm(stk::mesh::BulkData& bulkData, stk::mesh::Selector * selector,TurboOption turboOpt=TURBO_NONE,  bool is_surface_norm = false) :
        Base(bulkData, selector, turboOpt, is_surface_norm) {}
#endif

      virtual void operator()(Function& integrand, Function& result)
      {
        EXCEPTWATCH;
        m_is_norm_being_evaluated=true;

        /** contract:
         *
         * 1. @param integrand : maps input_phy_points([C],[P],[D]) to output_values([C],[P],[DOF])
         * 2. this function then copies output_values to the @param result (assumes it is a ConstantFunction)
         *
         * algorithm:
         * 1. loop over elements (later buckets) and fire @param integrand on each quadrature point/cell collection
         * 2. further process the output_values by the LN_Op function, which simply raises the results to the N'th power
         * 3. (note: later we can pass in an Op instead of this being hard-coded to norm_LN)
         * 4. fire a "finalOp" on the resulting parallel accumulation buffer (supplied by LN_Op)
         */

        //std::cout << "type= " << typeid(integrand).name() << " " << typeid(FieldFunction).name() << std::endl;
        if (typeid(integrand) == typeid(FieldFunction))
        {
          VERIFY_OP_ON(m_turboOpt, !=, TURBO_NONE, "Can't use TURBO_NONE with Norm of FieldFunction");
        }

          {
            // FIXME - make all percept code const-correct
            const stk::mesh::MetaData& meta = m_bulkData.mesh_meta_data();
            PerceptMesh eMesh(&meta, &m_bulkData);
            int spatialDim = eMesh.get_spatial_dim();
            H1_NormOp H1_op(integrand, spatialDim);
            //CompositeFunction H1_of_integrand("H1_of_integrand", integrand, H1_op);
            IntegratedOp integrated_H1_op(H1_op, m_turboOpt);
            integrated_H1_op.setCubDegree(m_cubDegree);

            const stk::mesh::Part& locally_owned_part = meta.locally_owned_part();
            stk::mesh::Selector selector(*m_selector & locally_owned_part);
            //eMesh.print_info("Norm");
            if (m_turboOpt == TURBO_NONE || m_turboOpt == TURBO_ELEMENT)
              {
                elementOpLoop(m_bulkData, integrated_H1_op, 0, &selector, m_is_surface_norm);
              }
            else if (m_turboOpt == TURBO_BUCKET)
              {
                bucketOpLoop(m_bulkData, integrated_H1_op, 0, &selector, m_is_surface_norm);
              }

            unsigned vec_sz = integrated_H1_op.getValue().size();
            std::vector<double> local = integrated_H1_op.getValue();

            //unsigned p_rank = m_bulkData.parallel_rank();
            //std::cout  << "P["<<p_rank<<"] value = " << local << std::endl;

            std::vector<double> global_sum(vec_sz, 0.0);
            stk::all_reduce_sum(m_bulkData.parallel(), &local[0], &global_sum[0], vec_sz);

            std::vector<double> result1(vec_sz);
            H1_op.finalOp(global_sum, result1);

            if (typeid(result) == typeid(ConstantFunction))
              {
                ConstantFunction& cf = *dynamic_cast<ConstantFunction *>(&result);
                cf.setValue(result1[0]);  // FIXME for multiple values
              }
            else if (typeid(result) == typeid(ConstantFunctionVec))
              {
                ConstantFunctionVec& cf = *dynamic_cast<ConstantFunctionVec *>(&result);
                cf.setValue(result1);
              }
          }

        m_is_norm_being_evaluated=false;
      }

      //void operator()(FieldFunction& input);
    };


  }

#endif
