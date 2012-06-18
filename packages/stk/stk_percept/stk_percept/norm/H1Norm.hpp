/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef stk_percept_H1Norm_hpp
#define stk_percept_H1Norm_hpp

#include <stk_percept/norm/Norm.hpp>

namespace stk
{
  namespace percept
  {

    /// compute the H1 norm or semi-norm
    class H1Norm : public typename Norm<2>
    {
    public:

      H1Norm(mesh::BulkData& bulkData, std::string partName, TurboOption turboOpt=TURBO_NONE, bool is_surface_norm=false) :
        Norm(bulkData, partName, turboOpt, is_surface_norm) {}

      H1Norm(mesh::BulkData& bulkData, MDArrayString& partNames, TurboOption turboOpt=TURBO_NONE, bool is_surface_norm=false) :
        Norm(bulkData, partNames, turboOpt, is_surface_norm) {}

      H1Norm(mesh::BulkData& bulkData, mesh::Part *part = 0, TurboOption turboOpt=TURBO_NONE, bool is_surface_norm=false) :
        Norm(bulkData, part, turboOpt, is_surface_norm) {}

#ifndef SWIG
      H1Norm(mesh::BulkData& bulkData, mesh::Selector * selector,TurboOption turboOpt=TURBO_NONE,  bool is_surface_norm = false) :
        Norm(bulkData, selector, turboOpt, is_surface_norm) {}
#endif

      virtual void operator()(Function& integrand, Function& result)
      {
        EXCEPTWATCH;

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

          {
            // FIXME - make all stk::percept code const-correct
            PerceptMesh eMesh(&mesh::fem::FEMMetaData::get(m_bulkData), &m_bulkData);
            LN_NormOp<2> LN_op(integrand);
            CompositeFunction LN_of_integrand("LN_of_integrand", integrand, LN_op);
            IntegratedOp integrated_LN_op(LN_of_integrand, m_turboOpt);
            integrated_LN_op.setCubDegree(m_cubDegree);

            const stk::mesh::Part& locally_owned_part = mesh::fem::FEMMetaData::get(m_bulkData).locally_owned_part();
            stk::mesh::Selector selector(*m_selector & locally_owned_part);
            //eMesh.print_info("Norm");
            if (m_turboOpt == TURBO_NONE || m_turboOpt == TURBO_ELEMENT)
              {
                eMesh.elementOpLoop(integrated_LN_op, 0, &selector, m_is_surface_norm);
              }
            else if (m_turboOpt == TURBO_BUCKET)
              {
                eMesh.bucketOpLoop(integrated_LN_op, 0, &selector, m_is_surface_norm);
              }

            unsigned vec_sz = integrated_LN_op.getValue().size();
            std::vector<double> local = integrated_LN_op.getValue();


            //unsigned p_rank = m_bulkData.parallel_rank();
            //std::cout  << "P["<<p_rank<<"] value = " << local << std::endl;

            std::vector<double> global_sum(vec_sz, 0.0);
            stk::all_reduce_sum(m_bulkData.parallel(), &local[0], &global_sum[0], vec_sz);

            std::vector<double> result1(vec_sz);
            LN_op.finalOp(global_sum, result1);

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

      }

      //void operator()(FieldFunction& input);
    };


  }
}


#endif
