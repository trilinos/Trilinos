// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_Norm_hpp
#define percept_Norm_hpp

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

#include <math.h>

#include <percept/PerceptMesh.hpp>
#include <percept/function/MDArray.hpp>

#include <percept/function/FunctionOperator.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/Function.hpp>
#include <percept/function/CompositeFunction.hpp>
#include <percept/function/StringFunction.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/function/MultipleFieldFunction.hpp>

#include <percept/function/internal/HasValue.hpp>
#include <percept/function/internal/IntegratedOp.hpp>

#include <percept/ExceptionWatch.hpp>

#include <percept/norm/IntrepidManager.hpp>
#include <percept/util/Loops.hpp>

  namespace percept
  {


    template<class ValueType>
    class HasFinalOp
    {
    public:
      virtual void finalOp(const ValueType& vin, ValueType& vout)=0;
    };

    /// [DEPRECATED]
    inline double square(double x) { return x*x; }

    /// for Power = -1, compute the inf-norm
    template<int Power=2>
    class LN_NormOp  : public Function, public HasFinalOp<std::vector<double> >
    {
    public:

      /// integrand tells what fields Intrepid2 should compute, etc.
      /// Note: this function is intended to be used to wrap a Function using CompositeFunction and thus its domain and codomain
      /// are the same as the wrapped function's codomain
      LN_NormOp(Function& integrand) : Function("LN_NormOp",integrand.getCodomainDimensions(), integrand.getCodomainDimensions()) {}

      void operator()(MDArray& integrand_values, MDArray& output_values, double time_value_optional=0.0)
      {
        VERIFY_OP(integrand_values.size(), ==, output_values.size(), "LN_NormOp::operator() bad sizes");
        for (size_t i = 0; i < integrand_values.size(); i++)
          {
            output_values[i] = std::pow(std::fabs(integrand_values[i]), double(std::fabs(Power)) );
          }
      }
      virtual void operator()(MDArray& domain, MDArray& codomain, const stk::mesh::Entity element, const MDArray& parametric_coords, double time_value_optional=0.0)
      {
        (*this)(domain, codomain, time_value_optional);
      }
      virtual void operator()(MDArray& domain, MDArray& codomain, const stk::mesh::Bucket& element, const MDArray& parametric_coords, double time_value_optional=0.0)
      {
        (*this)(domain, codomain, time_value_optional);
      }

      void finalOp(const std::vector<double>& vin, std::vector<double>& vout)
      {
        for (unsigned i = 0; i < vin.size(); i++)
          vout[i] = std::pow(vin[i], 1./(double(std::fabs(Power))) );
      }
    };

    class MaxOfNodeValues : public Function
    {
    public:

      MaxOfNodeValues(int spatialDim, Function& integrand) : Function("MaxOfNodeValues", Dimensions(spatialDim), integrand.getCodomainDimensions()), m_integrand(integrand), maxVal(0) {}
      Function& m_integrand;
      std::vector<double> maxVal;

      void operator()(MDArray& coords, MDArray& output_values, double time_value_optional)
      {
        //VERIFY_OP(coords.size(), ==, output_values.size(), "MaxOfNodeValues::operator() bad sizes");
        if (maxVal.size()==0) maxVal.resize(m_integrand.getCodomainDimensions()[0]);
        MDArray out("out",maxVal.size());
        m_integrand(coords, out);

        for (unsigned i = 0; i < maxVal.size(); i++)
          {
            maxVal[i] = std::max(maxVal[i], out(i));
          }
      }

      void operator()(MDArray& in, MDArray& out, const stk::mesh::Entity element, const MDArray& parametric_coords, double time_value_optional=0.0)
      {
        EXCEPTWATCH;
        throw std::runtime_error("Not implemented");
      }

      void operator()(MDArray& in, MDArray& out, const stk::mesh::Bucket& bucket, const MDArray& parametric_coords, double time_value_optional=0.0)
      {
        EXCEPTWATCH;
        throw std::runtime_error("Not implemented");
      }
    };

    struct NormBase
    {
      static bool m_is_norm_being_evaluated;
    };

    /// for Power = -1, compute the inf-norm
    template<int Power=2>
    class Norm : public FunctionOperator, public NormBase
    {
    protected:
      bool m_is_surface_norm;
      TurboOption m_turboOpt;
      unsigned m_cubDegree;
      stk::mesh::FieldBase *m_norm_field;
    public:

      Norm(stk::mesh::BulkData& bulkData, std::string partName, TurboOption turboOpt=TURBO_BUCKET, bool is_surface_norm=false) :
        FunctionOperator(bulkData, (stk::mesh::Part*)0), m_is_surface_norm(is_surface_norm), m_turboOpt(turboOpt), m_cubDegree(2), m_norm_field(0)
      {
        stk::mesh::Part * part = bulkData.mesh_meta_data().get_part(partName);
        if (!part) throw std::runtime_error(std::string("No part named ") +partName);
        init(part);
        error_check_is_surface_norm();
      }

      Norm(stk::mesh::BulkData& bulkData, MDArrayString& partNames, TurboOption turboOpt=TURBO_BUCKET, bool is_surface_norm=false) :
        FunctionOperator(bulkData, (stk::mesh::Part*)0), m_is_surface_norm(is_surface_norm), m_turboOpt(turboOpt), m_cubDegree(2), m_norm_field(0)
      {
        if (partNames.rank() != 1) throw std::runtime_error("Input array of strings should be rank 1 multi-d array (numpy array if from Python)");
        VERIFY_OP_ON(m_own_selector , ==, true, "logic error 1");
        VERIFY_OP_ON(m_selector, !=, 0, "logic error 2");
        delete m_selector;
        m_selector = new stk::mesh::Selector;
        for (size_t i = 0; i < partNames.extent_int(0); i++)
          {
            stk::mesh::Part * part = bulkData.mesh_meta_data().get_part(partNames(i));
            if (!part) throw std::runtime_error(std::string("No part named ") +partNames(i));
            *m_selector = (*m_selector) | (*part);
          }
        error_check_is_surface_norm();
      }

      Norm(stk::mesh::BulkData& bulkData, stk::mesh::Part *part = 0, TurboOption turboOpt=TURBO_BUCKET, bool is_surface_norm=false) :
        FunctionOperator(bulkData, part), m_is_surface_norm(is_surface_norm), m_turboOpt(turboOpt), m_cubDegree(2), m_norm_field(0)
      {
        // only select parts that are of the same rank
        if (!part)
          {
            stk::mesh::EntityRank rank = stk::topology::ELEMENT_RANK;
            if (m_is_surface_norm)
              {
                rank = bulkData.mesh_meta_data().side_rank();
              }
            if (!m_selector)
              m_selector = new stk::mesh::Selector();
            *m_selector = PerceptMesh::get_selector_of_rank(bulkData.mesh_meta_data(), rank);
          }
        // avoid inactive elements
        *m_selector &= !PerceptMesh::select_inactive_elements(bulkData);
        error_check_is_surface_norm();
      }

#ifndef SWIG
      Norm(stk::mesh::BulkData& bulkData, stk::mesh::Selector * selector,TurboOption turboOpt=TURBO_BUCKET,  bool is_surface_norm = false) :
        FunctionOperator(bulkData, selector), m_is_surface_norm(is_surface_norm), m_turboOpt(turboOpt), m_cubDegree(2), m_norm_field(0)
     {
        error_check_is_surface_norm();
     }
#endif

      virtual ~Norm() {}

      void setNormField(stk::mesh::FieldBase *norm_field) {m_norm_field=norm_field;}
      stk::mesh::FieldBase *getNormField() {return m_norm_field;}

      void setCubDegree(unsigned cubDegree) { m_cubDegree= cubDegree; }
      unsigned getCubDegree() { return m_cubDegree; }

      void set_is_surface_norm(bool is_surface_norm) {
        m_is_surface_norm = is_surface_norm;
        error_check_is_surface_norm();
      }
      bool get_is_surface_norm() { return m_is_surface_norm; }

      double evaluate(Function& integrand)
      {
        ConstantFunction sfx_res(0.0, "sfx_res");
        (*this)(integrand, sfx_res);
        return sfx_res.getValue();
      }

      /// if a Selector is specified with part(s) that are not auto-declared, make sure all parts are
      ///   of the same rank, and that m_is_surface_norm is set correctly (if not, warn...)
      void error_check_is_surface_norm()
      {
        stk::mesh::EntityRank element_rank = stk::topology::ELEMENT_RANK;
        stk::mesh::EntityRank side_rank    = m_bulkData.mesh_meta_data().side_rank();
        const stk::mesh::PartVector& parts = m_bulkData.mesh_meta_data().get_parts();
        stk::mesh::EntityRank all_ranks = stk::topology::NODE_RANK;
        unsigned nparts = parts.size();
        for (unsigned ipart=0; ipart < nparts; ipart++)
          {
            stk::mesh::Part& part = *parts[ipart];
            bool auto_part = 0 != part.attribute<AutoPart>();
            if (stk::mesh::is_auto_declared_part(part) || auto_part)
              continue;

            bool in_selector = (*m_selector)(part);
            stk::mesh::EntityRank rank = part.primary_entity_rank();
            //std::cout << "tmp srk Part= " << part.name() << " rank= " << rank << " in_selector= " << in_selector << std::endl;
            if (in_selector)
              {
                if (rank == element_rank || rank == side_rank)
                  {
                    if (all_ranks == 0)
                      all_ranks = rank;
                    //std::cout << "all_ranks= " << all_ranks << " rank= " << rank << std::endl;
                    if (rank != all_ranks)
                      {
                        throw std::runtime_error("all parts in the Selector must be the same rank");
                      }
                  }
              }
          }

        if ((all_ranks == 0 || all_ranks ==3) && m_is_surface_norm)
          {
            m_is_surface_norm = false;
            std::cout << "WARNING: Norm was constructed with is_surface_norm, but no surface Parts given. Turning is_surface_norm off." << std::endl;
          }
        if (all_ranks == 2 && !m_is_surface_norm)
          {
            m_is_surface_norm = true;
            std::cout << "WARNING: Norm was constructed with is_surface_norm=false, but surface Parts were given. Turning is_surface_norm on." << std::endl;
          }
      }

      virtual void operator()(Function& integrand, Function& result)
      {
        EXCEPTWATCH;
        m_is_norm_being_evaluated = true;
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
          VERIFY_OP_ON(m_turboOpt, !=, TURBO_NONE, "Can't use TURBO_NONE with Norm of FieldFunction");

        if (1 || typeid(integrand) == typeid(FieldFunction))
          {
            // FIXME - make all percept code const-correct
            LN_NormOp<Power> LN_op(integrand);
            CompositeFunction LN_of_integrand("LN_of_integrand", integrand, LN_op);
            IntegratedOp integrated_LN_op(LN_of_integrand, m_turboOpt);
            integrated_LN_op.setCubDegree(m_cubDegree);
            if (Power == -1) {
              integrated_LN_op.setAccumulationType(IntegratedOp::ACCUMULATE_MAX);
            }

            const stk::mesh::Part& locally_owned_part = m_bulkData.mesh_meta_data().locally_owned_part();
            stk::mesh::Selector selector(*m_selector & locally_owned_part);
            if (m_turboOpt == TURBO_NONE || m_turboOpt == TURBO_ELEMENT)
              {
                elementOpLoop(m_bulkData, integrated_LN_op, m_norm_field, &selector, m_is_surface_norm);
              }
            else if (m_turboOpt == TURBO_BUCKET)
              {
                bucketOpLoop(m_bulkData, integrated_LN_op, m_norm_field, &selector, m_is_surface_norm);
              }

            unsigned vec_sz = integrated_LN_op.getValue().size();
            std::vector<double> local = integrated_LN_op.getValue();

            if (Power == -1)
              {
                MaxOfNodeValues maxOfNodeValues(m_bulkData.mesh_meta_data().spatial_dimension(), integrand);
                nodalOpLoop(m_bulkData, maxOfNodeValues, 0, &selector);
		for (unsigned iDim = 0; iDim < local.size(); iDim++)
		  local[iDim] = std::max(local[iDim], maxOfNodeValues.maxVal[iDim]);
              }


            //unsigned p_rank = m_bulkData.parallel_rank();
            //std::cout  << "P["<<p_rank<<"] value = " << local << std::endl;

            std::vector<double> global_sum(vec_sz, 0.0);
            if (Power != -1)
              {
                stk::all_reduce_sum(m_bulkData.parallel(), &local[0], &global_sum[0], vec_sz);
              }
            else
              {
                for (unsigned iDim = 0; iDim < local.size(); iDim++)
                  {
                    double global = local[iDim];
                    stk::all_reduce( m_bulkData.parallel(), stk::ReduceMax<1>( &global ) );
                    global_sum[iDim] = global;
                  }
              }

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

        m_is_norm_being_evaluated = false;
      }

      //void operator()(FieldFunction& input);
    };


  }

#endif
