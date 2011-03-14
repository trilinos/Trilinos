/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef stk_percept_Norm_hpp
#define stk_percept_Norm_hpp

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

#include <math.h>

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/function/MDArray.hpp>

#include <stk_percept/function/FunctionOperator.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/Function.hpp>
#include <stk_percept/function/CompositeFunction.hpp>
#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/function/ElementOp.hpp>
#include <stk_percept/function/BucketOp.hpp>
#include <stk_percept/function/MultipleFieldFunction.hpp>

#include <stk_percept/function/internal/HasValue.hpp>
#include <stk_percept/function/internal/IntegratedOp.hpp>

#include <stk_percept/ExceptionWatch.hpp>

#include <stk_percept/norm/IntrepidManager.hpp>

namespace stk
{
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

#if 0
    class l2NormOp  : public Function, public HasFinalOp<double>
    {
    public:

      /// integrand tells what fields Intrepid should compute, etc.
      l2NormOp() : Function() {}

      void operator()(MDArray& integrand_values, MDArray& output_values, double time_value_optional=0.0)
      {
        VERIFY_OP(integrand_values.size(), ==, output_values.size(), "l2NormOp::operator() bad sizes");
        for (int i = 0; i < integrand_values.size(); i++)
          {
            output_values[i] = square(integrand_values[i]);
          }
      }
      using Function::operator();
      virtual void operator()(MDArray& domain, MDArray& codomain, const stk::mesh::Entity& element, const MDArray& parametric_coords, double time_value_optional=0.0)
      {
        (*this)(domain, codomain, time_value_optional);
      }

      void finalOp(const double& vin, double& vout)
      {
        vout = std::sqrt(vin);
      }
    };
#endif

    template<int Power=2>
    class LN_NormOp  : public Function, public HasFinalOp<std::vector<double> >
    {
    public:

      /// integrand tells what fields Intrepid should compute, etc.
      /// Note: this function is intended to be used to wrap a Function using CompositeFunction and thus its domain and codomain
      /// are the same as the wrapped function's codomain
      LN_NormOp(Function& integrand) : Function("LN_NormOp",integrand.getCodomainDimensions(), integrand.getCodomainDimensions()) {}

      void operator()(MDArray& integrand_values, MDArray& output_values, double time_value_optional=0.0)
      {
        VERIFY_OP(integrand_values.size(), ==, output_values.size(), "LN_NormOp::operator() bad sizes");
        for (int i = 0; i < integrand_values.size(); i++)
          {
            output_values[i] = std::pow(std::fabs(integrand_values[i]), double(Power) );
          }
      }
      virtual void operator()(MDArray& domain, MDArray& codomain, const stk::mesh::Entity& element, const MDArray& parametric_coords, double time_value_optional=0.0)
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
          vout[i] = std::pow(vin[i], 1./(double(Power)) );
      }
    };


    template<int Power=2>
    class Norm : public FunctionOperator
    {
      TurboOption m_turboOpt;
    public:
      Norm(mesh::BulkData& bulkData, mesh::Part *part = 0, TurboOption turboOpt=TURBO_NONE) :
        FunctionOperator(bulkData, part), m_turboOpt(turboOpt)
     {}

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

        if (1 || typeid(integrand) == typeid(FieldFunction))
          {
            // FIXME - make all stk::percept code const-correct
            PerceptMesh eMesh(&MetaData::get(m_bulkData), &m_bulkData);
            LN_NormOp<Power> LN_op(integrand);
            CompositeFunction LN_of_integrand("LN_of_integrand", integrand, LN_op);
            IntegratedOp integrated_LN_op(LN_of_integrand, m_turboOpt);
            eMesh.printInfo("Norm");
            if (m_turboOpt == TURBO_NONE || m_turboOpt == TURBO_ELEMENT)
              {
                eMesh.elementOpLoop(integrated_LN_op, 0, &MetaData::get(m_bulkData).locally_owned_part());
              }
            else if (m_turboOpt == TURBO_BUCKET)
              {
                eMesh.bucketOpLoop(integrated_LN_op, 0, &MetaData::get(m_bulkData).locally_owned_part());
              }

            unsigned vec_sz = integrated_LN_op.getValue().size();
            std::vector<double>& local = integrated_LN_op.getValue();

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
