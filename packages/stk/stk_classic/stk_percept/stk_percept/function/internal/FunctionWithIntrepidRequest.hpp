#ifndef stk_percept_FunctionWithIntrepidRequest_hpp
#define stk_percept_FunctionWithIntrepidRequest_hpp

#include <cmath>
#include <math.h>
#include <map>

#include <typeinfo>

#include <stk_percept/function/MDArray.hpp>

#include <stk_percept/function/FunctionOperator.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/Function.hpp>
#include <stk_percept/function/internal/HasValue.hpp>
#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/ElementOp.hpp>
#include <stk_percept/function/BucketOp.hpp>

#include <stk_percept/PerceptMesh.hpp>

#include <stk_percept/norm/IntrepidManager.hpp>

namespace stk_classic
{
  namespace percept
  {

    class FunctionWithIntrepidRequest : public Function
    {
      Function *m_function;
    public:
      Function *getFunction() { return m_function; }
      typedef FieldBase *FieldType;
      typedef MDArray *ContainerType;
      typedef std::map<FieldType, ContainerType> Request;
      typedef Request::iterator RequestIterator;
      Request m_order_values;
      Request m_order_gradient;
      Request m_order_higherDerivs;

      FunctionWithIntrepidRequest() : Function("noname") {}
      FunctionWithIntrepidRequest(Function *func,
                                  Request values = Request(),
                                  Request gradient = Request(),
                                  Request higherDerivs = Request()) : Function(*func),
                                                                      m_function(func), m_order_values(values), m_order_gradient(gradient), m_order_higherDerivs(higherDerivs) {}
    };

    class ExampleFunctionWithIntrepidRequest : public FunctionWithIntrepidRequest
    {
      MDArray m_pressure_array;
      MDArray m_pressure_grad_array;
    public:

      ExampleFunctionWithIntrepidRequest(BulkData& bulkData) : FunctionWithIntrepidRequest()
      {
        mesh::fem::FEMMetaData& metaData = mesh::fem::FEMMetaData::get(bulkData);
        FieldBase *pressure =  metaData.get_field<FieldBase>("pressure");
        m_pressure_array.resize(10, 1);
        m_pressure_grad_array.resize(10, 1, 3);
        m_order_values[pressure] = &m_pressure_array;
        m_order_gradient[pressure] = &m_pressure_grad_array;
      }

      /// input_phy_points: ([C], [P], [D])
      /// output_integrand_values: ([C], [P], [DOF])
      virtual void operator()(MDArray& input_phy_points, MDArray& output_integrand_values, double time_value_optional=0.0)
      {
        // FIXME check args
        int nCells = input_phy_points.dimension(0);
        int nPoints = input_phy_points.dimension(1);
        int nDof = output_integrand_values.dimension(2);
        int start = 0;
//         if (which_cell >= 0)
//           {
//             nCells = 1;
//             start = which_cell;
//           }
        int end = start+nCells;

        for (int iCell = start; iCell < end; iCell++)
          {
            for (int iPoints = 0; iPoints < nPoints; iPoints++)
              {
                for (int iDof = 0; iDof < nDof; iDof++)
                  {
                    output_integrand_values(iCell, iPoints, iDof) = std::sin(m_pressure_array(iCell, iPoints, iDof));
                    // or
                    output_integrand_values(iCell, iPoints, iDof) = square(m_pressure_array(iCell, iPoints, iDof));
                  }
              }
          }
      }
    };

    // This one uses the basis functions passed in from Intrepid
    class Example2FunctionWithIntrepidRequest : public FunctionWithIntrepidRequest
    {
      MDArray m_pressure_array;
      MDArray m_pressure_grad_array;
    public:

      Example2FunctionWithIntrepidRequest(BulkData& bulkData) : FunctionWithIntrepidRequest()
      {
        mesh::fem::FEMMetaData& metaData = mesh::fem::FEMMetaData::get(bulkData);
        FieldBase *pressure = metaData.get_field<FieldBase>("pressure");
        m_pressure_array.resize(10, 1);
        m_pressure_grad_array.resize(10, 1, 3);
        m_order_values[pressure] = &m_pressure_array;
        m_order_gradient[pressure] = &m_pressure_grad_array;
      }

      virtual void operator()(MDArray& input_phy_points, MDArray& output_integrand_values, double time_value_optional=0.0)
      {
        int npts = input_phy_points.dimension(0);
        int start = 0;
//         if (which_point >= 0)
//           {
//             npts = 1;
//             start = which_point;
//           }
        int end = start+npts;

        //MDArray& basis_values = *m_order_values[0]; // yes, "0" - by convention, the basis functions are associated with a null FieldBase
        MDArray& basis_gradients = *m_order_gradient[0]; // yes, "0" - by convention, the basis functions are associated with a null FieldBase

        int nBasis = 0; // FIXME
        int nSpaceDim = 0; // FIXME
        for (int iBasis = 0; iBasis < nBasis; iBasis++)
          {
            for (int ipts = start; ipts < end; ipts++)
              {
                double sum=0;
                // compute grad_N_{iBasis}.grad(pressure)
                for (int iSpaceDim = 0; iSpaceDim < nSpaceDim; iSpaceDim++)
                  {
                    sum += basis_gradients(ipts, iBasis, iSpaceDim) * m_pressure_grad_array(ipts, 0, iSpaceDim);
                  }
                output_integrand_values(ipts, iBasis) = sum;
              }
          }
      }
    };


    class l2NormOpScalar : public FunctionWithIntrepidRequest
    {
      MDArray m_field_array;
      //MDArray m_field_grad_array;
    public:

      l2NormOpScalar( FieldBase *field) : FunctionWithIntrepidRequest()
      {
        //FEMMetaData& metaData = FEMMetaData::get(BulkData);

#if 0
        int npts=0; // FIXME
        m_field_array.resize(npts, 1);

        m_order_values[field] = &m_field_array;  // if null, calling code will supply this array (temp only - don't hold on to it)
#else
        m_order_values[field] = 0;  // if null, calling code will supply this array (temp only - don't hold on to it)
#endif
        //m_order_gradient[field] = &m_field_grad_array;
      }

      /// input_phy_points: ([C], [P], [D])
      /// output_integrand_values: ([C], [P], [DOF])
      virtual void operator()(MDArray& input_phy_points, MDArray& output_integrand_values, double time_value_optional=0.0)
      {
        // FIXME check args
        int nCells = input_phy_points.dimension(0);
        int nPoints = input_phy_points.dimension(1);
        int nDof = output_integrand_values.dimension(2);
        int start = 0;
//         if (which_cell >= 0)
//           {
//             nCells = 1;
//             start = which_cell;
//           }
        int end = start+nCells;

        for (int iCell = start; iCell < end; iCell++)
          {
            for (int iPoints = 0; iPoints < nPoints; iPoints++)
              {
                for (int iDof = 0; iDof < nDof; iDof++)
                  {
                    output_integrand_values(iCell, iPoints, iDof) = square(m_field_array(iCell, iPoints, iDof));
                  }
              }
          }
      }

    };

    class l2NormOpScalarFunction : public FunctionWithIntrepidRequest
    {
      Function& m_f;
    public:

      l2NormOpScalarFunction(Function& f ) : FunctionWithIntrepidRequest(), m_f(f)
      {
        // request nothing to be computed by Intrepid, just integrate this function
        m_order_values.clear();
        //m_order_gradient.resize(0);
      }

      /// input_phy_points: ([C], [P], [D])
      /// output_integrand_values: ([C], [P], [DOF])
      virtual void operator()(MDArray& input_phy_points, MDArray& output_integrand_values, double time_value_optional=0.0)
      {
        m_f(input_phy_points, output_integrand_values, time);

        // FIXME check args
        int nCells = input_phy_points.dimension(0);
        int nPoints = input_phy_points.dimension(1);
        int nDof = output_integrand_values.dimension(2);
        int start = 0;
//         if (which_cell >= 0)
//           {
//             nCells = 1;
//             start = which_cell;
//           }
        int end = start+nCells;

        for (int iCell = start; iCell < end; iCell++)
          {
            for (int iPoints = 0; iPoints < nPoints; iPoints++)
              {
                for (int iDof = 0; iDof < nDof; iDof++)
                  {
                    output_integrand_values(iCell, iPoints, iDof) = square(output_integrand_values(iCell, iPoints, iDof));
                  }
              }
          }
      }

    };


  }
}

#endif
