// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_FunctionWithIntrepidRequest_hpp
#define percept_FunctionWithIntrepidRequest_hpp

#include <cmath>
#include <math.h>
#include <map>

#include <typeinfo>

#include <percept/function/MDArray.hpp>

#include <percept/function/FunctionOperator.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/Function.hpp>
#include <percept/function/internal/HasValue.hpp>
#include <percept/function/StringFunction.hpp>
#include <percept/function/ElementOp.hpp>
#include <percept/function/BucketOp.hpp>

#include <percept/PerceptMesh.hpp>

#include <percept/norm/IntrepidManager.hpp>

  namespace percept
  {

    class FunctionWithIntrepid2Request : public Function
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

      FunctionWithIntrepid2Request() : Function("noname") {}
      FunctionWithIntrepid2Request(Function *func,
                                  Request values = Request(),
                                  Request gradient = Request(),
                                  Request higherDerivs = Request()) : Function(*func),
                                                                      m_function(func), m_order_values(values), m_order_gradient(gradient), m_order_higherDerivs(higherDerivs) {}
    };

    class ExampleFunctionWithIntrepidRequest : public FunctionWithIntrepid2Request
    {
      MDArray m_pressure_array;
      MDArray m_pressure_grad_array;
    public:

      ExampleFunctionWithIntrepid2Request(BulkData& bulkData) : FunctionWithIntrepid2Request()
      {
        stk::mesh::MetaData& metaData = stk::mesh::MetaData::get(bulkData);
        FieldBase *pressure =  metaData.get_field(stk::topology::NODE_RANK, "pressure");
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

    // This one uses the basis functions passed in from Intrepid2
    class Example2FunctionWithIntrepid2Request : public FunctionWithIntrepid2Request
    {
      MDArray m_pressure_array;
      MDArray m_pressure_grad_array;
    public:

      Example2FunctionWithIntrepid2Request(BulkData& bulkData) : FunctionWithIntrepid2Request()
      {
        stk::mesh::MetaData& metaData = stk::mesh::MetaData::get(bulkData);
        FieldBase *pressure = metaData.get_field(stk::topology::NODE_RANK, "pressure");
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


    class l2NormOpScalar : public FunctionWithIntrepid2Request
    {
      MDArray m_field_array;
      //MDArray m_field_grad_array;
    public:

      l2NormOpScalar( FieldBase *field) : FunctionWithIntrepid2Request()
      {
        //MetaData& metaData = MetaData::get(BulkData);

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

    class l2NormOpScalarFunction : public FunctionWithIntrepid2Request
    {
      Function& m_f;
    public:

      l2NormOpScalarFunction(Function& f ) : FunctionWithIntrepid2Request(), m_f(f)
      {
        // request nothing to be computed by Intrepid2, just integrate this function
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

#endif
