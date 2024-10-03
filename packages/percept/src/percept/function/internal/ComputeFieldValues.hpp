// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_ComputeFieldValues_hpp
#define percept_ComputeFieldValues_hpp

#include <cmath>
#include <math.h>

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
    class ComputeFieldValues
    {
    public:

      // FIXME input_field_data_values: ([C],[F],[DOF])

      // transformed_basis_values: ([C],[F],[P]), or ([C],[F],[P],[D]) for GRAD
      // output_field_values: ([C],[P],[DOF])
      void get_fieldValues(const stk::mesh::BulkData& bulk, const stk::mesh::Entity element, MDArray& transformed_basis_values, stk::mesh::FieldBase* field, MDArray& output_field_values)
      {
        VERIFY_OP(output_field_values.rank(), ==, 3, "FieldValuesComputer::get_fieldValues output_field_values bad rank");
        VERIFY_OP(transformed_basis_values.rank(), ==, 3, "FieldValuesComputer::get_fieldValues transformed_basis_values bad rank");
        VERIFY_OP(output_field_values.extent_int(0), ==, transformed_basis_values.extent_int(0),
                  "FieldValuesComputer::get_fieldValues output_field_values.dim(0) doesn't match transformed_basis_values.dim(0)");
        VERIFY_OP(output_field_values.extent_int(1), ==, transformed_basis_values.extent_int(1),
                  "FieldValuesComputer::get_fieldValues output_field_values.dim(1) doesn't match transformed_basis_values.dim(1)");

        // [P] = num integration points
        int numInterpPoints = transformed_basis_values.extent_int(2);

        unsigned stride = 0;
        //double * fdata_bucket = stk::mesh::field_data( m_my_field , bucket, &stride);
        // intentionally ignoring return value to get around compiler warning
        //stk::mesh::field_data( field , bucket, &stride);
        unsigned nDOF = stride;

#ifndef NDEBUG
        int nOutDim = output_field_values.extent_int(2); // FIXME for tensor
        VERIFY_OP((int)nDOF, == , nOutDim,
                  "FieldValuesComputer::get_fieldValues: invalid dimensions nDof, m_codomain_dimensions[0]= ");
#endif

        int numCells = transformed_basis_values.extent_int(0); // FIXME for multiple cells

//         shards::CellTopology topo(bucket_cell_topo_data);
//         int numNodes = topo.getNodeCount();
//         int cellDim  = topo.getDimension();

//         if (0)
//           {
//             MDArray cellWorkset(numCells, numNodes, cellDim);
//             if (0) cellWorkset(0,0,0) = 0.0;
//           }

        int numBases = transformed_basis_values.extent_int(1);
        int numNodes = numBases;  // FIXME

        // ([C],[F],[P]), or ([C],[F],[P],[D]) for GRAD
        //MDArray transformed_basis_values(numCells, numBases, numInterpPoints);

        // FIXME - it appears that Intrepid2 only supports the evaluation of scalar-valued fields, so we have
        //   to copy the field one DOF at a time into a local array, evaluate, then copy back
        // ([C],[F])
        MDArray field_data_values("field_data_values", numCells, numBases);

        stk::mesh::Entity const* elem_nodes = bulk.begin_nodes(element);

        // ([P],[D])  [P] points in [D] dimensions

        // ([C],[P]) - place for results of evaluation
        MDArray loc_output_field_values("loc_output_field_values", numCells, numInterpPoints);

        // gather
        for (unsigned iDOF = 0; iDOF < nDOF; iDOF++)
          {
            for (int iCell = 0; iCell < numCells; iCell++)
              {
                for (int iNode = 0; iNode < numNodes; iNode++)
                  {
                    stk::mesh::Entity node = elem_nodes[iNode];
                    double * fdata = (double*)stk::mesh::field_data( *field , node);
                    field_data_values(iCell, iNode) = fdata[iDOF];
                  }
              }

            /// NOTE: this is needed since FunctionSpaceTools::evaluate method assumes the output array is initialized to 0
            Kokkos::deep_copy(loc_output_field_values,0.0);
            Intrepid2::FunctionSpaceTools<Kokkos::HostSpace>::evaluate(loc_output_field_values, field_data_values, transformed_basis_values);

            for (int iCell = 0; iCell < numCells; iCell++)
              {
                for (int iPoint = 0; iPoint < numInterpPoints; iPoint++)
                  {
                    output_field_values(iCell, iPoint, iDOF) = loc_output_field_values(iCell, iPoint);
                  }
              }
          }
      }
    };


  }

#endif
