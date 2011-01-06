#ifndef stk_percept_ComputeFieldValues_hpp
#define stk_percept_ComputeFieldValues_hpp

#include <cmath>
#include <math.h>

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

namespace stk
{
  namespace percept
  {
    class ComputeFieldValues
    {
    public:

      // FIXME input_field_data_values: ([C],[F],[DOF])

      // transformed_basis_values: ([C],[F],[P]), or ([C],[F],[P],[D]) for GRAD
      // output_field_values: ([C],[P],[DOF])
      void getFieldValues(const stk::mesh::Entity& element, MDArray& transformed_basis_values, FieldBase* field, MDArray& output_field_values)
      {
        VERIFY_OP(output_field_values.rank(), ==, 3, "FieldValuesComputer::getFieldValues output_field_values bad rank");
        VERIFY_OP(transformed_basis_values.rank(), ==, 3, "FieldValuesComputer::getFieldValues transformed_basis_values bad rank");
        VERIFY_OP(output_field_values.dimension(0), ==, transformed_basis_values.dimension(0), 
                  "FieldValuesComputer::getFieldValues output_field_values.dim(0) doesn't match transformed_basis_values.dim(0)");
        VERIFY_OP(output_field_values.dimension(1), ==, transformed_basis_values.dimension(1), 
                  "FieldValuesComputer::getFieldValues output_field_values.dim(1) doesn't match transformed_basis_values.dim(1)");

        // [P] = num integration points
        int numInterpPoints = transformed_basis_values.dimension(2);

        unsigned stride = 0;
        //double * fdata_bucket = PerceptMesh::field_data( m_my_field , bucket, &stride);
        // intentionally ignoring return value to get around compiler warning
        //PerceptMesh::field_data( field , bucket, &stride);
        unsigned nDOF = stride;
        int nOutDim = output_field_values.dimension(2); // FIXME for tensor
        VERIFY_OP((int)nDOF, == , nOutDim,
                  "FieldValuesComputer::getFieldValues: invalid dimensions nDof, m_codomain_dimensions[0]= ");

        int numCells = transformed_basis_values.dimension(0); // FIXME for multiple cells

//         shards::CellTopology topo(bucket_cell_topo_data);
//         int numNodes = topo.getNodeCount();
//         int cellDim  = topo.getDimension();

//         if (0)
//           {
//             MDArray cellWorkset(numCells, numNodes, cellDim);
//             if (0) cellWorkset(0,0,0) = 0.0;
//           }

        int numBases = transformed_basis_values.dimension(1);
        int numNodes = numBases;  // FIXME

        // ([C],[F],[P]), or ([C],[F],[P],[D]) for GRAD
        //MDArray transformed_basis_values(numCells, numBases, numInterpPoints); 

        // FIXME - it appears that Intrepid only supports the evaluation of scalar-valued fields, so we have
        //   to copy the field one DOF at a time into a local array, evaluate, then copy back
        // ([C],[F])
        MDArray field_data_values(numCells, numBases);

        const mesh::PairIterRelation elem_nodes = element.relations( mesh::Node );

        // ([P],[D])  [P] points in [D] dimensions

        // ([C],[P]) - place for results of evaluation
        MDArray loc_output_field_values(numCells, numInterpPoints);

        unsigned stride_node = 0;

        // gather
        for (unsigned iDOF = 0; iDOF < nDOF; iDOF++)
          {
            for (int iCell = 0; iCell < numCells; iCell++)
              {
                for (int iNode = 0; iNode < numNodes; iNode++)
                  {
                    mesh::Entity& node = *elem_nodes[iNode].entity();
                    double * fdata = PerceptMesh::field_data( field , node, &stride_node);
                    field_data_values(iCell, iNode) = fdata[iDOF];
                  }
              }

            /// NOTE: this is needed since FunctionSpaceTools::evaluate method assumes the output array is initialized to 0
            loc_output_field_values.initialize(0.0);
            FunctionSpaceTools::evaluate<double>(loc_output_field_values, field_data_values, transformed_basis_values);

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
}
#endif
