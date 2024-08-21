// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_ComputeBases_hpp
#define percept_ComputeBases_hpp

#include <cmath>
#include <math.h>

#include <typeinfo>

#include "Intrepid2_FunctionSpaceTools.hpp"

#include <percept/function/MDArray.hpp>

#include <percept/function/FunctionOperator.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/Function.hpp>
#include <percept/function/internal/HasValue.hpp>
#include <percept/function/StringFunction.hpp>

#include <percept/PerceptMesh.hpp>

#if defined(STK_PERCEPT_USE_INTREPID)
#include <percept/element/intrepid/BasisTable.hpp>
#endif

  namespace percept
  {

    class ComputeBases 
    {
      unsigned m_cached_topo_key;
      BasisTable::BasisType* m_cached_basis;
    public:
      ComputeBases() : m_cached_topo_key(0), m_cached_basis(0) {}

      // [P] = 1 (or numInterpPoints)
      // [C] = numCells
      // [B] = numBases
      // [D] = spatialDim

      /// ([P],[D])
      // parametric_coordinates: ([P],[D])
      // transformed_basis_values: ([C],[B],[P]), or ([C],[B],[P],[D]) for GRAD
      void getBases(const stk::mesh::Bucket &bucket, const MDArray& parametric_coordinates, MDArray& transformed_basis_values, int which_cell = -1)
      {
        VERIFY_OP(parametric_coordinates.rank(), ==, 2, "ComputeBases::operator() parametric_coordinates bad rank");
        VERIFY_OP(transformed_basis_values.rank(), ==, 3, "ComputeBases::operator() transformed_basis_values bad rank");
        VERIFY_OP(transformed_basis_values.extent_int(2), ==, parametric_coordinates.extent_int(0), 
                  "ComputeBases::operator() transformed_basis_values.dim(2) != parametric_coordinates.dim(0)");
        //VERIFY_OP(output_field_values.rank(), ==, 2, "FieldFunction::operator() output_field_values bad rank");

        // [P] = number of integration points
        int numInterpPoints = parametric_coordinates.extent_int(0);

        //const stk::mesh::Bucket & bucket = element->bucket();
        const CellTopologyData * const bucket_cell_topo_data = stk::mesh::get_cell_topology(bucket.topology()).getCellTopologyData();

        //unsigned stride = 0;
        //double * fdata_bucket = stk::mesh::field_data( m_my_field , bucket, &stride);
        // intentionally ignoring return value to get around compiler warning
        //stk::mesh::field_data( m_my_field , bucket, &stride);
        //unsigned nDOF = stride;
        //int nOutDim = m_codomain_dimensions.back(); // FIXME for tensor
        // FIXME 
        //VERIFY_OP((int)nDOF, == , nOutDim,
        //"FieldFunction::operator(): invalid dimensions nDof, m_codomain_dimensions[0]= ");

        int numCells = bucket.size(); // FIXME for multiple cells
        if (which_cell >= 0)
          numCells = 1;

        shards::CellTopology topo(bucket_cell_topo_data);
        int numNodes = topo.getNodeCount();
        int cellDim  = topo.getDimension();
        if (0)
          {
            MDArray cellWorkset("cellWorkset", numCells, numNodes, cellDim);
            if (0) cellWorkset(0,0,0) = 0.0;
          }

        // map cell topology to a basis
        BasisTable::BasisType *basis;
        //PerceptMesh::BasisTypeRCP basis;
        if (m_cached_topo_key != bucket_cell_topo_data->key)
          {
            BasisTable::BasisTypeRCP basisRCP = BasisTable::getInstance()->getBasis(topo);
            basis = basisRCP.get();
            m_cached_basis = basis;
            m_cached_topo_key = bucket_cell_topo_data->key;
          }
        else
          {
            basis = m_cached_basis;
          }

        int numBases = basis->getCardinality();
        if (numBases != numNodes)
          {
            throw std::runtime_error(" (numBases != numNodes) ");
          }

        // ([B],[P]), or ([B],[P],[D]) for GRAD
        MDArray basis_values("basis_values", numBases, numInterpPoints); 

        // ([C],[B],[P]), or ([C],[B],[P],[D]) for GRAD
        //MDArray transformed_basis_values(numCells, numBases, numInterpPoints); 

        basis->getValues(basis_values, parametric_coordinates, Intrepid2::OPERATOR_VALUE);

        // this function just spreads (copies) the values of the basis to all elements in the workset (numCells)
        Intrepid2::FunctionSpaceTools<Kokkos::HostSpace>::HGRADtransformVALUE(transformed_basis_values, basis_values);
      }

      void getBases(const stk::mesh::BulkData& bulk, const stk::mesh::Entity element, const MDArray& parametric_coordinates, MDArray& transformed_basis_values)
      {
        // FIXME this will need to be changed when the bases have gradients and thus depend on which element they're assoc with
        const stk::mesh::Bucket & bucket = bulk.bucket(element);
        getBases(bucket, parametric_coordinates, transformed_basis_values, 0);
      }

    };
  }

#endif
