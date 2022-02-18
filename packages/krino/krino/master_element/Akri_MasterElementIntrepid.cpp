// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MasterElementIntrepid.hpp>
#include <Akri_MasterElement.hpp>
#include "Akri_MasterElementCalc.hpp"

#ifdef __INTEL_COMPILER
#include <Intrepid_FieldContainer.hpp>
#else
//FieldContainer has shadowed variables
//this disables the checking on GCC only
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <Intrepid_FieldContainer.hpp>
#pragma GCC diagnostic pop
#endif

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>
#include <Intrepid_Basis.hpp>
#include <Intrepid_DefaultCubatureFactory.hpp>
#include <Intrepid_Cubature.hpp>
#include <Intrepid_HGRAD_LINE_Cn_FEM.hpp>
#include <Shards_CellTopology.hpp>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/MetaData.hpp> // for get_cell_topology

namespace krino {

const MasterElementIntrepid &
MasterElementIntrepid::getMasterElement(stk::topology t, const unsigned spatial_dimension)
{
  static std::vector<std::unique_ptr<MasterElementIntrepid>> all_master_elems(stk::topology::BEGIN_TOPOLOGY + stk::topology::NUM_TOPOLOGIES);
  std::unique_ptr<MasterElementIntrepid> & master_elem = all_master_elems[t()];
  if (nullptr == master_elem.get())
  {
    std::unique_ptr<Intrepid::Basis<double, Intrepid::FieldContainer<double>>> basis;
    switch(t())
    {
    case stk::topology::LINE_2:
        basis = std::make_unique<Intrepid::Basis_HGRAD_LINE_C1_FEM<double, Intrepid::FieldContainer<double>>>();
        break;
    case stk::topology::LINE_3:
        basis = std::make_unique<Intrepid::Basis_HGRAD_LINE_Cn_FEM<double, Intrepid::FieldContainer<double>>>(2, Intrepid::POINTTYPE_SPECTRAL);
        break;
    case stk::topology::TRI_3:
        basis = std::make_unique<Intrepid::Basis_HGRAD_TRI_C1_FEM<double, Intrepid::FieldContainer<double>>>();
        break;
    case stk::topology::TRI_6:
        basis = std::make_unique<Intrepid::Basis_HGRAD_TRI_C2_FEM<double, Intrepid::FieldContainer<double>>>();
        break;
    case stk::topology::QUAD_4:
        basis = std::make_unique<Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, Intrepid::FieldContainer<double>>>();
        break;
    case stk::topology::QUAD_9:
        basis = std::make_unique<Intrepid::Basis_HGRAD_QUAD_C2_FEM<double, Intrepid::FieldContainer<double>>>();
        break;
    case stk::topology::TRI_3_2D:
        basis = std::make_unique<Intrepid::Basis_HGRAD_TRI_C1_FEM<double, Intrepid::FieldContainer<double>>>();
        break;
    case stk::topology::TRI_6_2D:
        basis = std::make_unique<Intrepid::Basis_HGRAD_TRI_C2_FEM<double, Intrepid::FieldContainer<double>>>();
        break;
    case stk::topology::QUAD_4_2D:
        basis = std::make_unique<Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, Intrepid::FieldContainer<double>>>();
        break;
    case stk::topology::QUAD_9_2D:
        basis = std::make_unique<Intrepid::Basis_HGRAD_QUAD_C2_FEM<double, Intrepid::FieldContainer<double>>>();
        break;
    case stk::topology::TET_4:
        basis = std::make_unique<Intrepid::Basis_HGRAD_TET_C1_FEM<double, Intrepid::FieldContainer<double>>>();
        break;
    case stk::topology::TET_10:
        basis = std::make_unique<Intrepid::Basis_HGRAD_TET_C2_FEM<double, Intrepid::FieldContainer<double>>>();
        break;
    case stk::topology::HEX_8:
        basis = std::make_unique<Intrepid::Basis_HGRAD_HEX_C1_FEM<double, Intrepid::FieldContainer<double>>>();
        break;
    case stk::topology::HEX_27:
        basis = std::make_unique<Intrepid::Basis_HGRAD_HEX_C2_FEM<double, Intrepid::FieldContainer<double>>>();
        break;
    default:
        throw std::runtime_error("Element topology not found in MasterElementIntrepid::build: " + t.name());
        break;
    }
    master_elem = std::make_unique<MasterElementIntrepid>(t, std::move(basis), spatial_dimension);
  }
  return *master_elem;
}

MasterElementIntrepid::MasterElementIntrepid(
    stk::topology topology,
    std::unique_ptr<Intrepid::Basis<double, Intrepid::FieldContainer<double> > > basis,
    unsigned spatial_dimension)
: m_topology(topology),
  m_numCoordDims(spatial_dimension),
  m_intrepidBasis(std::move(basis))
{
  shards::CellTopology cellType = stk::mesh::get_cell_topology(topology);

  // set the cubature
  Intrepid::DefaultCubatureFactory<double> cubatureFactory;
  Teuchos::RCP<Intrepid::Cubature<double>> intrepidCubature = cubatureFactory.create(cellType, 2*m_intrepidBasis->getDegree());
  m_numIntgPts = intrepidCubature->getNumPoints();

  m_numNodes = topology.num_nodes();
  m_numElemDims  = intrepidCubature->getDimension();

  // Allocate reference data
  m_shapeFuncs.resize(m_numIntgPts*m_numNodes);
  m_pointGrads.resize(m_numIntgPts*m_numNodes*m_numElemDims );
  m_refPoints.resize(m_numIntgPts*m_numElemDims);
  m_refWeights.resize(m_numIntgPts);
  m_refCoords.resize(m_numNodes*m_numElemDims);

  // retrieve the cubature points and weights
  std::vector<int> refPointsDims = {m_numIntgPts, m_numElemDims};
  Intrepid::FieldContainer<double> refPointsFC( refPointsDims, m_refPoints.data() );
  std::vector<int> refWeightsDims = {m_numIntgPts};
  Intrepid::FieldContainer<double> refWeightsFC( refWeightsDims, m_refWeights.data() );
  intrepidCubature->getCubature(refPointsFC, refWeightsFC);

  // compute the refernce values and gradients at the integration points
  Intrepid::FieldContainer<double> pointVals(m_numNodes, m_numIntgPts);
  Intrepid::FieldContainer<double> pointGrads(m_numNodes, m_numIntgPts, m_numElemDims);
  m_intrepidBasis->getValues(pointVals, refPointsFC, Intrepid::OPERATOR_VALUE);
  m_intrepidBasis->getValues(pointGrads, refPointsFC, Intrepid::OPERATOR_GRAD);

  // re-order shape functions for consistency with other master elements
  for ( int ip(0); ip < m_numIntgPts; ++ip ) {
    for ( int node(0); node < m_numNodes; ++node ) {
      m_shapeFuncs[ip*m_numNodes + node] = pointVals(node, ip);
    }
  }
  for ( int ip(0); ip < m_numIntgPts; ++ip ) {
    for ( int node(0); node < m_numNodes; ++node ) {
      for ( int dim(0); dim < m_numElemDims; ++dim ) {
        m_pointGrads[(ip*m_numNodes + node)*m_numElemDims + dim] = pointGrads(node, ip, dim);
      }
    }
  }

  for ( int node(0); node < m_numNodes; ++node ) {
    const double * node_coords = Intrepid::CellTools<double>::getReferenceNode( cellType, node );
    for ( int dim(0); dim < m_numElemDims; ++dim ) {
      m_refCoords[node*m_numElemDims + dim] = node_coords[dim];
    }
  }
}

void
MasterElementIntrepid::determinant(
    const int nelem,
    const double* coords,  // (nvec,npe,nelem)
    double* det_J,         // (nelem,nint)
    double* error ) const  // (nelem)
{
  determinant(m_numIntgPts, m_numNodes, m_pointGrads.data(), nelem, coords, det_J, error);
}

void
MasterElementIntrepid::shape_fcn(
    const int nint,  // returns array(npe,nint)
    const double* p_coords,
    double* result) const
{
  Intrepid::FieldContainer<double> pointVals(m_numNodes, nint);

  // create the pcoordVec FC from the p_coords ptr and it's dimensions
  std::vector<int> pcoordDims = { nint, m_numElemDims };
  Intrepid::FieldContainer<double> pcoordVec( pcoordDims, const_cast< double * >(p_coords) );

  // compute the shape function values at the integration points
  m_intrepidBasis->getValues(pointVals, pcoordVec, Intrepid::OPERATOR_VALUE);

  // re-order shape functions for consistency with other master elements
  for ( int ip(0); ip < nint; ++ip ) {
    for ( int node(0); node < m_numNodes; ++node ) {
      result[ip*m_numNodes + node] = pointVals(node,ip);
    }
  }
}

void
MasterElementIntrepid::shape_fcn_deriv(
    const int nint,
    const double* p_coords,
    double* result ) const
{
  Intrepid::FieldContainer<double> resultVec(m_numNodes, nint, m_numElemDims);

  // create the pcoordVec FC from the p_coords ptr and it's dimensions
  std::vector<int> pcoordDims = { nint, m_numElemDims };
  Intrepid::FieldContainer<double> pcoordVec( pcoordDims, const_cast< double * >(p_coords) );

  m_intrepidBasis->getValues(resultVec, pcoordVec, Intrepid::OPERATOR_GRAD);

  // re-order shape functions for consistency with other master elements
  for ( int ip(0); ip < nint; ++ip ) {
    for ( int node(0); node < m_numNodes; ++node ) {
      for ( int dim(0); dim < m_numElemDims; ++dim ) {
        result[(ip*m_numNodes + node)*m_numElemDims + dim] = resultVec(node,ip,dim);
      }
    }
  }
}

void
MasterElementIntrepid::interpolate_point(
    const int  npar_coord,
    const double * par_coord,      // (npar_coord)
    const int  ncomp_field,
    const double * field,          // (num_nodes,ncomp_field)
    double * result ) const        // (ncomp_field)
{
  std::vector<double> shape(m_numNodes);
  shape_fcn(1, par_coord, shape.data());
  for ( int comp(0); comp < ncomp_field; ++comp ) {
    result[ comp ] = 0.0;
    for ( int node(0); node < m_numNodes; ++node ) {
      result[ comp ] += shape[node] * field[ m_numNodes * comp + node ];
    }
  }
}

void
MasterElementIntrepid::scalar_gradient(
    const int   nelem,      //:  number of elements to process
    const double* gradop,     //: (nvec,npe,nelem,nint)
    const double* det_J,      //: (nelem,nint)
    const double* sfield,     //: (npe,nelem)
    double* vector ) const    //: (nvec,nelem,nint)
{
  MasterElementCalc::scalar_gradient(m_numIntgPts, nelem, m_numElemDims, m_numNodes, gradop, det_J, sfield, vector);
}

void
MasterElementIntrepid::scalar_gradient(
        const int   nint,       //:  number of intg points
        const int   nelem,      //:  number of elements to process
        const double* gradop,     //: (nvec,npe,nelem,nint)
        const double* det_J,      //: (nelem,nint)
        const double* sfield,     //: (npe,nelem)
        double* vector ) const    //: (nvec,nelem,nint)
{
  MasterElementCalc::scalar_gradient(nint, nelem, m_numElemDims, m_numNodes, gradop, det_J, sfield, vector);
}

void
MasterElementIntrepid::determinant(
      const int  nint,
      const int npe_g,
      const double* deriv_g,  // (nvec,npe_g,nint)
      const int   nelem,
      const double* coords,   // (nvec,npe,nelem)
      double* det_J,          // (nelem,nint)
      double* error ) const   // (nelem)
{
  MasterElementCalc::determinant(m_numElemDims, m_numCoordDims, nint, npe_g, deriv_g, nelem, coords, det_J, error);
}

void
MasterElementIntrepid::gradient_operator(
        const int nint,
        const int npe_g,
        const double* deriv_g,   // (nvec,npe_g,nint)
        const int npe_f,
        const double* deriv_f,   // (nvec,npe_f,nint)
        const int   nelem,
        const double* coords,    // (nvec,npe,nelem)
        double* gradop,          // (nvec,npe,nelem,nint)
        double* det_J,           // (nelem,nint)
        double* error) const
{
  ThrowRequireMsg(m_numElemDims == m_numCoordDims, "MasterElementHybrid::gradient_operator does not support lower rank elements in higher dimensions (e.g. BAR,QUAD,TRI in 3D).");
  MasterElementCalc::gradient_operator(m_numElemDims, nint, npe_g, deriv_g, npe_f, deriv_f, nelem, coords, gradop, det_J, error);
}

} // namespace krino
