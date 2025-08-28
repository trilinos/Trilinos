// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MasterElementIntrepid.hpp>
#include <Akri_MasterElementCalc.hpp>

#include <Kokkos_DynRankView.hpp>
#include <Intrepid2_Basis.hpp>
#include <Intrepid2_Cubature.hpp>
#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_DefaultCubatureFactory.hpp>
#include <Intrepid2_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid2_HGRAD_HEX_C2_FEM.hpp>
#include <Intrepid2_HGRAD_LINE_C1_FEM.hpp>
#include <Intrepid2_HGRAD_LINE_C2_FEM.hpp>
#include <Intrepid2_HGRAD_QUAD_C1_FEM.hpp>
#include <Intrepid2_HGRAD_QUAD_C2_FEM.hpp>
#include <Intrepid2_HGRAD_TET_C1_FEM.hpp>
#include <Intrepid2_HGRAD_TET_C2_FEM.hpp>
#include <Intrepid2_HGRAD_TRI_C1_FEM.hpp>
#include <Intrepid2_HGRAD_TRI_C2_FEM.hpp>
#include <Shards_CellTopology.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino {

static std::unique_ptr<MasterElementIntrepid::IntrepidBasis> get_basis_for_topology(stk::topology t)
{
  using PointType = double;
  using OutputType = double;
  using ExecutionSpace = MasterElementIntrepid::ExecutionSpace;

  switch(t())
  {
  case stk::topology::LINE_2:
      return std::make_unique<Intrepid2::Basis_HGRAD_LINE_C1_FEM<ExecutionSpace, OutputType, PointType>>();
  case stk::topology::LINE_3:
      return std::make_unique<Intrepid2::Basis_HGRAD_LINE_C2_FEM<ExecutionSpace, OutputType, PointType>>();
  case stk::topology::TRI_3:
  case stk::topology::TRI_3_2D:
      return std::make_unique<Intrepid2::Basis_HGRAD_TRI_C1_FEM<ExecutionSpace, OutputType, PointType>>();
  case stk::topology::TRI_6:
  case stk::topology::TRI_6_2D:
      return std::make_unique<Intrepid2::Basis_HGRAD_TRI_C2_FEM<ExecutionSpace, OutputType, PointType>>();
  case stk::topology::QUAD_4:
  case stk::topology::QUAD_4_2D:
      return std::make_unique<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<ExecutionSpace, OutputType, PointType>>();
  case stk::topology::QUAD_9:
  case stk::topology::QUAD_9_2D:
      return std::make_unique<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<ExecutionSpace, OutputType, PointType>>();
  case stk::topology::TET_4:
      return std::make_unique<Intrepid2::Basis_HGRAD_TET_C1_FEM<ExecutionSpace, OutputType, PointType>>();
  case stk::topology::TET_10:
      return std::make_unique<Intrepid2::Basis_HGRAD_TET_C2_FEM<ExecutionSpace, OutputType, PointType>>();
  case stk::topology::HEX_8:
      return std::make_unique<Intrepid2::Basis_HGRAD_HEX_C1_FEM<ExecutionSpace, OutputType, PointType>>();
  case stk::topology::HEX_27:
      return std::make_unique<Intrepid2::Basis_HGRAD_HEX_C2_FEM<ExecutionSpace, OutputType, PointType>>();
  default:
      throw std::runtime_error("Element topology not found in get_basis_for_topology(): " + t.name());
  }
}

static double parametric_volume_for_topology(stk::topology topology)
{
  switch(topology())
    {
    case stk::topology::LINE_2:
    case stk::topology::LINE_3:
        return 2.;
    case stk::topology::TRI_3:
    case stk::topology::TRI_6:
    case stk::topology::TRI_3_2D:
    case stk::topology::TRI_6_2D:
        return 0.5;
    case stk::topology::QUAD_4:
    case stk::topology::QUAD_9:
    case stk::topology::QUAD_4_2D:
    case stk::topology::QUAD_9_2D:
        return 4.;
    case stk::topology::TET_4:
    case stk::topology::TET_10:
        return 1./6.;
    case stk::topology::HEX_8:
    case stk::topology::HEX_27:
        return 8.;
    default:
        throw std::runtime_error("Element topology not found in parametric_volume: " + topology.name());
        break;
    }
}

MasterElementIntrepid::MasterElementIntrepid(stk::topology topology)
: m_topology(topology)
{
  m_intrepidBasis = get_basis_for_topology(topology);
  using WeightType = double;
  shards::CellTopology cellType = stk::mesh::get_cell_topology(topology);
  auto intrepidCubature = Intrepid2::DefaultCubatureFactory().create<ExecutionSpace, PointType, WeightType>(cellType, 2*m_intrepidBasis->getDegree());

  m_numIntgPts = intrepidCubature->getNumPoints();
  m_numNodes = topology.num_nodes();
  m_numElemDims  = intrepidCubature->getDimension();

  // Allocate reference data
  m_shapeFuncs.resize(m_numIntgPts*m_numNodes);
  m_pointGrads.resize(m_numIntgPts*m_numNodes*m_numElemDims );
  m_refPoints.resize(m_numIntgPts*m_numElemDims);
  m_refWeights.resize(m_numIntgPts);
  m_refCoords.resize(m_numNodes*m_numElemDims);
  m_refVolume = parametric_volume_for_topology(topology);

  Kokkos::DynRankView<double, ExecutionSpace> refPoints(m_refPoints.data(), m_numIntgPts, m_numElemDims);
  Kokkos::DynRankView<double, ExecutionSpace> refWeights(m_refWeights.data(), m_numIntgPts);
  intrepidCubature->getCubature(refPoints, refWeights);

  Kokkos::DynRankView<double, ExecutionSpace> pointVals("pointVals", m_numNodes, m_numIntgPts);
  Kokkos::DynRankView<double, ExecutionSpace> pointGrads("pointGrads", m_numNodes, m_numIntgPts, m_numElemDims);

  m_intrepidBasis->getValues(pointVals, refPoints, Intrepid2::OPERATOR_VALUE);
  m_intrepidBasis->getValues(pointGrads, refPoints, Intrepid2::OPERATOR_GRAD);

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

  Kokkos::DynRankView<double, ExecutionSpace> nodeCoords("nodeCoords", m_numElemDims);
  for ( int node(0); node < m_numNodes; ++node ) {
    Intrepid2::CellTools<ExecutionSpace>::getReferenceNode( nodeCoords, cellType, node );
    for ( int dim(0); dim < m_numElemDims; ++dim ) {
      m_refCoords[node*m_numElemDims + dim] = nodeCoords[dim];
    }
  }
}

void
MasterElementIntrepid::determinant(
   const int numCoordDims,
   const int nelem,
    const double* coords,  // (nvec,npe,nelem)
    double* det_J,         // (nelem,nint)
    double* error ) const  // (nelem)
{
  determinant(numCoordDims, m_numIntgPts, m_numNodes, m_pointGrads.data(), nelem, coords, det_J, error);
}

void
MasterElementIntrepid::shape_fcn(
    const int nint,  // returns array(npe,nint)
    const double* p_coords,
    double* result) const
{
  Kokkos::DynRankView<double, ExecutionSpace> pointVals("pointVals", m_numNodes, nint);
  const Kokkos::DynRankView<double, ExecutionSpace> pcoordVec(const_cast<double*>(p_coords), nint, m_numElemDims);

  m_intrepidBasis->getValues(pointVals, pcoordVec, Intrepid2::OPERATOR_VALUE);

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
  Kokkos::DynRankView<double, ExecutionSpace> pointGrads("pointGrads", m_numNodes, nint, m_numElemDims);
  const Kokkos::DynRankView<double, ExecutionSpace> pcoordVec(const_cast<double*>(p_coords), nint, m_numElemDims);

  m_intrepidBasis->getValues(pointGrads, pcoordVec, Intrepid2::OPERATOR_GRAD);

  // re-order shape functions for consistency with other master elements
  for ( int ip(0); ip < nint; ++ip ) {
    for ( int node(0); node < m_numNodes; ++node ) {
      for ( int dim(0); dim < m_numElemDims; ++dim ) {
        result[(ip*m_numNodes + node)*m_numElemDims + dim] = pointGrads(node,ip,dim);
      }
    }
  }
}

void
MasterElementIntrepid::interpolate_point(
    const int  /*npar_coord*/,
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
      result[ comp ] += shape[node] * field[ ncomp_field * node + comp ];
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
      const int numCoordDims,
      const int  nint,
      const int npe_g,
      const double* deriv_g,  // (nvec,npe_g,nint)
      const int   nelem,
      const double* coords,   // (nvec,npe,nelem)
      double* det_J,          // (nelem,nint)
      double* error ) const   // (nelem)
{
  MasterElementCalc::determinant(m_numElemDims, numCoordDims, nint, npe_g, deriv_g, nelem, coords, det_J, error);
}

void
MasterElementIntrepid::gradient_operator(
        const int numCoordDims,
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
  STK_ThrowRequireMsg(m_numElemDims == numCoordDims, "MasterElementHybrid::gradient_operator does not support lower rank elements in higher dimensions (e.g. BAR,QUAD,TRI in 3D).");
  MasterElementCalc::gradient_operator(m_numElemDims, nint, npe_g, deriv_g, npe_f, deriv_f, nelem, coords, gradop, det_J, error);
}

} // namespace krino
