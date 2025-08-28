// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MasterElementHybrid.hpp>
#include <Akri_MasterElementBasis.hpp>
#include <Akri_MasterElement.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/MetaData.hpp> // for get_cell_topology
#include <Akri_MasterElementCalc.hpp>

#include <Kokkos_DynRankView.hpp>
#include <Intrepid2_Cubature.hpp>
#include <Intrepid2_DefaultCubatureFactory.hpp>

namespace krino {

static std::unique_ptr<Basis> get_basis_for_topology(stk::topology t)
{
  switch(t())
  {
  case stk::topology::LINE_2:
      return std::make_unique<Basis_LINE_2>();
  case stk::topology::LINE_3:
      return std::make_unique<Basis_LINE_3>();
  case stk::topology::TRI_3:
  case stk::topology::TRI_3_2D:
      return std::make_unique<Basis_TRI_3>();
  case stk::topology::TRI_6:
  case stk::topology::TRI_6_2D:
      return std::make_unique<Basis_TRI_6>();
  case stk::topology::QUAD_4:
  case stk::topology::QUAD_4_2D:
      return std::make_unique<Basis_QUAD_4>();
  case stk::topology::QUAD_9:
  case stk::topology::QUAD_9_2D:
      return std::make_unique<Basis_QUAD_9>();
  case stk::topology::TET_4:
      return std::make_unique<Basis_TET_4>();
  case stk::topology::TET_10:
      return std::make_unique<Basis_TET_10>();
  case stk::topology::HEX_8:
      return std::make_unique<Basis_HEX_8>();
  case stk::topology::HEX_27:
      return std::make_unique<Basis_HEX_27>();
  case stk::topology::WEDGE_6:
      return std::make_unique<Basis_WEDGE_6>();
  default:
      throw std::runtime_error("Element topology not found in get_basis_for_topology(): " + t.name());
  }
}

static stk::topology working_topology(const stk::topology topology)
{
  if (!topology.is_shell())
    return topology;
  switch(topology())
  {
  case stk::topology::SHELL_SIDE_BEAM_2:
    return stk::topology::LINE_2;
  case stk::topology::SHELL_SIDE_BEAM_3:
    return stk::topology::LINE_3;
  case stk::topology::SHELL_TRI_3:
  case stk::topology::SHELL_TRI_3_ALL_FACE_SIDES:
      return stk::topology::TRIANGLE_3;
  case stk::topology::SHELL_TRI_6:
  case stk::topology::SHELL_TRI_6_ALL_FACE_SIDES:
      return stk::topology::TRIANGLE_6;
  default:
      throw std::runtime_error("Shell topology not found in working_topology_for_processing_sides(): " + topology.name());
  }
}

MasterElementHybrid::MasterElementHybrid(stk::topology topology)
: m_topology(working_topology(topology))
{
  m_Basis = get_basis_for_topology(m_topology);
  using PointType = double;
  using WeightType = double;
  using ExecutionSpace = Kokkos::DefaultHostExecutionSpace;
  auto intrepid2Cubature = Intrepid2::DefaultCubatureFactory().create<ExecutionSpace, PointType, WeightType>(stk::mesh::get_cell_topology(m_topology), 2*m_Basis->degree());

  m_numIntgPts = intrepid2Cubature->getNumPoints();

  m_numNodes = m_topology.num_nodes();
  m_numElemDims  = intrepid2Cubature->getDimension();

  // Allocate reference data
  m_shapeFuncs.resize(m_numIntgPts*m_numNodes);
  m_pointGrads.resize(m_numIntgPts*m_numNodes*m_numElemDims );
  m_refPoints.resize(m_numIntgPts*m_numElemDims);
  m_refWeights.resize(m_numIntgPts);
  m_refCoords.resize(m_numNodes*m_numElemDims);
  m_refVolume = m_Basis->parametric_volume();

  Kokkos::DynRankView<double, ExecutionSpace> refPoints(m_refPoints.data(), m_numIntgPts, m_numElemDims);
  Kokkos::DynRankView<double, ExecutionSpace> refWeights(m_refWeights.data(), m_numIntgPts);
  intrepid2Cubature->getCubature(refPoints, refWeights);

  // compute the reference values and gradients at the integration points
  m_Basis->shape_fcn(m_numIntgPts, m_refPoints.data(), m_shapeFuncs.data());
  m_Basis->shape_fcn_deriv(m_numIntgPts, m_refPoints.data(), m_pointGrads.data());
  m_Basis->nodal_parametric_coordinates(m_refCoords.data());

  m_centroidParCoords.resize(m_numElemDims, 0.);
  for(int n=0; n < m_numNodes; ++n)
  {
    for(int d=0; d < m_numElemDims; ++d)
    {
      m_centroidParCoords[d] += m_refCoords[n*m_numElemDims + d];
    }
  }
}

MasterElementHybrid::~MasterElementHybrid()
{
}

void
MasterElementHybrid::determinant(
    const int numCoordDims,
    const int nelem,
    const double* coords,  // (numCoordDims,npe,nelem)
    double* det_J,         // (nelem,nint)
    double* error ) const  // (nelem)
{
  MasterElementCalc::determinant(m_numElemDims, numCoordDims, m_numIntgPts, m_numNodes, m_pointGrads.data(), nelem, coords, det_J, error);
}

void
MasterElementHybrid::shape_fcn(
    const int nint,  // returns array(npe,nint)
    const double* p_coords,
    double* result) const
{
  m_Basis->shape_fcn(nint, p_coords, result);
}

void
MasterElementHybrid::shape_fcn_deriv(
    const int nint,
    const double* p_coords,
    double* result ) const
{
  m_Basis->shape_fcn_deriv(nint, p_coords, result);
}


void
MasterElementHybrid::interpolate_point(
    const int  /*npar_coord*/,
    const double * par_coord,      // (npar_coord)
    const int  ncomp_field,
    const double * field,          // (ncomp_field,num_nodes)
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
MasterElementHybrid::scalar_gradient(
    const int   nelem,      //:  number of elements to process
    const double* gradop,     //: (nvec,npe,nelem,nint)
    const double* det_J,      //: (nelem,nint)
    const double* sfield,     //: (npe,nelem)
    double* vector ) const    //: (nvec,nelem,nint)
{
  MasterElementCalc::scalar_gradient(m_numIntgPts, nelem, m_numElemDims, m_numNodes, gradop, det_J, sfield, vector);
}

void
MasterElementHybrid::scalar_gradient(
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
MasterElementHybrid::determinant(
      const int numCoordDims,
      const int  nint,
      const int npe_g,
      const double* deriv_g,  // (m_numElemDims,npe_g,nint)
      const int   nelem,
      const double* coords,   // (numCoordDims,npe,nelem)
      double* det_J,          // (nelem,nint)
      double* error ) const   // (nelem)
{
  MasterElementCalc::determinant(m_numElemDims, numCoordDims, nint, npe_g, deriv_g, nelem, coords, det_J, error);
}

void
MasterElementHybrid::gradient_operator(
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
