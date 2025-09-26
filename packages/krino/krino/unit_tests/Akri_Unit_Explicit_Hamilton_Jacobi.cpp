// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <stk_math/StkVector.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/Relation.hpp>

namespace krino {

struct ProblemFields
{
  stk::mesh::Field<double> * levelSetField = nullptr;
  stk::mesh::Field<double> * coordsField = nullptr;
  stk::mesh::Field<double> * RHS = nullptr;
  stk::mesh::Field<double> * RHSNorm = nullptr;
  stk::mesh::Field<double> * speedField = nullptr;
};

void associate_input_mesh(stk::io::StkMeshIoBroker & stkIo, const std::string & meshName)
{
  stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RIB"));
  stkIo.add_mesh_database(meshName, stk::io::READ_MESH);
  stkIo.create_input_mesh();
}

stk::mesh::BulkData & read_mesh(stk::io::StkMeshIoBroker & stkIo)
{
  stkIo.populate_bulk_data();
  return stkIo.bulk_data();
}

void declare_fields(stk::mesh::MetaData & meta, ProblemFields & fields)
{
  fields.levelSetField = &meta.declare_field<double>(stk::topology::NODE_RANK, "LevelSet", 2);
  stk::mesh::put_field_on_mesh(*fields.levelSetField, meta.universal_part(), nullptr);
  fields.RHS = &meta.declare_field<double>(stk::topology::NODE_RANK, "RHS", 1);
  stk::mesh::put_field_on_mesh(*fields.RHS, meta.universal_part(), nullptr);
  fields.RHSNorm = &meta.declare_field<double>(stk::topology::NODE_RANK, "RHSNorm", 1);
  stk::mesh::put_field_on_mesh(*fields.RHSNorm, meta.universal_part(), nullptr);
  auto constCoordsField = static_cast<const stk::mesh::Field<double>*>(meta.coordinate_field());
  fields.coordsField = const_cast<stk::mesh::Field<double>*>(constCoordsField);

  if (true)
  {
    fields.speedField = &meta.declare_field<double>(stk::topology::ELEMENT_RANK, "Speed", 1);
    stk::mesh::put_field_on_mesh(*fields.speedField, meta.universal_part(), nullptr);
  }
}

const int TetFaceTable[4][3] = { {0, 1, 2},
        {1, 2, 3},
        {0, 2, 3},
        {0, 1, 3} };

const int TetEdgeNodeOrder[6][2] = // [edge][edge_node]
  { {0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3} };

void initialize_level_set(const stk::mesh::BulkData & mesh, const ProblemFields & fields, std::function<double(const double *)> initial_level_set, const double multiplier = 1.0)
{
  for (auto & bptr : mesh.buckets(stk::topology::NODE_RANK))
  {
    for (auto & node : *bptr)
    {
      const double * x = stk::mesh::field_data(*fields.coordsField, node);
      double* LS = stk::mesh::field_data(fields.levelSetField->field_of_state(stk::mesh::StateNP1), node);
      *LS = initial_level_set(x) * multiplier;
    }
  }
}

double initialize_constant_speed(const stk::mesh::BulkData & /*mesh*/, const ProblemFields & fields, const double speed)
{
  stk::mesh::field_fill(speed, *fields.speedField);
  return speed;
}

double compute_tet_vol(const stk::mesh::Entity * elemNodes, const stk::mesh::Field<double> & coordsField)
{
  std::array<stk::math::Vector3d, 4> xNode;
  for (int n=0; n<4; ++n) xNode[n] = stk::math::Vector3d(stk::mesh::field_data(coordsField, elemNodes[n]));

  return Dot(xNode[1]-xNode[0],Cross(xNode[2]-xNode[0],xNode[3]-xNode[0]))/6.;
}

double compute_tri_vol(const stk::mesh::Entity * elemNodes, const stk::mesh::Field<double> & coordsField)
{
  std::array<stk::math::Vector3d, 3> xNode;
  for (int n=0; n<3; ++n) xNode[n] = stk::math::Vector3d(stk::mesh::field_data(coordsField, elemNodes[n]),2);

  return 0.5*Cross(xNode[1]-xNode[0],xNode[2]-xNode[0]).length();
}

double compute_vol(const stk::mesh::Entity * elemNodes, const unsigned numElemNodes, const stk::mesh::Field<double> & coordsField)
{
  if (numElemNodes == 3)
    return compute_tri_vol(elemNodes, coordsField);
  else
    return compute_tet_vol(elemNodes, coordsField);
}

void compute_tet_gradOP_and_vol(const stk::mesh::Entity * elemNodes, const stk::mesh::Field<double> & coordsField, std::vector<stk::math::Vector3d> & gradOP, double & vol)
{
  std::array<stk::math::Vector3d, 4> xNode;
  for (int n=0; n<4; ++n) xNode[n] = stk::math::Vector3d(stk::mesh::field_data(coordsField, elemNodes[n]));

  std::array<stk::math::Vector3d, 4> xFace;
  for (int f=0; f<4; ++f) xFace[f] = (xNode[TetFaceTable[f][0]]+xNode[TetFaceTable[f][1]]+xNode[TetFaceTable[f][2]])/3.;

  std::array<stk::math::Vector3d, 6> xEdge;
  for (int e=0; e<6; ++e) xEdge[e] = 0.5*(xNode[TetEdgeNodeOrder[e][0]]+xNode[TetEdgeNodeOrder[e][1]]);

  vol = Dot(xNode[1]-xNode[0],Cross(xNode[2]-xNode[0],xNode[3]-xNode[0]))/6.;
  const double norm = 0.5/vol;

  gradOP[0] = (Cross(xEdge[3]-xEdge[0],xNode[0]-xFace[3])+
                         Cross(xEdge[2]-xEdge[3],xNode[0]-xFace[2])+
                         Cross(xEdge[0]-xEdge[2],xNode[0]-xFace[0]))*norm;
  gradOP[1] = (Cross(xEdge[0]-xEdge[4],xNode[1]-xFace[3])+
                         Cross(xEdge[1]-xEdge[0],xNode[1]-xFace[0])+
                         Cross(xEdge[4]-xEdge[1],xNode[1]-xFace[1]))*norm;
  gradOP[2] = (Cross(xEdge[2]-xEdge[1],xNode[2]-xFace[0])+
                         Cross(xEdge[5]-xEdge[2],xNode[2]-xFace[2])+
                         Cross(xEdge[1]-xEdge[5],xNode[2]-xFace[1]))*norm;
  gradOP[3] = (Cross(xEdge[5]-xEdge[4],xNode[3]-xFace[1])+
                         Cross(xEdge[3]-xEdge[5],xNode[3]-xFace[2])+
                         Cross(xEdge[4]-xEdge[3],xNode[3]-xFace[3]))*norm;
}

void compute_tri_gradOP_and_vol(const stk::mesh::Entity * elemNodes, const stk::mesh::Field<double> & coordsField, std::vector<stk::math::Vector3d> & gradOP, double & vol)
{
  std::array<stk::math::Vector3d, 3> xNode;
  for (int n=0; n<3; ++n) xNode[n] = stk::math::Vector3d(stk::mesh::field_data(coordsField, elemNodes[n]),2);

  vol = 0.5*Cross(xNode[1]-xNode[0],xNode[2]-xNode[0]).length();
  const double norm = 0.5/vol;

  gradOP[0] = (crossZ(xNode[0]-xNode[2])+crossZ(xNode[1]-xNode[0]))*norm;
  gradOP[1] = (crossZ(xNode[1]-xNode[0])+crossZ(xNode[2]-xNode[1]))*norm;
  gradOP[2] = (crossZ(xNode[2]-xNode[1])+crossZ(xNode[0]-xNode[2]))*norm;
}

void compute_gradOP_and_vol(const stk::mesh::Entity * elemNodes, const stk::mesh::Field<double> & coordsField, std::vector<stk::math::Vector3d> & gradOP, double &vol)
{
  if (gradOP.size() == 3)
    compute_tri_gradOP_and_vol(elemNodes, coordsField, gradOP, vol);
  else
    compute_tet_gradOP_and_vol(elemNodes, coordsField, gradOP, vol);
}

stk::math::Vector3d compute_scalar_gradient(const std::vector<stk::math::Vector3d> & nodalAreaVectors, const stk::mesh::Entity * elemNodes, const stk::mesh::Field<double> & levelSetField)
{
  stk::math::Vector3d grad(stk::math::Vector3d::ZERO);
  for (unsigned i=0; i<nodalAreaVectors.size(); ++i)
  {
    const double ls = *stk::mesh::field_data(levelSetField, elemNodes[i]);
    grad += ls * nodalAreaVectors[i];
  }
  return grad;
}

double mesh_minimum_length_scale_using_minimum_volume(const stk::mesh::BulkData & mesh, const stk::mesh::Field<double> & coordsField)
{
  const unsigned nodesPerElem = mesh.mesh_meta_data().spatial_dimension() + 1;

  double minVol = std::numeric_limits<double>::max();
  for (auto & bptr : mesh.buckets(stk::topology::ELEMENT_RANK))
  {
    for (auto & elem : *bptr)
    {
      const stk::mesh::Entity * elemNodes = mesh.begin_nodes(elem);
      minVol = std::min(minVol, compute_vol(elemNodes, nodesPerElem, coordsField));
    }
  }

  return std::pow(minVol, 1./mesh.mesh_meta_data().spatial_dimension());
}

double mesh_minimum_length_scale(const stk::mesh::BulkData & mesh, const stk::mesh::Field<double> & coordsField)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::math::Vector3d> gradOP(dim+1);

  double vol = 0.;
  double maxGradLength = 0.;
  for (auto & bptr : mesh.buckets(stk::topology::ELEMENT_RANK))
  {
    for (auto & elem : *bptr)
    {
      const stk::mesh::Entity * elemNodes = mesh.begin_nodes(elem);
      compute_gradOP_and_vol(elemNodes, coordsField, gradOP, vol);
      for (unsigned n=0; n<dim+1; ++n)
      {
        maxGradLength = std::max(maxGradLength, gradOP[n].length());
      }
    }
  }

  return 1./maxGradLength;
}

void write_fields(stk::io::StkMeshIoBroker &outStkIo, const size_t outputFileIndex, const double time)
{
  outStkIo.begin_output_step(outputFileIndex, time);
  outStkIo.write_defined_output_fields(outputFileIndex);
  outStkIo.end_output_step(outputFileIndex);
}

void assemble_residual_for_element_speed(const stk::mesh::BulkData & mesh, const ProblemFields & fields, const double /*ep*/)
{
  // Implementation of Barth-Sethian positive coefficient scheme for element speed fields.

  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::math::Vector3d> gradOP(dim+1);

  stk::mesh::field_fill(0.0, *fields.RHS);
  stk::mesh::field_fill(0.0, *fields.RHSNorm);

  for (auto & bptr : mesh.buckets(stk::topology::ELEMENT_RANK))
  {
    for (auto & elem : *bptr)
    {
      const unsigned numElemNodes = mesh.num_nodes(elem);
      const stk::mesh::Entity * elemNodes = mesh.begin_nodes(elem);
      double vol = 0.;
      compute_gradOP_and_vol(elemNodes, *fields.coordsField, gradOP, vol);
      const stk::math::Vector3d normalDir = compute_scalar_gradient(gradOP, elemNodes, fields.levelSetField->field_of_state(stk::mesh::StateN)).unit_vector();

      const double elementSpeed = *stk::mesh::field_data(*fields.speedField, elem);

      std::vector<double> volHamiltonianCoeffs(dim+1); // K_i in Barth-Sethian
      double sumNegCoeffs = 0.;
      double sumPosCoeffs = 0.;
      double sumNegContrib = 0.;
      double volHamiltonian = 0.;
      for (unsigned n=0; n<numElemNodes; ++n)
      {
        stk::mesh::Entity node = elemNodes[n];
        const double LSOld = *stk::mesh::field_data(fields.levelSetField->field_of_state(stk::mesh::StateN), node);
        volHamiltonianCoeffs[n] = elementSpeed*vol*Dot(normalDir, gradOP[n]);
        volHamiltonian += volHamiltonianCoeffs[n] * LSOld;

        if (volHamiltonianCoeffs[n] < 0.)
        {
          sumNegCoeffs += volHamiltonianCoeffs[n];
          sumNegContrib += volHamiltonianCoeffs[n] * LSOld;
        }
        else
        {
          sumPosCoeffs += volHamiltonianCoeffs[n];
        }
      }

      std::vector<double> alpha(dim+1, 0.); // delta phi_i in Barth-Sethian
      double sumPosAlpha = 0.;
      for (unsigned n=0; n<numElemNodes; ++n)
      {
        stk::mesh::Entity node = elemNodes[n];
        const double LSOld = *stk::mesh::field_data(fields.levelSetField->field_of_state(stk::mesh::StateN), node);
        if (volHamiltonianCoeffs[n] > 0.)
        {
          alpha[n] = volHamiltonianCoeffs[n]/sumPosCoeffs*(sumNegContrib-sumNegCoeffs*LSOld)/volHamiltonian;
          if (alpha[n] > 0.) sumPosAlpha += alpha[n];
        }
      }

      for (unsigned n=0; n<numElemNodes; ++n)
      {
        stk::mesh::Entity node = elemNodes[n];
        double & residual = *stk::mesh::field_data(*fields.RHS, node);
        double & residualNorm = *stk::mesh::field_data(*fields.RHSNorm, node);

        if (alpha[n] > 0.)
        {
          const double wt = alpha[n]/sumPosAlpha;
          residual += wt * volHamiltonian;
          residualNorm += wt*vol;
        }
      }
    }
  }
}

void assemble_residual_for_nodal_speed(const stk::mesh::BulkData & mesh, const ProblemFields & fields, const double /*eps*/)
{
  // Not extensively tested.  This is an adaptation of Barth-Sethian positive coefficient scheme for nodal speed fields

  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::math::Vector3d> gradOP(dim+1);

  stk::mesh::field_fill(0.0, *fields.RHS);
  stk::mesh::field_fill(0.0, *fields.RHSNorm);

  for (auto & bptr : mesh.buckets(stk::topology::ELEMENT_RANK))
  {
    for (auto & elem : *bptr)
    {
      const unsigned numElemNodes = mesh.num_nodes(elem);
      const stk::mesh::Entity * elemNodes = mesh.begin_nodes(elem);
      double vol = 0.;
      compute_gradOP_and_vol(elemNodes, *fields.coordsField, gradOP, vol);
      const stk::math::Vector3d normalDir = compute_scalar_gradient(gradOP, elemNodes, fields.levelSetField->field_of_state(stk::mesh::StateN)).unit_vector();

      for (unsigned n=0; n<numElemNodes; ++n)
      {
        stk::mesh::Entity node = elemNodes[n];
        const double LSOld = *stk::mesh::field_data(fields.levelSetField->field_of_state(stk::mesh::StateN), elemNodes[n]);

        //const double nodalSpeed = *stk::mesh::field_data(fields.speed, node);
        const double nodalSpeed = 1.0;

        const double volHamiltonianCoeffs_n = vol*nodalSpeed*Dot(normalDir, gradOP[n]); // K_n in Barth-Sethian, probably should use nodal volume not elem volume

        if (volHamiltonianCoeffs_n > 0.)
        {
          double sumNegCoeffs = 0.;
          double sumPosCoeffs = 0.;
          double sumNegContrib = 0.;
          double volHamiltonian = 0.;
          for (unsigned j=0; j<numElemNodes; ++j)
          {
            const double LSOldj = *stk::mesh::field_data(fields.levelSetField->field_of_state(stk::mesh::StateN), elemNodes[j]);
            const double volHamiltonianCoeffs_j = vol*nodalSpeed*Dot(normalDir, gradOP[j]);
            volHamiltonian += volHamiltonianCoeffs_j * LSOldj;

            if (volHamiltonianCoeffs_j < 0.)
            {
              sumNegCoeffs += volHamiltonianCoeffs_j;
              sumNegContrib += volHamiltonianCoeffs_j * LSOldj;
            }
            else
            {
              sumPosCoeffs += volHamiltonianCoeffs_j;
            }
          }

          const double wt = volHamiltonianCoeffs_n/sumPosCoeffs*(sumNegContrib-sumNegCoeffs*LSOld)/volHamiltonian;
          if (wt > 0.)
          {
            double & residual = *stk::mesh::field_data(*fields.RHS, node);
            double & residualNorm = *stk::mesh::field_data(*fields.RHSNorm, node);
            residual += wt * volHamiltonian;
            residualNorm += wt*vol;
          }
        }
      }
    }
  }
}

double assemble_and_update_Eikonal(const stk::mesh::BulkData & mesh, const ProblemFields & fields, const double eps, const double dt, const bool computeArrivalTime)
{
  // This is a hybrid between the Barth-Sethian positive coefficient scheme and the Morgan-Waltz scheme for reinitialization.
  // Uses element based speed and nodal sign function to assemble nodal contributions for Hamiltonian.
  // The assembled nodal Hamiltonian is then used with nodal source term to explicitly update signed distance
  // (or arrival time for non-unit speed).
  // Unlike the elemental algorithm developed by Barth-Sethian, this algorithm converges to the exact solution
  // for the "Distance Function Test" described in Morgan-Waltz.  Unlike the Morgan-Waltz algorithm, this
  // form converges much faster and is tolerant of meshes with obtuse angles because it uses the positive coefficient
  // form in Barth-Sethian.
  assert(!computeArrivalTime || nullptr != fields.speedField);

  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::math::Vector3d> gradOP(dim+1);

  stk::mesh::field_fill(0.0, *fields.RHS);
  stk::mesh::field_fill(0.0, *fields.RHSNorm);

  for (auto & bptr : mesh.buckets(stk::topology::ELEMENT_RANK))
  {
    for (auto & elem : *bptr)
    {
      const unsigned numElemNodes = mesh.num_nodes(elem);
      const stk::mesh::Entity * elemNodes = mesh.begin_nodes(elem);
      double vol = 0.;
      compute_gradOP_and_vol(elemNodes, *fields.coordsField, gradOP, vol);
      const stk::math::Vector3d normalDir = compute_scalar_gradient(gradOP, elemNodes, fields.levelSetField->field_of_state(stk::mesh::StateN)).unit_vector();

      double elementSpeed = 1.0;
      if (computeArrivalTime)
      {
        elementSpeed = *stk::mesh::field_data(*fields.speedField, elem);
      }

      for (unsigned n=0; n<numElemNodes; ++n)
      {
        stk::mesh::Entity node = elemNodes[n];
        const double LSOld = *stk::mesh::field_data(fields.levelSetField->field_of_state(stk::mesh::StateN), elemNodes[n]);

        double nodalSpeed = 0.;

        const double sign = LSOld/sqrt(LSOld*LSOld + eps*eps);
        nodalSpeed = sign*elementSpeed;

        const double volHamiltonianCoeffs_n = vol*nodalSpeed*Dot(normalDir, gradOP[n]); // K_n in Barth-Sethian, probably should use nodal volume not elem volume

        if (volHamiltonianCoeffs_n > 0.)
        {
          double sumNegCoeffs = 0.;
          double sumPosCoeffs = 0.;
          double sumNegContrib = 0.;
          double volHamiltonian = 0.;
          for (unsigned j=0; j<numElemNodes; ++j)
          {
            const double LSOldj = *stk::mesh::field_data(fields.levelSetField->field_of_state(stk::mesh::StateN), elemNodes[j]);
            const double volHamiltonianCoeffs_j = vol*nodalSpeed*Dot(normalDir, gradOP[j]);
            volHamiltonian += volHamiltonianCoeffs_j * LSOldj;

            if (volHamiltonianCoeffs_j < 0.)
            {
              sumNegCoeffs += volHamiltonianCoeffs_j;
              sumNegContrib += volHamiltonianCoeffs_j * LSOldj;
            }
            else
            {
              sumPosCoeffs += volHamiltonianCoeffs_j;
            }
          }

          const double wt = volHamiltonianCoeffs_n/sumPosCoeffs*(sumNegContrib-sumNegCoeffs*LSOld)/volHamiltonian;
          if (wt > 0.)
          {
            double & residual = *stk::mesh::field_data(*fields.RHS, node);
            double & residualNorm = *stk::mesh::field_data(*fields.RHSNorm, node);
            residual += wt * volHamiltonian;
            residualNorm += wt*vol;
          }
        }
      }
    }
  }

  double sumSqrResid = 0.;
  size_t sumCount = 0;

  for (auto & bptr : mesh.buckets(stk::topology::NODE_RANK))
  {
    for (auto & node : *bptr)
    {
      const double residField = *stk::mesh::field_data(*fields.RHS, node);
      const double residNormField = *stk::mesh::field_data(*fields.RHSNorm, node);

      const double Hamiltonian = (residNormField > 0.) ? (residField/residNormField) : 0.;

      const double LSOld = *stk::mesh::field_data(fields.levelSetField->field_of_state(stk::mesh::StateN), node);
      double & LS = *stk::mesh::field_data(fields.levelSetField->field_of_state(stk::mesh::StateNP1), node);
      const double sign = LSOld/sqrt(LSOld*LSOld + eps*eps);

      LS = LSOld - dt * (Hamiltonian - sign);

      if (computeArrivalTime || std::abs(LSOld) < eps)
      {
        sumSqrResid += (Hamiltonian - sign)*(Hamiltonian - sign);
        sumCount++;
      }
    }
  }
  return std::sqrt(sumSqrResid/sumCount);
}


double apply_level_set_update(const stk::mesh::BulkData & mesh, const ProblemFields & fields, const double eps, const double dt)
{
  double sumResid = 0.;
  size_t sumCount = 0;

  for (auto & bptr : mesh.buckets(stk::topology::NODE_RANK))
  {
    for (auto & node : *bptr)
    {
      const double residField = *stk::mesh::field_data(*fields.RHS, node);
      const double residNormField = *stk::mesh::field_data(*fields.RHSNorm, node);

      const double resid = (residNormField > 0.) ? (residField/residNormField) : 0.;

      const double LSOld = *stk::mesh::field_data(fields.levelSetField->field_of_state(stk::mesh::StateN), node);
      double & LS = *stk::mesh::field_data(fields.levelSetField->field_of_state(stk::mesh::StateNP1), node);

      LS = LSOld - dt * resid;

      if (std::abs(LSOld) < eps)
      {
        sumResid += std::abs(resid);
        sumCount++;
      }
    }
  }
  return sumResid/sumCount;
}

bool domain_contains_interface(const stk::mesh::BulkData & mesh, const stk::mesh::Field<double> & levelSetField)
{
  bool hasNeg = false;
  bool hasPos = false;

  for (auto & bptr : mesh.buckets(stk::topology::NODE_RANK))
  {
    for (auto & node : *bptr)
    {
      const double LS = *stk::mesh::field_data(levelSetField, node);

      if (LS < 0.) hasNeg = true;
      if (LS > 0.) hasPos = true;

      if (hasNeg && hasPos) return true;
    }
  }
  return false;
}

void evolve_level_set(const stk::mesh::BulkData & mesh, const ProblemFields & fields, const double eps, const double dt)
{
  assemble_residual_for_element_speed(mesh, fields, eps);
  apply_level_set_update(mesh, fields, eps, dt);
}

void reinitialize_level_set(
    stk::mesh::BulkData & mesh,
    const ProblemFields & fields,
    const double eps,
    const double dtau,
    stk::io::StkMeshIoBroker * stkIo = nullptr,
    const double outputStartTime = 0.,
    const double outputStopTime = 0.,
    const size_t outputFileIndex = 0)
{
  if (!domain_contains_interface(mesh, fields.levelSetField->field_of_state(stk::mesh::StateNP1)))
  {
    return;
  }

  const double convergedTol = 0.01;
  bool converged = false;
  const int maxIters = 1000;
  const int printFreq = 50;
  const double dOutputTime = (outputStopTime-outputStartTime)/(maxIters+1);
  for (int iter = 0; iter<maxIters; ++iter)
  {
    mesh.update_field_data_states(); // Is there a way to do this just for levelSet?

    const double averageNodalResidual = assemble_and_update_Eikonal(mesh, fields, eps, dtau, false);

    if ((iter+1) % printFreq == 0)
       std::cout << "After " << iter+1 << " iterations of reinitialization, the relative error is " << averageNodalResidual << std::endl;

    if (stkIo != nullptr) write_fields(*stkIo, outputFileIndex, outputStartTime+(iter+1)*dOutputTime);

    if (averageNodalResidual < convergedTol)
    {
      converged = true;
      std::cout << "Reinitialization converged after " << iter + 1
                << " iterations  with relative error of " << averageNodalResidual << std::endl;
      break;
    }
  }
  if (!converged)
  {
    std::cout << "Reinitialization failed to converge after " << maxIters << " iterations."
              << std::endl;
  }
}

void compute_arrival_time(
    stk::mesh::BulkData & mesh,
    const ProblemFields & fields,
    const double eps,
    const double dtau,
    const double convergedTol = 0.01,
    stk::io::StkMeshIoBroker * stkIo = nullptr,
    const size_t outputFileIndex = 0)
{
  if (!domain_contains_interface(mesh, fields.levelSetField->field_of_state(stk::mesh::StateNP1)))
  {
    std::cout << "Error, input level set field does not contain zero level set for initializing arrival time calculation." << std::endl;
    return;
  }

  bool converged = false;
  const int maxIters = 5000;
  const int printFreq = 50;
  for (int iter = 0; iter<maxIters; ++iter)
  {
    mesh.update_field_data_states(); // Is there a way to do this just for levelSet?

    const double averageNodalResidual =  assemble_and_update_Eikonal(mesh, fields, eps, dtau, true);

    if ((iter+1) % printFreq == 0)
       std::cout << "After " << iter+1 << " iterations in compute_arrival_time, the relative error is " << averageNodalResidual << std::endl;

    if (stkIo != nullptr) write_fields(*stkIo, outputFileIndex, (iter+1)*dtau);

    if (averageNodalResidual < convergedTol)
    {
      converged = true;
      std::cout << "Compute_arrival_time converged after " << iter + 1
                << " iterations  with relative error of " << averageNodalResidual << std::endl;
      break;
    }
  }
  if (!converged)
  {
    std::cout << "Compute_arrival_time failed to converge after " << maxIters << " iterations."
              << std::endl;
  }
}

size_t create_mesh(const std::string& filename, stk::io::StkMeshIoBroker &outStkIo)
{
  size_t outputFileIndex = outStkIo.create_output_mesh(filename, stk::io::WRITE_RESULTS);

  const stk::mesh::FieldVector fields = outStkIo.bulk_data().mesh_meta_data().get_fields();
  for(stk::mesh::FieldBase* field : fields)
  {
      const Ioss::Field::RoleType* fieldRole = stk::io::get_field_role(*field);
      if(fieldRole == nullptr || *fieldRole == Ioss::Field::TRANSIENT)
          outStkIo.add_field(outputFileIndex, *field);
  }

  outStkIo.write_output_mesh(outputFileIndex);
  return outputFileIndex;
}

double Heaviside(const double signedDist, const double eps)
{
  if (signedDist < -eps) return 0.;
  else if (signedDist > eps) return 1.;

  static const double pi = std::atan(1)*4;
  return 0.5*(1+ signedDist/eps + std::sin(pi*signedDist/eps)/pi);
}

double delta(const double signedDist, const double eps)
{
  if (signedDist < -eps || signedDist > eps) return 0.;

  static const double pi = std::atan(1)*4;
  return 0.5/eps * (1. + std::cos(pi*signedDist/eps));
}

double compute_level_set_volume(stk::mesh::BulkData & mesh, const ProblemFields & fields, const double eps)
{
  double vol = 0.;
  for (auto & bptr : mesh.buckets(stk::topology::ELEMENT_RANK))
  {
    for (auto & elem : *bptr)
    {
      const unsigned numElemNodes = mesh.num_nodes(elem);
      const stk::mesh::Entity * elemNodes = mesh.begin_nodes(elem);

      double avgLS = 0.;
      for (unsigned n=0; n<numElemNodes; ++n)
      {
        const double  LS = *stk::mesh::field_data(fields.levelSetField->field_of_state(stk::mesh::StateNP1), elemNodes[n]);
        avgLS += LS/numElemNodes;
      }

      vol += Heaviside(avgLS, eps) * compute_vol(elemNodes, numElemNodes, *fields.coordsField);
    }
  }
  return vol;
}

double compute_unit_radius_error_norm(stk::mesh::BulkData & mesh, const ProblemFields & fields)
{
  double norm = 0.;
  unsigned normCount = 0;
  for (auto & bptr : mesh.buckets(stk::topology::NODE_RANK))
  {
    for (auto & node : *bptr)
    {
      double & LS = *stk::mesh::field_data(fields.levelSetField->field_of_state(stk::mesh::StateNP1), node);
      const stk::math::Vector3d x(stk::mesh::field_data(*fields.coordsField, node), mesh.mesh_meta_data().spatial_dimension());

      const double error = LS - (x.length() - 1.0);
      norm += std::abs(error);
      ++normCount;
    }
  }
  return norm/normCount;
}

double poor_initial_condition_for_unit_circle(const double * x)
{
  return ((x[0]-1.)*(x[0]-1.)+(x[1]-1.)*(x[1]-1.)+0.1)*(std::sqrt(x[0]*x[0]+x[1]*x[1])-1.);
}

double poor_initial_condition_for_unit_sphere(const double * x)
{
  return ((x[0]-1.)*(x[0]-1.)+(x[1]-1.)*(x[1]-1.)+(x[2]-1.)*(x[2]-1.)+0.1)*(std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])-1.);
}

double flower_2D(const double * x)
{
  const int Nlobes = 3;
  const double r = std::sqrt(x[0]*x[0]+x[1]*x[1]);
  const double theta = std::atan2(x[1],x[0]);
  const double rSurf = 0.2 + 0.1*std::sin(Nlobes*theta);
  return r-rSurf;
}

TEST(HamiltonJacobi, 2DPoorInitialCondition_ComputingArrivalTimeProducesLowErrorEverywhere)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  if (parallel_size > 1) return;

  stk::io::StkMeshIoBroker stkIo(pm);
  ProblemFields fields;

  associate_input_mesh(stkIo, "square.g");
  declare_fields(stkIo.meta_data(),fields);
  stk::mesh::BulkData & mesh = read_mesh(stkIo);

  initialize_level_set(mesh, fields, poor_initial_condition_for_unit_circle);
  initialize_constant_speed(mesh, fields, 1.0);

  const double dx = mesh_minimum_length_scale(mesh, *fields.coordsField);
  const double eps = 1.5*dx; // Should have same units as level set
  const double dtau = 0.2*dx; // Reinitialization time step, based on unit advection speed used in reinitialization

  compute_arrival_time(mesh, fields, eps, dtau);

  const double errorNorm = compute_unit_radius_error_norm(mesh, fields);
  std::cout << "Error norm " << errorNorm << std::endl;
  EXPECT_LT(errorNorm, 0.001);
}

TEST(HamiltonJacobi, 3DPoorInitialCondition_ComputingArrivalTimeProducesLowErrorEverywhere)
{
#ifdef NDEBUG
#else
  return; // Optimized only due to length
#endif
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  if (parallel_size > 1) return;

  stk::io::StkMeshIoBroker stkIo(pm);
  ProblemFields fields;

  associate_input_mesh(stkIo, "cube_coarse.g");
  declare_fields(stkIo.meta_data(),fields);
  stk::mesh::BulkData & mesh = read_mesh(stkIo);

  initialize_level_set(mesh, fields, poor_initial_condition_for_unit_sphere);
  initialize_constant_speed(mesh, fields, 1.0);

  const double dx = mesh_minimum_length_scale(mesh, *fields.coordsField);
  const double eps = 1.5*dx; // Should have same units as level set
  const double dtau = 0.2*dx; // Reinitialization time step, based on unit advection speed used in reinitialization

  compute_arrival_time(mesh, fields, eps, dtau);

  const double errorNorm = compute_unit_radius_error_norm(mesh, fields);
  std::cout << "Error norm " << errorNorm << std::endl;
  EXPECT_LT(errorNorm, 0.01);
}

void test_circle_with_flower(const double speed)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  if (parallel_size > 1) return;

  stk::io::StkMeshIoBroker stkIo(pm);
  ProblemFields fields;

  associate_input_mesh(stkIo, "circle.g");
  declare_fields(stkIo.meta_data(),fields);
  stk::mesh::BulkData & mesh = read_mesh(stkIo);

  const double maxSpeed = initialize_constant_speed(mesh, fields, speed);
  initialize_level_set(mesh, fields, flower_2D, 1./maxSpeed);

  const double dx = mesh_minimum_length_scale(mesh, *fields.coordsField);
  const double eps = 1.5*dx / maxSpeed; // Should have same units as level set
  const double dtau = 0.2*dx / maxSpeed; // Units of time

  compute_arrival_time(mesh, fields, eps, dtau);
}

TEST(HamiltonJacobi, CircleWithFlowerIC_ArrivalTimeConvergesForAnySpeed)
{
  // Probably would be better to test that these converge in exactly the same number of steps
  test_circle_with_flower(1.0);
  test_circle_with_flower(10.0);
}

TEST(HamiltonJacobi, CircleWithFlowerIC_ReinitializationThenEvolveRuns)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  if (parallel_size > 1) return;

  stk::io::StkMeshIoBroker stkIo(pm);
  ProblemFields fields;

  associate_input_mesh(stkIo, "circle.g");
  declare_fields(stkIo.meta_data(), fields);
  stk::mesh::BulkData & mesh = read_mesh(stkIo);

  std::string outputFileName = "circle_flower.e";
  auto outputFileIndex = create_mesh(outputFileName, stkIo);

  const double maxSpeed = initialize_constant_speed(mesh, fields, 5.0);
  initialize_level_set(mesh, fields, flower_2D);

  const double Courant = 0.25;
  const double Ttotal = 0.1;
  const double dx = mesh_minimum_length_scale(mesh, *fields.coordsField);
  const unsigned Nt = (Ttotal+0.5*(Courant*dx/maxSpeed))/(Courant*dx/maxSpeed);
  const double dt = Ttotal/Nt;
  const double eps = 1.5*dx;
  const double dtau = 0.2*dx; // Reinitialization time step, based on unit advection speed used in reinitialization

  write_fields(stkIo, outputFileIndex, 0.0);
  reinitialize_level_set(mesh, fields, eps, dtau, &stkIo, 0.0, dt, outputFileIndex);

  std::cout << "Evolving to time " << Ttotal << " with " << Nt << " steps with Courant number " << Courant << std::endl;

  double time = 0.0;
  for (unsigned n=0; n<Nt; ++n)
  {
    mesh.update_field_data_states();
    evolve_level_set(mesh, fields, eps, dt);
    time += dt;
    write_fields(stkIo, outputFileIndex, time);
    reinitialize_level_set(mesh, fields, eps, dtau, &stkIo, time, time+dt, outputFileIndex);

    std::cout << "Time, Level set vol = " << time << " " << compute_level_set_volume(mesh, fields, eps) << std::endl;
  }

  const double finalVol = compute_level_set_volume(mesh, fields, eps);
  EXPECT_EQ(0., finalVol);

  std::remove(outputFileName.c_str()); // Remove output file unless debugging
}

TEST(HamiltonJacobi, CylinderWithFlowerIC_ReinitializationThenEvolveRuns)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  if (parallel_size > 1) return;

  stk::io::StkMeshIoBroker stkIo(pm);
  ProblemFields fields;

  associate_input_mesh(stkIo, "cylinder_coarse.g");
  declare_fields(stkIo.meta_data(), fields);
  stk::mesh::BulkData & mesh = read_mesh(stkIo);

  std::string outputFileName = "cylinder_flower.e";
  auto outputFileIndex = create_mesh(outputFileName, stkIo);

  const double maxSpeed = initialize_constant_speed(mesh, fields, 5.0);
  initialize_level_set(mesh, fields, flower_2D);

  const double Courant = 0.25;
  const double Ttotal = 0.01;
  const double dx = mesh_minimum_length_scale(mesh, *fields.coordsField);
  const unsigned Nt = (Ttotal+0.5*(Courant*dx/maxSpeed))/(Courant*dx/maxSpeed);
  const double dt = Ttotal/Nt;
  const double eps = 1.5*dx;
  const double dtau = 0.2*dx; // Reinitialization time step, based on unit advection speed used in reinitialization

  write_fields(stkIo, outputFileIndex, 0.0);
  reinitialize_level_set(mesh, fields, eps, dtau, &stkIo, 0.0, dt, outputFileIndex);

  std::cout << "Evolving to time " << Ttotal << " with " << Nt << " steps with Courant number " << Courant << std::endl;

  double time = 0.0;
  for (unsigned n=0; n<Nt; ++n)
  {
    mesh.update_field_data_states();
    evolve_level_set(mesh, fields, eps, dt);
    time += dt;
    write_fields(stkIo, outputFileIndex, time);
    reinitialize_level_set(mesh, fields, eps, dtau, &stkIo, time, time+dt, outputFileIndex);

    std::cout << "Time, Level set vol = " << time << " " << compute_level_set_volume(mesh, fields, eps) << std::endl;
  }

  const double finalVol = compute_level_set_volume(mesh, fields, eps);
  EXPECT_NEAR(0.545715, finalVol, 0.01);

  std::remove(outputFileName.c_str()); // Remove output file unless debugging
}

}
