// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#include <Akri_Eikonal_Calc.hpp>
#include <Akri_FastIterativeMethod.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_ParallelErrorMessage.hpp>
#include <Akri_Sign.hpp>
#include <math.h>

namespace krino{

void check_supported_element_topologies(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & selector)
{
  ParallelErrorMessage err(mesh.parallel());
  const stk::mesh::BucketVector& buckets = selector.get_buckets(stk::topology::ELEMENT_RANK);
  for (auto && bucket : buckets)
  {
    if (selector(*bucket) &&
        bucket->topology() != stk::topology::TRIANGLE_3_2D &&
        bucket->topology() != stk::topology::TETRAHEDRON_4)
    {
      err << "Topology " << bucket->topology().name() << " is not supported in Fast_Marching.\n";
    }
  }
  auto globalError = err.gather_message();
  if (globalError.first)
  {
    krinolog << globalError.second << stk::diag::dendl;
  }
  STK_ThrowRequireMsg(!globalError.first, "Unsupported element topology in Fast Marching or Fast Iterative Method.");
}

FastIterativeMethod::FastIterativeMethod(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & selector,
    const FieldRef& coordinates,
    const FieldRef& distance,
    const std::function<double(ParallelErrorMessage& err, stk::mesh::Entity)> & get_interface_speed,
    stk::diag::Timer & parentTimer)
  : myMesh(mesh),
    mySelector(selector),
    myCoordinates(coordinates),
    myDistance(distance),
    my_get_interface_speed(get_interface_speed),
    myTimer("Fast Iterative Method", parentTimer)
{
  check_supported_element_topologies(mesh, selector);
}

static size_t get_global_size(const stk::mesh::BulkData & mesh, const std::set<stk::mesh::Entity> & workingSet)
{
  size_t workingSetSize = 0;
  for (auto && node : workingSet)
    if (mesh.bucket(node).owned())
      ++workingSetSize;
  const size_t localWorkingSetSize = workingSetSize;
  stk::all_reduce_sum(mesh.parallel(), &localWorkingSetSize, &workingSetSize, 1);

  return workingSetSize;
}

static size_t get_global_num_selected_nodes(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & selector)
{
  size_t numNodes = 0;

  for ( auto&& bucket : mesh.get_buckets(stk::topology::NODE_RANK, selector) )
    numNodes += bucket->size();

  const size_t localWorkingSetSize = numNodes;
  stk::all_reduce_sum(mesh.parallel(), &localWorkingSetSize, &numNodes, 1);

  return numNodes;
}

void FastIterativeMethod::redistance()
{
  stk::diag::TimeBlock timer__(myTimer);

  {
    // make sure field is parallel consistent to start out
    std::vector<const stk::mesh::FieldBase *> parallel_fields(1, &myDistance.field());
    stk::mesh::copy_owned_to_shared(mesh(), parallel_fields);
  }

  ParallelErrorMessage err(mesh().parallel());

  const std::set<stk::mesh::Entity> initialNodes = initialize(err);
  std::set<stk::mesh::Entity> workingSet = build_initial_working_set(err, initialNodes);

  size_t workingSetSize = get_global_size(mesh(), workingSet);

  const size_t iterMax = get_global_num_selected_nodes(mesh(), field_selector() & mesh().mesh_meta_data().locally_owned_part()); // Should be upper bound
  size_t iter = 0;
  while (workingSetSize>0 && iter++<iterMax)
  {
    krinolog << "Fast iterative method iteration " << iter << " working set size " << workingSetSize << stk::diag::dendl;
    advance(err, initialNodes, workingSet);
    workingSetSize = get_global_size(mesh(), workingSet);
  }

  FieldRef oldDistance = myDistance.field_state(stk::mesh::StateOld);
  stk::mesh::field_copy(myDistance, oldDistance);
}

std::vector<double> FastIterativeMethod::compute_updates_and_communicate(ParallelErrorMessage& err, const std::set<stk::mesh::Entity> & nodes) const
{
  std::vector<double> distances(nodes.size(),0.);

  size_t index = 0;
  for (auto && node : nodes)
    distances[index++] = updated_node_distance(err, node);

  parallel_communicate_node_distances(err, nodes, distances);

  return distances;
}

void FastIterativeMethod::advance(ParallelErrorMessage& err, const std::set<stk::mesh::Entity> & initialNodes, std::set<stk::mesh::Entity> & workingSet)
{
  std::set<stk::mesh::Entity> oldWorkingSet;
  workingSet.swap(oldWorkingSet);

  const std::vector<double> updateDistances = compute_updates_and_communicate(err, oldWorkingSet);

  // Update changed nodes and put them into workingSet and put neighbors of unchanged nodes into neighborCandidatesForWorkingSet
  std::vector<stk::mesh::Entity> nodeNbrs;
  std::set<stk::mesh::Entity> neighborCandidatesForWorkingSet;
  size_t index = 0;
  for (auto && node : oldWorkingSet)
  {
    double & nodeDistance = *field_data<double>(myDistance, node);
    const double updateDistance = updateDistances[index++];
    if (nodeDistance == updateDistance)
    {
      fill_noninitial_node_neighbors(initialNodes, node, nodeNbrs);
      neighborCandidatesForWorkingSet.insert(nodeNbrs.begin(), nodeNbrs.end());
    }
    else
    {
      workingSet.insert(node);
      nodeDistance = updateDistance;
    }
  }

  parallel_communicate_nodes(err, neighborCandidatesForWorkingSet);

  const std::vector<double> nbrUpdateDistances = compute_updates_and_communicate(err, neighborCandidatesForWorkingSet);

  // Update changed neighbors and put them into workingSet
  index = 0;
  for (auto && nbr : neighborCandidatesForWorkingSet)
  {
    double & nbrDistance = *field_data<double>(myDistance, nbr);
    const double nbrUpdateDistance = nbrUpdateDistances[index++];
    if (nbrDistance != nbrUpdateDistance)
    {
      workingSet.insert(nbr);
      nbrDistance = nbrUpdateDistance;
    }
  }
}

double
FastIterativeMethod::update_triangle(const stk::mesh::Entity * elemNodes, int nodeToUpdate, const double speed) const
{
  static constexpr double far = std::numeric_limits<double>::max();
  const std::array<int,3> lnn = get_oriented_nodes_triangle(nodeToUpdate);
  const std::array<stk::math::Vector3d,3> x{stk::math::Vector3d(field_data<double>(myCoordinates, elemNodes[lnn[0]]),2), stk::math::Vector3d(field_data<double>(myCoordinates, elemNodes[lnn[1]]),2), stk::math::Vector3d(field_data<double>(myCoordinates, elemNodes[lnn[2]]),2)};
  const std::array<double,2> d{*field_data<double>(myDistance, elemNodes[lnn[0]]), *field_data<double>(myDistance, elemNodes[lnn[1]])};
  const int distSign = sign(*field_data<double>(myDistance, elemNodes[lnn[2]]));
  return eikonal_solve_triangle(x, d, distSign, far, speed);
}

double
FastIterativeMethod::update_tetrahedron(const stk::mesh::Entity * elemNodes, int nodeToUpdate, const double speed) const
{
  static constexpr double far = std::numeric_limits<double>::max();
  const std::array<int,4> lnn = get_oriented_nodes_tetrahedron(nodeToUpdate);
  const std::array<stk::math::Vector3d,4> x{field_data<double>(myCoordinates, elemNodes[lnn[0]]), field_data<double>(myCoordinates, elemNodes[lnn[1]]), field_data<double>(myCoordinates, elemNodes[lnn[2]]), field_data<double>(myCoordinates, elemNodes[lnn[3]])};
  const std::array<double,3> d{*field_data<double>(myDistance, elemNodes[lnn[0]]), *field_data<double>(myDistance, elemNodes[lnn[1]]), *field_data<double>(myDistance, elemNodes[lnn[2]])};
  const int distSign = sign(*field_data<double>(myDistance, elemNodes[lnn[3]]));
  return eikonal_solve_tetrahedron(x, d, distSign, far, speed);
}

double FastIterativeMethod::element_signed_distance_for_node(ParallelErrorMessage& err, const stk::mesh::Entity & elem, const stk::mesh::Entity & node) const
{
  const stk::mesh::Entity * elemNodes = mesh().begin_nodes(elem);
  int nodeOfElement = std::distance(elemNodes, std::find(elemNodes, mesh().end_nodes(elem), node));

  const double elementSpeed = my_get_interface_speed ? my_get_interface_speed(err, elem) : 1.0;

  const int npe_dist = myMesh.num_nodes(elem);
  if (3 == npe_dist)
  {
    return update_triangle(elemNodes, nodeOfElement, elementSpeed);
  }
  else if (4 == npe_dist)
  {
    return update_tetrahedron(elemNodes, nodeOfElement, elementSpeed);
  }
  STK_ThrowRequireMsg(false, "Unexpected number of nodes per element: " << npe_dist);
  return 0.0;
}

double FastIterativeMethod::updated_node_distance(ParallelErrorMessage& err, const stk::mesh::Entity & node) const
{
  const int distSign = sign(*field_data<double>(myDistance, node));
  double updateDistance = distSign * std::numeric_limits<double>::max();

  stk::mesh::Selector fieldNotGhost = field_not_ghost_selector();
  for (auto elem : StkMeshEntities{myMesh.begin_elements(node), myMesh.end_elements(node)})
  {
    if (fieldNotGhost(mesh().bucket(elem)))
    {
      const double elemNodeDistance = element_signed_distance_for_node(err, elem, node);
      updateDistance = distSign * std::min(std::abs(updateDistance), std::abs(elemNodeDistance));
    }
  }

  return updateDistance;
}

stk::mesh::Selector FastIterativeMethod::field_not_ghost_selector() const
{
  stk::mesh::Selector selector = (mesh().mesh_meta_data().locally_owned_part() | mesh().mesh_meta_data().globally_shared_part()) & mySelector & stk::mesh::selectField(myDistance);
  return selector;
}

stk::mesh::Selector FastIterativeMethod::field_selector() const
{
  stk::mesh::Selector selector = mySelector & stk::mesh::selectField(myDistance);
  return selector;
}

bool
FastIterativeMethod::have_crossing(const stk::mesh::Entity & elem) const
{
  const unsigned npe = mesh().bucket(elem).topology().num_nodes();
  STK_ThrowAssert(npe > 0);

  const stk::mesh::Entity * elem_nodes = mesh().begin(elem, stk::topology::NODE_RANK);
  const double * dist0 = field_data<double>(myDistance, elem_nodes[0]);
  STK_ThrowAssert(nullptr != dist0);
  for (unsigned n=1; n<npe; ++n)
  {
    const double * dist = field_data<double>(myDistance, elem_nodes[n]);
    STK_ThrowAssert(nullptr != dist);
    if (sign_change(*dist0, *dist))
    {
      return true;
    }
  }
  return false;
}

void
FastIterativeMethod::initialize_element(const stk::mesh::Entity & elem, const double speed)
{
  // To start the nodes of elements that have interfaces will be redistanced.
  // I have tried a few different methods for this distance calculation with varying success.
  // We can:
  // 1. Use only local facets to redistance the nodes of the element.  This appears to not work too
  //    well because the closest facet to a node might be in an element that does not support the node.
  // 2. Use local facets, but use the normal distance instead of the facet distance.  This seems to work
  //    better than (1) usually.  However, it seems prone to pathalogical behavior when small facets run
  //    parallel to an element side, producing inaccurate distance measures on this side.
  // 3. Use the standard parallel redistance calculation used in "normal" redistancing. This is susceptible
  //    to seeing through walls if the walls are thinner than the distance of the nodes of the cut element to
  //    the surface.
  // 4. Start the redistancing from the subelement facets, progressing through the subelements.
  //    Somewhat surprisingly, this does not work too well.  I think this may be due to the obtuse angles
  //    in the subelement decomposition.  The results are very similar to that produced by local redistancing (#1).
  // 5. Rescale each cut element so that it has a unit gradient (or prescribed gradient) and set the nodal
  //    distance to the minimum (magnitude) from all of the supporting elements.  This seems to work pretty
  //    well in the test cases and is completely local, and is not susceptible to seeing through walls.

  // Initialize using method #5 (element rescaling)

  const stk::mesh::Entity * elemNodes = mesh().begin(elem, stk::topology::NODE_RANK);
  const int npe = mesh().bucket(elem).topology().num_nodes();

  FieldRef oldDistance = myDistance.field_state(stk::mesh::StateOld);
  const double gradientMagnitude = calculate_gradient_magnitude(npe, elemNodes, oldDistance, myCoordinates);

  for (auto node : StkMeshEntities{myMesh.begin_nodes(elem), myMesh.end_nodes(elem)})
  {
    const double elemNodeSignedDistance = *field_data<double>(oldDistance, node) / (gradientMagnitude * speed);
    double & nodeDistance = *field_data<double>(myDistance, node);
    const int distSign = sign(nodeDistance);
    nodeDistance = distSign * std::min(std::abs(nodeDistance), std::abs(elemNodeSignedDistance));
  }
}

static
void pack_shared_nodes_for_sharing_procs(const stk::mesh::BulkData & mesh,
    const std::set<stk::mesh::Entity> & workingSet,
    stk::CommSparse &commSparse)
{
  std::vector<int> nodeSharedProcs;
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (auto node : workingSet)
    {
      if (mesh.bucket(node).shared())
      {
        const stk::mesh::EntityId nodeId = mesh.identifier(node);
        mesh.comm_shared_procs(node, nodeSharedProcs);
        for (int procId : nodeSharedProcs)
          commSparse.send_buffer(procId).pack(nodeId);
      }
    }
  });
}

static
void unpack_shared_nodes(const stk::mesh::BulkData & mesh,
    std::set<stk::mesh::Entity> & workingSet,
    stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      stk::mesh::EntityId nodeId;
      commSparse.recv_buffer(procId).unpack(nodeId);

      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeId);
      STK_ThrowRequire(mesh.is_valid(node));

      workingSet.insert(node);
    }
  });
}

void FastIterativeMethod::parallel_communicate_nodes(ParallelErrorMessage& /*err*/, std::set<stk::mesh::Entity> & nodes) const
{
  if (mesh().parallel_size() > 1)
  {
    stk::CommSparse commSparse(mesh().parallel());
    pack_shared_nodes_for_sharing_procs(mesh(), nodes, commSparse);
    unpack_shared_nodes(mesh(), nodes, commSparse);
  }
}

static
void pack_shared_node_distances_for_sharing_procs(const stk::mesh::BulkData & mesh,
    const std::set<stk::mesh::Entity> & nodes,
    const std::vector<double> & distances,
    stk::CommSparse &commSparse)
{
  std::vector<int> nodeSharedProcs;
  stk::pack_and_communicate(commSparse,[&]()
  {
    size_t index=0;
    for (auto && node : nodes)
    {
      const double distance = distances[index++];
      if (mesh.bucket(node).shared())
      {
        const stk::mesh::EntityId nodeId = mesh.identifier(node);
        mesh.comm_shared_procs(node, nodeSharedProcs);
        for (int procId : nodeSharedProcs)
        {
          commSparse.send_buffer(procId).pack(nodeId);
          commSparse.send_buffer(procId).pack(distance);
        }
      }
    }
  });
}

static
void unpack_shared_node_distances(const stk::mesh::BulkData & mesh,
    const std::set<stk::mesh::Entity> & nodes,
    std::vector<double> & distances,
    stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      stk::mesh::EntityId nodeId;
      commSparse.recv_buffer(procId).unpack(nodeId);
      double updateDistance = 0.;
      commSparse.recv_buffer(procId).unpack(updateDistance);

      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeId);
      const auto nodeIter = nodes.find(node);
      STK_ThrowAssert(nodeIter != nodes.end());

      double & nodeDistance = distances[std::distance(nodes.begin(), nodeIter)];
      const int distSign = sign(nodeDistance);
      nodeDistance = distSign * std::min(std::abs(nodeDistance), std::abs(updateDistance));
    }
  });
}

void FastIterativeMethod::parallel_communicate_node_distances(ParallelErrorMessage& /*err*/, const std::set<stk::mesh::Entity> & nodes, std::vector<double> & distances) const
{
  if (mesh().parallel_size() > 1)
  {
    stk::CommSparse commSparse(mesh().parallel());
    pack_shared_node_distances_for_sharing_procs(mesh(), nodes, distances, commSparse);
    unpack_shared_node_distances(mesh(), nodes, distances, commSparse);
  }
}

void FastIterativeMethod::parallel_communicate_nodes_and_nodal_distance(ParallelErrorMessage& err, std::set<stk::mesh::Entity> & nodes)
{
  if (mesh().parallel_size() > 1)
  {
    parallel_communicate_nodes(err, nodes);

    std::vector<double> distances;
    distances.reserve(nodes.size());

    for (auto && node : nodes)
      distances.push_back(*field_data<double>(myDistance, node));

    parallel_communicate_node_distances(err, nodes, distances);

    size_t index = 0;
    for (auto && node : nodes)
      *field_data<double>(myDistance, node) = distances[index++];
  }
}

void FastIterativeMethod::fill_noninitial_node_neighbors(const std::set<stk::mesh::Entity> & initialNodes, stk::mesh::Entity node, std::vector<stk::mesh::Entity> & nodeNbrs) const
{
  std::set<stk::mesh::Entity> allNodeNbrs;
  stk::mesh::Selector fieldNotGhost = field_not_ghost_selector();
  for (auto elem : StkMeshEntities{myMesh.begin_elements(node), myMesh.end_elements(node)})
    if (fieldNotGhost(mesh().bucket(elem)))
      for (auto nbr : StkMeshEntities{myMesh.begin_nodes(elem), myMesh.end_nodes(elem)})
        if (nbr != node)
          allNodeNbrs.insert(nbr);

  nodeNbrs.clear();
  for (auto && nodeNbr : allNodeNbrs)
    if (initialNodes.find(nodeNbr) == initialNodes.end())
      nodeNbrs.push_back(nodeNbr);
}

std::set<stk::mesh::Entity> FastIterativeMethod::initialize(ParallelErrorMessage& err)
{
  stk::mesh::Selector fieldNotGhost = field_not_ghost_selector();

  STK_ThrowRequireMsg(myDistance.number_of_states() > 1, "Fast iterative method requires multiple states.");

  FieldRef oldDistance = myDistance.field_state(stk::mesh::StateOld);
  stk::mesh::field_copy(myDistance, oldDistance);

  for ( auto&& bucket : mesh().get_buckets(stk::topology::NODE_RANK, field_selector()) )
  {
    for (auto && node : *bucket)
    {
      const double nodeOldDistance = *field_data<double>(oldDistance, node);
      double & nodeDistance = *field_data<double>(myDistance, node);
      nodeDistance = sign(nodeOldDistance) * std::numeric_limits<double>::max();
    }
  }

  std::set<stk::mesh::Entity> initialNodes;
  for ( auto&& bucket : mesh().get_buckets(stk::topology::ELEMENT_RANK, fieldNotGhost) )
  {
    for (auto&& elem : *bucket)
    {
      if (have_crossing(elem))
      {
        const double speed = my_get_interface_speed ? my_get_interface_speed(err, elem) : 1.0;
        initialize_element(elem, speed);
        initialNodes.insert(myMesh.begin_nodes(elem), myMesh.end_nodes(elem));
      }
    }
  }

  parallel_communicate_nodes_and_nodal_distance(err, initialNodes);

  return initialNodes;
}

std::set<stk::mesh::Entity> FastIterativeMethod::build_initial_working_set(ParallelErrorMessage& err, const std::set<stk::mesh::Entity> & initialNodes)
{
  std::set<stk::mesh::Entity> workingSet;
  std::vector<stk::mesh::Entity> nodeNbrs;
  for (auto && initialNode : initialNodes)
  {
    fill_noninitial_node_neighbors(initialNodes, initialNode, nodeNbrs);
    workingSet.insert(nodeNbrs.begin(), nodeNbrs.end());
  }

  parallel_communicate_nodes(err, workingSet);

  const std::vector<double> nodeDistanceUpdates = compute_updates_and_communicate(err, workingSet);

  size_t index = 0;
  for (auto && node : workingSet)
    *field_data<double>(myDistance, node) = nodeDistanceUpdates[index++];

  return workingSet;
}

bool FastIterativeMethod::check_converged_solution() const
{
  ParallelErrorMessage err(mesh().parallel());

  stk::mesh::Selector fieldNotGhost = field_not_ghost_selector();

  krinolog << "Checking converged solution..." << stk::diag::dendl;

  // Get all non-intial nodes by getting the initial nodes, parallel sync, check against this set

  std::set<stk::mesh::Entity> initialNodes;

  for ( auto&& bucket : mesh().get_buckets(stk::topology::ELEMENT_RANK, fieldNotGhost) )
    for (auto && elem : *bucket)
      if (have_crossing(elem))
        initialNodes.insert(myMesh.begin_nodes(elem), myMesh.end_nodes(elem));

  parallel_communicate_nodes(err, initialNodes);

  std::set<stk::mesh::Entity> nonInitialNodes;
  for ( auto&& bucket : mesh().get_buckets(stk::topology::NODE_RANK, fieldNotGhost) )
    for (auto && node : *bucket)
      if (initialNodes.find(node) == initialNodes.end())
        nonInitialNodes.insert(node);

  const std::vector<double> updateDistances = compute_updates_and_communicate(err, nonInitialNodes);

  size_t index = 0;
  for (auto && node : nonInitialNodes)
  {
    const double convergedDistance = *field_data<double>(myDistance, node);
    const double updateDistance = updateDistances[index++];
    if (updateDistance != convergedDistance)
      err << "Distance changed for node " << mesh().identifier(node) << " from " << convergedDistance << " to  " << updateDistance << " with change " << updateDistance-convergedDistance << "\n";
  }

  auto errMsg = err.gather_message();
  krinolog << errMsg.second;
  return !errMsg.first;
}

}
