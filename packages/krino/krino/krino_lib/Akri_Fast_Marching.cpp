// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>

#include <Akri_Fast_Marching.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_ParallelErrorMessage.hpp>
#include <Akri_Sign.hpp>
#include <math.h>

#include "Akri_Eikonal_Calc.hpp"
namespace krino{

Fast_Marching::Fast_Marching(const stk::mesh::BulkData & mesh,
      const stk::mesh::Selector & activeElementSelector,
      const FieldRef& coordinates,
      const FieldRef& distance,
      const std::function<double(ParallelErrorMessage& err, stk::mesh::Entity)> & get_interface_speed,
      stk::diag::Timer & parentTimer)
  : myMesh(mesh),
    mySelector(activeElementSelector),
    myCoordinates(coordinates),
    myDistance(distance),
    my_get_interface_speed(get_interface_speed),
    myTimer("Fast Marching", parentTimer),
    my_fm_node_less(mesh),
    trial_nodes(my_fm_node_less)
{
  ParallelErrorMessage err(mesh.parallel());
  stk::mesh::Selector fieldSelector(myDistance);
  const stk::mesh::BucketVector& buckets = activeElementSelector.get_buckets(stk::topology::ELEMENT_RANK);
  for (auto && bucket : buckets)
  {
    if (fieldSelector(*bucket) &&
        bucket->topology() != stk::topology::TRIANGLE_3_2D &&
        bucket->topology() != stk::topology::QUADRILATERAL_4_2D &&
        bucket->topology() != stk::topology::TETRAHEDRON_4 &&
        bucket->topology() != stk::topology::HEXAHEDRON_8)
    {
      err << "Topology " << bucket->topology().name() << " is not supported in Fast_Marching.\n";
    }
  }
  check_error(err, "Checking topology");
}

double & Fast_Marching::get_node_distance(const stk::mesh::Entity node)
{
  return get_scalar_field(myMesh, myDistance, node);
}

double Fast_Marching::get_node_distance(const stk::mesh::Entity node) const
{
  return get_scalar_field(myMesh, myDistance, node);
}

double Fast_Marching::get_element_interface_speed(ParallelErrorMessage& err, const stk::mesh::Entity elem) const
{
  return my_get_interface_speed ? my_get_interface_speed(err, elem) : 1.0;
}

void
Fast_Marching::check_error(const ParallelErrorMessage& err, const std::string & context) const
{
  auto globalError = err.gather_message();
  if (globalError.first)
  {
    krinolog<< "Error in " << context << ":" << stk::diag::dendl;
    krinolog << globalError.second << stk::diag::dendl;
  }
  STK_ThrowRequireMsg(!globalError.first, "Error in " << context << ".");
}

void Fast_Marching::redistance(const std::vector<stk::mesh::Entity> & initialNodes)
{
  stk::diag::TimeBlock timer__(myTimer);

  initialize_algorithm();
  initialize_given_nodes(initialNodes);
  propagate_distance_to_convergence();
}

void Fast_Marching::initialize_given_nodes(const std::vector<stk::mesh::Entity> & initialNodes)
{
  for (auto initialNode : initialNodes)
  {
    Fast_Marching_Node * fm_node = get_fm_node(initialNode);
    STK_ThrowAssert(nullptr != fm_node && fm_node->status() != STATUS_UNUSED);
    const double nodeDist = get_node_distance(initialNode);
    fm_node->set_signed_dist(nodeDist);
    fm_node->set_status(STATUS_INITIAL);
  }
}

void Fast_Marching::redistance()
{
  stk::diag::TimeBlock timer__(myTimer);

  initialize_algorithm();
  initialize_nodes_of_crossed_elements();
  propagate_distance_to_convergence();
}

void Fast_Marching::initialize_algorithm()
{
  // make sure field is parallel consistent to start out
  {
    std::vector<const stk::mesh::FieldBase *> parallel_fields(1, &myDistance.field());
    stk::mesh::copy_owned_to_shared(myMesh, parallel_fields);
  }

  fm_nodes.clear();
  fm_nodes.resize(stk::mesh::count_selected_entities(myMesh.mesh_meta_data().universal_part(), myMesh.buckets(stk::topology::NODE_RANK)));

  std::vector<stk::mesh::Entity> field_nodes;
  stk::mesh::get_selected_entities( selected_with_field_not_ghost_selector(), myMesh.buckets( stk::topology::NODE_RANK ), field_nodes );

  for ( auto&& node : field_nodes )
    initialize_node(node);
}

void Fast_Marching::initialize_nodes_of_crossed_elements()
{
  ParallelErrorMessage err(myMesh.parallel());

  std::vector<stk::mesh::Entity> field_elems;
  stk::mesh::get_selected_entities( selected_with_field_not_ghost_selector(), myMesh.buckets( stk::topology::ELEMENT_RANK ), field_elems );
  for (auto&& elem : field_elems)
    if (have_crossing(elem))
      initialize_element(elem, err);

  if (myMesh.parallel_size() > 1)
  {
    stk::mesh::Selector globally_shared_selector = selected_with_field_selector() & myMesh.mesh_meta_data().globally_shared_part();
    std::vector< stk::mesh::Entity> shared_nodes;
    stk::mesh::get_selected_entities( globally_shared_selector, myMesh.buckets( stk::topology::NODE_RANK ), shared_nodes );

    for ( auto && shared_node : shared_nodes )
    {
      Fast_Marching_Node * fm_node = get_fm_node(shared_node);
      if (nullptr != fm_node && fm_node->status() != STATUS_UNUSED)
      {
        double & node_dist = get_node_distance(fm_node->node());
        node_dist = fm_node->signed_dist()*fm_node->sign();
      }
    }
    stk::mesh::parallel_min(myMesh, {&myDistance.field()});

    for ( auto && shared_node : shared_nodes )
    {
      Fast_Marching_Node * fm_node = get_fm_node(shared_node);
      if (nullptr != fm_node && fm_node->status() != STATUS_UNUSED)
      {
        stk::mesh::Entity node = fm_node->node();
        const double min_node_unsigned_dist = get_node_distance(node);
        const double fm_node_unsigned_dist = fm_node->signed_dist()*fm_node->sign();
        fm_node->set_signed_dist(min_node_unsigned_dist*fm_node->sign());
        if (min_node_unsigned_dist < fm_node_unsigned_dist)
        {
          fm_node->set_status(STATUS_INITIAL);
        }
      }
    }
  }

  check_error(err, "Fast Marching Initialization");
}

void Fast_Marching::propagate_distance_to_convergence()
{
  ParallelErrorMessage err(myMesh.parallel());

  // neighbors of initial nodes are trial nodes
  for ( auto && fm_node : fm_nodes )
  {
    if (fm_node.status() == STATUS_INITIAL)
    {
      update_neighbors(fm_node, err);
    }
  }

  stk::mesh::Selector globally_shared_selector = selected_with_field_selector() & myMesh.mesh_meta_data().globally_shared_part();
  std::vector< stk::mesh::Entity> shared_nodes;
  stk::mesh::get_selected_entities( globally_shared_selector, myMesh.buckets( stk::topology::NODE_RANK ), shared_nodes );

  bool done = false;

  const size_t local_num_nodes = fm_nodes.size();
  size_t max_num_nodes = 0;
  stk::all_reduce_max(myMesh.parallel(), &local_num_nodes, &max_num_nodes, 1);

  const unsigned max_outer_steps = max_num_nodes;  // should be overkill
  unsigned num_outer_steps = 0;
  while (!done)
  {
    if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << "num_trial nodes = " << trial_nodes.size() << "\n";

    if (num_outer_steps++ > max_outer_steps)
    {
      err << "Error: Outer loop of Fast_Marching::redistance did not converge!\n";
      break;
    }

    const unsigned max_inner_steps = 10*max_outer_steps;  // should be overkill
    unsigned num_inner_steps = 0;
    while(!trial_nodes.empty())
    {
      auto begin = trial_nodes.begin();
      Fast_Marching_Node & fm_node = **begin;
      trial_nodes.erase(begin);
      fm_node.set_status(STATUS_ACCEPTED);
      update_neighbors(fm_node, err);
      if (num_inner_steps++ > max_inner_steps)
      {
        err << "Error: Inner loop of Fast_Marching::redistance did not converge! Number of trial nodes = " << trial_nodes.size() << "\n";
      }
      if (err.have_local_error()) break;
    }

    check_error(err, "Fast Marching Iteration");

    unsigned num_locally_updated = 0;
    if (myMesh.parallel_size() > 1)
    {
      for ( auto && shared_node : shared_nodes )
      {
        Fast_Marching_Node * fm_node = get_fm_node(shared_node);
        if (nullptr != fm_node && fm_node->status() != STATUS_UNUSED)
        {
          double & node_dist = get_node_distance(fm_node->node());
          node_dist = fm_node->signed_dist()*fm_node->sign();
        }
      }
      stk::mesh::parallel_min(myMesh, {&myDistance.field()});

      for ( auto && shared_node : shared_nodes )
      {
        Fast_Marching_Node * fm_node = get_fm_node(shared_node);
        if (nullptr != fm_node && fm_node->status() != STATUS_UNUSED)
        {
          stk::mesh::Entity node = fm_node->node();
          const double min_node_unsigned_dist = get_node_distance(node);
          const double fm_node_unsigned_dist = fm_node->signed_dist()*fm_node->sign();
          if(min_node_unsigned_dist < fm_node_unsigned_dist)
          {
            STK_ThrowAssertMsg(fm_node->status() != STATUS_INITIAL || fm_node->status() != STATUS_TRIAL, "Unexpected node to have INITIAL or TRIAL status.");
            fm_node->set_signed_dist(min_node_unsigned_dist*fm_node->sign());
            add_trial_node(*fm_node);
            ++num_locally_updated;
          }
        }
      }
    }

    unsigned num_globally_updated = 0;
    stk::all_reduce_sum(myMesh.parallel(), &num_locally_updated, &num_globally_updated, 1);

    done = (num_globally_updated == 0);
  }

  for ( auto && fm_node : fm_nodes )
  {
    if (fm_node.status() == STATUS_TRIAL || fm_node.status() == STATUS_FAR)
    {
      err << "Node " << myMesh.identifier(fm_node.node()) << " with status " << fm_node.status() << " with distance " << fm_node.signed_dist() << " did not get updated!\n";
    }
    if (fm_node.status() != STATUS_UNUSED)
    {
      double & node_dist = get_node_distance(fm_node.node());
      node_dist = fm_node.signed_dist();
    }
  }

  check_error(err, "Fast Marching Update");
}

Fast_Marching_Node * Fast_Marching::get_fm_node(stk::mesh::Entity node)
{
  if (myMesh.is_valid(node) && myMesh.local_id(node) < fm_nodes.size())
  {
    return &fm_nodes[myMesh.local_id(node)];
  }
  else
  {
    return nullptr;
  }
}

stk::mesh::Selector Fast_Marching::selected_with_field_not_ghost_selector() const
{
  stk::mesh::Selector selector = mySelector & stk::mesh::selectField(myDistance) & (myMesh.mesh_meta_data().locally_owned_part() | myMesh.mesh_meta_data().globally_shared_part());
  return selector;
}

stk::mesh::Selector Fast_Marching::selected_with_field_selector() const
{
  stk::mesh::Selector selector = mySelector & stk::mesh::selectField(myDistance);
  return selector;
}

void Fast_Marching::initialize_node(const stk::mesh::Entity node)
{
  const int currentSign = sign(get_node_distance(node));
  const stk::math::Vector3d coords = get_vector_field(myMesh, myCoordinates, node, myMesh.mesh_meta_data().spatial_dimension());

  Fast_Marching_Node * fm_node = get_fm_node(node);
  STK_ThrowAssert(nullptr != fm_node);
  *fm_node = Fast_Marching_Node(node,STATUS_FAR, currentSign*std::numeric_limits<double>::max(), currentSign, coords);
}

bool
Fast_Marching::have_crossing(const StkMeshEntities & elemNodes) const
{
  STK_ThrowAssert(elemNodes.size() > 0);
  const double dist0 = get_node_distance(elemNodes[0]);
  for (unsigned n=1; n<elemNodes.size(); ++n)
  {
    const double dist = get_node_distance(elemNodes[n]);
    if (sign_change(dist0, dist))
    {
      return true;
    }
  }
  return false;
}

bool
Fast_Marching::have_crossing(const stk::mesh::Entity & elem) const
{
  const StkMeshEntities elemNodes{myMesh.begin_nodes(elem), myMesh.end_nodes(elem)};
  return have_crossing(elemNodes);
}

static std::function<const stk::math::Vector3d &(stk::mesh::Entity)> build_get_fm_node_coordinates(Fast_Marching * fm)
{
  return [fm](stk::mesh::Entity node) -> const stk::math::Vector3d &
    {
      Fast_Marching_Node * fm_node = fm->get_fm_node(node);
      STK_ThrowAssert(fm_node);
      return fm_node->coords();
    };
}

template<size_t NNODES>
std::array<stk::mesh::Entity, NNODES> get_simplex_nodes_from_element_nodes(const stk::mesh::Entity* elemNodes, const std::array<int, NNODES> & simplexNodeIndices)
{
  std::array<stk::mesh::Entity, NNODES> simplexNodes;
  for (size_t n=0; n<NNODES; ++n)
    simplexNodes[n] = elemNodes[simplexNodeIndices[n]];
  return simplexNodes;
}

template<size_t NNODES>
void Fast_Marching::initialize_from_simplex(const std::array<stk::mesh::Entity, NNODES> & simplexNodes,
    const double elemSpeed,
    const std::function<const stk::math::Vector3d &(stk::mesh::Entity)> & get_coordinates)
{
  const StkMeshEntities simplexNodeEntities{simplexNodes.data(), simplexNodes.data()+NNODES};
  if (have_crossing(simplexNodeEntities))
  {
    const double mag_grad = calculate_gradient_magnitude(NNODES, simplexNodes.data(), myDistance, get_coordinates);

    for (size_t inode=0; inode<NNODES; ++inode)
    {
      Fast_Marching_Node * fm_node = get_fm_node(simplexNodes[inode]);
      STK_ThrowAssert(nullptr != fm_node && fm_node->status() != STATUS_UNUSED);
      const double elem_node_dist = get_node_distance(simplexNodes[inode]) / (mag_grad * elemSpeed);
      const int sign = fm_node->sign();
      fm_node->set_signed_dist(sign * std::min(std::abs(fm_node->signed_dist()), std::abs(elem_node_dist)));
      fm_node->set_status(STATUS_INITIAL);
    }
  }
}

template<size_t NSIMPLICES, size_t NNODES>
void Fast_Marching::initialize_from_simplices(const stk::mesh::Entity* elemNodes, const std::array<std::array<int, NNODES>,NSIMPLICES> & simplices,
    const double elemSpeed,
    const std::function<const stk::math::Vector3d &(stk::mesh::Entity)> & get_coordinates)
{
  for (auto & simplex : simplices)
    initialize_from_simplex(get_simplex_nodes_from_element_nodes(elemNodes, simplex), elemSpeed, get_coordinates);
}

static constexpr std::array<std::array<int,3>,1> triangleSimplices = {{ {{0,1,2}} }};
static constexpr std::array<std::array<int,3>,4> quadrilateralSimplices = {{ {{0,1,2}},  {{0,2,3}}, {{0,1,3}}, {{1,2,3}} }};
static constexpr std::array<std::array<int,4>,1> tetrahedronSimplices = {{ {{0,1,2,3}} }};
static constexpr std::array<std::array<int,4>,10> hexahedronSimplices = {{
    {{0,1,3,4}}, {{1,2,3,6}}, {{3,4,6,7}}, {{1,6,4,5}}, {{1,3,4,6}},
    {{0,1,2,5}}, {{0,2,3,7}}, {{0,5,7,4}}, {{2,7,5,6}}, {{0,2,7,5}}
}};

void
Fast_Marching::initialize_element(const stk::mesh::Entity & elem, ParallelErrorMessage& err)
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

  auto get_coordinates = build_get_fm_node_coordinates(this);

  const stk::topology elemTopology = myMesh.bucket(elem).topology();

  const double elemSpeed = get_element_interface_speed(err, elem);

  const stk::mesh::Entity* elemNodes = myMesh.begin_nodes(elem);

  switch(elemTopology())
  {
    case stk::topology::TRIANGLE_3_2D:
        initialize_from_simplices(elemNodes, triangleSimplices, elemSpeed, get_coordinates);
        break;
    case stk::topology::QUADRILATERAL_4_2D:
        initialize_from_simplices(elemNodes, quadrilateralSimplices, elemSpeed, get_coordinates);
        break;
    case stk::topology::TETRAHEDRON_4:
        initialize_from_simplices(elemNodes, tetrahedronSimplices, elemSpeed, get_coordinates);
        break;
    case stk::topology::HEXAHEDRON_8:
        initialize_from_simplices(elemNodes, hexahedronSimplices, elemSpeed, get_coordinates);
        break;
    default:
        err << "Unsupported element topology " << elemTopology.name() << " in initialize_element.\n";
  }
}

static bool do_make_neighbor_a_trial_node(const Fast_Marching_Node & acceptedNode, const Fast_Marching_Node & nbrNode)
{
  if (nbrNode.status() == STATUS_FAR)
    return true;
  if (nbrNode.status() == STATUS_ACCEPTED)
  {
    const double accepted_node_unsigned_dist = acceptedNode.signed_dist()*acceptedNode.sign();
    const double nbr_unsigned_dist = nbrNode.signed_dist()*nbrNode.sign();
    if(nbr_unsigned_dist > accepted_node_unsigned_dist)
      return true;
  }
  return false;
}

template<size_t NNODES>
void Fast_Marching::update_neighbors_from_simplex(const Fast_Marching_Node & acceptedNode, const std::array<stk::mesh::Entity, NNODES> & simplexNodes, const double elemSpeed)
{
  int node_to_update = -1;
  int num_trial = 0;

  std::array<Fast_Marching_Node *, NNODES> simplexFmNodes;

  for ( size_t i = 0; i < NNODES; ++i )
  {
    Fast_Marching_Node * fm_nbr = get_fm_node(simplexNodes[i]);
    STK_ThrowAssertMsg(nullptr != fm_nbr && (STATUS_INITIAL == fm_nbr->status() || STATUS_ACCEPTED == fm_nbr->status() || STATUS_FAR == fm_nbr->status() || STATUS_TRIAL == fm_nbr->status()), "Unexpected node status.");
    simplexFmNodes[i] = fm_nbr;
    if (do_make_neighbor_a_trial_node(acceptedNode, *fm_nbr))
    {
      add_trial_node(*fm_nbr);
    }
    if (fm_nbr->status() == STATUS_TRIAL)
    {
      ++num_trial;
      node_to_update = i;
    }
  }

  if (1 == num_trial)
  {
    update_trial_node_from_simplex(simplexFmNodes, node_to_update, elemSpeed);
  }
}

template<size_t NSIMPLICES, size_t NNODES>
void Fast_Marching::update_neighbors_from_simplices(const Fast_Marching_Node & acceptedNode,
    const stk::mesh::Entity* elemNodes,
    const std::array<std::array<int, NNODES>,NSIMPLICES> & simplices,
    const double elemSpeed)
{
  for (auto & simplex : simplices)
    update_neighbors_from_simplex(acceptedNode, get_simplex_nodes_from_element_nodes(elemNodes, simplex), elemSpeed);
}

void
Fast_Marching::update_neighbors(Fast_Marching_Node & acceptedNode, ParallelErrorMessage& err)
{
  STK_ThrowAssertMsg(STATUS_ACCEPTED == acceptedNode.status() || STATUS_INITIAL == acceptedNode.status(), "Expected ACCEPTED OR INITIAL status");

  stk::mesh::Entity node = acceptedNode.node();

  const stk::mesh::Selector elemSelector = selected_with_field_not_ghost_selector();

  for (auto elem : StkMeshEntities{myMesh.begin_elements(node), myMesh.end_elements(node)})
  {
    if (myMesh.is_valid(elem) && elemSelector(myMesh.bucket(elem)))
    {
      const double elemSpeed = get_element_interface_speed(err, elem);
      const stk::topology elemTopology = myMesh.bucket(elem).topology();

      const stk::mesh::Entity* elemNodes = myMesh.begin_nodes(elem);

      switch(elemTopology())
      {
        case stk::topology::TRIANGLE_3_2D:
            update_neighbors_from_simplices(acceptedNode, elemNodes, triangleSimplices, elemSpeed);
            break;
        case stk::topology::QUADRILATERAL_4_2D:
            update_neighbors_from_simplices(acceptedNode, elemNodes, quadrilateralSimplices, elemSpeed);
            break;
        case stk::topology::TETRAHEDRON_4:
            update_neighbors_from_simplices(acceptedNode, elemNodes, tetrahedronSimplices, elemSpeed);
            break;
        case stk::topology::HEXAHEDRON_8:
            update_neighbors_from_simplices(acceptedNode, elemNodes,  hexahedronSimplices, elemSpeed);
            break;
        default:
            err << "Unsupported element topology " << elemTopology.name() << " in initialize_element.\n";
      }
    }
  }
}

void
Fast_Marching::add_trial_node(Fast_Marching_Node & trial_node)
{
  STK_ThrowAssertMsg(trial_node.status() == STATUS_ACCEPTED || trial_node.status() == STATUS_FAR, "Expected ACCEPTED or FAR when adding trial node");
  trial_nodes.insert(&trial_node);
  trial_node.set_status(STATUS_TRIAL);
}

void
Fast_Marching::update_trial_node(Fast_Marching_Node & trial_node, const double dist)
{
  STK_ThrowAssertMsg(trial_node.status() == STATUS_TRIAL, "Unexpected node status when updating trial node");

  auto it = trial_nodes.find(&trial_node);
  STK_ThrowAssertMsg(it != trial_nodes.end(), "Can't find trial node");

  trial_nodes.erase(it);

  trial_node.set_signed_dist(dist);
  trial_nodes.insert(&trial_node);
}

template <size_t NNODES>
void Fast_Marching::update_trial_node_from_simplex(const std::array<Fast_Marching_Node *, NNODES> & simplexNodes, const int nodeToUpdate, const double speed)
{
  static_assert(NNODES == 3 || NNODES == 4);
  double dist = std::numeric_limits<double>::max();

  if constexpr(4 == NNODES)
    dist = update_tetrahedron(simplexNodes, nodeToUpdate, speed);
  else
    dist = update_triangle(simplexNodes, nodeToUpdate, speed);

  Fast_Marching_Node & fmNode = *simplexNodes[nodeToUpdate];
  if (dist*fmNode.sign() < fmNode.signed_dist()*fmNode.sign())
  {
    update_trial_node(fmNode, dist);
  }
}

double
Fast_Marching::update_triangle(const std::array<Fast_Marching_Node *, 3> & elemNodes, const int nodeToUpdate, const double speed)
{
  static constexpr double far = std::numeric_limits<double>::max();

  const std::array<int,3> lnn = get_oriented_nodes_triangle(nodeToUpdate);
  const std::array<stk::math::Vector3d,3> x{elemNodes[lnn[0]]->coords(), elemNodes[lnn[1]]->coords(), elemNodes[lnn[2]]->coords()};
  const std::array<double,2> d{elemNodes[lnn[0]]->signed_dist(), elemNodes[lnn[1]]->signed_dist()};
  const int sign = elemNodes[lnn[2]]->sign();
  const double signedDist = eikonal_solve_triangle(x, d, sign, far, speed);
  // For fast marching, we strictly enforce that the distance is increasing in magnitude
  return sign * std::max({std::abs(d[0]),std::abs(d[1]),std::abs(signedDist)});
}

double
Fast_Marching::update_tetrahedron(const std::array<Fast_Marching_Node *, 4> & elemNodes, const int nodeToUpdate, const double speed)
{
  static constexpr double far = std::numeric_limits<double>::max();
  const std::array<int,4> lnn = get_oriented_nodes_tetrahedron(nodeToUpdate);
  const std::array<stk::math::Vector3d,4> x{elemNodes[lnn[0]]->coords(), elemNodes[lnn[1]]->coords(), elemNodes[lnn[2]]->coords(), elemNodes[lnn[3]]->coords()};
  const std::array<double,3> d{elemNodes[lnn[0]]->signed_dist(), elemNodes[lnn[1]]->signed_dist(), elemNodes[lnn[2]]->signed_dist()};
  const int sign = elemNodes[lnn[3]]->sign();
  const double signedDist = eikonal_solve_tetrahedron(x, d, sign, far, speed);
  // For fast marching, we strictly enforce that the distance is increasing in magnitude
  return sign * std::max({std::abs(d[0]),std::abs(d[1]),std::abs(d[2]),std::abs(signedDist)});
}

} // namespace krino
