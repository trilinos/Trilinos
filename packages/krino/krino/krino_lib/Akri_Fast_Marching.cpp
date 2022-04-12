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
#include <Akri_AuxMetaData.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_ParallelErrorMessage.hpp>
#include <math.h>

#include "Akri_Eikonal_Calc.hpp"
namespace krino{

Fast_Marching::Fast_Marching(LevelSet & ls, const stk::mesh::Selector & selector, stk::diag::Timer & parent_timer)
  : my_ls(ls),
    my_selector(selector),
    my_timer("Fast Marching", parent_timer),
    my_tri_timer("Update Triangle", my_timer),
    my_tet_timer("Update Tetrahedron", my_timer),
    my_fm_node_less(ls.mesh()),
    trial_nodes(my_fm_node_less)
{
  ParallelErrorMessage err(mesh().parallel());
  stk::mesh::Selector fieldSelector(my_ls.get_distance_field());
  const stk::mesh::BucketVector& buckets = selector.get_buckets(stk::topology::ELEMENT_RANK);
  for (auto && bucket : buckets)
  {
    if (fieldSelector(*bucket) &&
        bucket->topology() != stk::topology::TRIANGLE_3_2D &&
        bucket->topology() != stk::topology::TETRAHEDRON_4)
    {
      err << "Topology " << bucket->topology().name() << " is not supported in Fast_Marching.\n";
    }
  }
  check_error(err, "Checking topology");
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
  ThrowRequireMsg(!globalError.first, "Error in " << context << ".");
}

void Fast_Marching::redistance()
{
  stk::diag::TimeBlock timer__(my_timer);
  const FieldRef& dRef = my_ls.get_distance_field();
  const FieldRef& olddRef = my_ls.get_old_distance_field();

  // make sure field is parallel consistent to start out
  {
    std::vector<const stk::mesh::FieldBase *> parallel_fields(1, &dRef.field());
    stk::mesh::copy_owned_to_shared(mesh(), parallel_fields);
  }

  ParallelErrorMessage err(mesh().parallel());

  initialize(err);

  // neighbors of initial nodes are trial nodes
  for ( auto && fm_node : fm_nodes )
  {
    if (fm_node.status() == STATUS_INITIAL)
    {
      update_neighbors(fm_node, err);
    }
  }

  check_error(err, "Fast Marching Initialization");

  stk::mesh::Selector globally_shared_selector = my_selector & mesh().mesh_meta_data().globally_shared_part();
  std::vector< stk::mesh::Entity> shared_nodes;
  stk::mesh::get_selected_entities( globally_shared_selector, mesh().buckets( stk::topology::NODE_RANK ), shared_nodes );

  bool done = false;

  const size_t local_num_nodes = fm_nodes.size();
  size_t max_num_nodes = 0;
  stk::all_reduce_max(mesh().parallel(), &local_num_nodes, &max_num_nodes, 1);

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
    if (mesh().parallel_size() > 1)
    {
      for ( auto && shared_node : shared_nodes )
      {
        Fast_Marching_Node * fm_node = get_fm_node(shared_node);
        if (nullptr != fm_node && fm_node->status() != STATUS_UNUSED)
        {
          double & node_dist = *field_data<double>(dRef, fm_node->node());
          node_dist = fm_node->signed_dist()*fm_node->sign();
        }
      }
      stk::mesh::parallel_min(mesh(), {&dRef.field()});

      for ( auto && shared_node : shared_nodes )
      {
        Fast_Marching_Node * fm_node = get_fm_node(shared_node);
        if (nullptr != fm_node && fm_node->status() != STATUS_UNUSED)
        {
          stk::mesh::Entity node = fm_node->node();
          const double min_node_unsigned_dist = *field_data<double>(dRef, node);
          const double fm_node_unsigned_dist = fm_node->signed_dist()*fm_node->sign();
          if(min_node_unsigned_dist < fm_node_unsigned_dist)
          {
            ThrowAssertMsg(fm_node->status() != STATUS_INITIAL || fm_node->status() != STATUS_TRIAL, "Unexpected node to have INITIAL or TRIAL status.");
            fm_node->set_signed_dist(min_node_unsigned_dist*fm_node->sign());
            add_trial_node(*fm_node);
            ++num_locally_updated;
          }
        }
      }
    }

    unsigned num_globally_updated = 0;
    stk::all_reduce_sum(mesh().parallel(), &num_locally_updated, &num_globally_updated, 1);

    done = (num_globally_updated == 0);
  }

  for ( auto && fm_node : fm_nodes )
  {
    if (fm_node.status() == STATUS_TRIAL || fm_node.status() == STATUS_FAR)
    {
      err << "Node " << mesh().identifier(fm_node.node()) << " with status " << fm_node.status() << " with distance " << fm_node.signed_dist() << " did not get updated!\n";
    }
    if (fm_node.status() != STATUS_UNUSED)
    {
      double & node_dist = *field_data<double>(dRef, fm_node.node());
      node_dist = fm_node.signed_dist();
    }
  }

  check_error(err, "Fast Marching Update");
  stk::mesh::field_copy(dRef, olddRef);
}

Fast_Marching_Node * Fast_Marching::get_fm_node(stk::mesh::Entity node)
{
  if (mesh().is_valid(node) && mesh().local_id(node) < fm_nodes.size())
  {
    return &fm_nodes[mesh().local_id(node)];
  }
  else
  {
    return nullptr;
  }
}

void Fast_Marching::initialize(ParallelErrorMessage& err)
{
  stk::mesh::BulkData& stk_mesh = mesh();
  const FieldRef xRef = my_ls.get_coordinates_field();
  const FieldRef dRef = my_ls.get_distance_field();

  const int dim = mesh().mesh_meta_data().spatial_dimension();

  fm_nodes.clear();
  fm_nodes.resize(stk::mesh::count_selected_entities(stk_mesh.mesh_meta_data().universal_part(), stk_mesh.buckets(stk::topology::NODE_RANK)));

  std::vector<stk::mesh::Entity> field_nodes;
  stk::mesh::Selector field_not_ghost = aux_meta().active_not_ghost_selector() & my_selector & stk::mesh::selectField(dRef);
  stk::mesh::get_selected_entities( field_not_ghost, stk_mesh.buckets( stk::topology::NODE_RANK ), field_nodes );

  for ( auto&& node : field_nodes )
  {
    const double * curr_node_dist = field_data<double>(dRef, node);
    ThrowAssert(nullptr != curr_node_dist);

    Vector3d coords = Vector3d::ZERO;
    const double * xptr = field_data<double>(xRef, node);
    ThrowAssert(nullptr != xptr);
    for ( int d = 0; d < dim; d++ )
    {
      coords[d] = xptr[d];
    }

    Fast_Marching_Node * fm_node = get_fm_node(node);
    ThrowAssert(nullptr != fm_node);
    *fm_node = Fast_Marching_Node(node,STATUS_FAR,LevelSet::sign(*curr_node_dist) * std::numeric_limits<double>::max(),LevelSet::sign(*curr_node_dist),coords);
  }

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

  {
    // Initialize using method #5 (element rescaling)
    std::vector<stk::mesh::Entity> field_elems;
    stk::mesh::get_selected_entities( field_not_ghost, stk_mesh.buckets( stk::topology::ELEMENT_RANK ), field_elems );
    for (auto&& elem : field_elems)
    {
      if (have_crossing(elem))
      {
        const double speed = my_ls.get_time_of_arrival_speed(elem, err);
        initialize_element(elem, speed);
      }
    }

    if (mesh().parallel_size() > 1)
    {
      stk::mesh::Selector globally_shared_selector = my_selector & mesh().mesh_meta_data().globally_shared_part();
      std::vector< stk::mesh::Entity> shared_nodes;
      stk::mesh::get_selected_entities( globally_shared_selector, stk_mesh.buckets( stk::topology::NODE_RANK ), shared_nodes );

      for ( auto && shared_node : shared_nodes )
      {
        Fast_Marching_Node * fm_node = get_fm_node(shared_node);
        if (nullptr != fm_node && fm_node->status() != STATUS_UNUSED)
        {
          double & node_dist = *field_data<double>(dRef, fm_node->node());
          node_dist = fm_node->signed_dist()*fm_node->sign();
        }
      }
      stk::mesh::parallel_min(mesh(), {&dRef.field()});

      for ( auto && shared_node : shared_nodes )
      {
        Fast_Marching_Node * fm_node = get_fm_node(shared_node);
        if (nullptr != fm_node && fm_node->status() != STATUS_UNUSED)
        {
          stk::mesh::Entity node = fm_node->node();
          const double min_node_unsigned_dist = *field_data<double>(dRef, node);
          const double fm_node_unsigned_dist = fm_node->signed_dist()*fm_node->sign();
          fm_node->set_signed_dist(min_node_unsigned_dist*fm_node->sign());
          if (min_node_unsigned_dist < fm_node_unsigned_dist)
          {
            fm_node->set_status(STATUS_INITIAL);
          }
        }
      }
    }
  }
}

bool
Fast_Marching::have_crossing(const stk::mesh::Entity & elem) const
{
  const FieldRef dRef = my_ls.get_distance_field();
  const unsigned npe = mesh().bucket(elem).topology().num_nodes();
  ThrowAssert(npe > 0);

  const stk::mesh::Entity * elem_nodes = mesh().begin(elem, stk::topology::NODE_RANK);
  const double * dist0 = field_data<double>(dRef, elem_nodes[0]);
  ThrowAssert(nullptr != dist0);
  for (unsigned n=1; n<npe; ++n)
  {
    const double * dist = field_data<double>(dRef, elem_nodes[n]);
    ThrowAssert(nullptr != dist);
    if (LevelSet::sign_change(*dist0, *dist))
    {
      return true;
    }
  }
  return false;
}

static std::function<const Vector3d &(stk::mesh::Entity)> build_get_fm_node_coordinates(Fast_Marching * fm)
{
  return [fm](stk::mesh::Entity node) -> const Vector3d &
    {
      Fast_Marching_Node * fm_node = fm->get_fm_node(node);
      ThrowAssert(fm_node);
      return fm_node->coords();
    };
}

void
Fast_Marching::initialize_element(const stk::mesh::Entity & elem, const double speed)
{
  // Still another way to initialize fast marching.
  // Here we go to each cut element and find the current distance gradient.
  // By comparing this to the desired gradient, we rescale each element.  The nodal
  // distance is set to the minimum (magnitude) for each of the rescaled elements that
  // support the node.
  const FieldRef dRef = my_ls.get_distance_field();
  const stk::mesh::Entity * elem_nodes = mesh().begin(elem, stk::topology::NODE_RANK);
  const int npe = mesh().bucket(elem).topology().num_nodes();

  auto get_coordinates = build_get_fm_node_coordinates(this);

  const double mag_grad = calculate_gradient_magnitude(npe, elem_nodes, dRef, get_coordinates);

  for (int inode=0; inode<npe; ++inode)
  {
    Fast_Marching_Node * fm_node = get_fm_node(elem_nodes[inode]);
    ThrowAssert(nullptr != fm_node && fm_node->status() != STATUS_UNUSED);
    const double elem_node_dist = *field_data<double>(dRef, elem_nodes[inode]) / (mag_grad * speed);
    const int sign = LevelSet::sign(fm_node->signed_dist());
    fm_node->set_signed_dist(sign * std::min(std::abs(fm_node->signed_dist()), std::abs(elem_node_dist)));
    fm_node->set_status(STATUS_INITIAL);
    fm_node->set_sign(sign);
  }
}

void
Fast_Marching::update_neighbors(Fast_Marching_Node & accepted_node, ParallelErrorMessage& err)
{
  const FieldRef dRef = my_ls.get_distance_field();
  const stk::mesh::Selector active_field_not_aura = my_selector & aux_meta().active_not_ghost_selector() & stk::mesh::selectField(dRef);

  stk::mesh::Entity node = accepted_node.node();

  ThrowAssertMsg(STATUS_ACCEPTED == accepted_node.status() || STATUS_INITIAL == accepted_node.status(), "Expected ACCEPTED OR INITIAL status");

  const int dim = mesh().mesh_meta_data().spatial_dimension();
  ThrowAssert(2 == dim || 3 == dim);
  const int npe_dist = (2==dim) ? 3 : 4;
  std::vector<Fast_Marching_Node *> elem_nodes(npe_dist);

  const unsigned num_node_elems = mesh().num_elements(node);
  const stk::mesh::Entity* node_elems = mesh().begin_elements(node);
  for (unsigned node_elem_index=0; node_elem_index<num_node_elems; ++node_elem_index)
  {
    stk::mesh::Entity elem = node_elems[node_elem_index];
    if (!mesh().is_valid(elem) || !active_field_not_aura(mesh().bucket(elem)))
    {
      continue;
    }

    const double speed = my_ls.get_time_of_arrival_speed(elem, err);

    int node_to_update = -1;
    int num_trial = 0;

    const stk::mesh::Entity* nodes = mesh().begin_nodes(elem);
    for ( int i = 0; i < npe_dist; ++i )
    {
      Fast_Marching_Node * fm_nbr = get_fm_node(nodes[i]);
      ThrowAssertMsg(nullptr != fm_nbr && (STATUS_INITIAL == fm_nbr->status() || STATUS_ACCEPTED == fm_nbr->status() || STATUS_FAR == fm_nbr->status() || STATUS_TRIAL == fm_nbr->status()), "Unexpected node status.");
      elem_nodes[i] = fm_nbr;
      bool do_add_trial_node = fm_nbr->status() == STATUS_FAR;
      if (fm_nbr->status() == STATUS_ACCEPTED)
      {
        const double accepted_node_unsigned_dist = accepted_node.signed_dist()*accepted_node.sign();
        const double nbr_unsigned_dist = fm_nbr->signed_dist()*fm_nbr->sign();
        if(nbr_unsigned_dist > accepted_node_unsigned_dist)
        {
          do_add_trial_node = true;
        }
      }
      if (do_add_trial_node)
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
      update_node(elem_nodes, node_to_update, speed);
    }
  }
}

void
Fast_Marching::add_trial_node(Fast_Marching_Node & trial_node)
{
  ThrowAssertMsg(trial_node.status() == STATUS_ACCEPTED || trial_node.status() == STATUS_FAR, "Expected ACCEPTED or FAR when adding trial node");
  trial_nodes.insert(&trial_node);
  trial_node.set_status(STATUS_TRIAL);
}

void
Fast_Marching::update_trial_node(Fast_Marching_Node & trial_node, const double dist)
{
  ThrowAssertMsg(trial_node.status() == STATUS_TRIAL, "Unexpected node status when updating trial node");

  auto it = trial_nodes.find(&trial_node);
  ThrowAssertMsg(it != trial_nodes.end(), "Can't find trial node");

  trial_nodes.erase(it);

  trial_node.set_signed_dist(dist);
  trial_nodes.insert(&trial_node);
}

void
Fast_Marching::update_node(std::vector<Fast_Marching_Node *> & elem_nodes, int node_to_update, const double speed)
{
  // update distance
  double dist = std::numeric_limits<double>::max();
  const int npe_dist = elem_nodes.size();

  if (3 == npe_dist)
  {
    dist = update_triangle(elem_nodes, node_to_update, speed);
  }
  else if (4 == npe_dist)
  {
    dist = update_tetrahedron(elem_nodes, node_to_update, speed);
  }
  else
  {
    ThrowAssertMsg(false, "Unexpected number of nodes per element: " << npe_dist);
  }

  Fast_Marching_Node & fm_node = *elem_nodes[node_to_update];
  if (dist*fm_node.sign() < fm_node.signed_dist()*fm_node.sign())
  {
    update_trial_node(fm_node, dist);
  }
}

double
Fast_Marching::update_triangle(std::vector<Fast_Marching_Node *> & elemNodes, int nodeToUpdate, const double speed)
{
  stk::diag::TimeBlock timer__(my_tri_timer);
  const int dim = mesh().mesh_meta_data().spatial_dimension();
  ThrowAssert(2 == dim || 3 == dim);

  const std::array<int,3> lnn = get_oriented_nodes_triangle(nodeToUpdate);
  const std::array<Vector3d,3> x{elemNodes[lnn[0]]->coords(), elemNodes[lnn[1]]->coords(), elemNodes[lnn[2]]->coords()};
  const std::array<double,3> d{elemNodes[lnn[0]]->signed_dist(), elemNodes[lnn[1]]->signed_dist(), elemNodes[lnn[2]]->signed_dist()};
  const int sign = elemNodes[lnn[2]]->sign();
  return eikonal_solve_triangle(x, d, sign, dim, speed);
}

double
Fast_Marching::update_tetrahedron(std::vector<Fast_Marching_Node *> & elemNodes, int nodeToUpdate, const double speed)
{
  stk::diag::TimeBlock timer__(my_tet_timer);

  const std::array<int,4> lnn = get_oriented_nodes_tetrahedron(nodeToUpdate);
  const std::array<Vector3d,4> x{elemNodes[lnn[0]]->coords(), elemNodes[lnn[1]]->coords(), elemNodes[lnn[2]]->coords(), elemNodes[lnn[3]]->coords()};
  const std::array<double,4> d{elemNodes[lnn[0]]->signed_dist(), elemNodes[lnn[1]]->signed_dist(), elemNodes[lnn[2]]->signed_dist(), elemNodes[lnn[3]]->signed_dist()};
  const int sign = elemNodes[lnn[3]]->sign();
  return eikonal_solve_tetrahedron(x, d, sign, speed);
}

} // namespace krino
