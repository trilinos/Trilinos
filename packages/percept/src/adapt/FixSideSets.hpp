// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef FixSideSets_hpp
#define FixSideSets_hpp

#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <percept/Percept.hpp>

#include <adapt/Refiner.hpp>
#include <adapt/RefinerUtil.hpp>
#include <percept/PerceptMesh.hpp>

#include <percept/MeshUtil.hpp>
#include <adapt/AdaptedMeshVerifier.hpp>

#define DEBUG_FSS 0

namespace percept {


  class FixSideSets {
  public:
    Refiner *m_refiner;
    PerceptMesh& m_eMesh;
    stk::mesh::PartVector& m_excludeParts;
    SidePartMap& m_side_part_map;
    std::string m_geomFile;
    bool m_avoidFixSideSetChecks;
    RefinerSelector *m_buildSideSetSelector;
    bool m_doProgress;
    bool m_debug{false};

    FixSideSets(Refiner *ref, PerceptMesh& eMesh, stk::mesh::PartVector& excludeParts, SidePartMap& side_part_map, const std::string& geomFile, bool avoidFixSideSetChecks, RefinerSelector *sel = 0, bool doProgress=false);

    void fix_permutation(SetOfEntities& side_set);
    bool connect(stk::mesh::Entity side, bool& valid_side_part_map, SetOfEntities* avoid_elems, bool onlyPosPerm=false);

    // if the element (element) has a side that matches  the given side (side), connect them but first delete old connections
    std::pair<bool,bool> connectSidesForced(stk::mesh::Entity element, stk::mesh::Entity side, bool& valid_side_part_map, stk::mesh::ConnectivityOrdinal *k_element_side, bool onlyPosPerm = false);
    void disconnect_entity(stk::mesh::Entity entity);

    void doProgressPrint(PerceptMesh& eMesh, const std::string& msg);

    void delete_unattached_sides(SetOfEntities& side_set, SetOfEntities *avoid_sides);
    bool bucket_acceptable(stk::mesh::Bucket& bucket, stk::mesh::EntityRank rank);
    void build_side_set(SetOfEntities& side_set, bool only_roots = false);
    void reconnect_sides(SetOfEntities& side_set, SetOfEntities *avoid_elems, bool onlyPosPerm);
    void check_connect(SetOfEntities& side_set, SetOfEntities *avoid_elems);
    void end_begin(const std::string& msg="");

    void move_sides_to_correct_surfaces();
    void move_side_to_correct_surface(stk::mesh::Part& surface, stk::mesh::Entity side, stk::mesh::Entity volume);

    std::pair<std::string, bool> get_new_sideset_part_name(const std::string& surfaceName, stk::mesh::Entity side, stk::mesh::Entity volume);
    void fill_change_parts(stk::mesh::Part& surface,
                           stk::mesh::Entity side, stk::mesh::Entity volume,
                           std::vector<stk::mesh::Part*>& add_parts, std::vector<stk::mesh::Part*>& remove_parts);

    // fast reconnector
    void fix_side_sets_2(bool allow_not_found, SetOfEntities *avoid_elems, SetOfEntities *avoid_sides, const std::string& msg);
  };


}
#endif
