// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/*
 * GameofLife.hpp
 *
 *  Created on: Jul 25, 2015
 *      Author: Jonathan Chu
 */
#ifndef GAMEOFLIFE_HPP
#define GAMEOFLIFE_HPP

#include <stddef.h>                     // for size_t
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_topology/topology.hpp>    // for topology
#include <string>                       // for string
#include <unordered_map>                // for unordered_map
#include <unordered_set>                // for unordered_set
#include "stk_mesh/base/HashEntityAndEntityKey.hpp"            // for hash
#include "stk_mesh/base/Types.hpp"      // for EntityIdVector, EntityId, etc
class GameofLifeMesh;
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { struct Entity; } }




typedef stk::mesh::Field<int> ScalarIntField;

/*
 * how to use:
 * Make a MeshBuilder, and then pass the mesh and whatever name you want into the constructor.
 * Pass some IDs to activate, and run the game of life for however many steps. It's that easy.
 * Just remember to pass it a MeshBuilder that had an aura.
 */
class GameofLife
{
public:
    GameofLife(GameofLifeMesh& Mesh, std::string meshName);

    virtual ~GameofLife() {}

    //basics
    virtual void activate_these_ids(stk::mesh::EntityIdVector& elemIdsToActivate);

    void run_game_of_life(int numSteps);

    virtual void write_mesh();

    //test functions
    inline unsigned get_num_elems_on_proc() const;

    unsigned get_num_active_elements() const;

    bool are_these_ids_active(const stk::mesh::EntityIdVector& elemIds) const;

protected:
    //big stuff
    stk::mesh::MetaData& m_metaData;
    stk::mesh::BulkData& m_bulkData;
    stk::mesh::EntityVector m_elements;

    // housekeeping
    unsigned m_elemsOnProc;

    // game of life
    ScalarIntField& m_lifeField;

    //useful stuff
    void update_this_element(stk::mesh::Entity elem);
    bool element_is_local(stk::mesh::Entity elem) const;

private:
    // io
    stk::io::StkMeshIoBroker m_stkIo;
    size_t m_fileHandler;
    std::string m_meshName;
    int m_time;

    // members
    const stk::topology m_elemType;

    //neighbors
    ScalarIntField& m_activeNeighborField;

    std::unordered_map<stk::mesh::Entity, std::unordered_set<stk::mesh::Entity,
    std::hash<stk::mesh::Entity>>, std::hash<stk::mesh::Entity>>
    m_neighborSets;

    //constructor
    void get_elements();

    void write_output_mesh();

    void fill_neighbor_sets_of_elements();
    void fill_neighbor_set_of_this_element(stk::mesh::Entity elem);
    void add_this_nodes_elements_to_set(stk::mesh::Entity node, stk::mesh::Entity elem);

    //activate elements
    virtual void activate_each_element_id_in_vector(stk::mesh::EntityIdVector& elemIdsToActivate)=0;

    //game of life
    void run_step_of_game_of_life();
    virtual void communicate_data() = 0;

    void update_neighbor_fields();
    void update_neighbor_field_of_this_element(stk::mesh::Entity elem);
    int  get_num_active_neighbors_of_this_element(stk::mesh::Entity elem);
    virtual bool element_is_active(stk::mesh::Entity elem) const = 0;

    virtual void update_each_element() = 0;

    void update_tri(stk::mesh::Entity elem);
    void update_quad(stk::mesh::Entity elem);
    void update_hex(stk::mesh::Entity elem);
    virtual void activate_element(stk::mesh::Entity elem) = 0;
    virtual void deactivate_element(stk::mesh::Entity elem) = 0;

    void write_output_step();
    void put_all_nodes_in_nodeset();
};

inline unsigned GameofLife::get_num_elems_on_proc() const
{
   return m_elemsOnProc;
}

class PartGameofLife : public GameofLife
{
public:
    PartGameofLife(GameofLifeMesh& Mesh, std::string meshName);
    virtual ~PartGameofLife(){}
private:
    //members
    stk::mesh::Part& m_activePart;
    stk::mesh::PartVector m_active;
    stk::mesh::PartVector m_empty;

    //activate element ids
    void activate_each_element_id_in_vector(stk::mesh::EntityIdVector& elemIdsToActivate);
    void activate_this_element_id(stk::mesh::EntityId elemId);

    //game of life
    virtual void communicate_data();

    virtual bool element_is_active(stk::mesh::Entity elem) const;
    virtual void update_each_element();
    virtual void activate_element(stk::mesh::Entity elem);
    virtual void deactivate_element(stk::mesh::Entity elem);
};

class FieldGameofLife : public GameofLife
{
public:
    FieldGameofLife(GameofLifeMesh& Mesh, std::string meshName);
    virtual ~FieldGameofLife(){}

private:
    //activate element ids
    void activate_each_element_id_in_vector(stk::mesh::EntityIdVector& elemIdsToActivate);
    void activate_this_element_id(stk::mesh::EntityId elemId);

    //game of life
    virtual void communicate_data();

    virtual bool element_is_active(stk::mesh::Entity elem) const;
    virtual void update_each_element();
    virtual void activate_element(stk::mesh::Entity elem);
    virtual void deactivate_element(stk::mesh::Entity elem);
};

class TodstwdGameOfLife : public FieldGameofLife
{
public:
    TodstwdGameOfLife(GameofLifeMesh& mesh, std::string meshName)
    : FieldGameofLife(mesh, meshName) { }
    virtual ~TodstwdGameOfLife() { }
    virtual void write_mesh() override;
};

#endif /*GameofLife.hpp*/
