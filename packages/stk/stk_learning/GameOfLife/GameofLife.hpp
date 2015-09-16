/*
 * GameofLife.hpp
 *
 *  Created on: Jul 25, 2015
 *      Author: Jonathan Chu
 */
#ifndef GAMEOFLIFE_HPP
#define GAMEOFLIFE_HPP

#include <vector>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <stdlib.h>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include "EntityKeyHash.hpp"
#include "GameofLifeMesh.hpp"

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

#endif /*GameofLife.hpp*/
