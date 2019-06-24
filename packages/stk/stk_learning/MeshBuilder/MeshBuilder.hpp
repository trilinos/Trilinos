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
 * MeshBuilder.h
 *
 *  Created on: Aug 17, 2015
 *      Author: jonchu
 */

#ifndef _MESHBUILDER_HPP_
#define _MESHBUILDER_HPP_

#include <iostream>                     // for basic_ostream::operator<<
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include <unordered_map>                // for unordered_map
#include <unordered_set>                // for unordered_set
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for EntityId, EntityIdVector
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { struct Cartesian2d; } }
namespace stk { namespace mesh { struct Cartesian3d; } }
struct ElemCoordPair;
struct ElemCoordTriple;




typedef stk::mesh::Field<int> ScalarIntField;

class MeshBuilder
{
public:
    MeshBuilder(stk::ParallelMachine comm, std::string name, stk::topology::topology_t m_elemType,
                int spacialDim);

    virtual ~MeshBuilder() {}

    //public methods
    void commit_meta();

    void begin_modification();

    void end_modification();

    void write_mesh();

    //gol
    void create_life_and_neighbor_fields(ScalarIntField*& lifeField, ScalarIntField*& neighborField);

    //test functions
    inline stk::mesh::MetaData& meta_data();

    inline stk::mesh::BulkData& bulk_data();

    inline int num_procs() const;

    inline int proc_rank() const;

    inline int spacial_dim() const;

    unsigned num_elems() const;

protected:
    stk::mesh::Entity create_element(stk::mesh::EntityId elemId,
                                     const stk::mesh::EntityIdVector& nodeIds,
                                     int chosenProc);

    void remove_element(stk::mesh::EntityId elemId);

private:
   stk::mesh::MetaData m_metaData;
   stk::mesh::BulkData m_bulkData;
   const int m_spacialDim;
   const int m_numProcs;
   const int m_procRank;

   const std::string m_name;
   unsigned m_time;

   stk::mesh::Part* m_elemPart;

   std::unordered_map<stk::mesh::EntityId, std::unordered_set<int>> m_nodeIdToSharingProcs;

   std::unordered_set<stk::mesh::EntityId> m_usedElemIds;

   //constructor
   virtual void declare_coordinates() = 0;

   //create element
   stk::mesh::Entity generate_element(stk::mesh::EntityId elemId,
                                      const stk::mesh::EntityIdVector& nodeIds);

       bool elem_id_not_used(stk::mesh::EntityId elemId);

       void share_nodes(const stk::mesh::EntityIdVector& nodeIds);

   void share_shared_nodes(const stk::mesh::EntityIdVector& nodeIds, int chosenProc);

       void update_node_to_processor_map(const stk::mesh::EntityIdVector& nodeIds, int proc);

};

//accessor
inline stk::mesh::MetaData& MeshBuilder::meta_data()
{
    return m_metaData;
}
inline stk::mesh::BulkData& MeshBuilder::bulk_data()
{
    return m_bulkData;
}
inline int MeshBuilder::num_procs() const
{
    return m_numProcs;
}
inline int MeshBuilder::proc_rank() const
{
    return m_procRank;
}
inline int MeshBuilder::spacial_dim() const
{
    return m_spacialDim;
}

class QuadMeshBuilder : public MeshBuilder
{
public:
    QuadMeshBuilder(stk::ParallelMachine comm, std::string name);

    virtual ~QuadMeshBuilder() {}

    void create_element(unsigned xCoord, unsigned yCoord, int chosenProc = 0);

    void fill_area(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned yUpper);

    void fill_area_randomly(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned yUpper);

    void fill_area_on_proc(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned yUpper,
                           int chosenProc);

    void fill_area_with_layers(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned yUpper);

    void remove_element(unsigned xCoord, unsigned yCoord);

    void remove_area(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned yUpper);

    //test functions
    double node_x_coord(stk::mesh::Entity node) const;

    double node_y_coord(stk::mesh::Entity node) const;

private:
    stk::mesh::Field<double, stk::mesh::Cartesian2d>* m_coordinates;

    virtual void declare_coordinates();

    void label_node_coordinates(const ElemCoordPair& elemCoords);

        void label_coordinates(stk::mesh::Entity node, unsigned nodeX, unsigned nodeY);

};

class HexMeshBuilder : public MeshBuilder
{
public:
    HexMeshBuilder(stk::ParallelMachine comm, std::string name);

    virtual ~HexMeshBuilder(){}

    void create_element(unsigned xCoord, unsigned yCoord, unsigned zCoord, int chosenProc = 0);

    void fill_area(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned yUpper, unsigned
                   zLower, unsigned zUpper);

    void fill_area_randomly(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned yUpper,
                            unsigned zLower, unsigned zUpper);

    void fill_area_on_proc(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned yUpper,
                           unsigned zLower, unsigned zUpper, int chosenProc);

    void fill_area_with_layers(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned yUpper,
                               unsigned zLower, unsigned zUpper);

    void remove_element(unsigned xCoord, unsigned yCoord, unsigned zCoord);

    void remove_area(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned yUpper, unsigned
                     zLower, unsigned zUpper);

    //test functions
    double node_x_coord(stk::mesh::Entity node) const;

    double node_y_coord(stk::mesh::Entity node) const;

    double node_z_coord(stk::mesh::Entity node) const;

private:
    stk::mesh::Field<double, stk::mesh::Cartesian3d>* m_coordinates;

    virtual void declare_coordinates();

    void label_node_coordinates(const ElemCoordTriple& elemCoords);

        void label_coordinates(stk::mesh::Entity node, unsigned nodeX, unsigned nodeY, unsigned nodeZ);


};


#endif /* R_MESHPRINTER_HPP_ */
