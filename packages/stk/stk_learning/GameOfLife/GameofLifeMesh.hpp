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
 *  Created on: Jul 10, 2015
 *      Author: jonchu
 */

#ifndef MESHBUILDER_HPP_
#define MESHBUILDER_HPP_

#include <ostream>                      // for basic_ostream::operator<<
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_topology/topology.hpp>    // for topology
#include <vector>                       // for vector
#include "mpi.h"                        // for ompi_communicator_t
#include "stk_mesh/base/Types.hpp"      // for EntityIdVector, etc
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
namespace stk { namespace mesh { class Part; } }



typedef stk::mesh::Field<int> ScalarIntField;

/*
    How to use:
    Pass in the comm, width, height, and depth (if hex), and also stk::mesh::BulkData::NO_AUTO_AURA
    if you don't want an aura. Basically the only use for this class is to make GameofLife, but
    hey, maybe you can come up with something else.
 */

class GameofLifeMesh
{
public:
    GameofLifeMesh(stk::ParallelMachine comm, stk::topology elemType,
                   unsigned spacialDim,
                stk::mesh::BulkData::AutomaticAuraOption auraOption);

    virtual ~GameofLifeMesh() {}


    // accessor funcitions for test functions and GameofLife
    inline stk::mesh::BulkData& bulk_data();

    inline stk::mesh::MetaData& meta_data();

    inline stk::mesh::Part& active_part();

    inline ScalarIntField& neighbor_field();

    inline ScalarIntField& life_field() const;

    inline stk::topology element_type() const;

    inline stk::ParallelMachine comm() const;

protected:
    //Base class needs to label nodes and assign elem node ids
    stk::mesh::EntityVector m_nodes;
    std::vector<stk::mesh::EntityIdVector> m_elemNodeIds;

    //Declared by base constructor, needed for math
    const unsigned m_numProcs;
    const unsigned m_procRank;

    //declared by derived constructor, needed for math
    unsigned m_elemsOnProc;
    unsigned m_elemProcOffset;

    //basics
    void fill_mesh();

private:
    //basics
    stk::mesh::MetaData m_metaData;
    stk::mesh::BulkData m_bulkData;
    stk::ParallelMachine m_comm;

    //parts and fields
    stk::mesh::Part* m_elemPart;
    stk::mesh::Part* m_activePart;
    ScalarIntField* m_lifeField;
    ScalarIntField* m_activeNeighborField;

    stk::mesh::Part* m_nodeset1;
    stk::mesh::Part* m_nodeset2;
    stk::mesh::Part* m_nodeset3;
    stk::mesh::Part* m_nodeset4;
    stk::mesh::Field<double>* m_distFactors;

    //members
    stk::topology m_elemType;
    stk::mesh::EntityIdVector m_elemIds;

    //constructor
    void declare_fields();

    void put_fields_on_parts();

    void declare_parts();

    void put_parts_in_io();

    //fill_mesh
    virtual void declare_element_nodes_ids()=0;

    void declare_element_ids();

    void create_entities();
    virtual bool only_one_active_proc() = 0;
    void declare_entities();
    virtual void share_nodes_between_processors()=0;

    virtual void number_coordinate_field()=0;

};
inline stk::mesh::MetaData& GameofLifeMesh::meta_data()
{
    return m_metaData;
}
inline stk::mesh::BulkData& GameofLifeMesh::bulk_data()
{
    return m_bulkData;
}
inline ScalarIntField& GameofLifeMesh::neighbor_field()
{
    return *m_activeNeighborField;
}
inline stk::mesh::Part& GameofLifeMesh::active_part()
{
    return *m_activePart;
}
inline ScalarIntField& GameofLifeMesh::life_field() const
{
    return *m_lifeField;
}
inline stk::topology GameofLifeMesh::element_type() const
{
    return m_elemType;
}
inline stk::ParallelMachine GameofLifeMesh::comm() const
{
    return m_comm;
}

class TwoDimGameofLifeMesh : public GameofLifeMesh
{
public:
    TwoDimGameofLifeMesh(stk::ParallelMachine comm, stk::topology elemType,
                              unsigned width, unsigned height,
                              stk::mesh::BulkData::AutomaticAuraOption auraOption);
    virtual ~TwoDimGameofLifeMesh() {}
protected:
    // coordinate field
    stk::mesh::Field<double,stk::mesh::Cartesian2d>* m_nodeCoords;
    // base constructor
    const unsigned m_width;
    const unsigned m_height;
    const unsigned m_rowsPerProc;
    const unsigned m_nodesPerRow;
    //derived constructor
    unsigned m_elemsPerRow;
private:
    //constructor
    virtual void declare_coordinate_field();

    //fill mesh
    virtual void declare_element_nodes_ids()=0;

    virtual bool only_one_active_proc();

    virtual void share_nodes_between_processors();
    virtual void share_top_nodes();
    virtual void share_bottom_nodes();
    void share_node_with_this_id_to_this_processor(unsigned nodeId, unsigned procNum);

    virtual void number_coordinate_field();
    void number_coordinate_field_of_node(unsigned nodeIndex);

};

class TriGameofLifeMesh : public TwoDimGameofLifeMesh
{
public:
    TriGameofLifeMesh(stk::ParallelMachine comm, unsigned width, unsigned
                      rowsPerProc,
                        stk::mesh::BulkData::AutomaticAuraOption auraOption =
                                stk::mesh::BulkData::AUTO_AURA);

    virtual ~TriGameofLifeMesh() {}

private:
    //fill_mesh

    virtual void declare_element_nodes_ids();
    void declare_node_ids_of_two_elements(unsigned index, unsigned initialOffset);
    void declare_node_ids_of_two_first_row_elements(unsigned index, unsigned initialOffset);
    void declare_node_ids_of_two_sucessive_row_elements(unsigned index);
    void declare_node_ids_of_this_element(unsigned index);

};

class QuadGameofLifeMesh : public TwoDimGameofLifeMesh
{
public:
    QuadGameofLifeMesh(stk::ParallelMachine comm, unsigned width, unsigned
                       rowsPerProc,
                       stk::mesh::BulkData::AutomaticAuraOption auraOption =
                            stk::mesh::BulkData::AUTO_AURA);
    virtual ~QuadGameofLifeMesh() {}
private:
    //fill_mesh
    virtual void declare_element_nodes_ids();
    void declare_node_ids_of_element(unsigned index, unsigned initialOffset);
    void declare_first_row_element_nodes(unsigned index, unsigned offset);
    void declare_remaining_element_nodes(unsigned index);

};

class ThreeDimGameofLifeMesh : public GameofLifeMesh
{
public:
    ThreeDimGameofLifeMesh(stk::ParallelMachine comm, stk::topology elemType,
                           unsigned width, unsigned height, unsigned depth,
                           stk::mesh::BulkData::AutomaticAuraOption
                           auraOption);

    virtual ~ThreeDimGameofLifeMesh() {}

protected:
    // coordinate field
    stk::mesh::Field<double,stk::mesh::Cartesian>* m_nodeCoords;

    //base constructor
    const unsigned m_width;
    const unsigned m_height;
    const unsigned m_depth;
    const unsigned m_slicesPerProc;
    const unsigned m_nodeWidth;
    const unsigned m_nodeHeight;
    const unsigned m_nodesPerSlice;

    //derived constructor
    unsigned m_elemsPerSlice;

private:
    //constructor
    virtual void declare_coordinate_field();

    //fill mesh
    virtual void declare_element_nodes_ids()=0;

    virtual bool only_one_active_proc();

    virtual void share_nodes_between_processors();
    void share_nodes_in_back();
    void share_nodes_in_front();
    void share_node_with_this_id_to_this_processor(unsigned nodeId,
                                                   unsigned procNum);

    virtual void number_coordinate_field();
    void number_coordinate_field_for_node(unsigned nodeIndex);

};

class HexGameofLifeMesh : public ThreeDimGameofLifeMesh
{
public:
    HexGameofLifeMesh(stk::ParallelMachine comm, unsigned width, unsigned height,
                      unsigned depth,
                      stk::mesh::BulkData::AutomaticAuraOption auraOption =
                           stk::mesh::BulkData::AUTO_AURA);

    virtual ~HexGameofLifeMesh() {}

private:
    //fill_mesh
    virtual void declare_element_nodes_ids();
    void declare_node_ids_of_element(unsigned index, unsigned offset);
    void declare_first_slice_element_node_ids(stk::mesh::EntityIdVector& V,
                                              unsigned index,
                                              unsigned offset);
    void declare_remaining_element_node_ids(stk::mesh::EntityIdVector& newer,
                                            stk::mesh::EntityIdVector& older,
                                            unsigned index);
};
#endif /* MeshBuilder.hpp */
