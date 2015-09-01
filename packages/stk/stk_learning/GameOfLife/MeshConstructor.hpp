/*
 * MeshBuilder.h
 *
 *  Created on: Jul 10, 2015
 *      Author: jonchu
 */

#ifndef MESHBUILDER_HPP_
#define MESHBUILDER_HPP_

#include <vector>
#include <algorithm>
#include <stdlib.h>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

typedef stk::mesh::Field<int> ScalarIntField;

/*
    How to use:
    Pass in the comm, width, height, and depth (if hex), and also stk::mesh::BulkData::NO_AUTO_AURA
    if you don't want an aura. Basically the only use for this class is to make GameofLife, but
    hey, maybe you can come up with something else.
 */

class MeshConstructor
{
public:
    MeshConstructor(stk::ParallelMachine comm, stk::topology elemType, unsigned spacialDim,
                stk::mesh::BulkData::AutomaticAuraOption auraOption);

    virtual ~MeshConstructor() {}


    // accessor funcitions for test functions and GameofLife
    inline stk::mesh::BulkData* bulk_data();

    inline stk::mesh::MetaData* meta_data();

    inline stk::mesh::Part* active_part();

    inline ScalarIntField* neighbor_field();

    inline ScalarIntField* life_field() const;

    inline stk::mesh::EntityVector* elements();

    inline stk::mesh::EntityIdVector* element_ids();

    inline stk::topology element_type() const;

    inline unsigned num_elems_on_proc() const;

    inline unsigned proc_rank() const;

    inline unsigned num_procs() const;

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
inline stk::mesh::MetaData* MeshConstructor::meta_data()
{
    return &m_metaData;
}
inline stk::mesh::BulkData* MeshConstructor::bulk_data()
{
    return &m_bulkData;
}
inline ScalarIntField* MeshConstructor::neighbor_field()
{
    return m_activeNeighborField;
}
inline stk::mesh::Part* MeshConstructor::active_part()
{
    return m_activePart;
}
inline ScalarIntField* MeshConstructor::life_field() const
{
    return m_lifeField;
}
inline stk::topology MeshConstructor::element_type() const
{
    return m_elemType;
}
inline stk::ParallelMachine MeshConstructor::comm() const
{
    return m_comm;
}

class TwoDimensionalMeshConstructor : public MeshConstructor
{
public:
    TwoDimensionalMeshConstructor(stk::ParallelMachine comm, stk::topology elemType,
                              unsigned width, unsigned height,
                              stk::mesh::BulkData::AutomaticAuraOption auraOption);
    virtual ~TwoDimensionalMeshConstructor() {}
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

class TriangleMeshConstructor : public TwoDimensionalMeshConstructor
{
public:
    TriangleMeshConstructor(stk::ParallelMachine comm, unsigned width, unsigned rowsPerProc,
                        stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::AUTO_AURA);

    virtual ~TriangleMeshConstructor() {}

private:
    //fill_mesh

    virtual void declare_element_nodes_ids();
    void declare_node_ids_of_two_elements(unsigned index, unsigned initialOffset);
    void declare_node_ids_of_two_first_row_elements(unsigned index, unsigned initialOffset);
    void declare_node_ids_of_two_sucessive_row_elements(unsigned index);
    void declare_node_ids_of_this_element(unsigned index);

};

class QuadMeshConstructor : public TwoDimensionalMeshConstructor
{
public:
    QuadMeshConstructor(stk::ParallelMachine comm, unsigned width, unsigned rowsPerProc,
                    stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::AUTO_AURA);
    virtual ~QuadMeshConstructor() {}
private:
    //fill_mesh
    virtual void declare_element_nodes_ids();
    void declare_node_ids_of_element(unsigned index, unsigned initialOffset);
    void declare_first_row_element_nodes(unsigned index, unsigned offset);
    void declare_remaining_element_nodes(unsigned index);

};

class ThreeDimensionalMeshConstructor : public MeshConstructor
{
public:
    ThreeDimensionalMeshConstructor(stk::ParallelMachine comm, stk::topology elemType, unsigned width,
                                unsigned height, unsigned depth,
                                stk::mesh::BulkData::AutomaticAuraOption auraOption);

    virtual ~ThreeDimensionalMeshConstructor() {}

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
    void share_node_with_this_id_to_this_processor(unsigned nodeId, unsigned procNum);

    virtual void number_coordinate_field();
    void number_coordinate_field_for_node(unsigned nodeIndex);

};

class HexMeshConstructor : public ThreeDimensionalMeshConstructor
{
public:
    HexMeshConstructor(stk::ParallelMachine comm, unsigned width, unsigned height, unsigned depth,
                   stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::AUTO_AURA);

    virtual ~HexMeshConstructor() {}

private:
    //fill_mesh
    virtual void declare_element_nodes_ids();
    void declare_node_ids_of_element(unsigned index, unsigned offset);
    void declare_first_slice_element_node_ids(stk::mesh::EntityIdVector& V, unsigned index,
                                              unsigned offset);
    void declare_remaining_element_node_ids(stk::mesh::EntityIdVector& newer,
                                            stk::mesh::EntityIdVector& older, unsigned index);
};
#endif /* MeshBuilder.hpp */
