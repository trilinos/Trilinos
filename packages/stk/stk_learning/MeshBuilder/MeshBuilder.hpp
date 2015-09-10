/*
 * MeshBuilder.h
 *
 *  Created on: Aug 17, 2015
 *      Author: jonchu
 */

#ifndef _MESHBUILDER_HPP_
#define _MESHBUILDER_HPP_

#include <iostream>
#include <vector>
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
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <stk_io/StkMeshIoBroker.hpp>

#include "../../../TPLs_src/Trilinos/packages/stk/stk_learning/MeshBuilder/CoordinateSets.hpp"

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

    //test functions
    inline stk::mesh::MetaData& meta_data();

    inline stk::mesh::BulkData& bulk_data();

    inline int num_procs() const;

    inline int spacial_dim() const;

    unsigned num_elems() const;

protected:
    stk::mesh::Entity create_element(stk::mesh::EntityId elemId,
                                     const stk::mesh::EntityIdVector& nodeIds,
                                     int chosenProc);

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

   //cojknstructor
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

    void fill_area_on_proc(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned yUpper,
                           int chosenProc);

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
