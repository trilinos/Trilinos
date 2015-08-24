/*
 * MeshBuilder.h
 *
 *  Created on: Aug 17, 2015
 *      Author: jonchu
 */

#ifndef MESHPRINTER_MESHPRINTER_H_
#define MESHPRINTER_MESHPRINTER_H_

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

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include "CoordinateSets.h"

class MeshPrinter
{
public:
    MeshPrinter(stk::ParallelMachine comm, std::string name, stk::topology::topology_t m_elemType,
                int spacialDim);

    virtual ~MeshPrinter() {}

    //public methods
    void commit_meta();

    void begin_modification();

    void end_modification();

    void write_mesh();

    //accessor
    inline stk::mesh::MetaData& meta_data();

    inline stk::mesh::BulkData& bulk_data();

    inline int num_procs() const;

protected:
    stk::mesh::Entity create_element(stk::mesh::EntityId elemId,
                                     const stk::mesh::EntityIdVector& nodeIds,
                                     int chosenProc);

private:
   stk::mesh::MetaData m_metaData;
   stk::mesh::BulkData m_bulkData;
   const int m_numProcs;
   const int m_procRank;

   const std::string m_name;
   unsigned m_time;

   stk::mesh::Part* m_elemPart;

   std::unordered_map<stk::mesh::EntityId, std::unordered_set<int>> m_nodeIdToSharingProcs;

   std::unordered_set<stk::mesh::EntityId> m_usedElemIds;

   //create element
   stk::mesh::Entity generate_element(stk::mesh::EntityId elemId,
                                      const stk::mesh::EntityIdVector& nodeIds);

       bool elem_id_not_used(stk::mesh::EntityId elemId);

       void share_nodes(const stk::mesh::EntityIdVector& nodeIds);

   void share_shared_nodes(const stk::mesh::EntityIdVector& nodeIds, int chosenProc);

       void update_node_to_processor_map(const stk::mesh::EntityIdVector& nodeIds, int proc);

};

//accessor
inline stk::mesh::MetaData& MeshPrinter::meta_data()
{
    return m_metaData;
}

inline stk::mesh::BulkData& MeshPrinter::bulk_data()
{
    return m_bulkData;
}

inline int MeshPrinter::num_procs() const
{
    return m_numProcs;
}

class QuadMeshPrinter : public MeshPrinter
{
public:
    QuadMeshPrinter(stk::ParallelMachine comm, std::string name);

    virtual ~QuadMeshPrinter() {}

    void create_element(unsigned xCoord, unsigned yCoord, int chosenProc = 0);

    //test functions
    double node_x_coord(stk::mesh::Entity node) const;

    double node_y_coord(stk::mesh::Entity node) const;

private:
    stk::mesh::Field<double, stk::mesh::Cartesian2d>* m_coordinates;

    void label_node_coordinates(const ElemCoordPair& elemCoords);

        void label_coordinates(stk::mesh::Entity node, unsigned nodeX, unsigned nodeY);

};

class HexMeshPrinter : public MeshPrinter
{
public:
    HexMeshPrinter(stk::ParallelMachine comm, std::string name);

    virtual ~HexMeshPrinter(){}

    void create_element(unsigned xCoord, unsigned yCoord, unsigned zCoord, int chosenProc = 0);

    //test functions
    double node_x_coord(stk::mesh::Entity node) const;

    double node_y_coord(stk::mesh::Entity node) const;

    double node_z_coord(stk::mesh::Entity node) const;

private:
    stk::mesh::Field<double, stk::mesh::Cartesian3d>* m_coordinates;

    void label_node_coordinates(const ElemCoordTriple& elemCoords);

        void label_coordinates(stk::mesh::Entity node, unsigned nodeX, unsigned nodeY, unsigned nodeZ);


};


#endif /* MESHPRINTER_MESHPRINTER_H_ */
