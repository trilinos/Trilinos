// Copyright(C) 2014
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#ifndef IOSS_Iogn_DashSurfaceMesh_h
#define IOSS_Iogn_DashSurfaceMesh_h

#include <generated/Iogn_GeneratedMesh.h>  // for GeneratedMesh
#include <stddef.h>                     // for size_t
#include <stdint.h>                     // for int64_t
#include <exception>                    // for exception
#include <string>                       // for string
#include <utility>                      // for pair
#include <vector>                       // for vector

namespace Iogn
{

enum
{
    INVALID                 = -1,
    NUMBER_OF_SURFACES      =  2,
    SPATIAL_DIMENSION       =  3,
    NUM_NODES_PER_QUAD_FACE =  4
};

struct SharedNode
{
    SharedNode() : nodeId(-1), procId(-1) {}
    int nodeId;
    int procId;
};

enum Topology {
    Shell4 = 4,
    Hex8   = 8
};

inline std::string getTopologyName(Topology topology)
{
    switch(topology)
    {
        case Shell4:
            return std::string("shell4");
        case Hex8:
            return std::string("hex8");
    }
    throw std::exception();
}

struct ExodusData {
    std::vector<double> coordinates;
    const std::vector< std::vector<int> > elementBlockConnectivity;
    const std::vector<int> globalNumberOfElementsInBlock;
    const std::vector<int> localNumberOfElementsInBlock;
    const std::vector<Topology> blockTopologicalData;

    const int globalNumberOfNodes;

    const std::vector<int> globalIdsOfLocalElements;
    const std::vector<int> globalIdsOfLocalNodes;

    std::vector<SharedNode> *sharedNodes;

    // A sideset' is basically an exodus sideset.  A
    // sideset has a list of elements and a corresponding local
    // element side (1-based) The side id is: side_id =
    // 10*element_id + local_side_number This assumes that all
    // sides in a sideset are boundary sides.
    std::vector< std::vector<int> > sidesetConnectivity;
    std::vector< std::vector<std::string> > sidesetTouchingBlocks;

    ExodusData(const std::vector<double> coords, const std::vector< std::vector<int> > elemBlockConnectivity,
               const std::vector<int> globalNumOfElemsInBlock, const std::vector<int> localNumOfElemsInBlock,
               const std::vector<Topology> blockTopoData, int globalNumNodes, const std::vector<int> globalIdsOfLocalElems,
               const std::vector<int> globalIdsLocalNodes,
               const std::vector< std::vector<int> > & sidesetConn = std::vector< std::vector<int> >(),
               const std::vector< std::vector<std::string> > & sidesetBlocks = std::vector< std::vector<std::string> >())
    : coordinates(coords), elementBlockConnectivity(elemBlockConnectivity),
      globalNumberOfElementsInBlock(globalNumOfElemsInBlock), localNumberOfElementsInBlock(localNumOfElemsInBlock),
      blockTopologicalData(blockTopoData), globalNumberOfNodes(globalNumNodes), globalIdsOfLocalElements(globalIdsOfLocalElems),
      globalIdsOfLocalNodes(globalIdsLocalNodes), sharedNodes(0), sidesetConnectivity(sidesetConn),
      sidesetTouchingBlocks(sidesetBlocks)
    {

    }
};

struct DashSurfaceData
{
    const std::vector<double> &coordinates;
    const std::vector<int> &surfaceAConnectivity;
    const std::vector<int> &surfaceBConnectivity;

    int globalNumberOfNodes;
    int globalNumberOfElements;

    int globalNumberOfElementsSurface1;
    int globalNumberOfElementsSurface2;

    std::vector<int> globalIdsOfLocalElements;
    std::vector<int> globalIdsOfLocalNodes;

    std::vector<SharedNode> *sharedNodes;

    DashSurfaceData(const std::vector<double> &coords, const std::vector<int> &connectivity1, const std::vector<int> &connectivity2)
    : coordinates(coords), surfaceAConnectivity(connectivity1), surfaceBConnectivity(connectivity2),
      sharedNodes(0)
    {
        this->setSerialDefaults();
    }

  private:

    void setSerialDefaults()
    {
        globalNumberOfNodes = coordinates.size()/SPATIAL_DIMENSION;

        globalNumberOfElementsSurface1 = surfaceBConnectivity.size()/NUM_NODES_PER_QUAD_FACE;
        globalNumberOfElementsSurface2 = surfaceAConnectivity.size()/NUM_NODES_PER_QUAD_FACE;
        globalNumberOfElements = globalNumberOfElementsSurface1 + globalNumberOfElementsSurface2;

        globalIdsOfLocalElements.resize(globalNumberOfElements);
        globalIdsOfLocalNodes.resize(globalNumberOfNodes);

        for (size_t i=0; i<globalIdsOfLocalElements.size();i++)
        {
            globalIdsOfLocalElements[i] = i+1;
        }

        for (size_t i=0; i<globalIdsOfLocalNodes.size();i++)
        {
            globalIdsOfLocalNodes[i] = i+1;
        }
    }
};

class DashSurfaceMesh : public GeneratedMesh
{
public:

    explicit DashSurfaceMesh(DashSurfaceData &dashSurfaceData) : mDashSurfaceData(dashSurfaceData) { }

    virtual ~DashSurfaceMesh() { }

    virtual int64_t node_count() const;
    virtual int64_t node_count_proc() const;

    virtual int64_t element_count() const;
    virtual int64_t element_count(int64_t block_number) const;
    virtual int64_t element_count_proc() const;
    virtual int64_t element_count_proc(int64_t block_number) const;

    virtual int64_t block_count() const;

    virtual int64_t nodeset_count() const;
    virtual int64_t nodeset_node_count_proc(int64_t id) const;

    virtual int64_t sideset_count() const;
    virtual int64_t sideset_side_count_proc(int64_t id) const;

    virtual int64_t communication_node_count_proc() const;

    virtual void coordinates(double *coord) const;
    virtual void coordinates(std::vector<double> &coord) const;
    virtual void coordinates(int component, std::vector<double> &xyz) const;
    virtual void coordinates(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z) const;

    virtual void connectivity(int64_t block_number, int* connect) const;

    virtual std::pair<std::string, int> topology_type(int64_t block_number) const;

    virtual void sideset_elem_sides(int64_t setId, std::vector<int64_t> &elem_sides) const;

    virtual void nodeset_nodes(int64_t nset_id, std::vector<int64_t> &nodes) const;

    virtual void node_communication_map(std::vector<int64_t> &map, std::vector<int> &proc);

    virtual void node_map(std::vector<int> &map) const;
    virtual void node_map(std::vector<int64_t> &map) const;

    virtual void element_map(int64_t block_number, std::vector<int> &map) const;
    virtual void element_map(int64_t block_number, std::vector<int64_t> &map) const;
    virtual void element_map(std::vector<int64_t> &map) const;
    virtual void element_map(std::vector<int> &map) const;

private:
    DashSurfaceData & mDashSurfaceData;
};

class ExodusMesh : public GeneratedMesh
{
public:

    explicit ExodusMesh(const ExodusData &exodusData);

    virtual ~ExodusMesh() { }

    virtual int64_t node_count() const;
    virtual int64_t node_count_proc() const;

    virtual int64_t element_count() const;
    virtual int64_t element_count(int64_t block_number) const;
    virtual int64_t element_count_proc() const;
    virtual int64_t element_count_proc(int64_t block_number) const;

    virtual int64_t block_count() const;

    virtual int64_t nodeset_count() const;
    virtual int64_t nodeset_node_count_proc(int64_t id) const;

    virtual int64_t sideset_count() const;
    virtual int64_t sideset_side_count_proc(int64_t id) const;

    virtual int64_t communication_node_count_proc() const;

    virtual void coordinates(double *coord) const;
    virtual void coordinates(std::vector<double> &coord) const;
    virtual void coordinates(int component, std::vector<double> &xyz) const;
    virtual void coordinates(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z) const;

    virtual void connectivity(int64_t block_number, int* connect) const;

    virtual std::pair<std::string, int> topology_type(int64_t block_number) const;

    virtual void sideset_elem_sides(int64_t setId, std::vector<int64_t> &elem_sides) const;

    virtual std::vector<std::string> sideset_touching_blocks(int64_t set_id) const;

    virtual void nodeset_nodes(int64_t nset_id, std::vector<int64_t> &nodes) const;

    virtual void node_communication_map(std::vector<int64_t> &map, std::vector<int> &proc);

    virtual void node_map(std::vector<int> &map) const;
    virtual void node_map(std::vector<int64_t> &map) const;

    virtual void element_map(int64_t block_number, std::vector<int> &map) const;
    virtual void element_map(int64_t block_number, std::vector<int64_t> &map) const;
    virtual void element_map(std::vector<int64_t> &map) const;
    virtual void element_map(std::vector<int> &map) const;

private:
    int64_t mGlobalNumberOfElements;
    int64_t mLocalNumberOfElements;

    const ExodusData & mExodusData;
    std::vector<int64_t> mElementOffsetForBlock;
};

}

#endif
