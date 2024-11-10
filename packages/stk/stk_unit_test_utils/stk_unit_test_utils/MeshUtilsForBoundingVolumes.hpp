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
// 

#ifndef STK_SEARCH_MESHUTILSFORBOUNDINGVOLUMES_H_
#define STK_SEARCH_MESHUTILSFORBOUNDINGVOLUMES_H_

#include "stk_unit_test_utils/Search_UnitTestUtils.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/ExodusTranslator.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"
#include <exodusII.h>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <stk_util/parallel/ParallelComm.hpp>

using stk::unit_test_util::build_mesh;

inline void findBoundingBoxCoordinates(const std::vector<double> &coordinates, std::vector<double>& boxCoordinates)
{
    int spatialDim = 3;
    double *minCoordinates = boxCoordinates.data();
    double *maxCoordinates = &boxCoordinates[spatialDim];

    for (int j=0;j<spatialDim;j++)
    {
        minCoordinates[j] = coordinates[j];
        maxCoordinates[j] = coordinates[j];
    }

    int numNodesPerElement = coordinates.size()/3;

    for (int i=1;i<numNodesPerElement;i++)
    {
        for (int j=0;j<spatialDim;j++)
        {
           minCoordinates[j] = std::min(minCoordinates[j], coordinates[spatialDim*i+j]);
           maxCoordinates[j] = std::max(maxCoordinates[j], coordinates[spatialDim*i+j]);
        }
    }
    bool inflateBox = true;
    double percentInflation = 10;
    if ( inflateBox )
    {
        for (int i=0;i<spatialDim;i++)
        {
            double dist = maxCoordinates[i]-minCoordinates[i];
            if ( dist <= 1e-8 ) dist = 0.001;
            double inflation = dist*(0.5*percentInflation)/100.0;
            minCoordinates[i] -= inflation;
            maxCoordinates[i] += inflation;
        }
    }
}


inline void createBoundingBoxesForSidesInSidesets(const stk::mesh::BulkData& bulk, std::vector<FloatBox>& domainBoxes)
{
    stk::mesh::ExodusTranslator exoTranslator(bulk);
    size_t numberBoundingBoxes = 0;
    std::vector<int64_t> sidesetIds;
    exoTranslator.fill_side_set_ids(sidesetIds);

    for (size_t i=0;i<sidesetIds.size();i++)
    {
        numberBoundingBoxes += exoTranslator.get_local_num_entities_for_id(sidesetIds[i], bulk.mesh_meta_data().side_rank());
    }

    domainBoxes.resize(numberBoundingBoxes);

    stk::mesh::FieldBase const * coords = bulk.mesh_meta_data().coordinate_field();

    size_t boxCounter = 0;

    std::vector<double> boxCoordinates(6);
    for (size_t ssetCounter=0;ssetCounter<sidesetIds.size();ssetCounter++)
    {
        const stk::mesh::Part* sideset = exoTranslator.get_exodus_part_of_rank(sidesetIds[ssetCounter], bulk.mesh_meta_data().side_rank());
        stk::mesh::EntityVector sides;
        stk::mesh::Selector sel = bulk.mesh_meta_data().locally_owned_part() & *sideset;
        const bool sortById = true;
        stk::mesh::get_entities(bulk, bulk.mesh_meta_data().side_rank(), sel, sides, sortById);

        for(size_t j=0;j<sides.size();++j)
        {
            unsigned num_nodes_per_side = bulk.num_nodes(sides[j]);
            const stk::mesh::Entity* nodes = bulk.begin_nodes(sides[j]);
            std::vector<double> coordinates(3*num_nodes_per_side,0);
            for(unsigned k=0;k<num_nodes_per_side;++k)
            {
                double *data = static_cast<double*>(stk::mesh::field_data(*coords, nodes[k]));
                coordinates[3*k] = data[0];
                coordinates[3*k+1] = data[1];
                coordinates[3*k+2] = data[2];
            }
            findBoundingBoxCoordinates(coordinates, boxCoordinates);
            domainBoxes[boxCounter].set_box(boxCoordinates[0], boxCoordinates[1], boxCoordinates[2],
                                            boxCoordinates[3], boxCoordinates[4], boxCoordinates[5]);
            boxCounter++;
        }
    }

    STK_ThrowRequireMsg(boxCounter == numberBoundingBoxes, "Program error. Please contact sierra-help for support");
}

inline void fillBoxesUsingSidesetsFromFile(MPI_Comm comm, const std::string& filename, std::vector<FloatBox> &domainBoxes)
{
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3, comm);
    stk::io::fill_mesh(filename, *bulk);

    createBoundingBoxesForSidesInSidesets(*bulk, domainBoxes);
}

inline int openFileAndGetId(const int numBoxes, const int num_element_blocks, const std::string &filename)
{
    int CPU_word_size = sizeof(double);
    int IO_word_size = 8;
    int exoid = ex_create (filename.c_str(), EX_CLOBBER, &CPU_word_size, &IO_word_size);
    int num_dim = 3;
    int num_elements = numBoxes;
    int num_nodes_per_element = 8;
    int num_nodes = num_nodes_per_element*num_elements;

    int num_ns = 0, num_ss = 0;
    ex_put_init(exoid, "Boxes", num_dim, num_nodes, num_elements, num_element_blocks, num_ns, num_ss);
    return exoid;
}

inline void setHexCoordinates(const double &xmin, const double &ymin, const double &zmin,
                       const double &xmax, const double &ymax, const double &zmax,
                       double* hexCoordinates)
{
//    int ordering[8] = { 4, 3, 2, 1, 8, 7, 6, 5 }; // one based!
    int ordering[8] = { 3, 2, 1, 0, 7, 6, 5, 4 };

    hexCoordinates[3*ordering[0]+0] = xmin;
    hexCoordinates[3*ordering[0]+1] = ymin;
    hexCoordinates[3*ordering[0]+2] = zmin;

    hexCoordinates[3*ordering[1]+0] = xmax;
    hexCoordinates[3*ordering[1]+1] = ymin;
    hexCoordinates[3*ordering[1]+2] = zmin;

    hexCoordinates[3*ordering[2]+0] = xmax;
    hexCoordinates[3*ordering[2]+1] = ymin;
    hexCoordinates[3*ordering[2]+2] = zmax;

    hexCoordinates[3*ordering[3]+0] = xmin;
    hexCoordinates[3*ordering[3]+1] = ymin;
    hexCoordinates[3*ordering[3]+2] = zmax;

    hexCoordinates[3*ordering[4]+0] = xmin;
    hexCoordinates[3*ordering[4]+1] = ymax;
    hexCoordinates[3*ordering[4]+2] = zmin;

    hexCoordinates[3*ordering[5]+0] = xmax;
    hexCoordinates[3*ordering[5]+1] = ymax;
    hexCoordinates[3*ordering[5]+2] = zmin;

    hexCoordinates[3*ordering[6]+0] = xmax;
    hexCoordinates[3*ordering[6]+1] = ymax;
    hexCoordinates[3*ordering[6]+2] = zmax;

    hexCoordinates[3*ordering[7]+0] = xmin;
    hexCoordinates[3*ordering[7]+1] = ymax;
    hexCoordinates[3*ordering[7]+2] = zmax;
}

inline void putCoordinatesInFile(const int exoid, const std::vector<FloatBox>& boxes)
{
    const int num_nodes_per_element = 8;
    const int spatialDim = 3;

    std::vector<double> x(num_nodes_per_element*boxes.size());
    std::vector<double> y(num_nodes_per_element*boxes.size());
    std::vector<double> z(num_nodes_per_element*boxes.size());

    for (size_t i=0; i<boxes.size(); i++)
    {
        double xmin = boxes[i].get_x_min();
        double ymin = boxes[i].get_y_min();
        double zmin = boxes[i].get_z_min();

        double xmax = boxes[i].get_x_max();
        double ymax = boxes[i].get_y_max();
        double zmax = boxes[i].get_z_max();

        double hexCoordinates[24];
        setHexCoordinates(xmin, ymin, zmin, xmax, ymax, zmax, hexCoordinates);

        unsigned offset = i*num_nodes_per_element;
        for (int j=0; j<num_nodes_per_element; j++)
        {
            x[offset+j] = hexCoordinates[spatialDim*j+0];
            y[offset+j] = hexCoordinates[spatialDim*j+1];
            z[offset+j] = hexCoordinates[spatialDim*j+2];
        }
    }

    ex_put_coord(exoid, x.data(), y.data(), z.data());
}

inline void fillNumElementsPerBlock(const int num_elements, std::vector<int> &numElementsPerBlock)
{
    const int minNumElementsPer = 1;
    const int maxNumElementsPer = 1000;

    int numElementsPer = num_elements / 100;

    numElementsPer = std::max( numElementsPer, minNumElementsPer );
    numElementsPer = std::min( numElementsPer, maxNumElementsPer );

    for (int i = 0; i < num_elements; i += numElementsPer)
    {
        int numElementsThisBlock = (i+numElementsPer) < num_elements ? numElementsPer : num_elements-i;
        numElementsPerBlock.push_back(numElementsThisBlock);
    }
}

inline void writeExodusFileUsingBoxes(const std::vector<FloatBox>& boxes, const std::string &filename)
{
    if ( boxes.size() == 0 )
    {
        // std::cerr << "Skipping writing of file. No boxes to write.\n";
        return;
    }

    const int num_nodes_per_elem = 8;
    const int num_attr = 0;
    const unsigned num_elements = boxes.size();
    std::vector<int> numElementsPerBlock;
    fillNumElementsPerBlock(num_elements, numElementsPerBlock);
    const int num_blocks = numElementsPerBlock.size();
    const int exoid = openFileAndGetId(boxes.size(), num_blocks, filename);
    putCoordinatesInFile(exoid, boxes);

    std::vector<int> connect(numElementsPerBlock[0]*num_nodes_per_elem);
    int ordering[8] = { 4, 3, 2, 1, 8, 7, 6, 5 }; // one based!
    unsigned offset = 0;
    for (int blockId=1;blockId<=num_blocks;blockId++)
    {
        const int num_elements_this_block = numElementsPerBlock[blockId-1];
        ex_put_block(exoid, EX_ELEM_BLOCK, blockId, "HEX", num_elements_this_block, num_nodes_per_elem, num_attr, 0, 0);

        for (int j=0;j<num_nodes_per_elem*num_elements_this_block;j++)
        {
            connect[j] = ordering[j%num_nodes_per_elem]+offset+num_nodes_per_elem*(j/num_nodes_per_elem);
        }
        offset += num_elements_this_block*num_nodes_per_elem;

        ex_put_conn(exoid, EX_ELEM_BLOCK, blockId, connect.data(), nullptr, nullptr);
    }

    ex_close(exoid);
}

inline std::vector<FloatBox>
fillDomainBoxes(MPI_Comm comm)
{
    std::vector<FloatBox> domainBoxes;
    std::string filename = stk::unit_test_util::get_option("-i", "input.exo");
    fillBoxesUsingSidesetsFromFile(comm, filename, domainBoxes);

    std::string exodusFilename = stk::unit_test_util::get_option("-o", "boxes.exo");
    if ( exodusFilename != "skip" )
    {
        writeExodusFileUsingBoxes(domainBoxes, exodusFilename);
    }

    return domainBoxes;
}

inline void convertFloatBoxesToDoubleBoxes(const FloatBoxIdentProcVector &floatBoxes, StkBoxIdentProcVector& doubleBoxes)
{
  doubleBoxes.resize(floatBoxes.size());
  for (size_t i=0;i<floatBoxes.size();i++) {
    Point min(floatBoxes[i].first.get_x_min(), floatBoxes[i].first.get_y_min(), floatBoxes[i].first.get_z_min());
    Point max(floatBoxes[i].first.get_x_max(), floatBoxes[i].first.get_y_max(), floatBoxes[i].first.get_z_max());
    doubleBoxes[i] = std::make_pair(StkBox(min,max), floatBoxes[i].second);
  }
}

inline void fillStkBoxesUsingFloatBoxes(const std::vector<FloatBox> &domainBoxes, const int procId, StkBoxIdentProcVector& stkBoxes)
{
    for (size_t i=0;i<domainBoxes.size();i++)
    {
        Point min(domainBoxes[i].get_x_min(), domainBoxes[i].get_y_min(), domainBoxes[i].get_z_min());
        Point max(domainBoxes[i].get_x_max(), domainBoxes[i].get_y_max(), domainBoxes[i].get_z_max());
        IdentProc domainBoxId(i, procId);
        stkBoxes[i] = std::make_pair(StkBox(min,max), domainBoxId);
    }
}

template<typename BoxType, typename IdentProcType>
inline void createBoundingBoxesForEntities(const stk::mesh::BulkData &bulk,
                                           stk::mesh::EntityRank rank,
                                           std::vector<std::pair<BoxType,IdentProcType>>& boundingBoxes)
{
    stk::mesh::EntityVector entities;
    const bool sortById = true;
    stk::mesh::get_entities(bulk, rank, bulk.mesh_meta_data().locally_owned_part(), entities, sortById);

    size_t numberBoundingBoxes = entities.size();
    boundingBoxes.resize(numberBoundingBoxes);

    stk::mesh::FieldBase const * coords = bulk.mesh_meta_data().coordinate_field();

    std::vector<double> boxCoordinates(6);

    for(size_t i=0;i<entities.size();++i)
    {
        unsigned num_nodes = bulk.num_nodes(entities[i]);
        std::vector<double> coordinates(3*num_nodes,0);
        const stk::mesh::Entity* nodes = bulk.begin_nodes(entities[i]);
        for(unsigned j=0;j<num_nodes;++j)
        {
            double* data = static_cast<double*>(stk::mesh::field_data(*coords, nodes[j]));
            coordinates[3*j] = data[0];
            coordinates[3*j+1] = data[1];
            coordinates[3*j+2] = data[2];
        }
        findBoundingBoxCoordinates(coordinates, boxCoordinates);
        unsigned id = bulk.identifier(entities[i]);
        IdentProcType domainBoxId;
        if constexpr (std::is_same_v<IdentProcType, IdentProc>) {
          domainBoxId = IdentProc(id, bulk.parallel_rank());
        }
        else {
            domainBoxId = id;
        }
        boundingBoxes[i] = std::make_pair(BoxType(boxCoordinates[0], boxCoordinates[1], boxCoordinates[2],
                                                  boxCoordinates[3], boxCoordinates[4], boxCoordinates[5]),
                                                  domainBoxId);

    }
}

template<typename BoxIdentProcType>
inline Kokkos::View<BoxIdentProcType *>
createBoundingBoxesForEntities(const stk::mesh::BulkData &bulk,
                                     stk::mesh::EntityRank rank)
{

  using BoxType = typename BoxIdentProcType::box_type;  
  using IdentProcType = typename BoxIdentProcType::second_type;  
  stk::mesh::EntityVector entities;
  const bool sortById = true;
  stk::mesh::get_entities(bulk, rank, bulk.mesh_meta_data().locally_owned_part(), entities, sortById);

  size_t numberBoundingBoxes = entities.size();
  Kokkos::View<BoxIdentProcType *> boundingBoxes("Bounding Boxes", numberBoundingBoxes);
  auto boundingBoxesHost = Kokkos::create_mirror_view(boundingBoxes);

  stk::mesh::FieldBase const * coords = bulk.mesh_meta_data().coordinate_field();

  std::vector<double> boxCoordinates(6);

  for (size_t i = 0; i < entities.size(); ++i) {
    unsigned num_nodes = bulk.num_nodes(entities[i]);
    std::vector<double> coordinates(3*num_nodes,0);
    const stk::mesh::Entity* nodes = bulk.begin_nodes(entities[i]);
    for (unsigned j = 0; j < num_nodes; ++j) {
      double* data = static_cast<double*>(stk::mesh::field_data(*coords, nodes[j]));
      coordinates[3*j] = data[0];
      coordinates[3*j+1] = data[1];
      coordinates[3*j+2] = data[2];
    }
    findBoundingBoxCoordinates(coordinates, boxCoordinates);

    int id = bulk.identifier(entities[i]);
    IdentProcType domainBoxId;
    if constexpr (std::is_same_v<IdentProcType, IdentProc>) {
      domainBoxId = IdentProc(id, bulk.parallel_rank());
    }
    else {
        domainBoxId = id;
    }

    boundingBoxesHost(i) = {BoxType(boxCoordinates[0], boxCoordinates[1], boxCoordinates[2],
                                                 boxCoordinates[3], boxCoordinates[4], boxCoordinates[5]), domainBoxId};
  }

  Kokkos::deep_copy(boundingBoxes, boundingBoxesHost);

  return boundingBoxes;
}

inline void createBoundingBoxesForElementsInElementBlocks(const stk::mesh::BulkData &bulk, FloatBoxIdentProcVector& domainBoxes)
{
  createBoundingBoxesForEntities(bulk, stk::topology::ELEM_RANK, domainBoxes);
}

inline void fillBoxesUsingElementBlocksFromFile(MPI_Comm comm, const std::string& volumeFilename, FloatBoxIdentProcVector &domainBoxes)
{
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3, comm);
    stk::io::fill_mesh(volumeFilename, *bulk);

    createBoundingBoxesForElementsInElementBlocks(*bulk, domainBoxes);
}

inline void fillBoundingVolumesUsingNodesFromFile(
        MPI_Comm comm, const std::string& sphereFilename, std::vector< std::pair<Sphere, IdentProc> > &spheres)
{
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3, comm);
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    stk::io::fill_mesh(sphereFilename, *bulk);

    stk::mesh::EntityVector nodes;
    const bool sortById = true;
    stk::mesh::get_entities(*bulk, stk::topology::NODE_RANK, meta.locally_owned_part(), nodes, sortById);

    spheres.clear();
    spheres.resize(nodes.size());

    stk::mesh::FieldBase const * coords = meta.coordinate_field();

    for (size_t i=0;i<nodes.size();i++)
    {
        stk::mesh::Entity node = nodes[i];
        double *data = static_cast<double*>(stk::mesh::field_data(*coords, node));

        double x=data[0];
        double y=data[1];
        double z=data[2];

        double radius=1e-5;
        unsigned id = bulk->identifier(node);
        spheres[i] = std::make_pair(Sphere(Point(x,y,z), radius), IdentProc(id, bulk->parallel_rank()));
    }
}

inline void fillBoundingVolumesUsingNodesFromFile(MPI_Comm comm, const std::string& sphereFilename, FloatBoxIdentProcVector &spheres)
{
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3, comm);
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    stk::io::fill_mesh(sphereFilename, *bulk);

    stk::mesh::EntityVector nodes;
    const bool sortById = true;
    stk::mesh::get_entities(*bulk, stk::topology::NODE_RANK, meta.locally_owned_part(), nodes, sortById);

    spheres.clear();
    spheres.resize(nodes.size());

    stk::mesh::FieldBase const * coords = meta.coordinate_field();

    for (size_t i=0;i<nodes.size();i++)
    {
        stk::mesh::Entity node = nodes[i];
        const double *data = static_cast<const double*>(stk::mesh::field_data(*coords, node));

        const double x=data[0];
        const double y=data[1];
        const double z=data[2];

        constexpr double radius=1e-5;
        const unsigned id = bulk->identifier(node);
        FloatBox box(x-radius, y-radius, z-radius, x+radius, y+radius, z+radius);
        spheres[i] = std::make_pair(box, IdentProc(id, bulk->parallel_rank()));
        STK_ThrowRequire(spheres[i].first == box);
    }
}

namespace simple_fields {

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void findBoundingBoxCoordinates(const std::vector<double> &coordinates, std::vector<double>& boxCoordinates);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void createBoundingBoxesForSidesInSidesets(const stk::mesh::BulkData& bulk, std::vector<FloatBox>& domainBoxes);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void fillBoxesUsingSidesetsFromFile(MPI_Comm comm, const std::string& filename, std::vector<FloatBox> &domainBoxes);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
int openFileAndGetId(const int numBoxes, const int num_element_blocks, const std::string &filename);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void setHexCoordinates(const double &xmin, const double &ymin, const double &zmin,
                       const double &xmax, const double &ymax, const double &zmax,
                       double* hexCoordinates);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void putCoordinatesInFile(const int exoid, const std::vector<FloatBox>& boxes);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void fillNumElementsPerBlock(const int num_elements, std::vector<int> &numElementsPerBlock);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void writeExodusFileUsingBoxes(const std::vector<FloatBox>& boxes, const std::string &filename);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
std::vector<FloatBox> fillDomainBoxes(MPI_Comm comm);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void fillStkBoxesUsingFloatBoxes(const std::vector<FloatBox> &domainBoxes, const int procId, StkBoxIdentProcVector& stkBoxes);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void createBoundingBoxesForElementsInElementBlocks(const stk::mesh::BulkData &bulk, FloatBoxIdentProcVector& domainBoxes);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void fillBoxesUsingElementBlocksFromFile(MPI_Comm comm, const std::string& volumeFilename, FloatBoxIdentProcVector &domainBoxes);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void fillBoundingVolumesUsingNodesFromFile(
        MPI_Comm comm, const std::string& sphereFilename, std::vector< std::pair<Sphere, IdentProc> > &spheres);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void fillBoundingVolumesUsingNodesFromFile(MPI_Comm comm, const std::string& sphereFilename, FloatBoxIdentProcVector &spheres);

} // namespace simple_fields

#endif
