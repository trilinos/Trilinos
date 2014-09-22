// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
// 

#ifndef STK_SEARCH_MESHUTILSFORBOUNDINGVOLUMES_H_
#define STK_SEARCH_MESHUTILSFORBOUNDINGVOLUMES_H_

#include <exodusMeshInterface.h>
#include <gtest/gtest.h>
#include <optionParsing/getOption.h>

inline void createBoundingBoxForElement(const sierra::Mesh::LocalNodeId *connectivity, const int numNodesPerElement,
        const std::vector<double> &coordinates, std::vector<double>& boxCoordinates)
{
    int spatialDim = 3;
    double *minCoordinates = &boxCoordinates[0];
    double *maxCoordinates = &boxCoordinates[spatialDim];

    int firstNode=0;
    for (int j=0;j<spatialDim;j++)
    {
        minCoordinates[j] = coordinates[spatialDim*connectivity[firstNode]+j];
        maxCoordinates[j] = coordinates[spatialDim*connectivity[firstNode]+j];
    }

    for (int i=1;i<numNodesPerElement;i++)
    {
        sierra::Mesh::LocalNodeId nodeId = connectivity[i];
        for (int j=0;j<spatialDim;j++)
        {
           minCoordinates[j] = std::min(minCoordinates[j], coordinates[spatialDim*nodeId+j]);
           maxCoordinates[j] = std::max(maxCoordinates[j], coordinates[spatialDim*nodeId+j]);
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

inline void createBoundingBoxesForSidesInSidesets(const sierra::Mesh &mesh, const std::vector<double> &coordinates,
        std::vector<GtkBox>& domainBoxes)
{
    size_t numberBoundingBoxes = 0;
    sierra::Mesh::SideSetIdVector sidesetIds;
    mesh.fillSideSetIds(sidesetIds);

    for (size_t i=0;i<sidesetIds.size();i++)
    {
        numberBoundingBoxes += mesh.getNumberSidesInSideSet(sidesetIds[i]);
    }

    domainBoxes.resize(numberBoundingBoxes);

    size_t boxCounter = 0;

    std::vector<double> boxCoordinates(6);
    for (size_t ssetCounter=0;ssetCounter<sidesetIds.size();ssetCounter++)
    {
        sierra::Mesh::LocalNodeIdVector nodeIds;
        std::vector<int> numNodesPerFace;
        mesh.fillSideSetLocalNodeIds(sidesetIds[ssetCounter], nodeIds, numNodesPerFace);
        size_t offset=0;
        for (size_t i=0;i<numNodesPerFace.size();i++)
        {
            createBoundingBoxForElement(&nodeIds[offset], numNodesPerFace[i], coordinates, boxCoordinates);
            domainBoxes[boxCounter].set_box(boxCoordinates[0], boxCoordinates[1], boxCoordinates[2],
                                            boxCoordinates[3], boxCoordinates[4], boxCoordinates[5]);
            boxCounter++;
            offset += numNodesPerFace[i];
        }
    }

    ASSERT_EQ(boxCounter, numberBoundingBoxes);
}

inline void fillBoxesUsingSidesetsFromFile(MPI_Comm comm, const std::string& filename, std::vector<GtkBox> &domainBoxes)
{
    sierra::ExodusMeshInterface mesh(filename, comm);

    std::vector<double> coordinates;
    mesh.fillCoordinates(coordinates);

    createBoundingBoxesForSidesInSidesets(mesh, coordinates, domainBoxes);
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

inline void putCoordinatesInFile(const int exoid, const std::vector<GtkBox>& boxes)
{
    const int num_nodes_per_element = 8;
    const int spatialDim = 3;
    double *x = new double[num_nodes_per_element*boxes.size()];
    double *y = new double[num_nodes_per_element*boxes.size()];
    double *z = new double[num_nodes_per_element*boxes.size()];

    for (size_t i=0;i<boxes.size();i++)
    {
        double xmin=0, ymin=0, zmin=0;
        xmin = boxes[i].get_x_min();
        ymin = boxes[i].get_y_min();
        zmin = boxes[i].get_z_min();

        double xmax=0, ymax=0, zmax=0;
        xmax = boxes[i].get_x_max();
        ymax = boxes[i].get_y_max();
        zmax = boxes[i].get_z_max();

        double hexCoordinates[24];
        setHexCoordinates(xmin, ymin, zmin, xmax, ymax, zmax, &hexCoordinates[0]);
        unsigned offset = i*num_nodes_per_element;
        for (int j=0;j<num_nodes_per_element;j++)
        {
            x[offset+j] = hexCoordinates[spatialDim*j+0];
            y[offset+j] = hexCoordinates[spatialDim*j+1];
            z[offset+j] = hexCoordinates[spatialDim*j+2];
        }
    }

    ex_put_coord(exoid, x, y, z);

    delete [] z;
    delete [] y;
    delete [] x;
}

inline void fillNumElementsPerBlock(const int num_elements, std::vector<int> &numElementsPerBlock)
{
    int numElementsPer=1;
    if ( num_elements < 100 )
    {
        numElementsPer = 1;
    }
    else if ( num_elements < 1000 )
    {
        numElementsPer=10;
    }
    else if ( num_elements < 10000 )
    {
        numElementsPer=100;
    }
    else
    {
        numElementsPer=1000;
    }

    for (int i=0;i<num_elements;i+=numElementsPer)
    {
        int numElementsThisBlock = (i+numElementsPer) < num_elements ? numElementsPer : num_elements-i;
        numElementsPerBlock.push_back(numElementsThisBlock);
    }
}

inline void writeExodusFileUsingBoxes(const std::vector<GtkBox>& boxes, const std::string &filename)
{
    if ( boxes.size() == 0 )
    {
        std::cerr << "Skipping writing of file. No boxes to write.\n";
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
        ex_put_elem_block(exoid, blockId, "HEX", num_elements_this_block, num_nodes_per_elem, num_attr);

        for (int j=0;j<num_nodes_per_elem*num_elements_this_block;j++)
        {
            connect[j] = ordering[j%num_nodes_per_elem]+offset+num_nodes_per_elem*(j/num_nodes_per_elem);
        }
        offset += num_elements_this_block*num_nodes_per_elem;

        ex_put_elem_conn(exoid, blockId, &connect[0]);
    }

    ex_close(exoid);
}

inline void fillDomainBoxes(MPI_Comm comm, std::vector<GtkBox>& domainBoxes)
{
    std::string filename = unitTestUtils::getOption("-i", "input.exo");
    fillBoxesUsingSidesetsFromFile(comm, filename, domainBoxes);

    std::string exodusFilename = unitTestUtils::getOption("-o", "boxes.exo");
    if ( exodusFilename != "skip" )
    {
        writeExodusFileUsingBoxes(domainBoxes, exodusFilename);
    }
}

inline void fillStkBoxesUsingGtkBoxes(const std::vector<GtkBox> &domainBoxes, const int procId, StkBoxVector& stkBoxes)
{
    for (size_t i=0;i<domainBoxes.size();i++)
    {
        Point min(domainBoxes[i].get_x_min(), domainBoxes[i].get_y_min(), domainBoxes[i].get_z_min());
        Point max(domainBoxes[i].get_x_max(), domainBoxes[i].get_y_max(), domainBoxes[i].get_z_max());
        Ident domainBoxId(i, procId);
        stkBoxes[i] = std::make_pair(StkBox(min,max), domainBoxId);
    }
}

inline void createBoundingBoxesForElementsInElementBlocks(const int procId, const sierra::Mesh &mesh, const std::vector<double> &coordinates, GtkBoxVector& domainBoxes)
{
    size_t numberBoundingBoxes = mesh.getNumberLocalElements();
    domainBoxes.resize(numberBoundingBoxes);

    sierra::Mesh::BlockIdVector blockIds;
    mesh.fillElementBlockIds(blockIds);

    sierra::Mesh::AnalystElementIdVector analystElementIds;
    mesh.fillAnalystElementIds(analystElementIds);
    size_t boxCounter = 0;

    std::vector<double> boxCoordinates(6);

    for (size_t elemBlockNum=0;elemBlockNum<mesh.getNumberElementBlocks();elemBlockNum++)
    {
        sierra::Mesh::LocalNodeIdVector connectivity;
        mesh.fillElementToLocalNodeConnectivityForBlock(blockIds[elemBlockNum], connectivity);
        int numNodesPerElement = mesh.getNumberNodesPerElement(blockIds[elemBlockNum]);
        size_t numElementsThisBlock = mesh.getNumberLocalElementsInBlock(blockIds[elemBlockNum]);
        for (size_t elemCounter=0;elemCounter<numElementsThisBlock;elemCounter++)
        {
            createBoundingBoxForElement(&connectivity[numNodesPerElement*elemCounter], numNodesPerElement, coordinates, boxCoordinates);
            Ident domainBoxId(analystElementIds[boxCounter], procId);
            domainBoxes[boxCounter] = std::make_pair(GtkBox(boxCoordinates[0], boxCoordinates[1], boxCoordinates[2],
                                                            boxCoordinates[3], boxCoordinates[4], boxCoordinates[5]),
                                                            domainBoxId);
            boxCounter++;
        }
    }
}

inline void fillBoxesUsingElementBlocksFromFile(
        MPI_Comm comm, const std::string& volumeFilename, GtkBoxVector &domainBoxes)
{
    sierra::ExodusMeshInterface volumeMesh(volumeFilename, comm);

    std::vector<double> coordinates;
    volumeMesh.fillCoordinates(coordinates);

    int procId=-1;
    MPI_Comm_rank(comm, &procId);

    createBoundingBoxesForElementsInElementBlocks(procId, volumeMesh, coordinates, domainBoxes);
}

inline void fillBoundingVolumesUsingNodesFromFile(
        MPI_Comm comm, const std::string& sphereFilename, std::vector< std::pair<Sphere, Ident> > &spheres)
{
    int procId=-1;
    MPI_Comm_rank(comm, &procId);

    sierra::ExodusMeshInterface sphereMesh(sphereFilename, comm);

    std::vector<double> coordinates;
    sphereMesh.fillCoordinates(coordinates);

    const int spatialDim = 3;
    const size_t numSpheres = coordinates.size()/spatialDim;

    spheres.clear();
    spheres.resize(numSpheres);
    sierra::Mesh::AnalystNodeIdVector analystIds;
    sphereMesh.fillAnalystNodeIds(analystIds);

    for (size_t i=0;i<numSpheres;i++)
    {
        double x=coordinates[spatialDim*i];
        double y=coordinates[spatialDim*i+1];
        double z=coordinates[spatialDim*i+2];
        double radius=1e-5;
        spheres[i] = std::make_pair(Sphere(Point(x,y,z), radius), Ident(analystIds[i],procId));
    }
}

inline void fillBoundingVolumesUsingNodesFromFile(
        MPI_Comm comm, const std::string& sphereFilename, GtkBoxVector &spheres)
{
    int procId=-1;
    MPI_Comm_rank(comm, &procId);

    sierra::ExodusMeshInterface sphereMesh(sphereFilename, comm);

    std::vector<double> coordinates;
    sphereMesh.fillCoordinates(coordinates);

    const int spatialDim = 3;
    const size_t numSpheres = coordinates.size()/spatialDim;

    spheres.clear();
    spheres.resize(numSpheres);
    sierra::Mesh::AnalystNodeIdVector analystIds;
    sphereMesh.fillAnalystNodeIds(analystIds);

    for (size_t i=0;i<numSpheres;i++)
    {
        double x=coordinates[spatialDim*i];
        double y=coordinates[spatialDim*i+1];
        double z=coordinates[spatialDim*i+2];
        double radius=1e-5;
        double offset = radius;
        GtkBox box(x-offset, y-offset, z-offset, x+offset, y+offset, z+offset);
        spheres[i] = std::make_pair(box, Ident(analystIds[i],procId));
    }
}

inline void gtk_search(GtkBoxVector& local_domain, GtkBoxVector& local_range, MPI_Comm comm, SearchResults& searchResults)
{
    int num_procs = -1;
    int proc_id   = -1;
    MPI_Comm_rank(comm, &proc_id);
    MPI_Comm_size(comm, &num_procs);

    std::vector<gtk::AxisAlignedBB> rangeBoxes(local_range.size());
    std::vector<gtk::AxisAlignedBB> domainBoxes(local_domain.size());

    for (size_t i=0;i<local_domain.size();i++)
    {
        domainBoxes[i] = local_domain[i].first;
    }

    for (size_t i=0;i<local_range.size();i++)
    {
        rangeBoxes[i] = local_range[i].first;
    }

    std::vector<int> ghost_indices;
    std::vector<int> ghost_procs;
    ACME::BoxA_BoxB_Ghost(domainBoxes, rangeBoxes, comm, ghost_indices, ghost_procs);

    std::vector< std::vector<gtk::AxisAlignedBB> > send_list(num_procs);
    std::vector< std::vector<gtk::AxisAlignedBB> > recv_list(num_procs);

    // i am sending proc 'ghost proc[i]' my range box 'ghost_indices[i]'
    // ghost_indices.size() is total number of communications that need to occur with all procs

    std::vector< std::vector<int> > send_indices(num_procs);
    std::vector< std::vector<int> > recv_indices(num_procs);

    for (size_t i=0;i<ghost_indices.size();i++)
    {
        send_list[ghost_procs[i]].push_back(rangeBoxes[ghost_indices[i]]);
        int id = local_range[ghost_indices[i]].second.id();
        send_indices[ghost_procs[i]].push_back(id);
    }

    ACME::Parallel_Data_Exchange(send_indices, recv_indices, comm );
    ACME::Parallel_Data_Exchange(send_list, recv_list, comm);

    for (size_t i=0;i<recv_list.size();i++)
    {
        for (size_t j=0;j<recv_list[i].size();j++)
        {
            rangeBoxes.push_back(recv_list[i][j]);
            local_range.push_back(std::make_pair(recv_list[i][j], Ident(recv_indices[i][j], i)));
        }
    }

    std::vector<int> interaction_list;
    std::vector<int> first_interaction;
    std::vector<int> last_interaction;

    gtk::BoxA_BoxB_Search(domainBoxes, rangeBoxes, interaction_list, first_interaction, last_interaction);

    typedef std::vector <std::pair<Ident,Ident> > localJunk;
    typedef std::set <std::pair<Ident,Ident> > localJunkSet;
    localJunk localResults;
    localResults.reserve(domainBoxes.size());

    // Ident box1, Ident box2
    for (size_t i=0;i<domainBoxes.size();i++)
    {
        Ident box1_ident = local_domain[i].second;
        for (int j=first_interaction[i];j<last_interaction[i];j++)
        {
            //            Ident box2_ident = local_range[j].second;
            Ident box2_ident = local_range[interaction_list[j]].second;
            localResults.push_back(std::make_pair(box1_ident, box2_ident));
        }
    }

//    std::cerr << "LocalResults.size = " << localResults.size() << std::endl;
//    localJunkSet ljset;
//    std::copy(localResults.begin(), localResults.end(), inserter(ljset, ljset.end()));

    localJunk tmp;
    tmp.reserve(localResults.size());
    stk::search::communicateVector< std::pair<Ident,Ident>, std::pair<Ident,Ident> >(comm, localResults, tmp);

    searchResults=tmp;
//    std::copy(tmp.begin(), tmp.end(), std::back_inserter(searchResults));
    std::sort(searchResults.begin(), searchResults.end());
}

enum NewSearchMethod { BOOST_RTREE, OCTREE, GTK };
inline stk::search::SearchMethod mapSearchMethodToStk( NewSearchMethod method )
{
    if ( method == BOOST_RTREE )
    {
        return stk::search::BOOST_RTREE;
    }
    else if ( method == OCTREE )
    {
        return stk::search::OCTREE;
    }
    else
    {
        ThrowRequireMsg(false, "GTK method not implemented for this.");
    }
    return stk::search::BOOST_RTREE;
}


inline void coarse_search_new(GtkBoxVector& local_domain, GtkBoxVector& local_range, NewSearchMethod algorithm, MPI_Comm comm, SearchResults& searchResults)
{
    if ( algorithm == GTK )
    {
        gtk_search(local_domain, local_range, comm, searchResults);
    }
    else if ( algorithm == OCTREE )
    {
        stk::search::coarse_search(local_domain, local_range, stk::search::OCTREE, comm, searchResults);
    }
    else if ( algorithm == BOOST_RTREE )
    {
        stk::search::coarse_search(local_domain, local_range, stk::search::BOOST_RTREE, comm, searchResults);
    }
    else
    {
        throw("Invalid search algorithm: not supported.\n");
    }
}

#endif
