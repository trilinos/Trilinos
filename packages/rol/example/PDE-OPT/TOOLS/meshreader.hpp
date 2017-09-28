// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  meshreader.hpp
    \brief Defines the MeshReader class.
*/

#ifndef MESHREADER_HPP
#define MESHREADER_HPP

#include "Shards_CellTopology.hpp"

#include <iostream>
#include <fstream>
#include <cstring>
#include "meshmanager.hpp"


/** \class  MeshReader
    \brief  This class implements the pure virtual MeshManager interface,
            through methods that read in text (ASCII) files generated
            from Exodus mesh files using the 'ncdump' tool.
*/
template <class Real>
class MeshReader : public MeshManager<Real> {
private:

  Teuchos::ParameterList parlist_;

  int spaceDim_;
  int numNodes_;
  int numCells_;
  int numEdges_;
  int numFaces_;
  int numSideSets_;
  int numNodesPerCell_;
  int numEdgesPerCell_;
  int numFacesPerCell_;
  int numNodesPerFace_;

  Teuchos::RCP<shards::CellTopology> cellTopo_;

  Teuchos::RCP<Intrepid::FieldContainer<Real> > meshNodes_;
  Teuchos::RCP<Intrepid::FieldContainer<int> >  meshCellToNodeMap_;
  Teuchos::RCP<Intrepid::FieldContainer<int> >  meshCellToEdgeMap_;
  Teuchos::RCP<Intrepid::FieldContainer<int> >  meshCellToFaceMap_;

  Teuchos::RCP<std::vector<std::vector<std::vector<int> > > >  meshSideSets_;

public:

  /** \brief Constructor.  Parses the mesh file, and fills all data structures.
  */
  MeshReader(Teuchos::ParameterList & parlist) : parlist_(parlist),
      spaceDim_(0), numNodes_(0), numCells_(0), numEdges_(0), numFaces_(0),
      numSideSets_(0), numNodesPerCell_(0), numEdgesPerCell_(0), numFacesPerCell_(0), numNodesPerFace_(0) {
    std::string   filename = parlist.sublist("Mesh").get("File Name", "mesh.txt");
    std::ifstream inputfile;
    std::string   line;
    std::string   token;

    inputfile.open(filename);

    // Check if file readable.
    if (!inputfile.good()) {
      throw std::runtime_error("\nMeshReader: Could not open mesh file!\n");
    }

    // Parse file header.
    bool processHeader = true;
    while (processHeader) {

      std::getline(inputfile, line);    // consider: while (getline(inputfile, line).good())
      std::stringstream ssline(line);

      while (ssline >> token) {

        // Get space dimension.
        if (token == "num_dim") {
          ssline >> token;  // skip "="
          ssline >> spaceDim_;
          break;
        }

        // Get number of nodes.
        if (token == "num_nodes") {
          ssline >> token;  // skip "="
          ssline >> numNodes_;
          break;
        }

        // Get number of cells.
        if (token == "num_elem") {
          ssline >> token;  // skip "="
          ssline >> numCells_;
          break;
        }

        // Get number of sidesets.
        if (token == "num_side_sets") {
          ssline >> token;  // skip "="
          ssline >> numSideSets_;
          processHeader = false;
          break;
        }

      } // end tokenizing

    } // end parse header

    cellTopo_ = Teuchos::rcp(new shards::CellTopology( shards::getCellTopologyData<shards::Hexahedron<8> >() ));

    // Set up internal storage.
    numNodesPerCell_ = cellTopo_->getVertexCount();
    numFacesPerCell_ = cellTopo_->getFaceCount();
    numEdgesPerCell_ = cellTopo_->getEdgeCount();
    numNodesPerFace_ = cellTopo_->getVertexCount(2,0);
    meshNodes_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numNodes_, spaceDim_));
    Intrepid::FieldContainer<Real> &nodes = *meshNodes_;
    meshCellToNodeMap_ = Teuchos::rcp(new Intrepid::FieldContainer<int>(numCells_, numNodesPerCell_));
    Intrepid::FieldContainer<int> &ctn = *meshCellToNodeMap_;

    // Parse node coordinates.
    std::vector<int> coordCount(3, 0);
    std::vector<std::string> coordLabel(3);
    coordLabel[0] = "coordx";
    coordLabel[1] = "coordy";
    coordLabel[2] = "coordz";
    char semicolon = ';';

    for (int d=0; d<spaceDim_; ++d) {

      bool processCoordsHeader = true;
      while (processCoordsHeader) {

        std::getline(inputfile, line);
        std::stringstream ssline(line);

        while (ssline >> token) {

          // Get node coordinates.
          if (token == coordLabel[d]) {
            ssline >> token;  // skip "="
            while (std::getline(ssline, token, ',')) {
              if (token != " ") {
                nodes(coordCount[d]++, d) = atof(token.c_str());
              }
            }
            processCoordsHeader = false;
          }
        }

      }

      bool processCoords = true;
      if (token.find(semicolon) != std::string::npos) {
        processCoords = false;
      }
      while (processCoords) {

        std::getline(inputfile, line);
        std::stringstream ssline(line);

        while (std::getline(ssline, token, ',')) {
          // Get node coordinates.
          if (token.find(semicolon) == std::string::npos) {
            if (token != " ") {
              nodes(coordCount[d]++, d) = atof(token.c_str());
            }
          }
          else { // encountered "some-number ;"
            std::string onlyNumber = token.substr(0, token.size()-1);
            nodes(coordCount[d]++, d) = atof(onlyNumber.c_str());
            processCoords = false;
          }
        }
      } // end process coords

    } // end dimension loop


    // Find connectivity information.
    bool searchConnect = true;
    while (searchConnect) {

      std::getline(inputfile, line);
      std::stringstream ssline(line);

      while (ssline >> token) {
        if (token == "connect1") { // found connectivity, cell-to-node map
          searchConnect = false;
          break;
        }
      }
    }

    // Process connectivity information.
    bool processConnectivity = true;
    int  cellCount = 0;
    int  nodeId = 0;
    while (processConnectivity) {

      std::getline(inputfile, line);
      std::stringstream ssline(line);

      while (std::getline(ssline, token, ',')) {
        // Get node ids.
        if (token.find(semicolon) == std::string::npos) {
          if (token != " ") {
            ctn(cellCount, nodeId++) = atoi(token.c_str())-1;
            if (nodeId == numNodesPerCell_) {
              nodeId = 0;
              cellCount++;
            }
          }
        }
        else { // encountered "some-number ;"
          std::string onlyNumber = token.substr(0, token.size()-1);
          ctn(cellCount, numNodesPerCell_-1) = atoi(onlyNumber.c_str())-1;
          processConnectivity = false;
        }
      }
    } // end process connectivity


    // Parse side sets.
    meshSideSets_ = Teuchos::rcp(new std::vector<std::vector<std::vector<int> > >(numSideSets_));
    for (int ss=0; ss<numSideSets_; ++ss) {
      (*meshSideSets_)[ss].resize(numFacesPerCell_);
    }
    std::vector<int> ssCellIds;

    for (int ss=0; ss<numSideSets_; ++ss) {

      ssCellIds.clear();

      bool processSSHeaderCell = true;
      while (processSSHeaderCell) {

        std::getline(inputfile, line);
        std::stringstream ssline(line);

        while (ssline >> token) {

          // Get cell ids.
          if (token.substr(0,7) == "elem_ss") {
            ssline >> token;  // skip "="
            while (std::getline(ssline, token, ',')) {
              if (token != " ") {
                ssCellIds.push_back(atoi(token.c_str())-1);
              }
            }
            processSSHeaderCell = false;
          }
        }

      }

      bool processSSCell = true;
      if (token.find(semicolon) != std::string::npos) {
        processSSCell = false;
      }
      while (processSSCell) {

        std::getline(inputfile, line);
        std::stringstream ssline(line);

        while (std::getline(ssline, token, ',')) {
          // Get cell ids.
          if (token.find(semicolon) == std::string::npos) {
            if (token != " ") {
              ssCellIds.push_back(atoi(token.c_str())-1);
            }
          }
          else { // encountered "some-number ;"
            std::string onlyNumber = token.substr(0, token.size()-1);
            ssCellIds.push_back(atoi(token.c_str())-1);
            processSSCell = false;
          }
        }
      } // end process side set cells

      bool processSSHeaderSide = true;
      int cellIdx = 0;
      while (processSSHeaderSide) {

        std::getline(inputfile, line);
        std::stringstream ssline(line);

        while (ssline >> token) {

          // Get local side ids.
          if (token.substr(0,7) == "side_ss") {
            ssline >> token;  // skip "="
            while (std::getline(ssline, token, ',')) {
              if (token != " ") {
                (*meshSideSets_)[ss][atoi(token.c_str())-1].push_back(ssCellIds[cellIdx]);
                cellIdx++;
              }
            }
            processSSHeaderSide = false;
          }
        }

      }

      bool processSSSide = true;
      if (token.find(semicolon) != std::string::npos) {
        processSSSide = false;
      }
      while (processSSSide) {

        std::getline(inputfile, line);
        std::stringstream ssline(line);

        while (std::getline(ssline, token, ',')) {
          // Get cell ids.
          if (token.find(semicolon) == std::string::npos) {
            if (token != " ") {
              (*meshSideSets_)[ss][atoi(token.c_str())-1].push_back(ssCellIds[cellIdx]);
              cellIdx++;
            }
          }
          else { // encountered "some-number ;"
            std::string onlyNumber = token.substr(0, token.size()-1);
            (*meshSideSets_)[ss][atoi(token.c_str())-1].push_back(ssCellIds[cellIdx]);
            ssCellIds.clear();
            processSSSide = false;
          }
        }
      } // end process side set sides


    } // end dimension loop

    inputfile.close();

    computeCellToEdgeMap();

    computeCellToFaceMap();

  }


  Teuchos::RCP<Intrepid::FieldContainer<Real> > getNodes() const {
    return meshNodes_;
  }


  Teuchos::RCP<Intrepid::FieldContainer<int> > getCellToNodeMap() const {
    return meshCellToNodeMap_;
  }


  Teuchos::RCP<Intrepid::FieldContainer<int> > getCellToEdgeMap() const {
    return meshCellToEdgeMap_;
  }


  Teuchos::RCP<Intrepid::FieldContainer<int> > getCellToFaceMap() const {
    return meshCellToFaceMap_;
  }


  Teuchos::RCP<std::vector<std::vector<std::vector<int> > > > getSideSets (
      const bool verbose = false,
      std::ostream & outStream = std::cout) const {
    if (verbose) {
      for (int i=0; i<numSideSets_; ++i) {
        outStream << "\nSideset " << i << std::endl;
        for (int j=0; j<numFacesPerCell_; ++j) {
          outStream << "    Local side " << j << ":";
          for (int k=0; k<static_cast<int>((*meshSideSets_)[i][j].size()); ++k) {
            outStream << " " << (*meshSideSets_)[i][j][k];
          }
          outStream << std::endl;
        }
      }
    }
    return meshSideSets_;
  }


  int getNumCells() const {
    return numCells_;
  } // getNumCells


  int getNumNodes() const {
    return numNodes_;
  } // getNumNodes


  int getNumEdges() const {
    return numEdges_;
  } // getNumEdges

  int getNumFaces() const {
    return numFaces_;
  } // getNumEdges

  int getNumSideSets() const {
    return numSideSets_;
  } // getNumSideSets

  void computeCellToEdgeMap() {
    meshCellToEdgeMap_ = Teuchos::rcp(new Intrepid::FieldContainer<int>(numCells_, numEdgesPerCell_));
    Intrepid::FieldContainer<int> &cte = *meshCellToEdgeMap_;
    Intrepid::FieldContainer<int> &ctn = *meshCellToNodeMap_;
    std::set<int> edgenodes;
    std::map<std::set<int>,int> edgemap;
    int edgeCt = 0;

    /* Build edge map. */
    for (int i=0; i<numCells_; ++i) { // loop over cells
      for (int j=0; j<numEdgesPerCell_; ++j) { // loop over local edges
        edgenodes.insert(ctn(i, cellTopo_->getNodeMap(1, j, 0)));
        edgenodes.insert(ctn(i, cellTopo_->getNodeMap(1, j, 1)));
        //std::pair<std::map<std::set<int>,int>::iterator,bool> ret;
        auto ret = edgemap.insert(std::pair<std::set<int>,int>(edgenodes, edgeCt++));
        cte(i, j) = edgemap[edgenodes];
        if (ret.second == false) {
          edgeCt--;
        }
        edgenodes.clear();
      }
    }

    numEdges_ = edgemap.size();

    /* Print edges and cell-to-edge map, for debugging. */
    /*for (auto it_em=edgemap.begin(); it_em != edgemap.end(); ++it_em) {
      for (auto it_ns=(it_em->first).begin(); it_ns != (it_em->first).end(); ++it_ns) {
        std::cout << *it_ns << " , ";
      }
      std::cout << " => " << it_em->second << std::endl;
    }
    std::cout << cte;*/

  }

  void computeCellToFaceMap() {
    meshCellToFaceMap_ = Teuchos::rcp(new Intrepid::FieldContainer<int>(numCells_, numFacesPerCell_));
    Intrepid::FieldContainer<int> &ctf = *meshCellToFaceMap_;
    Intrepid::FieldContainer<int> &ctn = *meshCellToNodeMap_;
    std::set<int> facenodes;
    std::map<std::set<int>,int> facemap;
    int faceCt = 0;

    /* Build face map. */
    for (int i=0; i<numCells_; ++i) { // loop over cells
      for (int j=0; j<numFacesPerCell_; ++j) { // loop over local faces
        for (int k=0; k<numNodesPerFace_; ++k) { // loop over nodes to insert in sorted order
          facenodes.insert(ctn(i, cellTopo_->getNodeMap(2, j, k)));
        }
        //std::pair<std::map<std::set<int>,int>::iterator,bool> ret;
        auto ret = facemap.insert(std::pair<std::set<int>,int>(facenodes, faceCt++));
        ctf(i, j) = facemap[facenodes];
        if (ret.second == false) {
          faceCt--;
        }
        facenodes.clear();
      }
    }

    numFaces_ = facemap.size();

    /* Print faces and cell-to-face map, for debugging. */
    /*for (auto it_fm=facemap.begin(); it_fm != facemap.end(); ++it_fm) {
      for (auto it_ns=(it_fm->first).begin(); it_ns != (it_fm->first).end(); ++it_ns) {
        std::cout << *it_ns << " , ";
      }
      std::cout << " => " << it_fm->second << std::endl;
    }
    std::cout << ctf;*/

  }

}; // MeshReader

#endif
