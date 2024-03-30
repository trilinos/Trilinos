/*
 * N2EExoWriter.cpp
 *
 *  Created on: Oct 13, 2020
 *      Author: rjmoral
 */

#include "N2EExoWriter.h"
#include "vector_data.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <tuple>

namespace ExoModules {

  N2EExoWriter::~N2EExoWriter()
  {
    if (this->exoFileID > 0) {
      ex_close(this->exoFileID);
      this->exoFileID = 0;
    }
  }

  bool N2EExoWriter::createDB(const std::string &name)
  {
    this->exoFileID = ex_create(name.c_str(), EX_CLOBBER, &this->CPU_ws, &this->IO_ws);
    return this->exoFileID > 0;
  }

  bool N2EExoWriter::setNodes(const std::vector<gridType> &gridpts)
  {
    try {
      this->gridList.reserve(gridpts.capacity());
      this->gridList = gridpts;
    }
    catch (...) {
      return false;
    }

    return true;
  }

  bool N2EExoWriter::setElements(const std::vector<elementType> &elist)
  {
    try {
      this->elementList.reserve(elist.capacity());
      this->elementList = elist;
    }
    catch (...) {
      return false;
    }

    return true;
  }

  bool N2EExoWriter::setSections(const std::vector<sectionType> &sList)
  {
    try {
      this->sections.reserve(sList.capacity());
      this->sections = sList;
    }
    catch (...) {

      return false;
    }

    return true;
  }

  void N2EExoWriter::setModelTitle(const std::string &title) { this->modelTitle = title; }

  bool N2EExoWriter::writeFile()
  {
    // This poor variable is going to be over-used just to save 1-bit
    // Tsk-tsk
    bool result = true;

    this->writtenBlocks = 0;
    this->writtenNodes  = 0;
    this->writtenTets   = 0;
    this->writtenHexes  = 0;

    result &= !this->sections.empty();
    result &= this->gridList.size() >= 4; // One Tet???
    result &= !this->elementList.empty(); // At least one element

    result &= this->writeFileParams();
    result &= this->writeCoords();
    result &= this->writeElements();

    return result;
  }

  bool N2EExoWriter::writeCoords()
  {
    bool result{true};

    // Now we write the nodes
    size_t              num_nodes = this->gridList.size();
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    x.reserve(num_nodes);
    y.reserve(num_nodes);
    z.reserve(num_nodes);
    for (size_t i = 0; i < num_nodes; i++) {

      N2EPoint3D crd = std::get<1>(this->gridList[i]);
      x.push_back(crd.x[0]);
      y.push_back(crd.x[1]);
      z.push_back(crd.x[2]);
    }

    int ret = ex_put_coord(this->exoFileID, Data(x), Data(y), Data(z));

    if (ret != 0) {
      std::cerr << "Problem writing node coordinates in N2EExoWriter::writeFile(). punching out.\n";
      result = false;
    }
    else {
      this->writtenNodes = num_nodes;
      result             = true;
    }

    return result;
  }

  bool N2EExoWriter::writeFileParams()
  {
    bool result{true};

    auto tmp = this->modelTitle.substr(0, MAX_LINE_LENGTH - 1);
    int  ret = ex_put_init(this->exoFileID, tmp.c_str(), 3 /* 3D models only*/,
                           this->gridList.size(), this->elementList.size(), this->sections.size(), 0,
                           0); // Make your fancy pants nodes and side sets elsewherem, laddy.

    if (ret != 0) {
      std::cerr << "Problem initializing model params in N2EExoWriter::writeFile(). punching out\n";
      result = false;
    }
    else {
      result = true;
    }

    return result;
  }

  bool N2EExoWriter::writeElements()
  {
    bool result{true};

    for (const sectionType &sect : this->sections) {

      std::vector<elementType> thisBlock;
      auto                     block = std::get<0>(sect);

      int retvalue{0};

      for (const elementType &elem : this->elementList) {

        if (std::get<1>(elem) == block) {
          thisBlock.emplace_back(elem);
        }
      }

      auto nodes_per_elem = std::get<2>(thisBlock[0]);

      int n = nodes_per_elem == 4 ? 0 : 1;

      // This easy until we support 3 types of elements
      supportedElements thisElType = ExoElTypes[n];

      retvalue = ex_put_block(this->exoFileID, thisElType.elementType, block, thisElType.elemDesc,
                              thisBlock.size(), thisElType.numNodesPerElem,
                              thisElType.numEdgesPerElem, thisElType.numFacesPerElem, 1);
      if (retvalue != 0) {

        std::cerr << "WARNING:The block: " << block << "\nMay not be configured correctly.\n";
      }
      this->writtenBlocks++;

      std::vector<int> elemCon(nodes_per_elem * thisBlock.size());

      int64_t numNodesCopied{0};

      for (const elementType &elem : thisBlock) {

        const N2EModules::N2EGridPtList &pts{std::get<3>(elem)};
        std::copy(pts.v, pts.v + nodes_per_elem, Data(elemCon) + numNodesCopied);
        numNodesCopied += nodes_per_elem;
      }

      retvalue = ex_put_conn(this->exoFileID, thisElType.elementType, block, Data(elemCon), nullptr,
                             nullptr);

      switch (thisElType.numNodesPerElem) {

      case 4: this->writtenTets += thisBlock.size(); break;

      case 8: this->writtenHexes += thisBlock.size(); break;
      }

      if (retvalue != 0) {
        result = false;
      }
    }

    return result;
  }

} // namespace ExoModules
