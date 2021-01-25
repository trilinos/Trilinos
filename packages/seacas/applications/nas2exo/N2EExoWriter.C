/*
 * N2EExoWriter.cpp
 *
 *  Created on: Oct 13, 2020
 *      Author: rjmoral
 */

#include "N2EExoWriter.h"

#include <algorithm>
#include <iostream>

namespace ExoModules {

  N2EExoWriter::N2EExoWriter()
  {

    this->exoFileID = -10;

    this->writtenBlocks = 0;
    this->writtenNodes  = 0;
    this->writtenTets   = 0;
    this->writtenHexes  = 0;
  }

  N2EExoWriter::~N2EExoWriter()
  {

    if (this->exoFileID > 0) {
      ex_close(this->exoFileID);
      this->exoFileID = -10;
    }
  }

  bool N2EExoWriter::createDB(string name)
  {

    this->exoFileID = ex_create(name.c_str(), EX_CLOBBER, &this->CPU_ws, &this->IO_ws);

    return this->exoFileID > 0;
  }

  bool N2EExoWriter::setNodes(vector<gridType> gridpts)
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

  bool N2EExoWriter::setElements(vector<elementType> elist)
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

  bool N2EExoWriter::setSections(vector<sectionType> sList)
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

  void N2EExoWriter::setModelTitle(string title)
  {

    string tmp = title;
    if (title.length() >= MAX_LINE_LENGTH - 1) {
      tmp = title.substr(0, MAX_LINE_LENGTH - 1);
    }

    strncat(this->modelTitle, tmp.c_str(), MAX_LINE_LENGTH - 1);
  }

  bool N2EExoWriter::writeFile()
  {
    // This poor variable is going to be over-used just to save 1-bit
    // Tsk-tsk
    bool result = true;

    this->writtenBlocks = 0;
    this->writtenNodes  = 0;
    this->writtenTets   = 0;
    this->writtenHexes  = 0;

    result &= this->sections.size() > 0;
    result &= this->gridList.size() >= 4;   // One Tet???
    result &= this->elementList.size() > 0; // At least one element

    result &= this->writeFileParams();
    result &= this->writeCoords();
    result &= this->writeElements();

    return result;
  }

  bool N2EExoWriter::writeCoords()
  {

    bool result{true};

    // Now we write the nodes
    unsigned num_nodes = this->gridList.size();
    double * x         = new double[num_nodes]{0.0};
    double * y         = new double[num_nodes]{0.0};
    double * z         = new double[num_nodes]{0.0};
    for (unsigned i = 0; i < num_nodes; i++) {

      N2EPoint3D crd = get<1>(this->gridList[i]);
      x[i]           = crd.x[0];
      y[i]           = crd.x[1];
      z[i]           = crd.x[2];
    }

    int ret = ex_put_coord(this->exoFileID, x, y, z);

    if (ret != 0) {
      std::cerr << "Problem writing node coordinates ";
      std::cerr << "in N2EExoWriter::writeFile().";
      std::cerr << "punching out";
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

    int ret = ex_put_init(this->exoFileID, this->modelTitle, 3 /* 3D models only*/,
                          this->gridList.size(), this->elementList.size(), this->sections.size(), 0,
                          0); // Make your fancy pants nodes and side sets elsewherem, laddy.

    if (ret != 0) {
      std::cerr << "Problem initializing model params ";
      std::cerr << "in N2EExoWriter::writeFile().";
      std::cerr << "punching out";
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

    for (sectionType sect : this->sections) {

      vector<elementType> thisBlock;
      int64_t             block = (int)get<0>(sect);

      int64_t                nodes_per_elem{0};
      std::unique_ptr<int[]> elemCon;
      int                    retvalue{0};

      for (elementType elem : this->elementList) {

        if ((int)get<1>(elem) == block) {
          thisBlock.emplace_back(elem);
        }
      }

      nodes_per_elem = (int)get<2>(thisBlock[0]);

      int n = nodes_per_elem == 4 ? 0 : 1;

      // This easy until we support 3 types of elements
      supportedElements thisElType = ExoElTypes[n];

      try {
        retvalue = ex_put_block(this->exoFileID, thisElType.elementType, block, thisElType.elemDesc,
                                thisBlock.size(), thisElType.numNodesPerElem,
                                thisElType.numEdgesPerElem, thisElType.numFacesPerElem, 1);
      }
      catch (...) {

        std::cerr << "ERROR!: In the block: " << block << std::endl;
        std::cerr << "Could not be configured.";
        return false;
      }
      if (retvalue != 0) {

        std::cerr << "WARNING:The block: " << block << std::endl;
        std::cerr << "May not be configured correctly.";
      }
      this->writtenBlocks++;

      // Statement below is supported C++ 14 and up.
      // elemCon = std::make_unique<int[]>(nodes_per_elem*thisBlock.size());
      // C++ 11 support
      elemCon.reset(new int[nodes_per_elem * thisBlock.size()]());

      int64_t numNodesCopied{0};

      for (elementType elem : thisBlock) {

        N2EGridPtList pts{get<3>(elem)};
        std::copy(pts.v, pts.v + nodes_per_elem, elemCon.get() + numNodesCopied);
        numNodesCopied += nodes_per_elem;
      }

      retvalue =
          ex_put_conn(this->exoFileID, thisElType.elementType, block, elemCon.get(), NULL, NULL);

      switch (thisElType.numNodesPerElem) {

      case 4: this->writtenTets += thisBlock.size(); break;

      case 8: this->writtenHexes += thisBlock.size(); break;
      }

      if (retvalue != 0) {
        result &= false;
      }
    }

    return result;
  }

} // namespace ExoModules
