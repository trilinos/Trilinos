/*
 * N2EExoWriter.h
 *
 *  Created on: Oct 10, 2020
 *      Author: Ramon J. Moral(Contractor, STRA LLC)
 * 		John Niederhouse(ORG 1443, SNL, Coordinator)
 *  Copyright: Sandia National Labs, OCT-2022
 */
#pragma once

#include "N2EDataTypes.h"
#include <string>
#include <vector>

using namespace N2EModules;

namespace N2ETestingSpace {
  class ExoWriterTester;
}

namespace ExoModules {

  class N2EExoWriter
  {

  public:
    N2EExoWriter() = default;
    virtual ~N2EExoWriter();

    bool createDB(const std::string &name);
    bool setNodes(const std::vector<gridType> &gridpts);
    bool setElements(const std::vector<elementType> &elist);
    bool setSections(const std::vector<sectionType> &sList);

    bool writeFile();

    void setModelTitle(const std::string &title);

    inline size_t getBlocksOut() { return this->writtenBlocks; };
    inline size_t getNodesOut() { return this->writtenNodes; };
    inline size_t getTetsOut() { return this->writtenTets; };
    inline size_t getHexesOut() { return this->writtenHexes; };

  protected:
    std::vector<sectionType> sections{};
    std::vector<gridType>    gridList{};
    std::vector<elementType> elementList{};

    std::string modelTitle{};

    int exoFileID{0};

    int CPU_ws = sizeof(double);
    int IO_ws  = sizeof(double);

    // Data persisteance functions
    bool writeCoords();
    bool writeFileParams();
    bool writeElements();

    size_t writtenBlocks{0};
    size_t writtenNodes{0};
    size_t writtenTets{0};
    size_t writtenHexes{0};
    // A Friend for testing
    friend class N2ETestingSpace::ExoWriterTester;
  };

} // namespace ExoModules
