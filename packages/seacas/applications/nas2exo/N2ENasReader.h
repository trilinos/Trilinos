//============================================================================
// Name        : testnas2exo.cpp
// Author      : Ramon J. Moral (STRA LLC), John Niederhaus (Coordinator, SNL)
// Version     :
// Copyright   : (c) Sandia National Labs 2020, 2021, 2022, 2024
// Description : Testing nas2exo Library, C++ 14
//============================================================================

#pragma once

#include "N2EDataTypes.h"
#include <fstream>
#include <memory>
#include <string>
#include <vector>

using namespace N2EModules;

namespace NasModules {

  class N2ENasReader
  {

  public:
    N2ENasReader(std::string ifname = "");

    inline unsigned    getLineCount() { return this->lineCount; };
    inline size_t      getNumberGridPts() { return this->gridList.size(); };
    inline size_t      getNumberElems() { return this->elementList.size(); };
    inline size_t      getNumberSects() { return this->sections.size(); };
    inline std::string getPath() { return this->inFileName; };

    bool processFile();

    std::vector<sectionType> getSections() { return this->sections; };
    std::vector<gridType>    getGridPoints() { return this->gridList; };
    std::vector<elementType> getElementList() { return this->elementList; };

    std::string getModelTitle();
    void        setModelTitle(const std::string &title = "");

  protected:
    std::string                    inFileName{};
    std::unique_ptr<std::ifstream> inStream{};
    unsigned                       lineCount{0u};

    std::vector<sectionType> sections{};
    std::vector<gridType>    gridList{};
    std::vector<elementType> elementList{};

    std::string modelTitle{};

    bool doesFileExist(const std::string &fname);
    // Local buffer for reading faster
    char _readBuffer[4096]{'\0'};

  private:
    unsigned                 lineCounter();
    std::vector<std::string> csvLineToTokens(char buff[]);
  };

} // namespace NasModules
