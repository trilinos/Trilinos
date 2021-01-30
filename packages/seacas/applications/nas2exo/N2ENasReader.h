//============================================================================
// Name        : testnas2exo.cpp
// Author      : Ramon J. Moral (STRA LLC), John Niederhaus (Coordinator, SNL)
// Version     :
// Copyright   : (c) Sandia National Labs 2020
// Description : Testing nas2exo Library, C++ 14
//============================================================================

#ifndef INCLUDE_N2ENASREADER_H_
#define INCLUDE_N2ENASREADER_H_

#include "N2EDataTypes.h"
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>

using namespace std;
using namespace N2EModules;

namespace NasModules {

  class N2ENasReader
  {

  public:
    N2ENasReader(string ifname = "");

    virtual ~N2ENasReader();

    inline unsigned getLineCount() { return this->lineCount; };
    inline unsigned getNumberGridPts() { return this->gridList.size(); };
    inline unsigned getNumberElems() { return this->elementList.size(); };
    inline unsigned getNumberSects() { return this->sections.size(); };
    inline string   getPath() { return this->inFileName; };

    bool processFile();

    vector<sectionType> getSections() { return this->sections; };
    vector<gridType>    getGridPoints() { return this->gridList; };
    vector<elementType> getElementList() { return this->elementList; };

    string getModelTitle();
    void   setModelTitle(string title = "");

  protected:
    string               inFileName;
    unique_ptr<ifstream> inStream;
    unsigned             lineCount{0u};

    vector<sectionType> sections;
    vector<gridType>    gridList;
    vector<elementType> elementList;

    char modelTitle[72]{'\0'};

    bool doesFileExist(string fname);
    // Local buffer for reading faster
    char _readBuffer[4096]{'\0'};

  private:
    unsigned       lineCounter();
    vector<string> csvLineToTokens(char buff[]);
  };

} // namespace NasModules

#endif /* INCLUDE_N2ENASREADER_H_ */
