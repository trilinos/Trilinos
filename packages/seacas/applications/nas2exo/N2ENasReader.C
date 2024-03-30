//============================================================================
// Name        : testnas2exo.cpp
// Author      : Ramon J. Moral (STRA LLC), John Niederhaus (Coordinator, SNL)
// Version     :
// Copyright   : (c) Sandia National Labs 2020, 2021, 2022, 2023
// Description : Testing nas2exo Library, C++ 14
//============================================================================

#include "N2ENasReader.h"
#include <iostream>
#include <sstream>

namespace NasModules {

  N2ENasReader::N2ENasReader(std::string ifname) : inFileName(std::move(ifname))
  {
    // TODO Auto-generated constructor stub
    if (!doesFileExist(this->inFileName)) {
      std::string msg = "This file does not exist: " + this->inFileName;
      throw std::invalid_argument(msg);
    }

    // C++ 14 version
    // this->inStream = make_unique<std::ifstream>(this->inFileName);
    // C++ 11 version
    this->inStream = std::make_unique<std::ifstream>(this->inFileName);

    // Let's set the buffer
    this->inStream->rdbuf()->pubsetbuf(this->_readBuffer, sizeof(this->_readBuffer));

    this->lineCount = this->lineCounter();
  }

  bool N2ENasReader::doesFileExist(const std::string &fname)
  {
    std::ifstream f(fname);
    return f.good();
  }

  std::string N2ENasReader::getModelTitle() { return this->modelTitle; }

  void N2ENasReader::setModelTitle(const std::string &title)
  {
    // This gymnastics with std::string is b/c
    // a Nastran title is limited to 72 chars.
    this->modelTitle = title.substr(0, 71);
  }

  unsigned N2ENasReader::lineCounter()
  {

    this->inStream->clear();
    this->inStream->seekg(0);

    char     buff[512];
    unsigned count = 0;

    for (;;) {

      this->inStream->getline(buff, 511);
      if (this->inStream->eof()) {
        break;
      }
      ++count;
    }

    return count;
  }

  bool N2ENasReader::processFile()
  {

    bool result = false;

    // Make sure the file pointer is parked at 0
    this->inStream->clear();
    this->inStream->seekg(0);

    char     buff[512];
    unsigned utmp1;
    unsigned utmp2;

    // Blowout the data so fresh data can be added
    this->sections.clear();
    this->elementList.clear();
    this->gridList.clear();

    for (;;) {

      std::vector<std::string> tokens;

      this->inStream->getline(buff, 511);

      if (this->inStream->eof()) {
        result = true;
        break;
      }

      // Ignore comments
      if (buff[0] == '$')
        continue;

      tokens = this->csvLineToTokens(buff);

      unsigned card_id = 0;
      for (int i = 1; i <= 4; i++) {

        if (tokens[0] == N2EFileCues[i]) {
          card_id = i;
          break;
        }
      }

      if (card_id == 0)
        continue;

      sectionType   sType;
      gridType      gType;
      elementType   eType;
      N2EPoint3D    pt;
      N2EGridPtList gLst;

      switch (card_id) {
      case 1:
        // Grid
        std::istringstream(tokens[1]) >> utmp1;
        std::istringstream(tokens[3]) >> pt.x[0];
        std::istringstream(tokens[4]) >> pt.x[1];
        std::istringstream(tokens[5]) >> pt.x[2];
        gType = std::make_tuple(utmp1, pt);
        this->gridList.emplace_back(gType);
        break;
      case 2:
        // CTETRA
        std::istringstream(tokens[1]) >> utmp1;
        std::istringstream(tokens[2]) >> utmp2;
        std::istringstream(tokens[3]) >> gLst.v[0];
        std::istringstream(tokens[4]) >> gLst.v[1];
        std::istringstream(tokens[5]) >> gLst.v[2];
        std::istringstream(tokens[6]) >> gLst.v[3];
        gLst.v[7] = gLst.v[6] = gLst.v[5] = gLst.v[4] = 0;
        eType                                         = std::make_tuple(utmp1, utmp2, 4, gLst);
        this->elementList.emplace_back(eType);
        break;

      case 3:
        // CHEXA 8 Node only!!!!!
        std::istringstream(tokens[1]) >> utmp1;
        std::istringstream(tokens[2]) >> utmp2;
        std::istringstream(tokens[3]) >> gLst.v[0];
        std::istringstream(tokens[4]) >> gLst.v[1];
        std::istringstream(tokens[5]) >> gLst.v[2];
        std::istringstream(tokens[6]) >> gLst.v[3];
        std::istringstream(tokens[7]) >> gLst.v[4];
        std::istringstream(tokens[8]) >> gLst.v[5];
        // Hexa 8 needs two more nodes from next line
        this->inStream->getline(buff, 511);
        tokens = csvLineToTokens(buff);
        std::istringstream(tokens[0]) >> gLst.v[6];
        std::istringstream(tokens[1]) >> gLst.v[7];
        eType = std::make_tuple(utmp1, utmp2, 8, gLst);
        this->elementList.emplace_back(eType);
        break;
      case 4:
        // Section
        std::istringstream(tokens[1]) >> utmp1;
        std::istringstream(tokens[2]) >> utmp2;
        sType = std::make_tuple(utmp1, utmp2);
        this->sections.emplace_back(sType);
      }
    }

    return result;
  }

  std::vector<std::string> N2ENasReader::csvLineToTokens(char buff[])
  {
    std::vector<std::string> toks;
    std::stringstream        sStream(buff);
    std::string              tmp;

    while (getline(sStream, tmp, ',')) {
      toks.push_back(tmp);
    }

    return toks;
  }

} // namespace NasModules
