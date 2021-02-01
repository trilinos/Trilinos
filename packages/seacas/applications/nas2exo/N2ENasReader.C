//============================================================================
// Name        : testnas2exo.cpp
// Author      : Ramon J. Moral (STRA LLC), John Niederhaus (Coordinator, SNL)
// Version     :
// Copyright   : (c) Sandia National Labs 2020
// Description : Testing nas2exo Library, C++ 14
//============================================================================

#include "N2ENasReader.h"
#include <iostream>
#include <sstream>

namespace NasModules {

  N2ENasReader::N2ENasReader(string ifname)
  {
    // TODO Auto-generated constructor stub
    this->inFileName = ifname;

    if (!doesFileExist(this->inFileName)) {

      // Big ass Error
      string msg = "This file does not exist: " + ifname;
      throw std::invalid_argument(msg);
    }

    // C++ 14 version
    // this->inStream = make_unique<ifstream>(this->inFileName);
    // C++ 11 version
    this->inStream.reset(new ifstream(this->inFileName));

    // Let's set the buffer
    this->inStream->rdbuf()->pubsetbuf(this->_readBuffer, sizeof(this->_readBuffer));

    this->lineCount = this->lineCounter();
  }

  N2ENasReader::~N2ENasReader()
  {

    // TODO Auto-generated destructor stub
  }

  bool N2ENasReader::doesFileExist(string fname)
  {

    ifstream f(fname);

    return f.good();
  }

  string N2ENasReader::getModelTitle() { return this->modelTitle; }

  void N2ENasReader::setModelTitle(string title)
  {

    string stmp = title;

    if (title.length() >= 72) {
      stmp = title.substr(0, 71);
    }

    // This gymnastics with string is b/c
    // a Nastran title is limited to 72 chars.
    strncat(this->modelTitle, stmp.c_str(), 71);
  }

  unsigned N2ENasReader::lineCounter()
  {

    this->inStream->clear();
    this->inStream->seekg(0);

    char     buff[512];
    unsigned count = 0;

    for (;;) {

      this->inStream->getline(buff, 511);
      if (this->inStream->eof())
        break;
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

      vector<string> tokens;

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
        istringstream(tokens[1]) >> utmp1;
        istringstream(tokens[3]) >> pt.x[0];
        istringstream(tokens[4]) >> pt.x[1];
        istringstream(tokens[5]) >> pt.x[2];
        gType = make_tuple(utmp1, pt);
        this->gridList.emplace_back(gType);
        break;
      case 2:
        // CTETRA
        istringstream(tokens[1]) >> utmp1;
        istringstream(tokens[2]) >> utmp2;
        istringstream(tokens[3]) >> gLst.v[0];
        istringstream(tokens[4]) >> gLst.v[1];
        istringstream(tokens[5]) >> gLst.v[2];
        istringstream(tokens[6]) >> gLst.v[3];
        gLst.v[7] = gLst.v[6] = gLst.v[5] = gLst.v[4] = 0;
        eType                                         = make_tuple(utmp1, utmp2, 4, gLst);
        this->elementList.emplace_back(eType);
        break;

      case 3:
        // CHEXA 8 Node only!!!!!
        istringstream(tokens[1]) >> utmp1;
        istringstream(tokens[2]) >> utmp2;
        istringstream(tokens[3]) >> gLst.v[0];
        istringstream(tokens[4]) >> gLst.v[1];
        istringstream(tokens[5]) >> gLst.v[2];
        istringstream(tokens[6]) >> gLst.v[3];
        istringstream(tokens[7]) >> gLst.v[4];
        istringstream(tokens[8]) >> gLst.v[5];
        // Hexa 8 needs two more nodes from next line
        this->inStream->getline(buff, 511);
        tokens = csvLineToTokens(buff);
        istringstream(tokens[0]) >> gLst.v[6];
        istringstream(tokens[1]) >> gLst.v[7];
        eType = make_tuple(utmp1, utmp2, 8, gLst);
        this->elementList.emplace_back(eType);
        break;
      case 4:
        // Section
        istringstream(tokens[1]) >> utmp1;
        istringstream(tokens[2]) >> utmp2;
        sType = make_tuple(utmp1, utmp2);
        this->sections.emplace_back(sType);
      }
    }

    return result;
  }

  vector<string> N2ENasReader::csvLineToTokens(char buff[])
  {
    vector<string> toks;
    stringstream   sStream(buff);
    string         tmp;

    while (getline(sStream, tmp, ',')) {
      toks.push_back(tmp);
    }

    return toks;
  }

} // namespace NasModules
