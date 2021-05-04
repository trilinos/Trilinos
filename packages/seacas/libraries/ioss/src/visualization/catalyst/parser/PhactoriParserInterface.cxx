// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "PhactoriParserInterface.h"
#include <iostream>
#include <fstream>
#include <sstream>

namespace Iovs {

void parsePhactoriString(const std::string & s,
    PhactoriParserInterface::ParseInfo & pinfo);

void PhactoriParserInterface::parseFile(const std::string & filepath,
    ParseInfo & pinfo) {

    std::ifstream f(filepath);
    if (f) {
        std::ostringstream ss;
        ss << f.rdbuf();
        parsePhactoriString(ss.str(), pinfo);
    }
    else {
        std::cerr << "Unable to open file for parsing: " << filepath << "\n";
    }
}

void PhactoriParserInterface::parseString(const std::string & s,
    ParseInfo & pinfo) {
    parsePhactoriString(s, pinfo);
}

} // namespace Iovs
