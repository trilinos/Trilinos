// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "PhactoriParserInterface.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>

int main(int argc, char **argv) {

  if (argc != 2) {
    std::cerr << "Usage: parse_phactori" << " FILE\n\n";
    std::cerr << "DESCRIPTION\n\n";
    std::cerr << "Parses input FILE and writes JSON result to "
        << "standard output.\n\n";
    std::cerr << "FILE\n\n";
    std::cerr << "File containing Phactori input syntax.\n";
    exit(1);
  }

  Iovs::PhactoriParserInterface::ParseInfo pinfo;
  Iovs::PhactoriParserInterface::parseFile(argv[1], pinfo);
  if (pinfo.parseFailed) {
      std::cerr << "Parse of input file failed.\n";
  }
  else {
      std::cout << pinfo.jsonParseResult << "\n";
  }
  exit(pinfo.parseFailed);
}
