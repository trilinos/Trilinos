// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystParserInterface.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

int main(int argc, char **argv)
{
  if ((argc != 2) && (argc != 3)) {
    std::cerr << argv[0] << " <file path to file to parse> [file path for json output]\n";
    exit(1);
  }

  std::ifstream fh(argv[1]);
  std::string   input;
  std::string   separator = "_";

  CatalystParserInterface::var_map  ev;
  CatalystParserInterface::var_map  nv;
  CatalystParserInterface::var_map  gv;
  CatalystParserInterface::id_range er = std::make_pair(1, 10000000000);
  CatalystParserInterface::id_range nr = std::make_pair(1, 10000000000);
  if (fh.is_open()) {
    std::string line;
    bool        icb = false;
    int         obc = 0;
    while (std::getline(fh, line)) {
      std::istringstream iss(line);
      std::string        tok1;
      std::string        tok2;
      std::string        tok3;
      std::string        tok4;
      std::string        tok5;
      std::string        tok6;

      iss >> tok1;
      iss >> tok2;
      std::transform(tok1.begin(), tok1.end(), tok1.begin(), ::tolower);
      std::transform(tok2.begin(), tok2.end(), tok2.begin(), ::tolower);

      if ((tok1 == "begin") && (tok2 == "catalyst"))
        icb = true;

      if (icb && tok1 == "begin")
        obc += 1;

      if (icb && tok1 == "end")
        obc -= 1;

      if (icb) {
        iss >> tok3 >> tok4 >> tok5 >> tok6;
        std::string tok3l = tok3;
        std::string tok4l = tok4;
        std::transform(tok3l.begin(), tok3l.end(), tok3l.begin(), ::tolower);
        std::transform(tok4l.begin(), tok4l.end(), tok4l.begin(), ::tolower);

        if ((tok1 == "x" || tok1 == "y" || tok1 == "z") && tok2 == "axis" && tok3l == "label" &&
            tok4l == "name" && tok5 == "=")
          input += tok1 + " axis label name = \"" + tok6 + "\"\n";
        else if (tok1 == "image" && tok2 == "name" && tok3l == "addon" && tok4 == "=")
          input += "image name addon = \"" + tok5 + "\"\n";
        else if (tok1 == "image" && tok2 == "basedirectory" && tok3 == "=")
          input += "image basedirectory = \"" + tok4 + "\"\n";
        else if (tok1 == "image" && tok2 == "basename" && tok3 == "=")
          input += "image basename = \"" + tok4 + "\"\n";
        else if (tok1 == "plot" && tok2 == "basedirectory" && tok3 == "=")
          input += "plot basedirectory = \"" + tok4 + "\"\n";
        else if (tok1 == "plot" && tok2 == "basename" && tok3 == "=")
          input += "plot basename = \"" + tok4 + "\"\n";
        else if (tok1 == "function" && tok2 == "=")
          input += "function = \"" + tok3 + "\"\n";
        else
          input += line + "\n";
      }
      else {
        iss >> tok3 >> tok4 >> tok5 >> tok6;
        std::string tok3l = tok3;
        std::string tok4l = tok4;
        std::string tok5l = tok5;
        std::transform(tok3l.begin(), tok3l.end(), tok3l.begin(), ::tolower);
        std::transform(tok4l.begin(), tok4l.end(), tok4l.begin(), ::tolower);
        std::transform(tok5l.begin(), tok5l.end(), tok5l.begin(), ::tolower);
        if (tok1 == "component" && tok2 == "separator" && tok3l == "character" && tok4 == "=")
          separator = tok5;
        if (tok1 == "component" && tok2 == "separator" && tok3l == "character" && tok4 == "=" &&
            tok5l == "none")
          separator = "";
        if (tok1 == "element" && tok2 == "variables" && tok3 == "=" && tok5l == "as")
          ev[tok6] = CatalystParserInterface::ALLTYPE;
        if (tok1 == "element" && tok2 == "variables" && tok3 == "=")
          ev[tok4] = CatalystParserInterface::ALLTYPE;
        if (tok1 == "nodal" && tok2 == "variables" && tok3 == "=" && tok5l == "as")
          nv[tok6] = CatalystParserInterface::ALLTYPE;
        if (tok1 == "nodal" && tok2 == "variables" && tok3 == "=")
          nv[tok4] = CatalystParserInterface::ALLTYPE;
        if (tok1 == "global" && tok2 == "variables" && tok3 == "=" && tok5l == "as")
          gv[tok6] = CatalystParserInterface::ALLTYPE;
        if (tok1 == "global" && tok2 == "variables" && tok3 == "=")
          gv[tok4] = CatalystParserInterface::ALLTYPE;
      }

      if (icb && !obc)
        icb = false;
    }
  }
  else {
    std::cerr << argv[0] << " unable to open file for input: " + std::string(argv[1]) << "\n";
    exit(1);
  }

  fh.close();

  separator.erase(std::remove(separator.begin(), separator.end(), ' '), separator.end());

  int                                 ret;
  CatalystParserInterface::parse_info pinfo;
  pinfo.node_vars    = &nv;
  pinfo.element_vars = &ev;
  pinfo.global_vars  = &gv;
  pinfo.nodeIDs      = &nr;
  pinfo.elementIDs   = &er;
  pinfo.separator    = separator;
  ret                = CatalystParserInterface::parseString(input, pinfo);

  if (!ret && (argc == 3)) {
    std::ofstream file(argv[2]);
    if (file.is_open()) {
      file << pinfo.json_result;
      file.close();
    }
    else {
      std::cerr << argv[0] << " unable to open file for json output: " + std::string(argv[2])
                << "\n";
      exit(1);
    }
  }
  exit(ret);
}
