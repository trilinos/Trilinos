/*
 * N2EDatatTypes.h
 *
 *  Created on: Oct 10, 2020
 *      Author: Ramon J. Moral(Contractor, STRA LLC)
 * 		John Niederhouse(ORG 1443, SNL, Coordinator)
 *  Copyright: Sandia National Labs, 2020, 2021, 2022, 2023
 */
#pragma once

#include <cstring>
#include <exodusII.h>
#include <string>
#include <tuple>

namespace N2EModules {

  const std::string N2EFileCues[5]{"UNSP", "GRID", "CTETRA", "CHEXA", "PSOLID"};
  // These are the supported element types for this
  // application.

  struct supportedElements
  {

    ex_entity_type elementType{};
    char           elemDesc[MAX_STR_LENGTH]{'\0'};
    int64_t        numNodesPerElem{};
    int64_t        numEdgesPerElem{};
    int64_t        numFacesPerElem{};
    int64_t        numAttrPerElem{};

    supportedElements(ex_entity_type elType, const std::string &elDesc, int64_t nodesPer,
                      int64_t edgesPer, int64_t facesPer, int64_t attrPer)
    {

      elementType = elType;
      strncpy(elemDesc, elDesc.c_str(), MAX_STR_LENGTH - 1);
      elemDesc[MAX_STR_LENGTH - 1] = '\0';
      numNodesPerElem              = nodesPer;
      numEdgesPerElem              = edgesPer;
      numFacesPerElem              = facesPer;
      numAttrPerElem               = attrPer;
    }
  };

  const supportedElements ExoElTypes[] = {
      supportedElements(ex_entity_type::EX_ELEM_BLOCK, std::string("TET4"), 4, 6, 4, 1),
      supportedElements(ex_entity_type::EX_ELEM_BLOCK, std::string("HEX8"), 8, 12, 6, 1)};

  struct N2EPoint3D
  {
    double x[3];
  };
  struct N2EGridPtList
  {
    int v[8];
  };
  using sectionType = std::tuple<unsigned /*propID*/, unsigned /*MATID*/>;
  using gridType    = std::tuple<unsigned /*GRIDID*/, N2EPoint3D>;
  using elementType = std::tuple<unsigned /*ID*/, unsigned /*propId*/, unsigned /*numNodes*/,
                                 N2EGridPtList /*nodes*/>;

} // namespace N2EModules
