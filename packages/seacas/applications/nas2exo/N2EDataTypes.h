/*
 * N2EDatatTypes.h
 *
 *  Created on: Oct 10, 2020
 *      Author: Ramon J. Moral(Contractor, STRA LLC)
 * 		John Niederhouse(ORG 1443, SNL, Coordinator)
 *  Copyright: Sandia National Labs, OCT-2020
 */

#ifndef INCLUDE_N2EDATATYPES_H_
#define INCLUDE_N2EDATATYPES_H_

#include <cstring>
#include <exodusII.h>
#include <fstream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

using namespace std;

namespace N2EModules {

  const string N2EFileCues[5]{"UNSP", "GRID", "CTETRA", "CHEXA", "PSOLID"};
  // These are the supported element types for this
  // application.

  struct supportedElements
  {

    ex_entity_type elementType;
    char           elemDesc[MAX_STR_LENGTH]{'\0'};
    int64_t        numNodesPerElem;
    int64_t        numEdgesPerElem;
    int64_t        numFacesPerElem;
    int64_t        numAttrPerElem;

    supportedElements(ex_entity_type elType, string elDesc, int64_t nodesPer, int64_t edgesPer,
                      int64_t facesPer, int64_t attrPer)
    {

      elementType = elType;
      strncpy(elemDesc, elDesc.c_str(), MAX_STR_LENGTH - 1);
      numNodesPerElem = nodesPer;
      numEdgesPerElem = edgesPer;
      numFacesPerElem = facesPer;
      numAttrPerElem  = attrPer;
    }
  };

  const supportedElements ExoElTypes[] = {
      supportedElements(ex_entity_type::EX_ELEM_BLOCK, string("TET"), 4, 6, 4, 1),
      supportedElements(ex_entity_type::EX_ELEM_BLOCK, string("HEX"), 8, 12, 6, 1)};

  struct N2EPoint3D
  {
    double x[3];
  };
  struct N2EGridPtList
  {
    double v[8];
  };
  using sectionType = tuple<unsigned /*propID*/, unsigned /*MATID*/>;
  using gridType    = tuple<unsigned /*GRIDID*/, N2EPoint3D>;
  using elementType =
      tuple<unsigned /*ID*/, unsigned /*propId*/, unsigned /*numNodes*/, N2EGridPtList /*nodes*/>;

} // namespace N2EModules
#endif /* INCLUDE_N2EDATATYPES_H_ */
