/*
 * N2EExoWriter.h
 *
 *  Created on: Oct 10, 2020
 *      Author: Ramon J. Moral(Contractor, STRA LLC)
 * 		John Niederhouse(ORG 1443, SNL, Coordinator)
 *  Copyright: Sandia National Labs, OCT-2020
 */

#ifndef _EXOWRITER_H_
#define _EXOWRITER_H_

#include "N2EDataTypes.h"
#include "exodusII.h"
#include <memory>

using namespace N2EModules;

namespace N2ETestingSpace {
  class ExoWriterTester;
}

namespace ExoModules {

  class N2EExoWriter
  {

  public:
    N2EExoWriter();
    virtual ~N2EExoWriter();

    bool createDB(string name);
    bool setNodes(vector<gridType> gridpts);
    bool setElements(vector<elementType> elist);
    bool setSections(vector<sectionType> sList);

    bool writeFile();

    void setModelTitle(string title);

    inline unsigned getBlocksOut() { return this->writtenBlocks; };
    inline unsigned getNodesOut() { return this->writtenNodes; };
    inline unsigned getTetsOut() { return this->writtenTets; };
    inline unsigned getHexesOut() { return this->writtenHexes; };

  protected:
    vector<sectionType> sections;
    vector<gridType>    gridList;
    vector<elementType> elementList;

    char modelTitle[MAX_LINE_LENGTH]{'\0'};

    int exoFileID;

    int CPU_ws = sizeof(double);
    int IO_ws  = sizeof(double);

    // Data persisteance functions
    bool writeCoords();
    bool writeFileParams();
    bool writeElements();

    unsigned writtenBlocks;
    unsigned writtenNodes;
    unsigned writtenTets;
    unsigned writtenHexes;
    // A Friend for testing
    friend class N2ETestingSpace::ExoWriterTester;
  };

} // namespace ExoModules

#endif /* INCLUDE_N2EEXOWRITER_H_ */
