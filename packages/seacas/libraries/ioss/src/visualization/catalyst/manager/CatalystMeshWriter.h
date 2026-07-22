// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __CATALYST_MESH_WRITER_H
#define __CATALYST_MESH_WRITER_H

#include <string>

class vtkDataObject;

namespace Iovs {

  class CatalystMeshWriter
  {

  public:
    CatalystMeshWriter();
    ~CatalystMeshWriter();

    bool outputCatalystMeshOneFileON();
    void setOutputCatalystMeshOneFilePrefix(std::string &prefix);

    bool outputCatalystMeshFilePerProcON();
    void setOutputCatalystMeshFilePerProcPrefix(std::string &prefix);

    void writeCatalystMeshOneFile(vtkDataObject *dobj, int timeStep);
    void writeCatalystMeshFilePerProc(vtkDataObject *dobj, int timeStep);

  private:
    bool        catalystMeshOneFile;
    bool        catalystMeshFilePerProc;
    std::string catalystMeshOneFilePrefix;
    std::string catalystMeshFilePerProcPrefix;
  };

} // namespace Iovs

#endif /* __CATALYST_MESH_WRITER_H */
