// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_GRID_SERIALXML_H
#define GALERI_GRID_SERIALXML_H

#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_StrUtils.hpp"

#include "Galeri_core_Workspace.h"
#include "Galeri_grid_Loadable.h"

namespace Galeri {

namespace grid {

class SerialXML 
{
  public:
    SerialXML() {}

    ~SerialXML() {}

    static
    map<std::string, Galeri::grid::Loadable>
    read(const Epetra_Comm& comm, const std::string& XMLFileName)
    {
      // read only on processor 0
      
      FileInputSource fileSrc(XMLFileName);
      Teuchos::XMLObject fileXML(fileSrc.getObject());

      map<std::string, Galeri::grid::Loadable> patches;

      int NumDimensions = fileXML.getRequiredInt("NumDimensions");
      Galeri::core::Workspace::setNumDimensions(NumDimensions);

      for (int i = 0; i < fileXML.numChildren(); ++i)
      {
        const XMLObject& child = fileXML.getChild(i);
        std::string tag = child.getTag();

        if (tag == "Patch")
        {
          std::string Label = child.getRequired("Label");
          std::string ElementType = child.getRequired("ElementType");
          Galeri::grid::Loadable patch;

          for (int j = 0; j < child.numChildren(); ++j)
          {
            const XMLObject& newChild = child.getChild(j);
            std::string tag = newChild.getTag();
            if (tag == "Elements")
            {
              int rows = newChild.getRequiredInt("rows");
              int cols = newChild.getRequiredInt("cols");

              // assign all elements to processor 0
              if (comm.MyPID()) rows = 0;
              patch.initialize(comm, -1, rows, ElementType);
              patch.setLabel(Label);

              int count = 0;
              if (comm.MyPID() == 0)
              {
                for (int k = 0; k < newChild.numContentLines(); ++k)
                {
                  const std::string& line = newChild.getContentLine(k);
                  Array<std::string> tokens = Teuchos::StrUtils::stringTokenizer(line);
                  if (tokens.size() != cols) continue;
                  for (int kk = 0; kk < cols; ++kk)
                  {
                    patch.setGlobalConnectivity(count, kk, Teuchos::StrUtils::atoi(tokens[kk]));
                  }
                  ++count;
                }
              }

              patch.freezeConnectivity();
            }
            else if (tag == "Vertices")
            {
              int rows = newChild.getRequiredInt("rows");
              int cols = newChild.getRequiredInt("cols");

              if (comm.MyPID() == 0)
              {
                for (int k = 0; k < newChild.numContentLines(); ++k)
                {
                  const std::string& line = newChild.getContentLine(k);
                  Array<std::string> tokens = Teuchos::StrUtils::stringTokenizer(line);
                  if (tokens.size() != NumDimensions + 1) continue;
                  int GVID = Teuchos::StrUtils::atoi(tokens[0]);
                  for (int kk = 0; kk < NumDimensions; ++kk)
                  {
                    double coord = Teuchos::StrUtils::atof(tokens[kk + 1]);
                    patch.setGlobalCoordinates(GVID, kk, coord);
                  }
                }
              }

              patch.freezeCoordinates();
            }
          }
          patches[Label] = patch;
        }
      }

      return(patches);
    }

  protected:

  private:



}; // class SerialXML


};

}; // namespace Galeri

#endif
