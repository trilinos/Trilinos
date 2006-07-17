// @HEADER
// ************************************************************************
//
//                  Galeri Matrix Generation Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef PHX_GRID_SERIALXML_H
#define PHX_GRID_SERIALXML_H

#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_StrUtils.hpp"

#include "phx_core_Utils.h"
#include "phx_grid_Loadable.h"

namespace phx {

namespace grid {

class SerialXML 
{
  public:
    SerialXML() {}

    ~SerialXML() {}

    static
    map<string, phx::grid::Loadable>
    read(const Epetra_Comm& comm, const string& XMLFileName)
    {
      // read only on processor 0
      
      FileInputSource fileSrc(XMLFileName);
      Teuchos::XMLObject fileXML(fileSrc.getObject());

      map<string, phx::grid::Loadable> patches;

      int NumDimensions = fileXML.getRequiredInt("NumDimensions");
      phx::core::Utils::setNumDimensions(NumDimensions);

      for (int i = 0; i < fileXML.numChildren(); ++i)
      {
        const XMLObject& child = fileXML.getChild(i);
        string tag = child.getTag();

        if (tag == "Patch")
        {
          string Label = child.getRequired("Label");
          string ElementType = child.getRequired("ElementType");
          phx::grid::Loadable patch;

          for (int j = 0; j < child.numChildren(); ++j)
          {
            const XMLObject& newChild = child.getChild(j);
            string tag = newChild.getTag();
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
                  const string& line = newChild.getContentLine(k);
                  Array<string> tokens = Teuchos::StrUtils::stringTokenizer(line);
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
                  const string& line = newChild.getContentLine(k);
                  Array<string> tokens = Teuchos::StrUtils::stringTokenizer(line);
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

}; // namespace phx

#endif
