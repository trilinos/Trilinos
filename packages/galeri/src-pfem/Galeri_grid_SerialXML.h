// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
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
    map<string, Galeri::grid::Loadable>
    read(const Epetra_Comm& comm, const string& XMLFileName)
    {
      // read only on processor 0
      
      FileInputSource fileSrc(XMLFileName);
      Teuchos::XMLObject fileXML(fileSrc.getObject());

      map<string, Galeri::grid::Loadable> patches;

      int NumDimensions = fileXML.getRequiredInt("NumDimensions");
      Galeri::core::Workspace::setNumDimensions(NumDimensions);

      for (int i = 0; i < fileXML.numChildren(); ++i)
      {
        const XMLObject& child = fileXML.getChild(i);
        string tag = child.getTag();

        if (tag == "Patch")
        {
          string Label = child.getRequired("Label");
          string ElementType = child.getRequired("ElementType");
          Galeri::grid::Loadable patch;

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

}; // namespace Galeri

#endif
