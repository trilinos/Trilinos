#ifndef PHX_GRID_SERIALXML_H
#define PHX_GRID_SERIALXML_H

#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_StrUtils.hpp"

#include "phx_grid_Loadable.h"

namespace phx {

namespace grid {

class SerialXML 
{
  public:
    SerialXML() {}

    ~SerialXML() {}

    static
    map<string, RefCountPtr<Loadable> > 
    readFile(const Epetra_Comm& Comm, const string& XMLFileName)
    {
      // read on all processors
      
      FileInputSource fileSrc(XMLFileName);
      Teuchos::XMLObject fileXML(fileSrc.getObject());

      map<string, RefCountPtr<phx::grid::Loadable> > patches;

      int NumDimensions = fileXML.getRequiredInt("NumDimensions");

      for (int i = 0; i < fileXML.numChildren(); ++i)
      {
        const XMLObject& child = fileXML.getChild(i);
        string tag = child.getTag();

        if (tag == "Patch")
        {
          string Label = child.getRequired("Label");
          string ElementType = child.getRequired("ElementType");
          RefCountPtr<phx::grid::Loadable> patch;

          RefCountPtr<Epetra_Map> ElementMap;
          RefCountPtr<phx::grid::Element> GridElement;

          for (int j = 0; j < child.numChildren(); ++j)
          {
            const XMLObject& newChild = child.getChild(j);
            string tag = newChild.getTag();
            if (tag == "Elements")
            {
              int rows = newChild.getRequiredInt("rows");
              int cols = newChild.getRequiredInt("cols");

              ElementMap = rcp(new Epetra_Map(rows, 0, Comm));
              if (ElementType == "Triangle")
              {
                GridElement = rcp(new phx::grid::Triangle(NumDimensions));
              }
              else if (ElementType == "Segment")
              {
                GridElement = rcp(new phx::grid::Segment(NumDimensions));
              }
              else
                throw(-1);

              patch = rcp(new phx::grid::Loadable(ElementMap, GridElement));

              int count = 0;
              for (int k = 0; k < newChild.numContentLines(); ++k)
              {
                const string& line = newChild.getContentLine(k);
                Array<string> tokens = Teuchos::StrUtils::stringTokenizer(line);
                if (tokens.size() != cols) continue;
                for (int kk = 0; kk < cols; ++kk)
                {
                  int LID = ElementMap->LID(count);
                  if (LID != -1)
                    patch->ADJ(LID, kk) = Teuchos::StrUtils::atoi(tokens[kk]);
                }
                ++count;
              }

              patch->freezeConnectivity();
            }
            else if (tag == "Vertices")
            {
              int rows = newChild.getRequiredInt("rows");
              int cols = newChild.getRequiredInt("cols");

              if (Comm.MyPID() == 0)
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
                    patch->setGlobalCoordinates(kk, 1, &GVID, &coord);
                  }
                }
              }

              patch->freezeCoordinates();
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
