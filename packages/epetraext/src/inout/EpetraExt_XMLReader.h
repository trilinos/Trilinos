#ifndef EPETRAEXT_XMLREADER_H
#define EPETRAEXT_XMLREADER_H

#include "EpetraExt_ConfigDefs.h"
#include "Teuchos_RefCountPtr.hpp"
#include <fstream>

class Epetra_Map;
class Epetra_Comm;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
namespace Teuchos {
  class FileXML;
  class XMLObject;
  class ParameterList;
}

namespace EpetraExt 
{
/*! 
\brief class XMLReader: A class for reading Epetra objects stored in XML files.

Class EpetraExt::XMLReader allows to read several Trilinos objects stored in
XML files. The XML data format is specified in the documentation of class
EpetraExt::XMLWriter, which also contains a MATLAB script. A typical usage of
this class is reported in file epetraext/example/inout/XML_IO.cpp.

This class requires Teuchos to be configured with the option \c --enable-teuchos-expat.

Reading objects from a file requires the following steps. First, we define an XMLReader object,
\code
EpetraExt::XMLReader XMLReader(Comm, "data.xml");
\endcode
Then, we define a set of pointers,
\code
Epetra_Map* MyMap;
Epetra_CrsMatrix* MyMatrix;
Epetra_MultiVector* MyLHS;
Epetra_MultiVector* MyRHS;
Teuchos::ParameterList MyParameters;
vector<string> Author;
vector<string> Date;
vector<string> MyContent;
\endcode
Reading simply goes as follows:
\code
XMLReader.Read("Author", Author);
XMLReader.Read("Date", Date);
XMLReader.Read("MyMap", MyMap);
XMLReader.Read("MyMatrix", MyMatrix);
XMLReader.Read("MyLHS", MyLHS);
XMLReader.Read("MyRHS", MyRHS);
XMLReader.Read("MyContent", MyContent);
XMLReader.Read("MyParameters", MyParameters);
\endcode 
In distributed environments, Epetra_MultiVector, Epetra_CrsGraph and
Epetra_CrsMatrix objects have a linear distribution. Epetra_Map objects can be
read only when  using the same number of processors used for writing.

\warning All the created objects must be deleted from the user using \c delete.

\author Marzio Sala, D-INFK/ETHZ

\date Last updated on 10-May-06.

*/
class XMLReader
{
  public: 
    // @{ \name Constructor and destructor.
    //! ctor
    XMLReader(const Epetra_Comm& Comm, const string& FileName); 

    //! dtor
    ~XMLReader() {}

    // @}
    // @{ \name Read operations
    
    //! Reads the Epetra_Map stored with label \c Label.
    void Read(const string& Label, Epetra_Map*& Map);

    //! Reads the Epetra_CrsGraph stored with label \c Label.
    void Read(const string& Label, Epetra_CrsGraph*& Graph);

    //! Reads the Epetra_CrsMatrix stored with label \c Label.
    void Read(const string& Label, Epetra_CrsMatrix*& Matrix);

    //! Reads the Epetra_MultiVector stored with label \c Label.
    void Read(const string& Label, Epetra_MultiVector*& MultiVector);

    //! Reads a vector of strings with label \c Label.
    void Read(const string& Label, vector<string>& Content);

    //! Reads the Teuchos::ParameterList stored with label \c Label.
    void Read(const string& Label, Teuchos::ParameterList& List);

    // @}
  private:
    //! If \c true, then the file has been successfully opened.
    bool IsOpen_;
    //! Communicator object.
    const Epetra_Comm& Comm_;
    //! parsed XML object.
    Teuchos::RefCountPtr<Teuchos::XMLObject> fileXML_;
};

} // namespace EpetraExt

#endif
