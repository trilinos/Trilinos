/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef EPETRAEXT_XMLREADER_H
#define EPETRAEXT_XMLREADER_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

#include "EpetraExt_ConfigDefs.h"
#include "Epetra_ConfigDefs.h"
#include "Teuchos_RCP.hpp"
#include <fstream>
#include <vector>

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
std::vector<std::string> Author;
std::vector<std::string> Date;
std::vector<std::string> MyContent;
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
    XMLReader(const Epetra_Comm& Comm, const std::string& FileName);

    //! dtor
    ~XMLReader() {}

    // @}
    // @{ \name Read operations

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    //! Reads the Epetra_Map stored with label \c Label.
    void Read(const std::string& Label, Epetra_Map*& Map);
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    //! Reads the Epetra_Map stored with label \c Label. Long Long version.
    void Read64(const std::string& Label, Epetra_Map*& Map);
#endif

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    //! Reads the Epetra_CrsGraph stored with label \c Label.
    void Read(const std::string& Label, Epetra_CrsGraph*& Graph);
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    //! Reads the Epetra_CrsGraph stored with label \c Label. Long Long version.
    void Read64(const std::string& Label, Epetra_CrsGraph*& Graph);
#endif

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    //! Reads the Epetra_CrsMatrix stored with label \c Label.
    void Read(const std::string& Label, Epetra_CrsMatrix*& Matrix);
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    //! Reads the Epetra_CrsMatrix stored with label \c Label. Long Long version.
    void Read64(const std::string& Label, Epetra_CrsMatrix*& Matrix);
#endif

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
        //! Reads the Epetra_MultiVector stored with label \c Label.
    void Read(const std::string& Label, Epetra_MultiVector*& MultiVector);
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
        //! Reads the Epetra_MultiVector stored with label \c Label. Long Long version.
    void Read64(const std::string& Label, Epetra_MultiVector*& MultiVector);
#endif

    //! Reads a std::vector of strings with label \c Label.
    void Read(const std::string& Label, std::vector<std::string>& Content);

    //! Reads the Teuchos::ParameterList stored with label \c Label.
    void Read(const std::string& Label, Teuchos::ParameterList& List);

    // @}
  private:
    //! If \c true, then the file has been successfully opened.
    bool IsOpen_;
    //! Communicator object.
    const Epetra_Comm& Comm_;
    //! parsed XML object.
    Teuchos::RCP<Teuchos::XMLObject> fileXML_;
};

} // namespace EpetraExt

#endif
