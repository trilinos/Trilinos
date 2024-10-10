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

#ifndef EPETRAEXT_XMLWRITER_H
#define EPETRAEXT_XMLWRITER_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

#include "EpetraExt_ConfigDefs.h"
#include "Teuchos_RCP.hpp"
#include <fstream>
#include <vector>

class Epetra_Map;
class Epetra_Comm;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_RowMatrix;
namespace Teuchos {
  class FileXML;
  class XMLObject;
  class ParameterList;
}

namespace EpetraExt
{
/*!
\brief class XMLWriter: A class for writing Trilinos objects to XML files.

Class EpetraExt::XMLWriter writes several Trilinos objects in an XML-compatible format.
The list of supported objects contains:
- Epetra_Map;
- Epetra_MultiVector;
- Epetra_CrsGraph;
- Epetra_CrsMatrix;
- Epetra_RowMatrix;
- Teuchos::ParameterList.

All objects can be read and written, with the std::exception of Epetra_RowMatrix
objects, that can only be written to files.

An example of usage is reported in file epetraext/example/inout/XML_IO.cpp.

Writing objects goes as follows. Let \c Map, \c Matrix, \c LHS and \c RHS an
Epetra_Map, Epetra_CrsMatrix, and two Epetra_MultiVector's, respectively. First, we define an XMLWriter object
\code
EpetraExt::XMLWriter XMLWriter(Comm, "data.xml");
\endcode
and we open the file using \c MyProblem label:
\code
XMLWriter.Create("MyProblem");
\endcode
Writing objects simply goes as
\code
XMLWriter.Write("MyMap", Map);
XMLWriter.Write("MyMatrix", Matrix);
XMLWriter.Write("MyLHS", LHS);
XMLWriter.Write("MyRHS", RHS);
\endcode
A \c Teuchos::ParameterList (List), a \c std::string, and a \c std::vector<std::string> can be written as
\code
XMLWriter.Write("MyParameters", List);
XMLWriter.Write("Author", "myself and others");
XMLWriter.Write("Date", "May 2006");
\endcode
Finally, we close the file
\code
XMLWriter.Close();
\endcode
Note that only processor 0 writes the Teuchos::ParameterList, \c std::string, and \c std::vector<std::string>.

The written file is as follows:
\code
<ObjectCollection Label="MyProblem">
<Text Label="Author">
myself and others
</Text>
<Text Label="Date">
May 2006
</Text>
<Map Label="MyMap" NumElements="4" IndexBase="0" NumProc="1" ElementsOnProc0="4">
<Proc ID="0">
0
1
2
3
</Proc>
</Map>
<PointMatrix Label="MyMatrix" Rows="4" Columns="4" Nonzeros="4" Type="double" StartingIndex="0">
0 0 1
1 1 1
2 2 1
3 3 1
</PointMatrix>
<MultiVector Label="MyLHS" Length="4" NumVectors="2" Type="double">
-0.232996 -0.893077
0.0388327 0.0594004
0.661931 0.342299
-0.930856 -0.984604
</MultiVector>
<MultiVector Label="MyRHS" Length="4" NumVectors="2" Type="double">
0 0
0 0
0 0
0 0
</MultiVector>
<Text Label="MyContent">
This is an example of description
The description is as long as desired,
just put it in a std::vector of strings.
</Text>
<List Label="MyParameters">
<ParameterList>
<Parameter name="double parameter" type="double" value="10"/>
<Parameter name="int parameter" type="int" value="10"/>
<Parameter name="std::string parameter" type="std::string" value="std::string"/>
</ParameterList>
</List>
</ObjectCollection>
\endcode

This class requires Teuchos to be configured with the option \c --enable-teuchos-expat.

\author Marzio Sala, D-INFK/ETHZ

\date Last updated on 10-May-06.

*/
class XMLWriter
{
  public:
    // @{ \name Constructor and destructor.
    //! ctor
    XMLWriter(const Epetra_Comm& Comm, const std::string& FileName);

    //! dtor
    ~XMLWriter() {}

    //! Creates the file, giving \c Label to the whole object.
    void Create(const std::string& Label);

    //! Closes the file. No Write operations can follow.
    void Close();

    // @}
    // @{ \name Read operations

    //! Writes an Epetra_Map using label \c Label.
    void Write(const std::string& Label, const Epetra_Map& Map);

    //! Writes an Epetra_RowMatrix using label \c Label.
    void Write(const std::string& Label, const Epetra_RowMatrix& Matrix);

    //! Writes an Epetra_MultiVector using label \c Label.
    void Write(const std::string& Label, const Epetra_MultiVector& MultiVector);

    //! Writes the std::vector of std::string's using label \c Label.
    void Write(const std::string& Label, const std::vector<std::string>& Content);

    //! Writes input std::string using label \c Label.
    void Write(const std::string& Label, const std::string& Text)
    {
      std::vector<std::string> Content;
      Content.push_back(Text);
      Write(Label, Content);
    }

    //! Writes a Teuchos::ParameterList using label \c Label.
    void Write(const std::string& Label, Teuchos::ParameterList& List);

    // @}
  private:
    //! Epetra communicator.
    const Epetra_Comm& Comm_;
    //! Name of the file.
    std::string FileName_;
    //! If \c true, the file has been successfully opened.
    bool IsOpen_;
};

} // namespace EpetraExt

#endif
