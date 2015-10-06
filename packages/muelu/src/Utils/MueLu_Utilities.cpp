// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "MueLu_Utilities_def.hpp"

#include <string>

#ifdef HAVE_MUELU_EPETRAEXT
#include "EpetraExt_Transpose_RowMatrix.h"
#endif

namespace MueLu {

  RCP<Xpetra::Matrix<double, int, int> > Utils2<double, int, int>::Transpose(Matrix& Op, bool optimizeTranspose, const std::string & label) {
   typedef double                                           Scalar;
   typedef int                                              LocalOrdinal;
   typedef int                                              GlobalOrdinal;
   typedef KokkosClassic::DefaultNode::DefaultNodeType      Node;

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
    std::string TorE = "epetra";
#else
    std::string TorE = "tpetra";
#endif

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
    try {
      Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2EpetraCrs(Op);
    }
    catch (...) {
      TorE = "tpetra";
    }
#endif

#ifdef HAVE_MUELU_TPETRA
    if (TorE == "tpetra") {
      try {
        const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpetraOp = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraCrs(Op);

        // Compute the transpose A of the Tpetra matrix tpetraOp.
        RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A;
        Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(rcpFromRef(tpetraOp),label);
        A = transposer.createTranspose();
        RCP<Xpetra::TpetraCrsMatrix<SC> > AA   = rcp(new Xpetra::TpetraCrsMatrix<SC>(A));
        RCP<Xpetra::CrsMatrix<SC> >       AAA  = rcp_implicit_cast<Xpetra::CrsMatrix<SC> >(AA);
        RCP<Xpetra::CrsMatrixWrap<SC> >   AAAA = rcp( new Xpetra::CrsMatrixWrap<SC> (AAA));

        return AAAA;
      }
      catch (std::exception& e) {
        std::cout << "threw exception '" << e.what() << "'" << std::endl;
        throw Exceptions::RuntimeError("Utils::Transpose failed, perhaps because matrix is not a Crs matrix");
      }
    } //if
#endif

    if (TorE == "tpetra") {
#ifdef HAVE_MUELU_TPETRA
#else
      throw Exceptions::RuntimeError("Tpetra");
#endif // HAVE_MUELU_TPETRA

    } else {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
      Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("ZZ Entire Transpose"));
      // Epetra case
      Epetra_CrsMatrix& epetraOp = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstEpetraCrs(Op);
      EpetraExt::RowMatrix_Transpose transposer;
      Epetra_CrsMatrix * A = dynamic_cast<Epetra_CrsMatrix*>(&transposer(epetraOp));
      transposer.ReleaseTranspose(); // So we can keep A in Muelu...

      RCP<Epetra_CrsMatrix> rcpA(A);
      RCP<EpetraCrsMatrix>            AA   = rcp(new EpetraCrsMatrix(rcpA));
      RCP<Xpetra::CrsMatrix<SC> >     AAA  = rcp_implicit_cast<Xpetra::CrsMatrix<SC> >(AA);
      RCP<Xpetra::CrsMatrixWrap<SC> > AAAA = rcp( new Xpetra::CrsMatrixWrap<SC>(AAA));
      AAAA->fillComplete(Op.getRangeMap(), Op.getDomainMap());

      return AAAA;
#else
      throw Exceptions::RuntimeError("Epetra (Err. 2)");
#endif
    }
    return Teuchos::null;
  } //Transpose

  // -- ------------------------------------------------------- --

  void Utils2<double,int,int>::MyOldScaleMatrix_Epetra(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector,
                               bool doFillComplete,
                               bool doOptimizeStorage) {
#ifdef HAVE_MUELU_EPETRA
    try {
      const Epetra_CrsMatrix& epOp = Utils<double,int,int>::Op2NonConstEpetraCrs(Op);

      Epetra_Map const &rowMap = epOp.RowMap();
      int nnz;
      double *vals;
      int *cols;

      for (int i = 0; i < rowMap.NumMyElements(); ++i) {
        epOp.ExtractMyRowView(i, nnz, vals, cols);
        for (int j = 0; j < nnz; ++j)
          vals[j] *= scalingVector[i];
      }

    } catch (...){
      throw Exceptions::RuntimeError("Only Epetra_CrsMatrix types can be scaled");
    }
#else
    throw Exceptions::RuntimeError("Matrix scaling is not possible because Epetra has not been enabled.");
#endif // HAVE_MUELU_EPETRA
  } //Utils2::MyOldScaleMatrix_Epetra()

  // -- ------------------------------------------------------- --

  RCP<Xpetra::MultiVector<double,int,int> > Utils2<double,int,int>::ReadMultiVector(const std::string& fileName, const RCP<const Map>& map) {
    Xpetra::UnderlyingLib lib = map->lib();

    if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
      Epetra_MultiVector * MV;
      EpetraExt::MatrixMarketFileToMultiVector(fileName.c_str(), toEpetra(map), MV);
      return Xpetra::toXpetra<int>(rcp(MV));
#else
      throw Exceptions::RuntimeError("MueLu has not been compiled with Epetra and EpetraExt support.");
#endif
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      typedef Tpetra::CrsMatrix<SC,LO,GO,NO>                    sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type>  reader_type;
      typedef Tpetra::Map<LO,GO,NO>                             map_type;
      typedef Tpetra::MultiVector<SC,LO,GO,NO>                  multivector_type;

      RCP<const map_type>   temp = toTpetra(map);
      RCP<multivector_type> TMV  = reader_type::readDenseFile(fileName,map->getComm(),map->getNode(),temp);
      RCP<MultiVector>      rmv  = Xpetra::toXpetra(TMV);
      return rmv;
#else
      throw Exceptions::RuntimeError("MueLu has not been compiled with Tpetra support.");
#endif
    } else {
      throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
    }

    return Teuchos::null;
  }

  RCP<const Xpetra::Map<int,int> > Utils2<double,int,int>::ReadMap(const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int> >& comm) {
    if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
        Epetra_Map *eMap;
        int rv = EpetraExt::MatrixMarketFileToMap(fileName.c_str(), *(Xpetra::toEpetra(comm)), eMap);
        if (rv != 0)
          throw Exceptions::RuntimeError("Error reading matrix with EpetraExt::MatrixMarketToMap (returned " + toString(rv) + ")");

        RCP<Epetra_Map> eMap1 = rcp(new Epetra_Map(*eMap));
        return Xpetra::toXpetra<int>(*eMap1);
#else
        throw Exceptions::RuntimeError("MueLu has not been compiled with Epetra and EpetraExt support.");
#endif
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      typedef Tpetra::CrsMatrix<double,int,int,NO> sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type>                          reader_type;

      RCP<NO> node = rcp(new NO());

      RCP<const Tpetra::Map<int,int,NO> > tMap = reader_type::readMapFile(fileName, comm, node);
      if (tMap.is_null())
        throw Exceptions::RuntimeError("The Tpetra::Map returned from readSparseFile() is null.");

      return Xpetra::toXpetra(tMap);
#else
      throw Exceptions::RuntimeError("MueLu has not been compiled with Tpetra support.");
#endif
    } else {
      throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
    }
  }


  /* Removes the following non-serializable data (A,P,R,Nullspace,Coordinates)
     from level-specific sublists from inList
     and moves it to nonSerialList.  Everything else is copied to serialList.
     This function returns the level number of the highest level for which
     non-serializable data was provided.
  */
  long ExtractNonSerializableData(const Teuchos::ParameterList& inList, Teuchos::ParameterList& serialList, Teuchos::ParameterList& nonSerialList) {
    using Teuchos::ParameterList;

    ParameterList dummy;
    long maxLevel = 0;

    for (ParameterList::ConstIterator it = inList.begin(); it != inList.end(); it++) {
      const std::string& levelName = it->first;

      // Check for mach of the form "level X" where X is a positive integer
      if (inList.isSublist(levelName) && levelName.find("level ") == 0 && levelName.size() > 6) {
        int levelID = strtol(levelName.substr(6).c_str(), 0, 0);
        if (maxLevel < levelID)
          maxLevel = levelID;

        // Split the sublist
        const ParameterList& levelList = inList.sublist(levelName);
        for (ParameterList::ConstIterator it2 = levelList.begin(); it2 != levelList.end(); it2++) {
          const std::string& name = it2->first;
          if (name == "A" || name == "P" || name == "R" || name == "Nullspace" || name == "Coordinates")
            nonSerialList.sublist(levelName).setEntry(name, it2->second);
          #ifdef HAVE_MUELU_MATLAB
          else if(IsParamMuemexVariable(name))
          {
            nonSerialList.sublist(levelName).setEntry(name, it2->second);
          }
          #endif
          else
            serialList.sublist(levelName).setEntry(name, it2->second);
        }

      } else {
        serialList.setEntry(it->first, it->second);
      }
    }

    return maxLevel;
  }

  void TokenizeStringAndStripWhiteSpace(const std::string& stream, std::vector<std::string>& tokenList, const char* delimChars)
  {
    //note: default delimiter string is ","
    // Take a comma-separated list and tokenize it, stripping out leading & trailing whitespace.  Then add to tokenList
    char* buf = (char*) malloc(stream.size() + 1);
    strcpy(buf, stream.c_str());
    char* token = strtok(buf, delimChars);
    if(token == NULL)
    {
      free(buf);
      return;
    }
    while(token)
    {
      //token points to start of string to add to tokenList
      //remove front whitespace...
      char* tokStart = token;
      char* tokEnd = token + strlen(token) - 1;
      while(*tokStart == ' ' && tokStart < tokEnd)
        tokStart++;
      while(*tokEnd == ' ' && tokStart < tokEnd)
        tokEnd--;
      tokEnd++;
      if(tokStart < tokEnd)
      {
        std::string finishedToken(tokStart, tokEnd - tokStart); //use the constructor that takes a certain # of chars
        tokenList.push_back(finishedToken);
      }
      token = strtok(NULL, delimChars);
    }
    free(buf);
  }

  bool IsParamMuemexVariable(const std::string& name)
  {
    //see if paramName is exactly two "words" - like "OrdinalVector myNullspace" or something
    char* str = (char*) malloc(name.length() + 1);
    strcpy(str, name.c_str());
    //Strip leading and trailing whitespace
    char* firstWord = strtok(str, " ");
    if(!firstWord)
      return false;
    char* secondWord = strtok(NULL, " ");
    if(!secondWord)
      return false;
    char* thirdWord = strtok(NULL, " ");
    if(thirdWord)
      return false;
    //convert first word to all lowercase for case insensitive compare
    char* tolowerIt = firstWord;
    while(*tolowerIt)
    {
      *tolowerIt = (char) tolower(*tolowerIt);
      tolowerIt++;
    }
    //See if the first word is one of the custom variable names
    if(strstr(firstWord, "matrix") ||
       strstr(firstWord, "multivector") ||
       strstr(firstWord, "map") ||
       strstr(firstWord, "ordinalvector") ||
       strstr(firstWord, "int") ||
       strstr(firstWord, "scalar") ||
       strstr(firstWord, "double") ||
       strstr(firstWord, "complex") ||
       strstr(firstWord, "string"))
      //Add name to list of keys to remove
    {
      free(str);
      return true;
    }
    else
    {
      free(str);
      return false;
    }
  }

} // namespace MueLu
