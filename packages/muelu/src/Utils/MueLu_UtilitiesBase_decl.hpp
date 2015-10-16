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
#ifndef MUELU_UTILITIESBASE_DECL_HPP
#define MUELU_UTILITIESBASE_DECL_HPP

#include <unistd.h> //necessary for "sleep" function in debugging methods
#include <string>

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Xpetra_BlockedCrsMatrix_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_ExportFactory.hpp>

#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>


#include "MueLu_Exceptions.hpp"




namespace MueLu {

  /*!
    @class Utilities
    @brief MueLu utility class.

    This class provides a number of static helper methods. Some are temporary and will eventually
    go away, while others should be moved to Xpetra.
  */
  template <class Scalar,
            class LocalOrdinal  = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node          = KokkosClassic::DefaultNode::DefaultNodeType>
  class UtilitiesBase {
  public:
#undef MUELU_UTILITIESBASE_SHORT
//#include "MueLu_UseShortNames.hpp"
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;


    static RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >                Crs2Op(RCP<Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Op) {
      if (Op.is_null())
        return Teuchos::null;
      return rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Op));
    }

    /*! @brief Extract Matrix Diagonal

    Returns Matrix diagonal in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static Teuchos::ArrayRCP<Scalar> GetMatrixDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A) {
      size_t numRows = A.getRowMap()->getNodeNumElements();
      Teuchos::ArrayRCP<Scalar> diag(numRows);
      Teuchos::ArrayView<const LocalOrdinal> cols;
      Teuchos::ArrayView<const Scalar> vals;
      for (size_t i = 0; i < numRows; ++i) {
        A.getLocalRowView(i, cols, vals);
        LocalOrdinal j = 0;
        for (; j < cols.size(); ++j) {
          if (Teuchos::as<size_t>(cols[j]) == i) {
            diag[i] = vals[j];
            break;
          }
        }
        if (j == cols.size()) {
          // Diagonal entry is absent
          diag[i] = Teuchos::ScalarTraits<Scalar>::zero();
        }
      }
      return diag;
    }

    /*! @brief Extract Matrix Diagonal

    Returns inverse of the Matrix diagonal in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > GetMatrixDiagonalInverse(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps()*100) {
      RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap = A.getRowMap();
      RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > diag = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(rowMap);
      ArrayRCP<Scalar> diagVals = diag->getDataNonConst(0);
      size_t numRows = rowMap->getNodeNumElements();
      Teuchos::ArrayView<const LocalOrdinal> cols;
      Teuchos::ArrayView<const Scalar> vals;
      for (size_t i = 0; i < numRows; ++i) {
        A.getLocalRowView(i, cols, vals);
        LocalOrdinal j = 0;
        for (; j < cols.size(); ++j) {
          if (Teuchos::as<size_t>(cols[j]) == i) {
            if(Teuchos::ScalarTraits<Scalar>::magnitude(vals[j]) > tol)
              diagVals[i] = Teuchos::ScalarTraits<Scalar>::one() / vals[j];
            else
              diagVals[i]=Teuchos::ScalarTraits<Scalar>::zero();
            break;
          }
        }
        if (j == cols.size()) {
          // Diagonal entry is absent
          diagVals[i]=Teuchos::ScalarTraits<Scalar>::zero();
        }
      }
      diagVals=null;
      return diag;
    }



    /*! @brief Extract Matrix Diagonal of lumped matrix

    Returns Matrix diagonal of lumped matrix in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static Teuchos::ArrayRCP<Scalar> GetLumpedMatrixDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A) {
      size_t numRows = A.getRowMap()->getNodeNumElements();
      Teuchos::ArrayRCP<Scalar> diag(numRows);
      Teuchos::ArrayView<const LocalOrdinal> cols;
      Teuchos::ArrayView<const Scalar> vals;
      for (size_t i = 0; i < numRows; ++i) {
        A.getLocalRowView(i, cols, vals);
        diag[i] = Teuchos::ScalarTraits<Scalar>::zero();
        for (LocalOrdinal j = 0; j < cols.size(); ++j) {
          diag[i] += Teuchos::ScalarTraits<Scalar>::magnitude(vals[j]);
        }
      }
      return diag;
    }

    /*! @brief Extract Overlapped Matrix Diagonal

    Returns overlapped Matrix diagonal in ArrayRCP.

    The local overlapped diagonal has an entry for each index in A's column map.
    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > GetMatrixOverlappedDiagonal(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A) {
      RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap = A.getRowMap(), colMap = A.getColMap();
      RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > localDiag = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(rowMap);

      try {
         const Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>* crsOp = dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>*>(&A);
         if (crsOp == NULL) {
           throw Exceptions::RuntimeError("cast to CrsMatrixWrap failed");
         }
         Teuchos::ArrayRCP<size_t> offsets;
         crsOp->getLocalDiagOffsets(offsets);
         crsOp->getLocalDiagCopy(*localDiag,offsets());
      }
      catch (...) {
        ArrayRCP<Scalar>   localDiagVals = localDiag->getDataNonConst(0);
        Teuchos::ArrayRCP<Scalar> diagVals = GetMatrixDiagonal(A);
        for (LocalOrdinal i = 0; i < localDiagVals.size(); i++)
          localDiagVals[i] = diagVals[i];
        localDiagVals = diagVals = null;
      }

      RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > diagonal = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(colMap);
      RCP< const Xpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > importer;
      importer = A.getCrsGraph()->getImporter();
      if (importer == Teuchos::null) {
        importer = Xpetra::ImportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(rowMap, colMap);
      }
      diagonal->doImport(*localDiag, *(importer), Xpetra::INSERT);
      return diagonal;
    }


    // TODO: should NOT return an Array. Definition must be changed to:
    // - ArrayRCP<> ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS)
    // or
    // - void ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS, Array &)
    static Teuchos::Array<Magnitude> ResidualNorm(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS) {
      TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != RHS.getNumVectors(), Exceptions::RuntimeError, "Number of solution vectors != number of right-hand sides")
       const size_t numVecs = X.getNumVectors();
       RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > RES = Residual(Op, X, RHS);
       Teuchos::Array<Magnitude> norms(numVecs);
       RES->norm2(norms);
       return norms;
    }

    static RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Residual(const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Op, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& RHS) {
      TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != RHS.getNumVectors(), Exceptions::RuntimeError, "Number of solution vectors != number of right-hand sides")
        const size_t numVecs = X.getNumVectors();
        Scalar one = Teuchos::ScalarTraits<Scalar>::one(), negone = -one, zero = Teuchos::ScalarTraits<Scalar>::zero();
        RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > RES = Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(Op.getRangeMap(), numVecs, false); // no need to initialize to zero
        Op.apply(X, *RES, Teuchos::NO_TRANS, one, zero);
        RES->update(one, RHS, negone);
        return RES;
    }

#ifndef _WIN32
#include <unistd.h>
    static void PauseForDebugger() {
      RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
      int myPID = comm->getRank();
      int pid   = getpid();
      char hostname[80];
      for (int i = 0; i <comm->getSize(); i++) {
        if (i == myPID) {
          gethostname(hostname, sizeof(hostname));
          std::cout << "Host: " << hostname << "\tMPI rank: " << myPID << ",\tPID: " << pid << "\n\tattach " << pid << std::endl;
          sleep(1);
        }
      }
      if (myPID == 0) {
        std::cout << "** Enter a character to continue > " << std::endl;
        char go = ' ';
        int r = scanf("%c", &go);
        (void)r;
        assert(r > 0);
      }
      comm->barrier();
    }
#else
    static void PauseForDebugger() {
         throw(Exceptions::RuntimeError("MueLu Utils: PauseForDebugger not implemented on Windows."));
     }
#endif

    /*! @brief Simple transpose for Tpetra::CrsMatrix types

        Note:  This is very inefficient, as it inserts one entry at a time.
    */

    /*! @brief Power method.

    @param A matrix
    @param scaleByDiag if true, estimate the largest eigenvalue of \f$ D^; A \f$.
    @param niters maximum number of iterations
    @param tolerance stopping tolerance
    @verbose if true, print iteration information

    (Shamelessly grabbed from tpetra/examples.)
    */
    static Scalar PowerMethod(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, bool scaleByDiag = true,
                              LocalOrdinal niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123) {
      TEUCHOS_TEST_FOR_EXCEPTION(!(A.getRangeMap()->isSameAs(*(A.getDomainMap()))), Exceptions::Incompatible,
          "Utils::PowerMethod: operator must have domain and range maps that are equivalent.");

      // Create three vectors, fill z with random numbers
      RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > q = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(A.getDomainMap());
      RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > r = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(A.getRangeMap());
      RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > z = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(A.getRangeMap());

      z->setSeed(seed);  // seed random number generator
      z->randomize(true);// use Xpetra implementation: -> same results for Epetra and Tpetra

      Teuchos::Array<Magnitude> norms(1);

      typedef Teuchos::ScalarTraits<Scalar> STS;

      const Scalar zero = STS::zero(), one = STS::one();

      Scalar lambda = zero;
      Magnitude residual = STS::magnitude(zero);

      // power iteration
      RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > diagInvVec;
      if (scaleByDiag) {
        RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > diagVec = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(A.getRowMap());
        A.getLocalDiagCopy(*diagVec);
        diagInvVec = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(A.getRowMap());
        diagInvVec->reciprocal(*diagVec);
      }

      for (int iter = 0; iter < niters; ++iter) {
        z->norm2(norms);                                  // Compute 2-norm of z
        q->update(one/norms[0], *z, zero);                // Set q = z / normz
        A.apply(*q, *z);                                  // Compute z = A*q
        if (scaleByDiag)
          z->elementWiseMultiply(one, *diagInvVec, *z, zero);
        lambda = q->dot(*z);                              // Approximate maximum eigenvalue: lamba = dot(q,z)

        if (iter % 100 == 0 || iter + 1 == niters) {
          r->update(1.0, *z, -lambda, *q, zero);          // Compute A*q - lambda*q
          r->norm2(norms);
          residual = STS::magnitude(norms[0] / lambda);
          if (verbose) {
            std::cout << "Iter = " << iter
                      << "  Lambda = " << lambda
                      << "  Residual of A*q - lambda*q = " << residual
                      << std::endl;
          }
        }
        if (residual < tolerance)
          break;
      }
      return lambda;
    }



    static RCP<Teuchos::FancyOStream> MakeFancy(std::ostream& os) {
      RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(os));
      return fancy;
    }

    /*! @brief Squared distance between two rows in a multivector

       Used for coordinate vectors.
    */
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType Distance2(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& v, LocalOrdinal i0, LocalOrdinal i1) {
      size_t numVectors = v.getNumVectors();

      Scalar d = Teuchos::ScalarTraits<Scalar>::zero();
      for (size_t j = 0; j < numVectors; j++) {
        Teuchos::ArrayRCP<const Scalar> vv = v.getData(j);
        d += (vv[i0] - vv[i1])*(vv[i0] - vv[i1]);
      }
      return Teuchos::ScalarTraits<Scalar>::magnitude(d);
    }

    /*! @brief Detect Dirichlet rows

        @param[in] A matrix
        @param[in] tol If a row entry's magnitude is less than or equal to this tolerance, the entry is treated as zero.

        @return boolean array.  The ith entry is true iff row i is a Dirichlet row.
    */
    static Teuchos::ArrayRCP<const bool> DetectDirichletRows(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A, const Magnitude& tol = Teuchos::ScalarTraits<Scalar>::zero()) {
      LocalOrdinal numRows = A.getNodeNumRows();
      typedef Teuchos::ScalarTraits<Scalar> STS;
      ArrayRCP<bool> boundaryNodes(numRows, true);
      for (LocalOrdinal row = 0; row < numRows; row++) {
        ArrayView<const LocalOrdinal> indices;
        ArrayView<const Scalar> vals;
        A.getLocalRowView(row, indices, vals);
        size_t nnz = A.getNumEntriesInLocalRow(row);
        if (nnz > 1)
          for (size_t col = 0; col < nnz; col++)
            if ( (indices[col] != row) && STS::magnitude(vals[col]) > tol) {
              boundaryNodes[row] = false;
              break;
            }
      }
      return boundaryNodes;
    }

    /*! @brief Set seed for random number generator.

      Distribute the seeds evenly in [1,INT_MAX-1].  This guarantees nothing
      about where random number streams on difference processes will intersect.
      This does avoid overflow situations in parallel when multiplying by a PID.
      It also avoids the pathological case of having the *same* random number stream
      on each process.
    */

    static void SetRandomSeed(const Teuchos::Comm<int> &comm) {
      // Distribute the seeds evenly in [1,maxint-1].  This guarantees nothing
      // about where in random number stream we are, but avoids overflow situations
      // in parallel when multiplying by a PID.  It would be better to use
      // a good parallel random number generator.
      double one = 1.0;
      int maxint = INT_MAX; //= 2^31-1 = 2147483647 for 32-bit integers
      int mySeed = Teuchos::as<int>((maxint-1) * (one -(comm.getRank()+1)/(comm.getSize()+one)) );
      if (mySeed < 1 || mySeed == maxint) {
        std::ostringstream errStr;
        errStr << "Error detected with random seed = " << mySeed << ". It should be in the interval [1,2^31-2].";
        throw Exceptions::RuntimeError(errStr.str());
      }
      std::srand(mySeed);
      // For Tpetra, we could use Kokkos' random number generator here.
      Teuchos::ScalarTraits<Scalar>::seedrandom(mySeed);
      // Epetra
      //   MultiVector::Random() -> Epetra_Util::RandomDouble() -> Epetra_Utils::RandomInt()
      // Its own random number generator, based on Seed_. Seed_ is initialized in Epetra_Util constructor with std::rand()
      // So our setting std::srand() affects that too
    }


  }; // class Utils



  ///////////////////////////////////////////


#if 0 // not sure whether we want to have the MueMex routines here!
  /*! Removes the following non-serializable data (A,P,R,Nullspace,Coordinates) from level-specific sublists from inList
    and moves it to nonSerialList.  Everything else is copied to serialList.  This function returns the level number of the highest level
    for which non-serializable data was provided.
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

  /*! Tokenizes a (comma)-separated string, removing all leading and trailing whitespace
    WARNING: This routine is not threadsafe on most architectures
  */
  void TokenizeStringAndStripWhiteSpace(const std::string & stream, std::vector<std::string> & tokenList, const char* delimChars = ",") {
    //note: default delimiter string is ","
    // Take a comma-separated list and tokenize it, stripping out leading & trailing whitespace.  Then add to tokenList
    char* buf = (char*) malloc(stream.size() + 1);
    strcpy(buf, stream.c_str());
    char* token = strtok(buf, delimChars);
    if(token == NULL) {
      free(buf);
      return;
    }
    while(token) {
      //token points to start of string to add to tokenList
      //remove front whitespace...
      char* tokStart = token;
      char* tokEnd = token + strlen(token) - 1;
      while(*tokStart == ' ' && tokStart < tokEnd)
        tokStart++;
      while(*tokEnd == ' ' && tokStart < tokEnd)
        tokEnd--;
      tokEnd++;
      if(tokStart < tokEnd) {
        std::string finishedToken(tokStart, tokEnd - tokStart); //use the constructor that takes a certain # of chars
        tokenList.push_back(finishedToken);
      }
      token = strtok(NULL, delimChars);
    }
    free(buf);
  }

  /*! Returns true if a parameter name is a valid Muemex custom level variable, e.g. "MultiVector myArray"
  */
  bool IsParamMuemexVariable(const std::string& name) {
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
    while(*tolowerIt) {
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
        strstr(firstWord, "string")) {
      //Add name to list of keys to remove
      free(str);
      return true;
    }
    else {
      free(str);
      return false;
    }
  }

  //! Little helper function to convert non-string types to strings
  template<class T>
  std::string toString(const T& what) {
    std::ostringstream buf;
    buf << what;
    return buf.str();
  }
#endif



} //namespace MueLu

#define MUELU_UTILITIESBASE_SHORT
#endif // MUELU_UTILITIESBASE_DECL_HPP
