// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_Array.hpp>
#include "Xpetra_EpetraCrsMatrix.hpp"

namespace Xpetra {

  EpetraCrsMatrix::EpetraCrsMatrix(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : mtx_(Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra(rowMap), maxNumEntriesPerRow, toEpetra(pftype)))), isFillResumed_(false) { }

  EpetraCrsMatrix::EpetraCrsMatrix(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : isFillResumed_(false)
  {
    Teuchos::Array<int> numEntriesPerRowToAlloc(NumEntriesPerRowToAlloc.begin(), NumEntriesPerRowToAlloc.end()); // convert array of "size_t" to array of "int"
    mtx_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra(rowMap), numEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)));
  }

  EpetraCrsMatrix::EpetraCrsMatrix(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : mtx_(Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra(rowMap), toEpetra(colMap), maxNumEntriesPerRow, toEpetra(pftype)))), isFillResumed_(false) { }

  EpetraCrsMatrix::EpetraCrsMatrix(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : isFillResumed_(false)
  {
    Teuchos::Array<int> numEntriesPerRowToAlloc(NumEntriesPerRowToAlloc.begin(), NumEntriesPerRowToAlloc.end()); // convert array of "size_t" to array of "int"
    mtx_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra(rowMap), toEpetra(colMap), numEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)));
  }

  EpetraCrsMatrix::EpetraCrsMatrix(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node, LocalMatOps > > &graph, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : mtx_(Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra(graph)))), isFillResumed_(false) { }

#ifdef HAVE_XPETRA_EXPERIMENTAL
  EpetraCrsMatrix::EpetraCrsMatrix(const EpetraCrsMatrix& matrix)
    : mtx_(Teuchos::rcp(new Epetra_CrsMatrix(*(matrix.mtx_)))), isFillResumed_(false) { }
#endif


  EpetraCrsMatrix::EpetraCrsMatrix(const Teuchos::RCP<const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
				   const Import<LocalOrdinal,GlobalOrdinal,Node> &importer,
				   const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
				   const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
				   const Teuchos::RCP<Teuchos::ParameterList>& params):
    isFillResumed_(false)
  {
    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrix, *sourceMatrix, tSourceMatrix, "Xpetra::EpetraCrsMatrix constructor only accepts Xpetra::EpetraCrsMatrix as an input argument.");
    XPETRA_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Xpetra::EpetraCrsMatrix constructor only accepts Xpetra::EpetraImport as an input argument.");

    const Epetra_Map* myDomainMap = (domainMap!=Teuchos::null)? &toEpetra(domainMap): 0;
    const Epetra_Map* myRangeMap  = (rangeMap !=Teuchos::null)? &toEpetra(rangeMap) : 0;

    // Follows the Tpetra parameters
    bool restrictComm=false;
    if(!params.is_null()) restrictComm = params->get("Restrict Communicator",restrictComm);

    mtx_ = Teuchos::rcp(new Epetra_CrsMatrix(*tSourceMatrix.getEpetra_CrsMatrix(),*tImporter.getEpetra_Import(),myDomainMap,myRangeMap,restrictComm));
  }

  EpetraCrsMatrix::EpetraCrsMatrix(const Teuchos::RCP<const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
				   const Export<LocalOrdinal,GlobalOrdinal,Node> &exporter,
				   const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
				   const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
				   const Teuchos::RCP<Teuchos::ParameterList>& params):
    isFillResumed_(false)
  {
    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrix, *sourceMatrix, tSourceMatrix, "Xpetra::EpetraCrsMatrix constructor only accepts Xpetra::EpetraCrsMatrix as an input argument.");
    XPETRA_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Xpetra::EpetraCrsMatrix constructor only accepts Xpetra::EpetraExport as an input argument.");

    const Epetra_Map* myDomainMap = (domainMap!=Teuchos::null)? &toEpetra(domainMap): 0;
    const Epetra_Map* myRangeMap  = (rangeMap !=Teuchos::null)? &toEpetra(rangeMap) : 0;

    // Follows the Tpetra parameters
    bool restrictComm=false;
    if(!params.is_null()) restrictComm = params->get("Restrict Communicator",restrictComm);

    mtx_ = Teuchos::rcp(new Epetra_CrsMatrix(*tSourceMatrix.getEpetra_CrsMatrix(),*tExporter.getEpetra_Export(),myDomainMap,myRangeMap,restrictComm));
  }



  void EpetraCrsMatrix::insertGlobalValues(int globalRow, const ArrayView<const int> &cols, const ArrayView<const double> &vals) {
    XPETRA_MONITOR("EpetraCrsMatrix::insertGlobalValues");
    XPETRA_ERR_CHECK(mtx_->InsertGlobalValues(globalRow, vals.size(), vals.getRawPtr(), cols.getRawPtr()));
  }

  void EpetraCrsMatrix::insertLocalValues(int localRow, const ArrayView<const int> &cols, const ArrayView<const double> &vals) {
    XPETRA_MONITOR("EpetraCrsMatrix::insertLocalValues");
    XPETRA_ERR_CHECK(mtx_->InsertMyValues(localRow, vals.size(), vals.getRawPtr(), cols.getRawPtr()));
  }

  void EpetraCrsMatrix::replaceGlobalValues(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &indices, const ArrayView< const Scalar > &values) {
    XPETRA_MONITOR("EpetraCrsMatrix::replaceGlobalValues");

    {
      const std::string tfecfFuncName("replaceGlobalValues");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! isFillActive(), std::runtime_error,
                                            ": Fill must be active in order to call this method.  If you have already "
                                            "called fillComplete(), you need to call resumeFill() before you can "
                                            "replace values.");

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(values.size() != indices.size(),
                                            std::runtime_error, ": values.size() must equal indices.size().");
    }

    XPETRA_ERR_CHECK(mtx_->ReplaceGlobalValues(globalRow, indices.size(), values.getRawPtr(), indices.getRawPtr()));

  }

  void EpetraCrsMatrix::replaceLocalValues(LocalOrdinal localRow, const ArrayView< const LocalOrdinal > &indices, const ArrayView< const Scalar > &values) {
    XPETRA_MONITOR("EpetraCrsMatrix::replaceLocalValues");

    {
      const std::string tfecfFuncName("replaceGlobalValues");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! isFillActive(), std::runtime_error,
                                            ": Fill must be active in order to call this method.  If you have already "
                                            "called fillComplete(), you need to call resumeFill() before you can "
                                            "replace values.");

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(values.size() != indices.size(),
                                            std::runtime_error, ": values.size() must equal indices.size().");
    }

    XPETRA_ERR_CHECK(mtx_->ReplaceMyValues(localRow, indices.size(), values.getRawPtr(), indices.getRawPtr()));

  }

  void EpetraCrsMatrix::allocateAllValues(size_t numNonZeros,ArrayRCP<size_t> & rowptr, ArrayRCP<int> & colind, ArrayRCP<double> & values)
  {
    // NOTE: Tpetra's insistence on using size_t's necessitate not directly accessing the rowptr here... Thus is is the ONLY thing we copy in setAllValues.
     XPETRA_MONITOR("EpetraCrsMatrix::allocateAllValues");

    // Grab pointers
    Epetra_IntSerialDenseVector& myColind = mtx_->ExpertExtractIndices();
    double *& myValues = mtx_->ExpertExtractValues();

    // Resize
    myColind.Resize(numNonZeros);
    delete [] myValues; myValues = new double[numNonZeros];

    // Wrap in array RCPs w/o the memory control.
    rowptr.resize(numNonZeros);
    colind = Teuchos::arcp(myColind.Values(),numNonZeros,false);
    values = Teuchos::arcp(myValues,numNonZeros,false);
  }

  void EpetraCrsMatrix::setAllValues(const ArrayRCP<size_t> & rowptr, const ArrayRCP<int> & colind, const ArrayRCP<double> & values)
  {
    XPETRA_MONITOR("EpetraCrsMatrix::setAllValues");

    TEUCHOS_TEST_FOR_EXCEPTION((size_t)rowptr.size()!=getNodeNumRows()+1, Xpetra::Exceptions::RuntimeError, "An exception is thrown to let you know that the size of your rowptr array is incorrect.");
    TEUCHOS_TEST_FOR_EXCEPTION(colind.getRawPtr()!=mtx_->ExpertExtractIndices().Values(), Xpetra::Exceptions::RuntimeError, "An exception is thrown to let you know that you mismatched your pointers.");
    TEUCHOS_TEST_FOR_EXCEPTION(values.getRawPtr()!=mtx_->ExpertExtractValues(), Xpetra::Exceptions::RuntimeError, "An exception is thrown to let you know that you mismatched your pointers.");
    TEUCHOS_TEST_FOR_EXCEPTION(values.size()!=colind.size(), Xpetra::Exceptions::RuntimeError, "An exception is thrown to let you know that you mismatched your pointers.");

    // NOTE: Tpetra's insistence on using size_t's necessitate not directly accessing the rowptr here... Thus is is the ONLY thing we copy in setAllValues.
    Epetra_IntSerialDenseVector& myRowptr = mtx_->ExpertExtractIndexOffset();
    size_t N = getNodeNumRows();

    myRowptr.Resize(N+1);
    for(size_t i=0; i<N; i++)
      myRowptr[i] = Teuchos::as<int>(rowptr[i]);
  }


  void EpetraCrsMatrix::resumeFill(const RCP< ParameterList > &params) {
    XPETRA_MONITOR("EpetraCrsMatrix::resumeFill");

    // According to Tpetra documentation, resumeFill() may be called repeatedly.
    isFillResumed_ = true;
  }

  bool EpetraCrsMatrix::isFillComplete() const { XPETRA_MONITOR("EpetraCrsMatrix::isFillComplete"); if (isFillResumed_) return false; else return mtx_->Filled(); }

  bool EpetraCrsMatrix::isFillActive() const { XPETRA_MONITOR("EpetraCrsMatrix::isFillActive"); return !isFillComplete(); }

  bool EpetraCrsMatrix::supportsRowViews() const { XPETRA_MONITOR("EpetraCrsMatrix::supportsRowViews"); return true; }

  //TODO: throw same exception as Tpetra
  void EpetraCrsMatrix::getLocalRowCopy(int LocalRow, const ArrayView<int> &Indices, const ArrayView<double> &Values, size_t &NumEntries) const {
    XPETRA_MONITOR("EpetraCrsMatrix::getLocalRowCopy");

    int numEntries=-1;
    XPETRA_ERR_CHECK(mtx_->ExtractMyRowCopy(LocalRow, Indices.size(), numEntries, Values.getRawPtr(), Indices.getRawPtr()));
    NumEntries = numEntries;
  }

  void EpetraCrsMatrix::getGlobalRowView(int GlobalRow, ArrayView<const int> &indices, ArrayView<const double> &values) const {
    XPETRA_MONITOR("EpetraCrsMatrix::getGlobalRowView");

    int      numEntries;
    double * eValues;
    int    * eIndices;

    XPETRA_ERR_CHECK(mtx_->ExtractGlobalRowView(GlobalRow, numEntries, eValues, eIndices));
    if (numEntries == 0) { eValues = NULL; eIndices = NULL; } // Cf. TEUCHOS_TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    indices = ArrayView<const int>(eIndices, numEntries);
    values  = ArrayView<const double>(eValues, numEntries);
  }

  void EpetraCrsMatrix::getLocalRowView(int LocalRow, ArrayView<const int> &indices, ArrayView<const double> &values) const {
    XPETRA_MONITOR("EpetraCrsMatrix::getLocalRowView");

    int      numEntries;
    double * eValues;
    int    * eIndices;

    XPETRA_ERR_CHECK(mtx_->ExtractMyRowView(LocalRow, numEntries, eValues, eIndices));
    if (numEntries == 0) { eValues = NULL; eIndices = NULL; } // Cf. TEUCHOS_TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    indices = ArrayView<const int>(eIndices, numEntries);
    values  = ArrayView<const double>(eValues, numEntries);
  }

  void EpetraCrsMatrix::apply(const MultiVector<double,int,int> & X, MultiVector<double,int,int> &Y, Teuchos::ETransp mode, double alpha, double beta) const {
    XPETRA_MONITOR("EpetraCrsMatrix::apply");

    //TEUCHOS_TEST_FOR_EXCEPTION((alpha != 1) || (beta != 0), Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraCrsMatrix.multiply() only accept alpha==1 and beta==0");

    XPETRA_DYNAMIC_CAST(const EpetraMultiVector, X, eX, "Xpetra::EpetraCrsMatrix->apply() only accept Xpetra::EpetraMultiVector as input arguments.");
    XPETRA_DYNAMIC_CAST(      EpetraMultiVector, Y, eY, "Xpetra::EpetraCrsMatrix->apply() only accept Xpetra::EpetraMultiVector as input arguments.");

    TEUCHOS_TEST_FOR_EXCEPTION((mode != Teuchos::NO_TRANS) && (mode != Teuchos::TRANS), Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraCrsMatrix->apply() only accept mode == NO_TRANS or mode == TRANS");
    bool eTrans = toEpetra(mode);

    // /!\ UseTranspose value
    TEUCHOS_TEST_FOR_EXCEPTION(mtx_->UseTranspose(), Xpetra::Exceptions::NotImplemented, "An exception is throw to let you know that Xpetra::EpetraCrsMatrix->apply() do not take into account the UseTranspose() parameter of Epetra_CrsMatrix.");

    RCP<Epetra_MultiVector> epY = eY.getEpetra_MultiVector();

    // helper vector: tmp = A*x
    RCP<Epetra_MultiVector> tmp = Teuchos::rcp(new Epetra_MultiVector(*epY));
    tmp->PutScalar(0.0);
    XPETRA_ERR_CHECK(mtx_->Multiply(eTrans, *eX.getEpetra_MultiVector(), *tmp));

    // calculate alpha * A * x + beta * y
    XPETRA_ERR_CHECK(eY.getEpetra_MultiVector()->Update(alpha,*tmp,beta));
  }

  std::string EpetraCrsMatrix::description() const {
    XPETRA_MONITOR("EpetraCrsMatrix::description");

    // This implementation come from Tpetra_CrsMatrix_def.hpp (without modification)
    std::ostringstream oss;
    //TODO: oss << DistObject<char, LocalOrdinal,GlobalOrdinal>::description();
    if (isFillComplete()) {
      oss << "{status = fill complete"
          << ", global rows = " << getGlobalNumRows()
          << ", global cols = " << getGlobalNumCols()
          << ", global num entries = " << getGlobalNumEntries()
          << "}";
    }
    else {
      oss << "{status = fill not complete"
          << ", global rows = " << getGlobalNumRows()
          << "}";
    }
    return oss.str();

  }

  void EpetraCrsMatrix::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    XPETRA_MONITOR("EpetraCrsMatrix::describe");

    // This implementation come from Tpetra_CrsMatrix_def.hpp (without modification)
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;
    RCP<const Comm<int> > comm = this->getComm();
    const int myImageID = comm->getRank(),
      numImages = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t>(width,11) + 2;
    Teuchos::OSTab tab(out);
    //    none: print nothing
    //     low: print O(1) info from node 0
    //  medium: print O(P) info, num entries per node
    //    high: print O(N) info, num entries per row
    // extreme: print O(NNZ) info: print indices and values
    //
    // for medium and higher, print constituent objects at specified verbLevel
    if (vl != VERB_NONE) {
      if (myImageID == 0) out << this->description() << std::endl;
      // O(1) globals, minus what was already printed by description()
      if (isFillComplete() && myImageID == 0) {
        out << "Global number of diagonals = " << getGlobalNumDiags() << std::endl;
        out << "Global max number of entries = " << getGlobalMaxNumRowEntries() << std::endl;
      }
      // constituent objects
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        if (myImageID == 0) out << "\nRow map: " << std::endl;
        getRowMap()->describe(out,vl);
        //
        if (getColMap() != null) {
          if (getColMap() == getRowMap()) {
            if (myImageID == 0) out << "\nColumn map is row map.";
          }
          else {
            if (myImageID == 0) out << "\nColumn map: " << std::endl;
            getColMap()->describe(out,vl);
          }
        }
        if (getDomainMap() != null) {
          if (getDomainMap() == getRowMap()) {
            if (myImageID == 0) out << "\nDomain map is row map.";
          }
          else if (getDomainMap() == getColMap()) {
            if (myImageID == 0) out << "\nDomain map is row map.";
          }
          else {
            if (myImageID == 0) out << "\nDomain map: " << std::endl;
            getDomainMap()->describe(out,vl);
          }
        }
        if (getRangeMap() != null) {
          if (getRangeMap() == getDomainMap()) {
            if (myImageID == 0) out << "\nRange map is domain map." << std::endl;
          }
          else if (getRangeMap() == getRowMap()) {
            if (myImageID == 0) out << "\nRange map is row map." << std::endl;
          }
          else {
            if (myImageID == 0) out << "\nRange map: " << std::endl;
            getRangeMap()->describe(out,vl);
          }
        }
        if (myImageID == 0) out << std::endl;
      }
      // O(P) data
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << "Node ID = " << imageCtr << std::endl;
            // TODO: need a graph
            //               if (staticGraph_->indicesAreAllocated() == false) {
            //                 out << "Node not allocated" << std::endl;
            //               }
            //               else {
            //                 out << "Node number of allocated entries = " << staticGraph_->getNodeAllocationSize() << std::endl;
            //               }

            // TMP:
            //            const Epetra_CrsGraph & staticGraph_ = mtx_->Graph();
            // End of TMP

            out << "Node number of entries = " << getNodeNumEntries() << std::endl;
            if (isFillComplete()) {
              out << "Node number of diagonals = " << getNodeNumDiags() << std::endl;
            }
            out << "Node max number of entries = " << getNodeMaxNumRowEntries() << std::endl;
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
      // O(N) and O(NNZ) data
      if (vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << std::setw(width) << "Node ID"
                << std::setw(width) << "Global Row"
                << std::setw(width) << "Num Entries";
            if (vl == VERB_EXTREME) {
              out << std::setw(width) << "(Index,Value)";
            }
            out << std::endl;
            for (size_t r=0; r < getNodeNumRows(); ++r) {
              const size_t nE = getNumEntriesInLocalRow(r);
              GlobalOrdinal gid = getRowMap()->getGlobalElement(r);
              out << std::setw(width) << myImageID
                  << std::setw(width) << gid
                  << std::setw(width) << nE;
              if (vl == VERB_EXTREME) {
                if (isGloballyIndexed()) {
                  ArrayView<const GlobalOrdinal> rowinds;
                  ArrayView<const Scalar> rowvals;
                  getGlobalRowView(gid,rowinds,rowvals);
                  for (size_t j=0; j < nE; ++j) {
                    out << " (" << rowinds[j]
                        << ", " << rowvals[j]
                        << ") ";
                  }
                }
                else if (isLocallyIndexed()) {
                  ArrayView<const LocalOrdinal> rowinds;
                  ArrayView<const Scalar> rowvals;
                  getLocalRowView(r,rowinds,rowvals);
                  for (size_t j=0; j < nE; ++j) {
                    out << " (" << getColMap()->getGlobalElement(rowinds[j])
                        << ", " << rowvals[j]
                        << ") ";
                  }
                }
              }
              out << std::endl;
            }
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
    }

  }

  // TODO: use toEpetra()
  void EpetraCrsMatrix::doImport(const DistObject<char, int, int> &source,
                                 const Import<int, int> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsMatrix::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrix, source, tSource, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraCrsMatrix as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tSource.getEpetra_CrsMatrix();
    int err = mtx_->Import(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraCrsMatrix::doExport(const DistObject<char, int, int> &dest,
                                 const Import<int, int>& importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsMatrix::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrix, dest, tDest, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraCrsMatrix as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tDest.getEpetra_CrsMatrix();
    int err = mtx_->Export(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraCrsMatrix::doImport(const DistObject<char, int, int> &source,
                                 const Export<int, int>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsMatrix::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrix, source, tSource, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraCrsMatrix as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tSource.getEpetra_CrsMatrix();
    int err = mtx_->Import(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");

  }

  void EpetraCrsMatrix::doExport(const DistObject<char, int, int> &dest,
                                 const Export<int, int>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraCrsMatrix::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrix, dest, tDest, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraCrsMatrix as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tDest.getEpetra_CrsMatrix();
    int err = mtx_->Export(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

    void EpetraCrsMatrix::fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap,
                                       const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap,
                                       const RCP< ParameterList > &params) {
      XPETRA_MONITOR("EpetraCrsMatrix::fillComplete");

      // For Epetra matrices, resumeFill() is a fictive operation. There is no need for a fillComplete after some resumeFill() operations.
      if (isFillResumed_ == true) { isFillResumed_ = false; return; }

      bool doOptimizeStorage = true;
      if (params != null && params->get("Optimize Storage",true) == false) doOptimizeStorage = false;
      mtx_->FillComplete(toEpetra(domainMap), toEpetra(rangeMap), doOptimizeStorage);
    }

    void EpetraCrsMatrix::fillComplete(const RCP< ParameterList > &params) {
      XPETRA_MONITOR("EpetraCrsMatrix::fillComplete");

      // For Epetra matrices, resumeFill() is a fictive operation. There is no need for a fillComplete after some resumeFill() operations.
      if (isFillResumed_ == true) { isFillResumed_ = false; return; }

      bool doOptimizeStorage = true;
      if (params != null && params->get("Optimize Storage",true) == false) doOptimizeStorage = false;
      mtx_->FillComplete(doOptimizeStorage);
    }


  void EpetraCrsMatrix::replaceDomainMapAndImporter(const Teuchos::RCP< const  Map< LocalOrdinal, GlobalOrdinal, Node > >& newDomainMap, Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> >  & newImporter) {
      XPETRA_MONITOR("EpetraCrsMatrix::replaceDomainMapAndImporter");
      XPETRA_DYNAMIC_CAST(const EpetraImport, *newImporter, eImporter, "Xpetra::EpetraCrsMatrix::replaceDomainMapAndImporter only accepts Xpetra::EpetraImport.");

      const RCP<const Epetra_Import> & myImport = eImporter.getEpetra_Import();
      int rv=0;
      if(myImport==Teuchos::null)
	rv=mtx_->ReplaceDomainMapAndImporter( toEpetra(newDomainMap),0);
      else
	rv=mtx_->ReplaceDomainMapAndImporter( toEpetra(newDomainMap),&*myImport);
      TEUCHOS_TEST_FOR_EXCEPTION(rv != 0, std::runtime_error, "Xpetra::EpetraCrsMatrix::replaceDomainMapAndImporter FAILED!");
  }


  void EpetraCrsMatrix::expertStaticFillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap,
						 const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap,
						 const RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > &importer,
						 const RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > &exporter,
						 const RCP<ParameterList> & params)
  {
    XPETRA_MONITOR("EpetraCrsMatrix::expertStaticFillComplete");
    XPETRA_DYNAMIC_CAST(const EpetraImport, *importer, eImporter, "Xpetra::EpetraCrsMatrix::expertStaticFillComplete only accepts Xpetra::EpetraImport.");
    //    XPETRA_DYNAMIC_CAST(const EpetraExport, *exporter, eExporter, "Xpetra::EpetraCrsMatrix::expertStaticFillComplete only accepts Xpetra::EpetraImport.");

    const Epetra_Import * myimport = eImporter.getEpetra_Import().getRawPtr();
    TEUCHOS_TEST_FOR_EXCEPTION(exporter != Teuchos::null, std::runtime_error, "Xpetra::EpetraCrsMatrix::expertStaticFillComplete can't handle a non-null exporter.");
    //    const Epetra_Export * myexport = eExporter.getEpetra_Export().getRawPtr();

    int rv=mtx_->ExpertStaticFillComplete(toEpetra(domainMap), toEpetra(rangeMap), myimport);
    TEUCHOS_TEST_FOR_EXCEPTION(rv != 0, std::runtime_error, "Xpetra::EpetraCrsMatrix::expertStaticFillComplete FAILED!");
  }




}
