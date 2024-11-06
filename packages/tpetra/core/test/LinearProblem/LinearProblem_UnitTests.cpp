// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_LinearProblem.hpp"

#include "Tpetra_CrsMatrix.hpp"


namespace { // (anonymous)

  using Tpetra::TestingUtilities::getDefaultComm;
  using Tpetra::createContigMapWithNode;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::tuple;
  //using Teuchos::NO_TRANS;
  //using Teuchos::TRANS;
  //using Teuchos::CONJ_TRANS;
  using std::endl;
  using GST = Tpetra::global_size_t;
  

  /// \brief Print out pretty version of RowMatrix.
  template <typename Scalar, typename LO, typename GO, typename Node>
  void Display_CrsMatrix (std::string label, RCP<Tpetra::RowMatrix<Scalar, LO, GO, Node> > A, RCP< const Comm< int > > comm, Teuchos::FancyOStream& myOut)
  {
    using local_ordinal_type = typename Tpetra::Vector<Scalar, LO, GO, Node>::local_ordinal_type;
    
    using crs_matrix_type = typename Tpetra::CrsMatrix<Scalar, LO, GO, Node>;
  
    using crs_local_inds_host_view_type = typename crs_matrix_type::local_inds_host_view_type;
    using crs_values_host_view_type = typename crs_matrix_type::values_host_view_type;

    const local_ordinal_type INVALID = Teuchos::OrdinalTraits<local_ordinal_type>::invalid();
    
    // Get the number of rows and columns
    GST numRows = A->getGlobalNumRows();
    GST numCols = A->getGlobalNumCols();
  
    // Loop over all global rows
    for (GST globalRow = 0; globalRow < numRows; ++globalRow) {
        // Check if this row belongs to the current process
        if (A->getRowMap()->getLocalElement(globalRow) != INVALID) {
            myOut << "Row " << std::setw(2) << globalRow << " [ ";
  
            // Extract the row view
            crs_local_inds_host_view_type localIndices;
            crs_values_host_view_type     values;
            size_t localRow = A->getRowMap()->getLocalElement(globalRow);
            A->getLocalRowView(localRow, localIndices, values);
  
            // Initialize a vector to track printed entries
            std::vector<bool> printed(numCols, false);
  
            // Print the entries in the row
            size_t numEntries = A->getNumEntriesInLocalRow(localRow);
            for (size_t k = 0; k < numEntries; ++k) {
                // Convert local index to global index
                GST globalIndex = A->getColMap()->getGlobalElement(localIndices(k));
                printed[globalIndex] = true; // Mark the index as having a non-zero entry
            }
  
            // Print the values for each column
            for (GST j = 0; j < numCols; ++j) {
                if (printed[j]) {
                    // Find the corresponding value for the global index
                    for (size_t k = 0; k < numEntries; ++k) {
                        // Convert local index to global index
                        GST globalIndex = A->getColMap()->getGlobalElement(localIndices(k));
                        if (globalIndex == j) {
                            myOut << std::setw(8) << values(k) << " ";
                            break;
                        }
                    }
                } else {
                    myOut << std::setw(8) << 0 << " ";
                }
            }
            myOut << "]" << endl;
        }
        // Synchronize processes before printing
        comm->barrier();
    }
  }

  template <typename Scalar, typename LO, typename GO, typename Node>
  void Display_MultiVector (std::string label, Teuchos::RCP<Tpetra::MultiVector<Scalar, LO, GO, Node>> multivector,  Teuchos::RCP< const Teuchos::Comm< int > > comm, Teuchos::FancyOStream& myOut)
  {
    using local_ordinal_type = typename Tpetra::Vector<>::local_ordinal_type;
    const local_ordinal_type INVALID = Teuchos::OrdinalTraits<local_ordinal_type>::invalid();
  
    auto map = multivector->getMap();
    const size_t myImageID = comm->getRank();
  
    if (myImageID==0) {
      myOut << label << endl;
      myOut << std::setw(8) << "Rank" << std::setw(12) << "GID" << std::setw(20) << "Value(s)" << endl;
    }
    GST numRows = multivector->getGlobalLength();
    for (GST globalRow = 0; globalRow < numRows; ++globalRow) {
      // Check if this row belongs to the current process
      if (map->getLocalElement(globalRow) != INVALID) {
        size_t localElement = map->getLocalElement(globalRow);
        myOut << std::setw(8) << myImageID << std::setw(12) << globalRow << "          ";
        for (size_t j = 0; j < multivector->getNumVectors(); ++j) {
          myOut << std::setw(10) << multivector->getData(j)[localElement];
        }
        myOut << endl;
      }
      comm->barrier();
    }
  }

  template <typename Scalar, typename LO, typename GO, typename Node>
  void Display_Vector (std::string label, Teuchos::RCP<Tpetra::Vector<Scalar, LO, GO, Node>> vector, Teuchos::RCP< const Teuchos::Comm< int > > comm, Teuchos::FancyOStream& myOut)
  {
    auto multivector = Teuchos::rcp_dynamic_cast<Tpetra::MultiVector<Scalar, LO, GO, Node>> (vector);
    Display_MultiVector(label, multivector, comm, myOut);
  }


  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( LinearProblem, basic, LO, GO, Scalar, Node )
  {
    using map_type = Tpetra::Map<LO, GO, Node>;
    using ST = Teuchos::ScalarTraits<Scalar>;
    using mag_type = typename ST::magnitudeType;
    
    using MAT = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;
    using VT = Tpetra::Vector<Scalar,LO,GO,Node>;
    using MV = Tpetra::MultiVector<Scalar,LO,GO,Node>;
    using LPT = Tpetra::LinearProblem<Scalar,LO,GO,Node>;
    //using local_ordinal_type = typename Tpetra::Vector<Scalar, LO, GO, Node>::local_ordinal_type;
    using global_ordinal_type = typename Tpetra::Vector<Scalar, LO, GO, Node>::global_ordinal_type;

    
    const global_ordinal_type INVALID = Teuchos::OrdinalTraits<global_ordinal_type>::invalid();
    constexpr bool debug = true;

    RCP<Teuchos::FancyOStream> outPtr = debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;
    Teuchos::OSTab tab0 (myOut);

    myOut << "Test: LinearProblem, Constructors" << endl;

    RCP<const Comm<int> > comm = getDefaultComm();
    //const size_t numImages = comm->getSize();
    const size_t myImageID = comm->getRank();
    // create a Map
    const size_t numLocal = 10;
    const size_t numVecs  = 1;
    RCP<const map_type> map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    GO base = numLocal*myImageID;
    GST globalNumElements = map->getGlobalNumElements();
    RCP<Tpetra::RowMatrix<Scalar,LO,GO,Node> > A;
    {
      RCP<MAT> A_crs = rcp(new MAT(map,3));
      for (size_t i=0; i<numLocal; ++i) {
        //A_crs->insertGlobalValues(base+i,tuple<GO>(base+i),tuple<Scalar>(ST::one()));
        A_crs->insertGlobalValues(base + i, tuple<GO>(base + i), tuple<Scalar>(2.0)); // Diagonal entry

        GST globalIndex = base + i;
        // Insert the first subdiagonal entry if not in the first row
        if (globalIndex > 0) {
            A_crs->insertGlobalValues(base + i, tuple<GO>(base + i - 1), tuple<Scalar>(1.0)); // Subdiagonal entry
        }

        // Insert the first superdiagonal entry if not in the last row
        if (globalIndex < globalNumElements - 1) {
            A_crs->insertGlobalValues(base + i, tuple<GO>(base + i + 1), tuple<Scalar>(1.0)); // Superdiagonal entry
        }
      }
      A_crs->fillComplete();
      A = A_crs;
    }

    // create solution, rhs and scaling vector
    RCP<MV> X = rcp (new MV (map, numVecs));
    RCP<MV> B = rcp (new MV (map, numVecs));
    RCP<VT> S = rcp (new VT (map));

    // Assign values to the MultiVector based on the global index
    for (size_t j = 0; j < numVecs; ++j) { // Loop over each vector (column)
        for (GST i = 0; i < globalNumElements; ++i) {
            // Assign a value (for example, the global index plus the vector index)
            X->replaceGlobalValue(i, j, Teuchos::as<Scalar>(i + j + 1));
            B->replaceGlobalValue(i, j, Teuchos::as<Scalar>(i + j + 1));
        }
    }
    for (GST i = 0; i < globalNumElements; ++i) {
        S->replaceGlobalValue(i, Teuchos::as<Scalar>(i + 1));
    }

    RCP<LPT> linearProblem = rcp(new LPT());

    linearProblem->setOperator(A);
    linearProblem->setLHS(X);
    linearProblem->setRHS(B);

    linearProblem->checkInput();

    //if (myImageID==0) myOut << "Original LinearProblem" << endl;
    //Display_CrsMatrix("A", A, comm, myOut);
    //Display_MultiVector("Solution Vector", X, comm, myOut);
    //Display_MultiVector("RHS Vector", B, comm, myOut);
    //Display_Vector("Scaling Vector", S, comm, myOut);

    // Original LinearProblem
    GST N = globalNumElements;
    double normF = std::sqrt(6*N - 2);
    TEST_FLOATING_EQUALITY(linearProblem->getMatrix()->getFrobeniusNorm(),
                           Teuchos::as<Scalar>(normF), Teuchos::as<Scalar>(1.0e-14));
                           //Teuchos::as<Scalar>(7.615773105863909), Teuchos::as<Scalar>(1.0e-14));

    Array<mag_type> norms(numVecs);
    linearProblem->getLHS()->norm1(norms());
    size_t vector_sum = N*(N+1)/2;
    TEST_FLOATING_EQUALITY(norms[0], Teuchos::as<Scalar>(vector_sum), Teuchos::as<Scalar>(1.0e-14));
    linearProblem->getRHS()->norm1(norms());
    TEST_FLOATING_EQUALITY(norms[0], Teuchos::as<Scalar>(vector_sum), Teuchos::as<Scalar>(1.0e-14));
    
    // Left Scaling
    linearProblem->leftScale(S);

    size_t vector_sum_squared = N*(N+1)*(2*N+1)/6;
    normF = std::sqrt(6*vector_sum_squared - N*N - 1);
    TEST_FLOATING_EQUALITY(linearProblem->getMatrix()->getFrobeniusNorm(),
                           Teuchos::as<Scalar>(normF), Teuchos::as<Scalar>(1.0e-14));
    linearProblem->getLHS()->norm1(norms());
    TEST_FLOATING_EQUALITY(norms[0], Teuchos::as<Scalar>(vector_sum), Teuchos::as<Scalar>(1.0e-14));
    linearProblem->getRHS()->norm1(norms());
    TEST_FLOATING_EQUALITY(norms[0], Teuchos::as<Scalar>(vector_sum_squared), Teuchos::as<Scalar>(1.0e-14));
    
    //if (myImageID==0) myOut << "After Left Scaling" << endl;
    //Display_CrsMatrix("A", A, comm, myOut);
    //Display_MultiVector("Solution Vector", X, comm, myOut);
    //Display_MultiVector("RHS Vector", B, comm, myOut);
    //Display_Vector("Scaling Vector", S, comm, myOut);

    // Right Scaling
    linearProblem->rightScale(S);

    N = N-1;
    size_t off_diags = 2.0*((N * (N + 1) * (2 * N + 1) * (3 * N * N + 3 * N - 1)) / 30.0
                            + (N * N * (N + 1) * (N + 1)) / 2.0
                            + (N * (N + 1) * (2 * N + 1)) / 6.0);
    N = N+1;
    size_t diag = (2.0 * N * (N + 1) * (2 * N + 1) * (3 * N * N + 3 * N - 1)) / 15.0;
    normF = std::sqrt(diag + off_diags);
    TEST_FLOATING_EQUALITY(linearProblem->getMatrix()->getFrobeniusNorm(),
                           Teuchos::as<Scalar>(normF), Teuchos::as<Scalar>(1.0e-14));
    linearProblem->getLHS()->norm1(norms());
    TEST_FLOATING_EQUALITY(norms[0], Teuchos::as<Scalar>(N), Teuchos::as<Scalar>(1.0e-14));
    linearProblem->getRHS()->norm1(norms());
    TEST_FLOATING_EQUALITY(norms[0], Teuchos::as<Scalar>(vector_sum_squared), Teuchos::as<Scalar>(1.0e-14));
    
    //if (myImageID==0) myOut << "After Right Scaling" << endl;
    //Display_CrsMatrix("A", A, comm, myOut);
    //Display_MultiVector("Solution Vector", X, comm, myOut);
    //Display_MultiVector("RHS Vector", B, comm, myOut);
    //Display_Vector("Scaling Vector", S, comm, myOut);

    // Constructor with matrix
    {
      RCP<LPT> linearProblem_Matrix = rcp(new LPT(A,X,B));
      linearProblem_Matrix->checkInput();
    }
      
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( LinearProblem, basic,  LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}
