// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Belos_TpetraTestFramework_hpp
#define __Belos_TpetraTestFramework_hpp

/// \file BelosTpetraTestFramework.hpp
/// \brief A common test framework for Tpetra instantiations of Belos solvers.
/// \author Mark Hoemmen

#include "BelosTpetraAdapter.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Trilinos_Util_iohb.h"

namespace Belos {
  namespace Tpetra {

    /// \function convertData
    template<class scalar>
    void convertData( const double* dPtr, const int nnz, scalar* retVals );

    template<>
    void convertData( const double* dPtr, const int nnz, double* retVals )
    {
      for ( int i=0; i<nnz; ++i )
        retVals[i] = dPtr[i];
    }

    template<>
    void convertData( const double* dPtr, const int nnz, float* retVals )
    {
      for ( int i=0; i<nnz; ++i )
        retVals[i] = (float) dPtr[i];
    }


    template<>
    void convertData( const double* dPtr, const int nnz, std::complex<double>* retVals )
    {
      for ( int i=0; i<nnz; ++i )
        retVals[i] = std::complex<double>( dPtr[2*i], dPtr[2*i+1] );
    }

    /// \class HarwellBoeingReader
    /// \brief Read a Harwell-Boeing format file into a Tpetra::CrsMatrix.
    /// \author Mark Hoemmen
    ///
    /// \warning This reader is incomplete.  It currently only handles
    ///   double-precision floating-point data and some of the possible
    ///   representations to create float, real, and complex-valued matrices.
    template<class SparseMatrixType>
    class HarwellBoeingReader {

      typedef SparseMatrixType sparse_matrix_type;
      typedef typename SparseMatrixType::scalar_type scalar_type;
      typedef typename SparseMatrixType::local_ordinal_type local_ordinal_type;
      typedef typename SparseMatrixType::global_ordinal_type global_ordinal_type;
      typedef typename SparseMatrixType::node_type node_type;

    public:
      /// \typedef multivector_type
      /// \brief Specialization of \c Tpetra::MultiVector.
      ///
      /// This class deduces the \c Tpeta::MultiVector specialization
      /// from the \c Tpetra::CrsMatrix specialization.  This ensures
      /// compatibility of the sparse matrix type with the dense
      /// multivector type.
      typedef ::Tpetra::MultiVector<scalar_type,
                                    local_ordinal_type,
                                    global_ordinal_type,
                                    node_type> multivector_type;
      /// \typedef map_type
      /// \brief Specialization of \c Tpetra::Map.
      ///
      /// This class deduces the \c Tpeta::Map specialization from the
      /// \c Tpetra::CrsMatrix specialization.  This ensures
      /// compatibility of the sparse matrix type with the map type.
      typedef ::Tpetra::Map<local_ordinal_type,
                            global_ordinal_type,
                            node_type> map_type;

      /// \brief Constructor.
      ///
      /// \param comm [in] Communicator over which to distribute the
      ///   sparse matrices.
      HarwellBoeingReader (const Teuchos::RCP<const Teuchos::Comm<int> >& comm) :
        comm_ (comm)
      {}

      //! Read the sparse matrix from the file with the given name.
      Teuchos::RCP<sparse_matrix_type>
      readFromFile (const std::string& filename)
      {
        //
        // Get the data from the HB file and build the Map,Matrix
        //
        int MyPID = rank(*comm_);
        int info = 0;
        int dim,dim2,nnz,rnnzmax;
        double *dvals;
        int *colptr,*rowind;
        nnz = -1;
        if (MyPID == 0) {
          info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
          // find maximum NNZ over all rows
          std::vector<int> rnnz(dim,0);
          for (int *ri=rowind; ri<rowind+nnz; ++ri) {
            ++rnnz[*ri-1];
          }  
          rnnzmax = *std::max_element(rnnz.begin(),rnnz.end());
        }
        else {
          // address uninitialized data warnings
          dvals = NULL;
          colptr = NULL;
          rowind = NULL;
        }
        Teuchos::broadcast(*comm_,0,&info);
        Teuchos::broadcast(*comm_,0,&nnz);
        Teuchos::broadcast(*comm_,0,&dim);
        Teuchos::broadcast(*comm_,0,&rnnzmax);

        TEUCHOS_TEST_FOR_EXCEPTION((info == 0 || nnz < 0), std::runtime_error,
          "HarwellBoeingReader::readFromFile: Failed to read Harwell-Boeing "
          "sparse matrix file: readHB_newmat_double() returned INFO = "
          << info << ".");

        // create map
        Teuchos::RCP<const ::Tpetra::Map<> > map = Teuchos::rcp (new ::Tpetra::Map<> (dim, 0, comm_));
        Teuchos::RCP<sparse_matrix_type> A = Teuchos::rcp(new sparse_matrix_type(map, rnnzmax));
        if (MyPID == 0) {
          // HB format is compressed column. CrsMatrix is compressed row.
          scalar_type* svals = new scalar_type[ nnz ];
          convertData<scalar_type>( dvals, nnz, svals );
          const scalar_type *sptr = svals;
          const int *rptr = rowind;
          for (int c=0; c<dim; ++c) {
            for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
              A->insertGlobalValues (static_cast<global_ordinal_type> (*rptr++ - 1), 
                                     Teuchos::tuple<global_ordinal_type> (c), 
                                     Teuchos::tuple (sptr[0]));
              sptr++;
            }
          }

          // Clean up.
          free( svals );
          free( dvals );
          free( colptr );
          free( rowind );
        }

        A->fillComplete();
        return A;
      }

    private:
      //! Communicator over which to distribute the sparse matrices.
      Teuchos::RCP<const Teuchos::Comm<int> > comm_;

      /// \brief Convert the Harwell-Boeing data to a Tpetra::CrsMatrix.
      ///
      /// \note We allow the values and indices read from the file to
      ///   have different types than the types of the sparse matrix's
      ///   values resp. global indices.  This is why this method is
      ///   templated on ValueType and IndexType, which are the types
      ///   of values resp. indices found in the file.
      template<class ValueType>
      void
      insertValues (Teuchos::RCP<sparse_matrix_type>& A,
                    const double* dptr, const int *rind, const int *cptr );

    };

    /// \class ProblemMaker
    /// \brief Make a test linear problem using Tpetra objects.
    /// \author Mark Hoemmen
    ///
    /// \tparam SparseMatrixType A specialization of Tpetra::CrsMatrix.
    ///
    /// This class is useful for testing the Tpetra specialization of
    /// Belos solvers.  It is meant as a rudimentary version of Galeri
    /// (which currently only creates Epetra objects) and avoids
    /// duplicated code in the Tpetra tests.
    ///
    /// This class is templated on a specialization of \c
    /// Tpetra::CrsMatrix.  This saves users the trouble of needing to
    /// know all five template arguments of their Tpetra::CrsMatrix
    /// specialization in order to get the right sparse matrix type.
    /// This class uses the typedefs of \c Tpetra::CrsMatrix to deduce
    /// the \c Tpetra::MultiVector specialization to use.
    ///
    /// \note To Belos developers: using objects in the Tpetra
    ///   namespace in this class requires prefixing the Tpetra
    ///   namespace with ::, as in \c ::Tpetra::Map for example.
    template<class SparseMatrixType>
    class ProblemMaker : public Teuchos::ParameterListAcceptorDefaultBase {
    public:
      typedef SparseMatrixType sparse_matrix_type;
      typedef typename SparseMatrixType::scalar_type scalar_type;
      typedef typename SparseMatrixType::local_ordinal_type local_ordinal_type;
      typedef typename SparseMatrixType::global_ordinal_type global_ordinal_type;
      typedef typename SparseMatrixType::node_type node_type;
      /// \typedef multivector_type
      /// \brief Specialization of \c Tpetra::MultiVector.
      ///
      /// This class deduces the \c Tpeta::MultiVector specialization
      /// from the \c Tpetra::CrsMatrix specialization.  This ensures
      /// compatibility of the sparse matrix type with the dense
      /// multivector type.
      typedef ::Tpetra::MultiVector<scalar_type,
                                    local_ordinal_type,
                                    global_ordinal_type,
                                    node_type> multivector_type;

    private:
      typedef Teuchos::ScalarTraits<scalar_type> STS;
      typedef ::Tpetra::Map<local_ordinal_type,
                            global_ordinal_type,
                            node_type> map_type;

    public:
      /// \brief Constructor.
      ///
      /// \param comm [in] Communicator over which to distribute the
      ///   sparse matrices and dense vectors created by this class
      ///   instance.
      /// \param node [in] Kokkos Node instance with which to create
      ///   Tpetra objects.
      /// \param out [out] Output stream that prints only on Rank 0
      ///   of the given communicator.
      /// \param plist [in/out] On input: nonnull list of parameters.
      ///   On output, missing parameters will be filled in with their
      ///   default values.
      ProblemMaker (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                    const Teuchos::RCP<Teuchos::FancyOStream>& out,
                    const Teuchos::RCP<Teuchos::ParameterList> plist) :
        comm_ (comm), out_ (out), tolerant_ (false), debug_ (false)
      {
        setParameterList (plist);
      }

      /// \brief Constructor.
      ///
      /// \param comm [in] Communicator over which to distribute the
      ///   sparse matrices and dense vectors created by this class
      ///   instance.
      /// \param node [in] Kokkos Node instance with which to create
      ///   Tpetra objects.
      /// \param out [out] Output stream that prints only on Rank 0
      ///   of the given communicator.
      /// \param tolerant [in] Whether to parse files tolerantly.
      /// \param debug [in] Whether to print copious debugging output.
      ProblemMaker (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                    const Teuchos::RCP<Teuchos::FancyOStream>& out,
                    const bool tolerant,
                    const bool debug) :
        comm_ (comm),
        out_ (out),
        tolerant_ (tolerant),
        debug_ (debug)
      {
        using Teuchos::ParameterList;
        using Teuchos::parameterList;
        using Teuchos::RCP;

        RCP<ParameterList> plist = parameterList ();
        plist->set ("Parse files tolerantly", tolerant);
        plist->set ("Debug", debug);

        this->setMyParamList (plist);
      }

      //! @name Implementation of Teuchos::ParameterListAcceptor
      //@{

      Teuchos::RCP<const Teuchos::ParameterList>
      getValidParameters() const
      {
        using Teuchos::ParameterList;
        using Teuchos::parameterList;
        using Teuchos::RCP;

        RCP<ParameterList> validParams = parameterList ();
        const bool tolerant = false;
        const bool debug = false;
        validParams->set ("Parse files tolerantly", tolerant);
        validParams->set ("Debug", debug);

        return validParams;
      }

      void
      setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist)
      {
        using Teuchos::ParameterList;
        using Teuchos::RCP;

        RCP<const ParameterList> defaultParams = getValidParameters ();
        tolerant_ = plist->get<bool> ("Parse files tolerantly");
        debug_ = plist->get<bool> ("Debug");
        this->setMyParamList (plist);
      }

      //@}

      /// \brief Generate a sparse matrix using the given parameters.
      ///
      /// This function accepts the following parameters:
      /// - "Input matrix filename" (std::string): If provided and not
      ///   "", then read the sparse matrix in from the file with the
      ///   given name.  Currently, we assume that the file is in
      ///   Matrix Market format.
      /// - "Output matrix filename" (std::string): If provided and
      ///   not "", write the sparse matrix (whether generated or
      ///   read) out in Matrix Market format to the file with the
      ///   given name.
      /// - "Global number of rows" (global_ordinal_type): If provided
      ///   and if generating the matrix (rather than reading it from
      ///   a file), use this as the global number of rows (and
      ///   columns) in the matrix.  Otherwise ignored.  The type of
      ///   this parameter must be the same as \c
      ///   SparseMatrixType::global_ordinal_type.
      /// - "Problem type" (std::string): If generating the matrix,
      ///   the type of matrix to generate.  You should provide a
      ///   global number of rows if you provide this parameter.
      ///   "Symmetric" means generate a symmetric matrix.  "GMRES
      ///   unfriendly" means generate a matrix for which GMRES
      ///   converges slowly.  The default is "Nonsymmetric", which
      ///   generates a nonsymmetric matrix for which GMRES
      ///   convergence is not too slow.
      Teuchos::RCP<sparse_matrix_type>
      makeMatrix (const Teuchos::RCP<Teuchos::ParameterList>& params)
      {
        const std::string inputMatrixFilename =
          params->isType<std::string> ("Input matrix filename") ?
          params->get<std::string> ("Input matrix filename") : "";
        const std::string outputMatrixFilename =
          params->isType<std::string> ("Output matrix filename") ?
          params->get<std::string> ("Output matrix filename") : "";
        const global_ordinal_type globalNumRows =
          params->isType<global_ordinal_type> ("Global number of rows") ?
          params->get<global_ordinal_type> ("Global number of rows") :
          comm_->getSize() * 100; // A sensible default
        const bool generated = params->isType<bool> ("Generated") ?
          params->get<bool> ("Generated") : inputMatrixFilename == "";
        bool symmetric = false;
        bool gmresUnfriendly = false;
        if (params->isParameter ("Problem type")) {
          const std::string probType = params->get<std::string> ("Problem type");
          if (probType == "Symmetric") {
            symmetric = true;
          }
          else if (probType == "GMRES unfriendly") {
            gmresUnfriendly = true;
          }
        }

        Teuchos::RCP<sparse_matrix_type> A;
        if (inputMatrixFilename != "") {
          A = readSparseMatrixFromFile (inputMatrixFilename);
        } else if (generated) {
          if (gmresUnfriendly) {
            A = generateGmresUnfriendlyMatrix (globalNumRows);
          } else {
            A = generateTestMatrix (globalNumRows, symmetric);
          }
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(! (inputMatrixFilename == "" || generated),
            std::logic_error, "Should never get here!");
        }

        if (outputMatrixFilename != "") {
          typedef ::Tpetra::MatrixMarket::Writer<sparse_matrix_type> writer_type;
          writer_type::writeSparseFile (outputMatrixFilename, A, "", "", debug_);
        }
        return A;
      }

      /// \brief Read in a Tpetra sparse matrix from a Matrix Market file.
      ///
      /// \param matrixFilename [in] Name of the Matrix Market format
      ///   sparse matrix file to read (on Rank 0 of the communicator
      ///   only).
      ///
      /// \return Sparse matrix (represented as a Tpetra::CrsMatrix).
      ///
      /// You can access this functionality from \c makeMatrix().  We
      /// merely provide it as a commonly used special case, so you
      /// don't have to construct a parameter list just to read a
      /// matrix from a file.
      Teuchos::RCP<sparse_matrix_type>
      readSparseMatrixFromFile (const std::string& filename)
      {
        using Teuchos::RCP;
        using std::endl;
        typedef ::Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;

        RCP<Teuchos::FancyOStream> out = getDebugStream ();
        Teuchos::OSTab tab (out); // Indent for debug output

        // Read the sparse matrix A from the file.
        *out << "Reading sparse matrix A from file \"" << filename
             << "\"" << endl;
        const bool callFillComplete = true;
        RCP<sparse_matrix_type> A =
          reader_type::readSparseFile (filename, comm_, callFillComplete,
                                       tolerant_, debug_);
        return A;
      }

      /// \brief Given a sparse matrix A, make corresponding linear system vectors.
      ///
      /// \param A [in] The sparse matrix.  We require it as input
      ///   because we use its domain and range maps to create
      ///   right-hand side, exact solution, and initial guess
      ///   vectors.
      /// \param X_guess [out] The initial guess for the linear system AX=B.
      /// \param X_exact [out] The exact solution X of the linear system AX=B.
      /// \param B [out] The right-hand side of the linear system AX=B.
      /// \param params [in/out] List of parameters for generating the
      ///   vectors.  This method accepts the following parameters:
      ///   - "Input RHS filename" (std::string): If provided and not
      ///     "", then read the right-hand side (RHS) in from the file
      ///     with the given name.  Currently, we assume that the file
      ///     is in Matrix Market dense matrix format.
      ///   - "Output RHS filename" (std::string): If provided and not
      ///     "", write the right-hand side(s) (RHS) out in Matrix
      ///     Market format to the file with the given name.
      ///   - "Number of RHS" (size_t): If provided, and if no input
      ///     filename was specified, B, X_guess, and X_exact will all
      ///     have this number of columns.  Defaults to 1.
      ///   - "Random solution" (bool): If true, pick a random (via
      ///     MultiVector::randomize()) exact solution, and compute
      ///     the right-hand side(s) accordingly.
      ///   - "Problem type" (std::string): Optional: The type of
      ///     problem to generate.  Currently, the only value for
      ///     which this matters is "GMRES unfriendly", which means
      ///     generate a right-hand side corresponding to the linear
      ///     system for which GMRES converges slowly.
      void
      makeVectors (const Teuchos::RCP<const sparse_matrix_type>& A,
                   Teuchos::RCP<multivector_type>& X_guess,
                   Teuchos::RCP<multivector_type>& X_exact,
                   Teuchos::RCP<multivector_type>& B,
                   const Teuchos::RCP<Teuchos::ParameterList>& params)
      {
        using Teuchos::as;
        using Teuchos::RCP;
        using Teuchos::rcp;
        using std::endl;

        typedef global_ordinal_type GO;
        typedef multivector_type MV;

        RCP<Teuchos::FancyOStream> out = getDebugStream ();
        Teuchos::OSTab tab (out); // Indent for verbose output

        const std::string inRhsFilename =
          params->isType<std::string> ("Input RHS filename") ?
          params->get<std::string> ("Input RHS filename") : "";
        const std::string outRhsFilename =
          params->isType<std::string> ("Output RHS filename") ?
          params->get<std::string> ("Output RHS filename") : "";

        bool gmresUnfriendly = false;
        if (params->isParameter ("Problem type")) {
          const std::string probType = params->get<std::string> ("Problem type");
          if (probType == "GMRES unfriendly") {
            gmresUnfriendly = true;
          }
        }

        // This applies only if generating the right-hand side(s), and
        // if the problem type is not "GMRES unfriendly".
        const bool randomSolution = ! gmresUnfriendly && inRhsFilename == "" &&
          params->isParameter ("Random solution") &&
          params->get<bool> ("Random solution");

        // Number of right-hand side(s).
        size_t numRHS = 1;
        if (inRhsFilename == "") {
          // Be liberal about the accepted type of this parameter.
          if (params->isType<size_t> ("Number of RHS")) {
            numRHS = params->get<size_t> ("Number of RHS");
          }
          else if (params->isType<long> ("Number of RHS")) {
            numRHS = as<size_t> (params->get<long> ("Number of RHS"));
          }
          else if (params->isType<int> ("Number of RHS")) {
            numRHS = as<size_t> (params->get<int> ("Number of RHS"));
          }
          // If reading the right-hand side(s) from a file, numRHS
          // will be (re)set accordingly.
        }

        RCP<const map_type> domainMap = A->getDomainMap();
        RCP<const map_type> rangeMap = A->getRangeMap();

        if (inRhsFilename != "") {
          // If reading the right-hand side(s) from a file, don't set
          // the exact solution(s).  Later we could try to read those
          // from a file too.
          *out << "Reading right-hand side(s) B from Matrix Market file" << endl;
          Teuchos::OSTab tab2 (out);
          typedef ::Tpetra::MatrixMarket::Reader<SparseMatrixType> reader_type;
          B = reader_type::readDenseFile (inRhsFilename, comm_, rangeMap,
                                          tolerant_, debug_);
          TEUCHOS_TEST_FOR_EXCEPTION(B.is_null (), std::runtime_error, "Failed "
            "to read right-hand side(s) B from Matrix Market file \""
            << inRhsFilename << "\".");
          numRHS = B->getNumVectors ();
        }

        // Construct X_guess (the initial guess for the solution of
        // AX=B) from the domain of the matrix A, and fill it with
        // zeros.
        {
          *out << "Constructing initial guess vector X_guess" << endl;
          X_guess = rcp (new MV (domainMap, numRHS));
          X_guess->putScalar (STS::zero());
        }

        // Compute the exact solution vector(s).
        X_exact = rcp (new MV (domainMap, numRHS));
        if (inRhsFilename != "") {
          // FIXME (mfh 03 Apr 2012) We don't compute the exact
          // solution vector(s) in this case.  We could always fire up
          // Amesos2 and do a direct solve, but instead, we just fill
          // with zeros and be done with it.
          X_exact->putScalar (STS::zero ());
        }
        else { // if (inRhsFilename == "")
          // Our choice of exact solution and right-hand side depend
          // on the test problem.  If we generated the
          // GMRES-unfriendly example, we need B = e_{globalNumRows}
          // and therefore X_exact = e_{globalNumRows} as well.
          // Otherwise, we pick X_exact first and compute B via the
          // product B = A * X_exact.
          X_exact = rcp (new MV (domainMap, numRHS));

          // Construct the right-hand side B from the range of the
          // matrix A.  Don't just clone X_guess, since the range may
          // differ from the domain.
          B = rcp (new MV (rangeMap, numRHS));

          if (gmresUnfriendly) {
            *out << "Constructing B and X_exact for canonical \"GMRES-"
              "unfriendly\" example" << endl;
            X_exact->putScalar (STS::zero ());

            // Only the owning process is allowed to call
            // MV::replaceGlobalValue().  First check via the domain
            // (for X_exact) resp. range (for B) Maps whether the
            // calling process owns that row, and if so, make the
            // change there.  (Domain and range Maps of a CrsMatrix
            // should never be overlapping, so we don't have to worry
            // about multiple processes owning the row.)
            const GO rowToChange = as<GO> (A->getGlobalNumRows() - 1);
            if (domainMap->isNodeGlobalElement (rowToChange)) {
              X_exact->replaceGlobalValue (rowToChange, 0, STS::one ());
            }
            B->putScalar (STS::zero());
            if (rangeMap->isNodeGlobalElement (rowToChange)) {
              B->replaceGlobalValue (rowToChange, 0, STS::one ());
            }
          } else {
            if (randomSolution) {
              *out << "Randomizing X_exact" << endl;
              X_exact->randomize ();
            }
            else {
              // Construct the exact solution vector and fill it with
              // all ones.  Tacky, but deterministic.  Not good if we
              // expect A to be singular with rigid body modes.
              *out << "Setting X_exact = [1; ...; 1]" << endl;
              X_exact->putScalar (STS::one ());
            }

            // Compute the right-hand side B := A*X_exact.
            *out << "Computing B := A*X_exact" << endl;
            A->apply (*X_exact, *B);
          }
        }
      }

      /// \brief Same as \c makeMatrix() followed by \c makeVectors().
      ///
      /// This method first uses \c makeMatrix() with the given
      /// parameters to construct a sparse matrix A.  It then uses \c
      /// makeVectors() with the same parameter list to construct the
      /// right-hand side, initial guess, and exact solution vectors.
      ///
      /// See the documentation of \c makeMatrix() and \c
      /// makeVectors() for the parameters that this method accepts.
      void
      makeProblem (Teuchos::RCP<sparse_matrix_type>& A,
                   Teuchos::RCP<multivector_type>& X_guess,
                   Teuchos::RCP<multivector_type>& X_exact,
                   Teuchos::RCP<multivector_type>& B,
                   const Teuchos::RCP<Teuchos::ParameterList>& params)
      {
        A = makeMatrix (params);
        makeVectors (A, X_guess, X_exact, B, params);
      }

    private:
      //! Communicator over which to distribute sparse matrices and dense vectors.
      Teuchos::RCP<const Teuchos::Comm<int> > comm_;
      //! Output stream for indented verbose output.
      Teuchos::RCP<Teuchos::FancyOStream> out_;
      //! Whether to parse files tolerantly.
      bool tolerant_;
      //! Whether to print debug output.
      bool debug_;

      /// \brief Generate, distribute, and return a test sparse matrix.
      ///
      /// \param globalNumRows [in] Global (over all MPI process(es))
      ///   number of rows in the sparse matrix.
      ///
      /// \param symmetric [in] Whether to generate a symmetric test problem.
      ///   Storage is nonsymmetric regardless; symmetry here only applies to
      ///   the entries' locations and numerical values.
      ///
      /// \return The sparse matrix (global, distributed).
      Teuchos::RCP<sparse_matrix_type>
      generateTestMatrix (const global_ordinal_type globalNumRows,
                          const bool symmetric)
      {
        using Teuchos::RCP;
        using Teuchos::rcp;
        using Teuchos::tuple;
        using ::Tpetra::createUniformContigMapWithNode;
        using std::endl;
        typedef local_ordinal_type LO;
        typedef global_ordinal_type GO;
        typedef node_type NT;

        // For a square matrix, we only need a Map for the range of the matrix.
        RCP<const map_type> pRangeMap =
          createUniformContigMapWithNode<LO, GO, NT> (globalNumRows, comm_);
        // The sparse matrix object to fill.
        RCP<sparse_matrix_type> pMat = rcp (new sparse_matrix_type (pRangeMap, 0));

        const int myRank = comm_->getRank();
        if (myRank == 0) {
          const scalar_type leftVal = -STS::one();
          const scalar_type rightVal = symmetric ? -STS::one() : +2*STS::one();
          // Boost the diagonal if nonsymmetric, hopefully to make
          // convergence a little faster.  The point isn't to make a
          // difficult matrix, just an initial test of a solver.
          const scalar_type centerVal = symmetric ? 2*STS::one() : 8*STS::one();

          for (GO curRow = 0; curRow < globalNumRows; ++curRow) {
            if (curRow > 0) {
              pMat->insertGlobalValues (curRow, tuple(curRow-1), tuple(leftVal));
            }
            pMat->insertGlobalValues (curRow, tuple(curRow), tuple(centerVal));
            if (curRow < globalNumRows-1) {
              pMat->insertGlobalValues (curRow, tuple(curRow+1), tuple(rightVal));
            }
          }
        }
        // Make sure Rank 0 is done filling in the matrix.
        Teuchos::barrier (*comm_);
        // fillComplete() doesn't need any arguments if the domain and
        // range maps are the same.
        pMat->fillComplete ();
        return pMat;
      }

      /// \brief Generate a sparse matrix for which GMRES converges slowly.
      ///
      /// \param globalNumRows [in] Global number of rows (and
      ///   columns) in the sparse matrix.
      Teuchos::RCP<sparse_matrix_type>
      generateGmresUnfriendlyMatrix (const global_ordinal_type globalNumRows)
      {
        using Teuchos::RCP;
        using Teuchos::rcp;
        using Teuchos::tuple;
        using ::Tpetra::createUniformContigMapWithNode;
        using std::endl;

        typedef typename SparseMatrixType::local_ordinal_type LO;
        typedef typename SparseMatrixType::global_ordinal_type GO;
        typedef typename SparseMatrixType::node_type NT;

        // For a square matrix, we only need a Map for the range of the matrix.
        RCP<const map_type> pRangeMap =
          createUniformContigMapWithNode<LO, GO, NT> (globalNumRows, comm_);

        // The sparse matrix object to fill.
        RCP<sparse_matrix_type> pMat = rcp (new sparse_matrix_type (pRangeMap, 0));

        const int myRank = comm_->getRank();
        if (myRank == 0) {
          const scalar_type val = STS::one();
          for (GO curRow = 0; curRow < globalNumRows; ++curRow) {
            const GO curCol = (curRow == 0) ? (globalNumRows-1) : (curRow-1);
            pMat->insertGlobalValues (curRow, tuple(curCol), tuple(val));
          }
        }
        // Make sure Rank 0 is done filling in the matrix.
        Teuchos::barrier (*comm_);
        // fillComplete() doesn't need any arguments if the domain and
        // range maps are the same.
        pMat->fillComplete ();
        return pMat;
      }

      //! Return an output stream that only prints in debug mode.
      Teuchos::RCP<Teuchos::FancyOStream>
      getDebugStream () const
      {
        using Teuchos::getFancyOStream;

        if (debug_) {
          return getFancyOStream (out_);
        } else {
          return getFancyOStream (Teuchos::rcp (new Teuchos::oblackholestream));
        }
      }
    };

  } // namespace Tpetra
} // namespace Belos

#endif // __Belos_TpetraTestFramework_hpp
