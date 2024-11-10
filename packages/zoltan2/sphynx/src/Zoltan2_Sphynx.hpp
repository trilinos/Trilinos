// @HEADER
// *****************************************************************************
//                            Sphynx
//
// Copyright 2020 NTESS and the Sphynx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_SPHYNXALGORITHM_HPP_
#define _ZOLTAN2_SPHYNXALGORITHM_HPP_


////////////////////////////////////////////////////////////////////////////////
// This file contains the Sphynx algorithm.
//
// Sphynx is a graph partitioning algorithm that is based on a spectral method.
// It has three major steps:
// 1) compute the Laplacian matrix of the input graph,
// 2) compute logK+1 eigenvectors on the Laplacian matrix,
// 3) use eigenvector entries as coordinates and compute a K-way partition on
//    them using a geometric method.
//
// Step1: Sphynx provides three eigenvalue problems and hence Laplacian matrix:
//        i) combinatorial (Lx = \lambdax, where L = A - D)
//        ii) generalized (Lx = \lambdaDx, where L = A - D)
//        iii) normalized  (L_nx, \lambdax, where Ln = D^{-1/2}LD^{-1/2}
//             and L = A - D)
//
// Step2: Sphynx calls the LOBPCG algorithm provided in Anasazi to obtain
//        logK+1 eigenvectors. 
//        (Alternately, uses the experimental Randomized eigensolver.)
// Step3: Sphynx calls the MJ algorithm provided in Zoltan2Core to compute the
//        partition.
////////////////////////////////////////////////////////////////////////////////

#include "Zoltan2Sphynx_config.h"

#include "Zoltan2_XpetraCrsGraphAdapter.hpp"
#include "Zoltan2_XpetraMultiVectorAdapter.hpp"

#include "Zoltan2_CoordinateModel.hpp"
#include "Zoltan2_AlgMultiJagged.hpp"

#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziGeneralizedDavidsonSolMgr.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziRandomizedSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziTpetraAdapter.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosTpetraOperator.hpp"

#include "Ifpack2_Factory.hpp"

#include "Teuchos_TimeMonitor.hpp"

#ifdef HAVE_ZOLTAN2SPHYNX_MUELU
#include "MueLu_CreateTpetraPreconditioner.hpp"
#endif

namespace Zoltan2 {

  template <typename Adapter>
    class Sphynx : public Algorithm<Adapter>
  {

    public:

      using scalar_t = double; // Sphynx with scalar_t=double obtains better cutsize
      using lno_t = typename Adapter::lno_t;
      using gno_t = typename Adapter::gno_t;
      using node_t = typename Adapter::node_t;
      using offset_t = typename Adapter::offset_t;
      using part_t = typename Adapter::part_t;
      using weight_t = typename Adapter::scalar_t;

      using graph_t = Tpetra::CrsGraph<lno_t, gno_t, node_t>;
      using matrix_t = Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t>;
      using mvector_t = Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t>;
      using op_t = Tpetra::Operator<scalar_t, lno_t, gno_t, node_t>;

      enum problemType {COMBINATORIAL, GENERALIZED, NORMALIZED};
      typedef Anasazi::MultiVecTraits<scalar_t,mvector_t> MVT;

      ///////////////////////////////////////////////////////////////////////////
      ///////////////////////// CONSTRUCTORS ////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////

      // Takes the graph from the input adapter and computes the Laplacian matrix
      Sphynx(const RCP<const Environment> &env,
          const RCP<Teuchos::ParameterList> &params,
          const RCP<Teuchos::ParameterList> &sphynxParams,
          const RCP<const Comm<int> > &comm,
          const RCP<const XpetraCrsGraphAdapter<graph_t> > &adapter) :
        env_(env), params_(params), sphynxParams_(sphynxParams), comm_(comm), adapter_(adapter)
    {

      // Don't compute the Laplacian if the number of requested parts is 1
      const Teuchos::ParameterEntry *pe = params_->getEntryPtr("num_global_parts");
      numGlobalParts_ = pe->getValue<int>(&numGlobalParts_);
      if(numGlobalParts_ > 1){

        Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Sphynx::Laplacian"));

        // The verbosity is common throughout the algorithm, better to check and set now
        pe = sphynxParams_->getEntryPtr("sphynx_verbosity");
        if (pe)
          verbosity_ = pe->getValue<int>(&verbosity_);

        // Do we need to pre-process the input?
        pe = sphynxParams_->getEntryPtr("sphynx_skip_preprocessing");
        if (pe)
          skipPrep_ = pe->getValue<bool>(&skipPrep_);

        // Get the graph from XpetraCrsGraphAdapter if skipPrep_ is true
        // We assume the graph has all of the symmetric and diagonal entries
        if(skipPrep_)
          graph_ = adapter_->getUserGraph();
        else {
          throw std::runtime_error("\nSphynx Error: Preprocessing has not been implemented yet.\n");
        }

        // Check if the graph is regular
        determineRegularity();

        // Set Sphynx defaults: preconditioner, problem type, tolerance, initial vectors.
        setDefaults();

        // Compute the Laplacian matrix
        computeLaplacian();

        if(problemType_ == GENERALIZED)
          computeDegreeMatrix();

      }
    }

      ///////////////////////////////////////////////////////////////////////////
      ///////////////////// FORWARD DECLARATIONS  ///////////////////////////////
      ///////////////////////////////////////////////////////////////////////////

      void partition(const Teuchos::RCP<PartitioningSolution<Adapter> > &solution);

      int AnasaziWrapper(const int numEigenVectors);

      template<typename problem_t>
        void setPreconditioner(Teuchos::RCP<problem_t> &problem);

      template<typename problem_t>
        void setMueLuPreconditioner(Teuchos::RCP<problem_t> &problem);

      template<typename problem_t>
        void setJacobiPreconditioner(Teuchos::RCP<problem_t> &problem);

      template<typename problem_t>
        void setPolynomialPreconditioner(Teuchos::RCP<problem_t> &problem);

      void eigenvecsToCoords(Teuchos::RCP<mvector_t> &eigenVectors,
          int computedNumEv,
          Teuchos::RCP<mvector_t> &coordinates);

      void computeWeights(std::vector<const weight_t *> vecweights,
          std::vector<int> strides);

      void MJwrapper(const Teuchos::RCP<const mvector_t> &coordinates,
          std::vector<const weight_t *> weights,
          std::vector<int> strides,
          const Teuchos::RCP<PartitioningSolution<Adapter>> &solution);

      void setUserEigenvectors(const Teuchos::RCP<mvector_t> &userEvects);
      ///////////////////////////////////////////////////////////////////////////
      ///////////////////// MEMBER FUNCTIONS - Laplacian-related ones ///////////
      ///////////////////////////////////////////////////////////////////////////

      ///////////////////////////////////////////////////////////////////////////
      // Determine input graph's regularity = maxDegree/AvgDegree < 10.
      // If Laplacian type is not specified, then use combinatorial for regular
      // and generalized for irregular.
      // MueLu settings will differ depending on the regularity, too.
      void determineRegularity()
      {
        // Get the row pointers in the host
        auto rowOffsets = graph_->getLocalGraphDevice().row_map;
        auto rowOffsets_h = Kokkos::create_mirror_view(rowOffsets);
        Kokkos::deep_copy(rowOffsets_h, rowOffsets);

        // Get size information
        const size_t numGlobalEntries = graph_->getGlobalNumEntries();
        const size_t numLocalRows = graph_->getLocalNumRows();
        const size_t numGlobalRows = graph_->getGlobalNumRows();

        // Compute local maximum degree
        size_t localMax = 0;
        for(size_t i = 0; i < numLocalRows; i++){
          if(rowOffsets_h(i+1) - rowOffsets_h(i) - 1 > localMax)
            localMax = rowOffsets_h(i+1) - rowOffsets_h(i) - 1;
        }

        // Compute global maximum degree
        size_t globalMax = 0;
        Teuchos::reduceAll<int, size_t> (*comm_, Teuchos::REDUCE_MAX, 1, &localMax, &globalMax);

        double avg = static_cast<double>(numGlobalEntries-numGlobalRows)/numGlobalRows;
        double maxOverAvg = static_cast<double>(globalMax)/avg;

        // Use generalized Laplacian on irregular graphs
        irregular_ = false;
        if(maxOverAvg > 10) {
          irregular_ = true;
        }

        // Let the user know about what we determined if verbose
        if(verbosity_ > 0) {
          if(comm_->getRank() == 0) {
            std::cout << "Regularity of Graph ----------------" << std::endl;  
            std::cout << "  Maximum degree: " << globalMax << std::endl;
            std::cout << "  Average degree: " << avg << std::endl;
            std::cout << "  Max/Avg:        " << globalMax/avg << std::endl;
            std::cout << "  Regular graph?: " << !irregular_ << std::endl;
          }
        }
      }

      ///////////////////////////////////////////////////////////////////////////
      // If preconditioner type is not specified:
      //    use muelu if it is enabled, and jacobi otherwise.
      // If eigenvalue problem type is not specified:
      //    use combinatorial for regular and
      //        normalized for irregular with polynomial preconditioner,
      //        generalized for irregular with other preconditioners.
      // If convergence tolerance is not specified:
      //    use 1e-3 for regular with jacobi and polynomial, and 1e-2 otherwise.
      // If how to decide the initial vectors is not specified:
      //    use random for regular and constant for irregular
      void setDefaults()
      {

        // Set the default preconditioner to muelu if it is enabled, jacobi otherwise.
        precType_ = "jacobi";
#ifdef HAVE_ZOLTAN2SPHYNX_MUELU
        precType_ = "muelu";
#endif

        // Override the preconditioner type with the user's preference
        const Teuchos::ParameterEntry *pe = sphynxParams_->getEntryPtr("sphynx_preconditioner_type");
        if (pe) {
          precType_ = pe->getValue<std::string>(&precType_);
          if(precType_ != "muelu" && precType_ != "jacobi" && precType_ != "polynomial")
            throw std::runtime_error("\nSphynx Error: " + precType_ + " is not recognized as a preconditioner.\n"
                + "              Possible values: muelu (if enabled), jacobi, and polynomial\n");
        }

        solverType_ = sphynxParams_->get("sphynx_eigensolver","LOBPCG");
        TEUCHOS_TEST_FOR_EXCEPTION(!(solverType_ == "LOBPCG" || solverType_ == "randomized" || solverType_ == "BlockDavidson" || solverType_ == "GeneralizedDavidson" || solverType_ == "BlockKrylovSchur" ), 
            std::invalid_argument, "Sphynx: sphynx_eigensolver must be set to LOBPCG, randomized, BlockDavidson, GeneralizedDavidson, or BlockKrylovSchur.");

        // Set the default problem type
        problemType_ = COMBINATORIAL;
        if(irregular_) {
          if(precType_ == "polynomial")
            problemType_ = NORMALIZED;
          else
            problemType_ = GENERALIZED;
        }

        // Override the problem type with the user's preference
        pe = sphynxParams_->getEntryPtr("sphynx_problem_type");
        if (pe) {
          std::string probType = "";
          probType = pe->getValue<std::string>(&probType);

          if(probType == "combinatorial")
            problemType_ = COMBINATORIAL;
          else if(probType == "generalized")
            problemType_ = GENERALIZED;
          else if(probType == "normalized")
            problemType_ = NORMALIZED;
          else
            throw std::runtime_error("\nSphynx Error: " + probType + " is not recognized as a problem type.\n"
                + "              Possible values: combinatorial, generalized, and normalized.\n");
        }


        // Set the default for tolerance
        tolerance_ = 1.0e-2;
        if(!irregular_ && (precType_ == "jacobi" || precType_ == "polynomial"))
          tolerance_ = 1.0e-3;


        // Override the tolerance with the user's preference
        pe = sphynxParams_->getEntryPtr("sphynx_tolerance");
        if (pe)
          tolerance_ = pe->getValue<scalar_t>(&tolerance_);


        // Set the default for initial vectors
        randomInit_ = true;
        if(irregular_)
          randomInit_ = false;

        // Override the initialization method with the user's preference
        pe = sphynxParams_->getEntryPtr("sphynx_initial_guess");
        if (pe) {
          std::string initialGuess = "";
          initialGuess = pe->getValue<std::string>(&initialGuess);

          if(initialGuess == "random")
            randomInit_ = true;
          else if(initialGuess == "constants")
            randomInit_ = false;
          else
            throw std::runtime_error("\nSphynx Error: " + initialGuess + " is not recognized as an option for initial guess.\n"
                + "              Possible values: random and constants.\n");
        }

      }

      ///////////////////////////////////////////////////////////////////////////
      // Compute the Laplacian matrix depending on the eigenvalue problem type.
      // There are 3 options for the type: combinatorial, generalized, and normalized.
      // Combinatorial and generalized share the same Laplacian but generalized
      // also needs a degree matrix.
      void computeLaplacian()
      {
        if(solverType_ == "randomized")
          laplacian_ = computeNormalizedLaplacian(true);
        else if(problemType_ == NORMALIZED)
          laplacian_ = computeNormalizedLaplacian();
        else
          laplacian_ = computeCombinatorialLaplacian();
      }

      ///////////////////////////////////////////////////////////////////////////
      // Compute a diagonal matrix with the vertex degrees in the input graph
      void computeDegreeMatrix()
      {

        // Get the row pointers in the host
        auto rowOffsets = graph_->getLocalGraphHost().row_map;

        // Create the degree matrix with max row size set to 1
        Teuchos::RCP<matrix_t> degMat(new matrix_t (graph_->getRowMap(),
              graph_->getRowMap(),
              1));

        scalar_t *val = new scalar_t[1];
        lno_t *ind = new lno_t[1];
        lno_t numRows = static_cast<lno_t>(graph_->getLocalNumRows());

        // Insert the diagonal entries as the degrees
        for (lno_t i = 0; i < numRows; ++i) {
          val[0] = rowOffsets(i+1) - rowOffsets(i) - 1;
          ind[0] = i;
          degMat->insertLocalValues(i, 1, val, ind);
        }
        delete [] val;
        delete [] ind;

        degMat->fillComplete(graph_->getDomainMap(), graph_->getRangeMap());

        degMatrix_ = degMat;
      }

      ///////////////////////////////////////////////////////////////////////////
      // Compute the combinatorial Laplacian: L = D - A.
      // l_ij = degree of vertex i if i = j
      // l_ij = -1 if i != j and a_ij != 0
      // l_ij = 0 if i != j and a_ij = 0
      Teuchos::RCP<matrix_t> computeCombinatorialLaplacian()
      {
        using range_policy = Kokkos::RangePolicy<
          typename node_t::device_type::execution_space, Kokkos::IndexType<lno_t>>;
        using values_view_t = Kokkos::View<scalar_t*, typename node_t::device_type>;
        using offset_view_t = Kokkos::View<size_t*, typename node_t::device_type>;

        const size_t numEnt = graph_->getLocalNumEntries();
        const size_t numRows = graph_->getLocalNumRows();

        // Create new values for the Laplacian, initialize to -1
        values_view_t newVal (Kokkos::view_alloc("CombLapl::val", Kokkos::WithoutInitializing), numEnt);
        Kokkos::deep_copy(newVal, -1);

        // Get the diagonal offsets
        offset_view_t diagOffsets(Kokkos::view_alloc("Diag Offsets", Kokkos::WithoutInitializing), numRows);
        graph_->getLocalDiagOffsets(diagOffsets);

        // Get the row pointers in the host
        auto rowOffsets = graph_->getLocalGraphDevice().row_map;

        // Compute the diagonal entries as the vertex degrees
        Kokkos::parallel_for("Combinatorial Laplacian", range_policy(0, numRows),
            KOKKOS_LAMBDA(const lno_t i){
            newVal(rowOffsets(i) + diagOffsets(i)) = rowOffsets(i+1) - rowOffsets(i) - 1;
            }
            );
        Kokkos::fence ();

        // Create the Laplacian matrix using the input graph and with the new values
        Teuchos::RCP<matrix_t> laplacian (new matrix_t(graph_, newVal));
        laplacian->fillComplete (graph_->getDomainMap(), graph_->getRangeMap());

        // Create the Laplacian maatrix using the input graph and with the new values
        return laplacian;
      }

      ///////////////////////////////////////////////////////////////////////////
      // For AHat = false:
      // Compute the normalized Laplacian: L_N = D^{-1/2} L D^{-1/2}, where L = D - A.
      // l_ij = 1
      // l_ij = -1/(sqrt(deg(v_i))sqrt(deg(v_j)) if i != j and a_ij != 0
      // l_ij = 0 if i != j and a_ij = 0
      //
      // For AHat = true:
      // AHat is turned to true if (and only if) using the randomized Eigensolver.
      // For the randomized Eigensolver, we find eigenvalues of the matrix
      // AHat =  2*I - L_N, and this is the matrix computed by this function.
      Teuchos::RCP<matrix_t> computeNormalizedLaplacian(bool AHat = false)
      {
        using range_policy = Kokkos::RangePolicy<
          typename node_t::device_type::execution_space, Kokkos::IndexType<lno_t>>;
        using values_view_t = Kokkos::View<scalar_t*, typename node_t::device_type>;
        using offset_view_t = Kokkos::View<size_t*, typename node_t::device_type>;
        using vector_t = Tpetra::Vector<scalar_t, lno_t, gno_t, node_t>;
        using dual_view_t = typename vector_t::dual_view_type;
        using KAT = Kokkos::ArithTraits<scalar_t>;

        const size_t numEnt = graph_->getLocalNumEntries();
        const size_t numRows = graph_->getLocalNumRows();

        // Create new values for the Laplacian, initialize to -1
        values_view_t newVal (Kokkos::view_alloc("NormLapl::val", Kokkos::WithoutInitializing), numEnt);
        if(AHat){
          Kokkos::deep_copy(newVal, 1);
        }
        else{
          Kokkos::deep_copy(newVal, -1);
        }

        // D^{-1/2}
        dual_view_t dv ("MV::DualView", numRows, 1);
        auto deginvsqrt = dv.d_view;

        // Get the diagonal offsets
        offset_view_t diagOffsets(Kokkos::view_alloc("Diag Offsets", Kokkos::WithoutInitializing), numRows);
        graph_->getLocalDiagOffsets(diagOffsets);

        // Get the row pointers
        auto rowOffsets = graph_->getLocalGraphDevice().row_map;

        //      if(!AHat){
        // Compute the diagonal entries as the vertex degrees
        Kokkos::parallel_for("Combinatorial Laplacian", range_policy(0, numRows),
            KOKKOS_LAMBDA(const lno_t i){
            newVal(rowOffsets(i) + diagOffsets(i)) = rowOffsets(i+1) - rowOffsets(i) - 1;
            deginvsqrt(i,0) = 1.0/KAT::sqrt(rowOffsets(i+1) - rowOffsets(i) - 1);
            }
            );
        Kokkos::fence ();

        // Create the Laplacian graph using the same graph structure with the new values
        Teuchos::RCP<matrix_t> laplacian (new matrix_t(graph_, newVal));
        laplacian->fillComplete (graph_->getDomainMap(), graph_->getRangeMap());

        // Create the vector for D^{-1/2}
        vector_t degInvSqrt(graph_->getRowMap(), dv);

        // Normalize the laplacian matrix as D^{-1/2} L D^{-1/2}
        laplacian->leftScale(degInvSqrt);
        laplacian->rightScale(degInvSqrt);

        return laplacian;
      }

      ///////////////////////////////////////////////////////////////////////////
      ///////////////////////////// DATA MEMBERS ////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////

    private:

      // User-provided members
      const Teuchos::RCP<const Environment> env_;
      Teuchos::RCP<Teuchos::ParameterList> params_;
      Teuchos::RCP<Teuchos::ParameterList> sphynxParams_;
      const Teuchos::RCP<const Teuchos::Comm<int> > comm_;
      const Teuchos::RCP<const Adapter> adapter_;

      // Internal members
      int numGlobalParts_;
      Teuchos::RCP<const graph_t> graph_;
      Teuchos::RCP<matrix_t> laplacian_;
      Teuchos::RCP<const matrix_t> degMatrix_;
      Teuchos::RCP<mvector_t> eigenVectors_;

      bool irregular_;           // decided internally
      std::string  precType_;    // obtained from user params or decided internally
      std::string  solverType_;  // obtained from user params or decided internally
      problemType problemType_;  // obtained from user params or decided internally
      double tolerance_;         // obtained from user params or decided internally
      bool randomInit_;          // obtained from user params or decided internally
      int verbosity_;            // obtained from user params
      bool skipPrep_;            // obtained from user params
  };

  ///////////////////////////////////////////////////////////////////////////
  /////////////////////// MORE MEMBER FUNCTIONS  ////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Allows the user to manually set eigenvectors for the Sphynx partitioner
  // to use rather than solving for them with Anasazi. Mainly intended 
  // for debugging purposes.
  ///////////////////////////////////////////////////////////////////////////
  template <typename Adapter>
    void Sphynx<Adapter>::setUserEigenvectors(const Teuchos::RCP<mvector_t> &userEvects)
    {
      eigenVectors_ = userEvects; 
    }

  ///////////////////////////////////////////////////////////////////////////
  // Compute a partition using the Laplacian matrix (and possibly the degree
  // matrix as well). First call LOBPCG (or Randomized solver) 
  // to compute logK+1 eigenvectors, then
  // transform the eigenvectors to coordinates, and finally call MJ to compute
  // a partition on the coordinates.
  template <typename Adapter>
    void Sphynx<Adapter>::partition(const Teuchos::RCP<PartitioningSolution<Adapter>> &solution)
    {
      // Return a trivial solution if only one part is requested
      if(numGlobalParts_ == 1) {

        size_t numRows =adapter_->getUserGraph()->getLocalNumRows();
        Teuchos::ArrayRCP<part_t> parts(numRows,0);
        solution->setParts(parts);

        return;

      }

      // The number of eigenvectors to be computed
      int numEigenVectors = (int) log2(numGlobalParts_)+1;
      int computedNumEv;

      if(eigenVectors_ == Teuchos::null){
        // Compute the eigenvectors using LOBPCG
        // or Randomized eigensolver
        computedNumEv = Sphynx::AnasaziWrapper(numEigenVectors);

        if(computedNumEv <= 1 && (solverType_ == "LOBPCG" || solverType_ == "GeneralizedDavidson" || solverType_ == "BlockDavidson" || solverType_ == "BlockKrylovSchur")) { 
        throw
            std::runtime_error("\nAnasazi Error: LOBPCGSolMgr::solve() returned unconverged.\n"
                "Sphynx Error:  LOBPCG could not compute any eigenvectors.\n"
                "               Increase either max iters or tolerance.\n");
        }
      }
      else{
        // Use eigenvectors provided by user. 
        computedNumEv = (int) eigenVectors_->getNumVectors();
        if(computedNumEv <= numEigenVectors) {
          throw 
            std::runtime_error("\nSphynx Error: Number of eigenvectors given by user\n"
                " is less than number of Eigenvectors needed for partition." );
        }
      }

      // Transform the eigenvectors into coordinates 
      Teuchos::RCP<mvector_t> coordinates;

      Sphynx::eigenvecsToCoords(eigenVectors_, computedNumEv, coordinates);

      // Get the weights from the adapter
      std::vector<const weight_t *> weights;
      std::vector<int> wstrides;
      Sphynx::computeWeights(weights, wstrides);


      // Compute the partition using MJ on coordinates
      Sphynx::MJwrapper(coordinates, weights, wstrides, solution);

    }

  ///////////////////////////////////////////////////////////////////////////
  // Call LOBPCG on the Laplacian matrix.
  // Or use the randomized eigensolver.
  template <typename Adapter>
    int Sphynx<Adapter>::AnasaziWrapper(const int numEigenVectors)
    {

      Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Sphynx::Anasazi"));

      // Set defaults for the parameters
      // and get user-set values.
      std::string which = (solverType_ == "randomized" ? "LM" : "SR");
      std::string ortho = "SVQB";
      bool relConvTol = false;
      int maxIterations = sphynxParams_->get("sphynx_max_iterations",1000);
      int blockSize = sphynxParams_->get("sphynx_block_size",numEigenVectors);
      int orthoFreq = sphynxParams_->get("sphynx_ortho_freq", 0);
      int resFreq = sphynxParams_->get("sphynx_res_freq", 0);
      bool isHermitian = true;
      bool relLockTol = false;
      bool lock = false;
      bool useFullOrtho = sphynxParams_->get("sphynx_use_full_ortho",true);

      // Information to output in a verbose run
      int numfailed = 0;
      int iter = 0;

      // Set Anasazi verbosity level
      int anasaziVerbosity = Anasazi::Errors + Anasazi::Warnings;
      if (verbosity_ >= 1)  // low
        anasaziVerbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
      if (verbosity_ >= 2)  // medium
        anasaziVerbosity += Anasazi::IterationDetails;
      if (verbosity_ >= 3)  // high
        anasaziVerbosity += Anasazi::StatusTestDetails + Anasazi::OrthoDetails
          + Anasazi::Debug;

      // Create the parameter list to pass into solver
      Teuchos::ParameterList anasaziParams;
      anasaziParams.set("Verbosity", anasaziVerbosity);
      anasaziParams.set("Which", which);
      anasaziParams.set("Convergence Tolerance", tolerance_);
      anasaziParams.set("Maximum Iterations", maxIterations);
      anasaziParams.set("Block Size", blockSize);
      anasaziParams.set("Relative Convergence Tolerance", relConvTol);
      anasaziParams.set("Orthogonalization", ortho);
      anasaziParams.set("Use Locking", lock);
      anasaziParams.set("Relative Locking Tolerance", relLockTol);
      anasaziParams.set("Full Ortho", useFullOrtho);
      anasaziParams.set("Orthogonalization Frequency", orthoFreq);
      anasaziParams.set("Residual Frequency", resFreq);
     
      if(solverType_ == "GeneralizedDavidson" || solverType_ == "BlockKrylovSchur" || solverType_ == "BlockDavidson"){ 
        anasaziParams.set( "Num Blocks", maxIterations+1 );
        anasaziParams.set( "Maximum Restarts", 0 );
        anasaziParams.set( "Maximum Subspace Dimension", (maxIterations+1)*blockSize );
      }

      // Create and set initial vectors
      auto map = laplacian_->getRangeMap();
      Teuchos::RCP<mvector_t> ivec( new mvector_t(map, numEigenVectors));

      if (randomInit_) {
        // 0-th vector constant 1, rest random
        Anasazi::MultiVecTraits<scalar_t, mvector_t>::MvRandom(*ivec);
        ivec->getVectorNonConst(0)->putScalar(1.);
      }
      else { // This implies we will use constant initial vectors.
        // 0-th vector constant 1, other vectors constant per block
        // Create numEigenVectors blocks, but only use numEigenVectors-1 of them.
        // This assures orthogonality.
        ivec->getVectorNonConst(0)->putScalar(1.);
        for (int j = 1; j < numEigenVectors; j++)
          ivec->getVectorNonConst(j)->putScalar(0.);

        gno_t blkSize = map->getGlobalNumElements() / numEigenVectors;
        TEUCHOS_TEST_FOR_EXCEPTION(blkSize <= 0, std::runtime_error, "Blocksize too small for \"constants\" initial guess. Try \"random\".");

        for (size_t lid = 0; lid < ivec->getLocalLength(); lid++) {
          gno_t gid = map->getGlobalElement(lid);
          for (int j = 1; j < numEigenVectors; j++){
            if (((j-1)*blkSize <= gid) && (j*blkSize > gid))
              ivec->replaceLocalValue(lid,j,1.);
          }
        }
      }

      // Create the eigenproblem to be solved
      using problem_t = Anasazi::BasicEigenproblem<scalar_t, mvector_t, op_t>;
      Teuchos::RCP<problem_t> problem (new problem_t(laplacian_, ivec));
      problem->setHermitian(isHermitian);
      problem->setNEV(numEigenVectors);

      if(solverType_ != "randomized"){
        // Set preconditioner
        Sphynx::setPreconditioner(problem);
        if(problemType_ == Sphynx::GENERALIZED) problem->setM(degMatrix_);
      }

      // Inform the eigenproblem that you are finished passing it information
      bool boolret = problem->setProblem();
      if (boolret != true) {
        throw std::runtime_error("\nAnasazi::BasicEigenproblem::setProblem() returned with error.\n");
      }
      // Set Eigensolver
      Teuchos::RCP<Anasazi::SolverManager<scalar_t, mvector_t, op_t>> solver;

      if(solverType_ == "LOBPCG"){
        solver = Teuchos::rcp(new Anasazi::LOBPCGSolMgr<scalar_t, mvector_t, op_t>(problem, anasaziParams));
      }
      else if(solverType_ == "BlockDavidson"){
        solver = Teuchos::rcp(new Anasazi::BlockDavidsonSolMgr<scalar_t, mvector_t, op_t>(problem, anasaziParams));
      }
      else if(solverType_ == "GeneralizedDavidson"){
        solver = Teuchos::rcp(new Anasazi::GeneralizedDavidsonSolMgr<scalar_t, mvector_t, op_t>(problem, anasaziParams));
      }
      else if(solverType_ == "BlockKrylovSchur"){
        solver = Teuchos::rcp(new Anasazi::BlockKrylovSchurSolMgr<scalar_t, mvector_t, op_t>(problem, anasaziParams));
      }
      else{
        solver = Teuchos::rcp(new Anasazi::Experimental::RandomizedSolMgr<scalar_t, mvector_t, op_t>(problem, anasaziParams));
        //numFailed = solver->getNumFailed();
      }

      if (verbosity_ > 0 && comm_->getRank() == 0)
        anasaziParams.print(std::cout);

      // Solve the problem
      Anasazi::ReturnType returnCode = solver->solve();

      // Check convergence, niters, and solvetime
      iter = solver->getNumIters();
      //solvetime = (solver->getTimers()[0])->totalElapsedTime();

      // Retrieve the solution
      using solution_t = Anasazi::Eigensolution<scalar_t, mvector_t>;
      solution_t sol = problem->getSolution();
      eigenVectors_ = sol.Evecs;
      int numev = sol.numVecs;
      std::vector<Anasazi::Value<double>> eigenValues_ = sol.Evals;

      // Summarize iteration counts and solve time
      if (verbosity_ > 0 && comm_->getRank() == 0) {
        std::cout << std::endl;
        std::cout << "ANASAZI SUMMARY" << std::endl;
        std::cout << "Failed to converge:    " << numfailed << std::endl;
        std::cout << "No of iterations :     " << iter << std::endl;
        std::cout << "Solve time:            " << std::endl; // solvetime << std::endl;
        std::cout << "No of comp. vecs. :    " << numev << std::endl;
      }

      std::cout << "Solver type: " << solverType_ << std::endl;
      // Compute residuals (LOBPCG does this internally)
      if(solverType_ == "randomized") {
        std::vector<double> normR(numev);
        mvector_t Aevec (map, numev);

        if (numev > 0) { 
          Teuchos::SerialDenseMatrix<int,double> T (numev, numev);
          T.putScalar(0.0);
          for (int i = 0; i < numev; ++i) T(i,i) = eigenValues_[i].realpart;
          laplacian_->apply(*eigenVectors_, Aevec);
          MVT::MvTimesMatAddMv(-1.0, *eigenVectors_, T, 1.0, Aevec);
          MVT::MvNorm(Aevec, normR);
        } // End residual computation

        if(comm_->getRank() == 0 && verbosity_ > 0) {
          std::cout << std::endl << "Solver manager returned "
            << (returnCode == Anasazi::Converged ? "converged." : "unconverged.")
            << std::endl << std::endl
            << std::setw(16) << "Eigenvalue"
            << std::setw(18) << "Direct Residual"
            << std::endl
            << "------------------------------------------------------" << std::endl;

          for (int i = 0; i < numev; ++i) {
            std::cout << std::setw(16) << 2.0-eigenValues_[i].realpart
              << std::setw(18) << normR[i] / eigenValues_[i].realpart
              << std::endl;
          }
          std::cout << "------------------------------------------------------" << std::endl;
        }
      }

      return numev;

    }

  ///////////////////////////////////////////////////////////////////////////
  template <typename Adapter>
    template <typename problem_t>
    void Sphynx<Adapter>::setPreconditioner(Teuchos::RCP<problem_t> &problem)
    {
      if(comm_->getRank() == 0) std::cout << "precType_ is: " << precType_ << std::endl;
      // Set the preconditioner
      if(precType_ == "muelu") {
        Sphynx<Adapter>::setMueLuPreconditioner(problem);
      }
      else if(precType_ == "polynomial") {
        Sphynx<Adapter>::setPolynomialPreconditioner(problem);
      }
      else if(precType_ == "jacobi") {
        Sphynx<Adapter>::setJacobiPreconditioner(problem);
      }
      // else: do we want a case where no preconditioning is applied?
    }

  ///////////////////////////////////////////////////////////////////////////
  template <typename Adapter>
    template <typename problem_t>
    void Sphynx<Adapter>::setMueLuPreconditioner(Teuchos::RCP<problem_t> &problem)
    {
#ifdef HAVE_ZOLTAN2SPHYNX_MUELU
      Teuchos::ParameterList paramList;
      if(verbosity_ == 0)
        paramList.set("verbosity", "none");
      else if(verbosity_ == 1)
        paramList.set("verbosity", "low");
      else if(verbosity_ == 2)
        paramList.set("verbosity", "medium");
      else if(verbosity_ == 3)
        paramList.set("verbosity", "high");
      else if(verbosity_ >= 4)
        paramList.set("verbosity", "extreme");

      paramList.set("smoother: type", "CHEBYSHEV");
      Teuchos::ParameterList smootherParamList;
      smootherParamList.set("chebyshev: degree", 3);
      smootherParamList.set("chebyshev: ratio eigenvalue", 7.0);
      smootherParamList.set("chebyshev: eigenvalue max iterations", irregular_ ? 100 : 10);
      paramList.set("smoother: params", smootherParamList);
      paramList.set("use kokkos refactor", true);
      paramList.set("transpose: use implicit", true);

      if(irregular_) {

        paramList.set("multigrid algorithm", "unsmoothed");

        paramList.set("coarse: type", "CHEBYSHEV");
        Teuchos::ParameterList coarseParamList;
        coarseParamList.set("chebyshev: degree", 3);
        coarseParamList.set("chebyshev: ratio eigenvalue", 7.0);
        paramList.set("coarse: params", coarseParamList);

        paramList.set("max levels", 5);
        paramList.set("aggregation: drop tol", 0.40);

      }
      using prec_t = MueLu::TpetraOperator<scalar_t, lno_t, gno_t, node_t>;
      Teuchos::RCP<prec_t> prec = MueLu::CreateTpetraPreconditioner<
        scalar_t, lno_t, gno_t, node_t>(laplacian_, paramList);

      problem->setPrec(prec);

#else
      throw std::runtime_error("\nSphynx Error: MueLu requested but not compiled into Trilinos.\n");
#endif
    }

  ///////////////////////////////////////////////////////////////////////////
  template <typename Adapter>
    template <typename problem_t>
    void Sphynx<Adapter>::setPolynomialPreconditioner(Teuchos::RCP<problem_t> &problem)
    {
      int verbosity2 = Belos::Errors;
      if(verbosity_ == 1)
        verbosity2 += Belos::Warnings;
      else if(verbosity_ == 2)
        verbosity2 += Belos::Warnings + Belos::FinalSummary;
      else if(verbosity_ == 3)
        verbosity2 += Belos::Warnings + Belos::FinalSummary + Belos::TimingDetails;
      else if(verbosity_ >= 4)
        verbosity2 += Belos::Warnings + Belos::FinalSummary + Belos::TimingDetails
          + Belos::StatusTestDetails;

      Teuchos::ParameterList paramList;
      paramList.set("Polynomial Type", "Roots");
      paramList.set("Orthogonalization","ICGS");
      paramList.set("Maximum Degree", laplacian_->getGlobalNumRows() > 100 ? 25 : 5);
      paramList.set("Polynomial Tolerance", 1.0e-6 );
      paramList.set("Verbosity", verbosity2 );
      paramList.set("Random RHS", false );
      paramList.set("Outer Solver", "");
      paramList.set("Timer Label", "Belos Polynomial Solve" );

      // Construct a linear problem for the polynomial solver manager
      using lproblem_t = Belos::LinearProblem<scalar_t, mvector_t, op_t>;
      Teuchos::RCP<lproblem_t> innerPolyProblem(new lproblem_t());
      innerPolyProblem->setOperator(laplacian_);

      using btop_t = Belos::TpetraOperator<scalar_t, lno_t, gno_t, node_t>;
      Teuchos::RCP<btop_t> polySolver(new btop_t(innerPolyProblem,
            Teuchos::rcpFromRef(paramList),
            "GmresPoly", true));
      problem->setPrec(polySolver);
    }

  ///////////////////////////////////////////////////////////////////////////
  template <typename Adapter>
    template <typename problem_t>
    void Sphynx<Adapter>::setJacobiPreconditioner(Teuchos::RCP<problem_t> &problem)
    {

      Teuchos::RCP<Ifpack2::Preconditioner<scalar_t, lno_t, gno_t, node_t>> prec;
      std::string precType = "RELAXATION";

      prec = Ifpack2::Factory::create<matrix_t> (precType, laplacian_);

      Teuchos::ParameterList precParams;
      precParams.set("relaxation: type", "Jacobi");
      precParams.set("relaxation: fix tiny diagonal entries", true);
      precParams.set("relaxation: min diagonal value", 1.0e-8);

      prec->setParameters(precParams);
      prec->initialize();
      prec->compute();

      problem->setPrec(prec);

    }

  ///////////////////////////////////////////////////////////////////////////
  // Transform the computed eigenvectors into coordinates
  template <typename Adapter>
    void Sphynx<Adapter>::eigenvecsToCoords(Teuchos::RCP<mvector_t> &eigenVectors,
        int computedNumEv,
        Teuchos::RCP<mvector_t> &coordinates)
    {
      // Extract the meaningful eigenvectors by getting rid of the first one
      Teuchos::Array<size_t> columns (computedNumEv-1);
      for (int j = 0; j < computedNumEv-1; ++j) {
        columns[j] = j+1;
      }
      coordinates = eigenVectors->subCopy (columns());
      coordinates->setCopyOrView (Teuchos::View);
    }


  ///////////////////////////////////////////////////////////////////////////
  // If user provided some weights, use them by getting them from the adapter.
  // If user didn't provide weights but told us to use degree as weight, do so.
  // If user neither provided weights nor told us what to do, use degree as weight.
  template <typename Adapter>
    void Sphynx<Adapter>::computeWeights(std::vector<const weight_t *> vecweights,
        std::vector<int> strides)
    {

      int numWeights = adapter_->getNumWeightsPerVertex();
      int numConstraints = numWeights > 0 ? numWeights : 1;

      size_t myNumVertices = adapter_->getLocalNumVertices();
      weight_t ** weights = new weight_t*[numConstraints];
      for(int j = 0; j < numConstraints; j++)
        weights[j] = new weight_t[myNumVertices];

      // Will be needed if we use degree as weight
      const offset_t *offset;
      const gno_t *columnId;

      // If user hasn't set any weights, use vertex degrees as weights
      if(numWeights == 0) {

        // Compute the weight of vertex i as the number of nonzeros in row i.
        adapter_->getEdgesView(offset, columnId);
        for (size_t i = 0; i < myNumVertices; i++)
          weights[0][i] = offset[i+1] - offset[i] - 1;

        vecweights.push_back(weights[0]);
        strides.push_back(1);
      }
      else {

        // Use the weights if there are any already set in the adapter
        for(int j = 0; j < numConstraints; j++) {

          if(adapter_->useDegreeAsVertexWeight(j)) {
            // Compute the weight of vertex i as the number of nonzeros in row i.
            adapter_->getEdgesView(offset, columnId);
            for (size_t i = 0; i < myNumVertices; i++)
              weights[j][i] = offset[i+1] - offset[i];
          }
          else{
            int stride;
            const weight_t *wgt = NULL;
            adapter_->getVertexWeightsView(wgt, stride, j);

            for (size_t i = 0; i < myNumVertices; i++)
              weights[j][i] = wgt[i];
          }

          vecweights.push_back(weights[j]);
          strides.push_back(1);

        }
      }

    }


  ///////////////////////////////////////////////////////////////////////////
  // Compute a partition by calling MJ on given coordinates with given weights
  template <typename Adapter>
    void Sphynx<Adapter>::MJwrapper(const Teuchos::RCP<const mvector_t> &coordinates,
        std::vector<const weight_t *> weights,
        std::vector<int> strides,
        const Teuchos::RCP<PartitioningSolution<Adapter>> &solution)
    {

      Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Sphynx::MJ"));

      using mvector_adapter_t = Zoltan2::XpetraMultiVectorAdapter<mvector_t>;
      using base_adapter_t = typename mvector_adapter_t::base_adapter_t;
      using mj_t = Zoltan2::Zoltan2_AlgMJ<mvector_adapter_t>;
      using solution_t = Zoltan2::PartitioningSolution<mvector_adapter_t>;

      size_t myNumVertices = coordinates->getLocalLength();

      // Create the base adapter for the multivector adapter
      Teuchos::RCP<mvector_adapter_t> adapcoordinates(new mvector_adapter_t(coordinates,
            weights,
            strides));
      Teuchos::RCP<const base_adapter_t> baseAdapter =
        Teuchos::rcp(dynamic_cast<const base_adapter_t *>(adapcoordinates.get()), false);

      // Create the coordinate model using the base multivector adapter
      Zoltan2::modelFlag_t flags;

      // Create the MJ object
      Teuchos::RCP<const Comm<int>> comm2 = comm_;
      Teuchos::RCP<mj_t> mj(new mj_t(env_, comm2, baseAdapter));

      // Partition with MJ
      Teuchos::RCP<solution_t> vectorsolution( new solution_t(env_, comm2, 1, mj));
      mj->partition(vectorsolution);

      // Transform the solution
      Teuchos::ArrayRCP<part_t> parts(myNumVertices);
      for(size_t i = 0; i < myNumVertices; i++) parts[i] = vectorsolution->getPartListView()[i];
      solution->setParts(parts);
    }

} // namespace Zoltan2

#endif
