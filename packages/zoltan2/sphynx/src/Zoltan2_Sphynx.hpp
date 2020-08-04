// @HEADER
//
// ***********************************************************************
//
//                            Sphynx
//           Copyright 2020 National Technology & Engineering
//                  Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Seher Acer        (sacer@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//                    Karen Devine      (kddevin@sandia.gov)
//
// ***********************************************************************
//
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
// Step3: Sphynx calls the MJ algorithm provided in Zoltan2Core to compute the 
//        partition. 
////////////////////////////////////////////////////////////////////////////////

#include "Zoltan2Sphynx_config.h"

#include "Zoltan2_XpetraCrsGraphAdapter.hpp"
#include "Zoltan2_XpetraMultiVectorAdapter.hpp"

#include "Zoltan2_CoordinateModel.hpp"
#include "Zoltan2_AlgMultiJagged.hpp"

#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziTpetraAdapter.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosTpetraOperator.hpp"

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
    
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////// CONSTRUCTORS ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    // Takes the graph from the input adapter and computes the Laplacian matrix
    Sphynx(const RCP<const Environment> &env,
	      const RCP<Teuchos::ParameterList> &params,
	      const RCP<const Comm<int> > &comm,
	      const RCP<const XpetraCrsGraphAdapter<graph_t> > &adapter) :
      env_(env), params_(params), comm_(comm), adapter_(adapter)
    { 

      // Don't compute the Laplacian if the number of requested parts is 1
      const Teuchos::ParameterEntry *pe = params_->getEntryPtr("num_global_parts");
      numGlobalParts_ = pe->getValue<int>(&numGlobalParts_);
      if(numGlobalParts_ > 1){

	Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Sphynx::Laplacian"));

	// The verbosity is common throughout the algorithm, better to check and set now
	pe = params_->getEntryPtr("sphynx_verbosity");
	if (pe)
	  verbosity_ = pe->getValue<int>(&verbosity_);

	// Do we need to pre-process the input?
	pe = params_->getEntryPtr("sphynx_skip_preprocessing");
	if (pe)
	  skipPrep_ = pe->getValue<bool>(&skipPrep_);

	// Get the graph from XpetraCrsGraphAdapter if skipPrep_ is true
	// We assume the graph has all of the symmetric and diagonal entries
	if(skipPrep_) 
	  graph_ = adapter_->getUserGraph();
	else {
	  throw std::runtime_error("\nSphynx Error: Preprocessing has not been implemented yet.\n");
	} 
	
	// Check if the graph is regular and set the problem type accordingly
	determineRegularity();

	// Compute the Laplacian matrix
	computeLaplacian();

	if(problemType_ == GENERALIZED)
	  computeDegreeMatrix();

      }
    }

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////// FORWARD DECLARATIONS  ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    void partition(const RCP<PartitioningSolution<Adapter> > &solution);

 
    int LOBPCGwrapper(const int numEigenVectors);

    template<typename problem_t>
    void setPreconditioner(Teuchos::RCP<problem_t> &problem);


    void eigenvecsToCoords(Teuchos::RCP<mvector_t> &eigenVectors,
			   int computedNumEv,
			   Teuchos::RCP<mvector_t> &coordinates);
    

    void computeWeights(std::vector<const weight_t *> vecweights,
			std::vector<int> strides);
		       
  
    void MJwrapper(const Teuchos::RCP<const mvector_t> &coordinates,
		   std::vector<const weight_t *> weights,
		   std::vector<int> strides,
		   const Teuchos::RCP<PartitioningSolution<Adapter>> &solution);


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
      auto rowOffsets = graph_->getLocalGraph().row_map;
      auto rowOffsets_h = Kokkos::create_mirror_view(rowOffsets);
      Kokkos::deep_copy(rowOffsets_h, rowOffsets);

      // Get size information 
      const size_t numGlobalEntries = graph_->getGlobalNumEntries();
      const size_t numLocalRows = graph_->getNodeNumRows();
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

      // Get the user preference on which laplacian to use
      // The default is the combinatorial Laplacian
      const Teuchos::ParameterEntry *pe = params_->getEntryPtr("sphynx_problem_type");
      std::string probType = "";
      bool userSetProbType = false;
      if (pe)
	probType = pe->getValue<std::string>(&probType);

      if(probType == "combinatorial"){
	problemType_ = COMBINATORIAL;
	userSetProbType = true;
      }
      else if(probType == "generalized"){
	problemType_ = GENERALIZED;
	userSetProbType = true;
      }
      else if(probType == "normalized"){
	problemType_ = NORMALIZED;
	userSetProbType = true;
      }

      
      if(!userSetProbType){
	if(irregular_)
	  problemType_ = GENERALIZED;
	else
	  problemType_ = COMBINATORIAL;
      }

      // Let the user know about what we determined if verbose  
      if(verbosity_ > 0) {
	if(comm_->getRank() == 0) {
	  std::cout << "Determining Regularity --  max degree: " << globalMax 
		    << ", avg degree: " << avg << ", max/avg: " << globalMax/avg << std::endl
		    << "Determined Regularity --  regular: " << !irregular_ << std::endl; 
	  
	}
      }
    }
    

    ///////////////////////////////////////////////////////////////////////////
    // Compute the Laplacian matrix depending on the eigenvalue problem type.
    // There are 3 options for the type: combinatorial, generalized, and normalized.
    // Combinatorial and generalized share the same Laplacian but generalized 
    // also needs a degree matrix.
    void computeLaplacian()
    {
  
      if(problemType_ == NORMALIZED)
	laplacian_ = computeNormalizedLaplacian();
      else 
	laplacian_ = computeCombinatorialLaplacian();
    }

    ///////////////////////////////////////////////////////////////////////////
    // Compute a diagonal matrix with the vertex degrees in the input graph
    void computeDegreeMatrix()
    {

	// Get the row pointers in the host
	auto rowOffsets = graph_->getLocalGraph().row_map;
	auto rowOffsets_h = Kokkos::create_mirror_view(rowOffsets);
	Kokkos::deep_copy(rowOffsets_h, rowOffsets);

	// Create the degree matrix with max row size set to 1
 	Teuchos::RCP<matrix_t> degMat(new matrix_t (graph_->getRowMap(), 
						    graph_->getRowMap(), 
						    1, Tpetra::StaticProfile));

	scalar_t *val = new scalar_t[1];
	lno_t *ind = new lno_t[1];
	lno_t numRows = static_cast<lno_t>(graph_->getNodeNumRows());

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
  
      const size_t numEnt = graph_->getNodeNumEntries();
      const size_t numRows = graph_->getNodeNumRows();

      // Create new values for the Laplacian, initialize to -1 
      values_view_t newVal (Kokkos::view_alloc("CombLapl::val", Kokkos::WithoutInitializing), numEnt);
      Kokkos::deep_copy(newVal, -1);

      // Get the diagonal offsets
      offset_view_t diagOffsets(Kokkos::view_alloc("Diag Offsets", Kokkos::WithoutInitializing), numRows);
      graph_->getLocalDiagOffsets(diagOffsets);

      // Get the row pointers in the host
      auto rowOffsets = graph_->getLocalGraph().row_map;

      // Compute the diagonal entries as the vertex degrees
      Kokkos::parallel_for("Combinatorial Laplacian", range_policy(0, numRows),
			   KOKKOS_LAMBDA(const lno_t i){
			     newVal(rowOffsets(i) + diagOffsets(i)) = rowOffsets(i+1) - rowOffsets(i) - 1;
			   }
			   );
      Kokkos::fence ();

      // Create the Laplacian maatrix using the input graph and with the new values
      Teuchos::RCP<matrix_t> laplacian (new matrix_t(graph_, newVal));
      laplacian->fillComplete (graph_->getDomainMap(), graph_->getRangeMap());

      return laplacian;

    }


    ///////////////////////////////////////////////////////////////////////////
    // Compute the normalized Laplacian: L_N = D^{-1/2} L D^{-1/2}, where L = D - A.
    // l_ij = 1
    // l_ij = -1/(sqrt(deg(v_i))sqrt(deg(v_j)) if i != j and a_ij != 0
    // l_ij = 0 if i != j and a_ij = 0 
    Teuchos::RCP<matrix_t> computeNormalizedLaplacian()
    {
      using range_policy = Kokkos::RangePolicy<
	typename node_t::device_type::execution_space, Kokkos::IndexType<lno_t>>;
      using values_view_t = Kokkos::View<scalar_t*, typename node_t::device_type>;
      using offset_view_t = Kokkos::View<size_t*, typename node_t::device_type>;
      using vector_t = Tpetra::Vector<scalar_t, lno_t, gno_t, node_t>;
      using dual_view_t = typename vector_t::dual_view_type;
      using KAT = Kokkos::Details::ArithTraits<scalar_t>;

      const size_t numEnt = graph_->getNodeNumEntries();
      const size_t numRows = graph_->getNodeNumRows();

      // Create new values for the Laplacian, initialize to -1 
      values_view_t newVal (Kokkos::view_alloc("NormLapl::val", Kokkos::WithoutInitializing), numEnt);
      Kokkos::deep_copy(newVal, -1);

      // D^{-1/2}
      dual_view_t dv ("MV::DualView", numRows, 1);
      auto deginvsqrt = dv.d_view;

      // Get the diagonal offsets
      offset_view_t diagOffsets(Kokkos::view_alloc("Diag Offsets", Kokkos::WithoutInitializing), numRows);
      graph_->getLocalDiagOffsets(diagOffsets);

      // Get the row pointers
      auto rowOffsets = graph_->getLocalGraph().row_map;

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
    const Teuchos::RCP<const Teuchos::Comm<int> > comm_;
    const Teuchos::RCP<const Adapter> adapter_;
    
    // Internal members
    int numGlobalParts_; 
    Teuchos::RCP<const graph_t> graph_;
    Teuchos::RCP<matrix_t> laplacian_;
    Teuchos::RCP<const matrix_t> degMatrix_;
    Teuchos::RCP<mvector_t> eigenVectors_;
    problemType problemType_;  // obtained from user params or decided internally
    bool irregular_;           // decided internally
    int verbosity_;            // obtained from user params
    bool skipPrep_;            // obtained from user params
  };



  ///////////////////////////////////////////////////////////////////////////
  /////////////////////// MORE MEMBER FUNCTIONS  ////////////////////////////
  ///////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////
  // Compute a partition using the Laplacian matrix (and possibly the degree
  // matrix as well). First call LOBPCG to compute logK+1 eigenvectors, then 
  // transform the eigenvectors to coordinates, and finally call MJ to compute
  // a partition on the coordinates.
  template <typename Adapter>
  void Sphynx<Adapter>::partition(const RCP<PartitioningSolution<Adapter>> &solution)
  {
    // Return a trivial solution if only one part is requested
    if(numGlobalParts_ == 1) {

      size_t numRows =adapter_->getUserGraph()->getNodeNumRows();
      Teuchos::ArrayRCP<part_t> parts(numRows,0);
      solution->setParts(parts);
      
      return;

    }

    // The number of eigenvectors to be computed
    int numEigenVectors = (int) log2(numGlobalParts_)+1;

    // Compute the eigenvectors using LOBPCG
    int computedNumEv = Sphynx::LOBPCGwrapper(numEigenVectors);

    if(computedNumEv <= 1) {
      throw 
	std::runtime_error("\nAnasazi Error: LOBPCGSolMgr::solve() returned unconverged.\n"
			   "Sphynx Error:  LOBPCG could not compute any eigenvectors.\n"
			   "               Increase either max iters or tolerance.\n");
    
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
  template <typename Adapter>
  int Sphynx<Adapter>::LOBPCGwrapper(const int numEigenVectors)
  {

    Teuchos::TimeMonitor t(*Teuchos::TimeMonitor::getNewTimer("Sphynx::LOBPCG"));

    // Set defaults for the parameters
    std::string which = "SR";
    std::string ortho = "SVBQ";
    double tolerance = 1.0e-2;
    bool relConvTol = false;
    int maxIterations = 1000;
    bool isHermitian = true;
    bool relLockTol = false;
    bool lock = false;

    // Information to output in a verbose run
    int numfailed = 0;
    int iter = 0;
    double solvetime = 0;
  
    // Get the user values for the parameters
    const Teuchos::ParameterEntry *pe;

    pe = params_->getEntryPtr("sphynx_tolerance");
    if (pe)
      tolerance = pe->getValue<scalar_t>(&tolerance);

    pe = params_->getEntryPtr("sphynx_maxiterations");
    if (pe)
      maxIterations = pe->getValue<int>(&maxIterations);

  
    // Set Anasazi verbosity level
    int anasaziVerbosity = Anasazi::Errors + Anasazi::Warnings;
    if (verbosity_ >= 1)  // low
      anasaziVerbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
    if (verbosity_ >= 2)  // medium
      anasaziVerbosity += Anasazi::IterationDetails;
    if (verbosity_ >= 3)  // high
      anasaziVerbosity += Anasazi::StatusTestDetails
	+ Anasazi::OrthoDetails
	+ Anasazi::Debug;

  
    // Create the parameter list to pass into solver
    Teuchos::ParameterList anasaziParams;
    anasaziParams.set("Verbosity", anasaziVerbosity);
    anasaziParams.set("Which", which);
    anasaziParams.set("Convergence Tolerance", tolerance);
    anasaziParams.set("Maximum Iterations", maxIterations);
    anasaziParams.set("Block Size", numEigenVectors);
    anasaziParams.set("Relative Convergence Tolerance", relConvTol);
    anasaziParams.set("Orthogonalization", ortho);
    anasaziParams.set("Use Locking", lock);
    anasaziParams.set("Relative Locking Tolerance", relLockTol);


    // Create and set initial vectors
    RCP<mvector_t> ivec( new mvector_t(laplacian_->getRangeMap(), numEigenVectors));
    Anasazi::MultiVecTraits<scalar_t, mvector_t>::MvRandom(*ivec);
    for (size_t i = 0; i < ivec->getLocalLength(); i++)
      ivec->replaceLocalValue(i,0,1.);
    

    // Create the eigenproblem to be solved
    using problem_t = Anasazi::BasicEigenproblem<scalar_t, mvector_t, op_t>;
    Teuchos::RCP<problem_t> problem (new problem_t(laplacian_, ivec));
    problem->setHermitian(isHermitian);
    problem->setNEV(numEigenVectors);


    // Set preconditioner
    Sphynx::setPreconditioner(problem);

    if(problemType_ == Sphynx::GENERALIZED)
      problem->setM(degMatrix_);

    // Inform the eigenproblem that you are finished passing it information
    bool boolret = problem->setProblem();
    if (boolret != true) {
      throw std::runtime_error("\nAnasazi::BasicEigenproblem::setProblem() returned with error.\n");
    }

    // Set LOBPCG
    using solver_t = Anasazi::LOBPCGSolMgr<scalar_t, mvector_t, op_t>;
    solver_t solver(problem, anasaziParams);

    if (verbosity_ > 0 && comm_->getRank() == 0) 
      anasaziParams.print(std::cout);
  
    // Solve the problem
    if (verbosity_ > 0 && comm_->getRank() == 0) 
      std::cout << "Beginning the LOBPCG solve..." << std::endl;
    Anasazi::ReturnType returnCode = solver.solve();

    // Check convergence, niters, and solvetime
    if (returnCode != Anasazi::Converged) {
      ++numfailed;
    }
    iter = solver.getNumIters();
    solvetime = (solver.getTimers()[0])->totalElapsedTime();
  
 
    // Retrieve the solution
    using solution_t = Anasazi::Eigensolution<scalar_t, mvector_t>;
    solution_t sol = problem->getSolution();
    eigenVectors_ = sol.Evecs;
    int numev = sol.numVecs;

    // Summarize iteration counts and solve time
    if (verbosity_ > 0 && comm_->getRank() == 0) {
      std::cout << std::endl;
      std::cout << "LOBPCG SUMMARY" << std::endl;
      std::cout << "Failed to converge:    " << numfailed << std::endl;
      std::cout << "No of iterations :     " << iter << std::endl;
      std::cout << "Solve time:            " << solvetime << std::endl;
      std::cout << "No of comp. vecs. :    " << numev << std::endl;
    }
  
    return numev;
  
  }

  ///////////////////////////////////////////////////////////////////////////
  // Determine which preconditioner to use and set it in the given problem.
  // There are two options: polynomial preconditioner (from Belos) and Muelu.
  // Since MueLu is an optional dependency, we use polynomial when MueLu is 
  // not enabled. When MueLu is enabled, using MueLu is the default setting
  // but the user can set otherwise.
  template <typename Adapter>
  template <typename problem_t>
  void Sphynx<Adapter>::setPreconditioner(Teuchos::RCP<problem_t> &problem)
  {
    bool hasSet = false;
    bool usePoly = false;
  
    //Get the user values for the parameters
    const Teuchos::ParameterEntry *pe = params_->getEntryPtr("sphynx_preconditioner_poly");
    if (pe)
      usePoly = pe->getValue<bool>(&usePoly);

  
#ifdef HAVE_ZOLTAN2SPHYNX_MUELU

    if(!usePoly) {

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
      smootherParamList.set("chebyshev: eigenvalue max iterations", 100);
      paramList.set("smoother: params", smootherParamList);
      paramList.set("use kokkos refactor", false); 


      if(irregular_) {
	
	paramList.set("multigrid algorithm", "unsmoothed");

	paramList.set("coarse: type", "CHEBYSHEV");
	Teuchos::ParameterList coarseParamList;
	coarseParamList.set("chebyshev: degree", 3);
	coarseParamList.set("chebyshev: ratio eigenvalue", 7.0);
	coarseParamList.set("chebyshev: eigenvalue max iterations", 100);
	paramList.set("coarse: params", coarseParamList);

	paramList.set("max levels", 5);
	paramList.set("aggregation: drop tol", 0.40);
      
      }

      using prec_t = MueLu::TpetraOperator<scalar_t, lno_t, gno_t, node_t>;
      Teuchos::RCP<prec_t> prec = MueLu::CreateTpetraPreconditioner<
	scalar_t, lno_t, gno_t, node_t>(laplacian_, paramList);
  
      problem->setPrec(prec);
  
      hasSet = true;
    }
#endif

    if(!hasSet) {

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
      paramList.set("Maximum Degree", 25);
      paramList.set("Polynomial Tolerance", 1.0e-6 );
      paramList.set("Verbosity", verbosity2 );
      paramList.set("Random RHS", false );
      paramList.set("Outer Solver", ""); 
      paramList.set("Timer Label", "Belos Polynomial Solve" );

      // Construct a linear problem for the polynomial solver manager
      using lproblem_t = Belos::LinearProblem<scalar_t, mvector_t, op_t>;
      Teuchos::RCP<lproblem_t> innerPolyProblem(new lproblem_t());
      innerPolyProblem->setOperator(laplacian_);

      using btop_t = Belos::TpetraOperator<scalar_t>;
      Teuchos::RCP<btop_t> polySolver(new btop_t(innerPolyProblem, 
						 Teuchos::rcpFromRef(paramList), 
						 "GmresPoly", true));
      problem->setPrec(polySolver);
    }
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
    using cmodel_t = Zoltan2::CoordinateModel<base_adapter_t>;
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
    Teuchos::RCP<const cmodel_t> coordModel (new cmodel_t(baseAdapter, env_, comm_, flags));

    // Create the MJ object
    Teuchos::RCP<const Comm<int>> comm2 = comm_;
    Teuchos::RCP<mj_t> mj(new mj_t(env_, comm2, coordModel));

    // Partition with MJ 
    Teuchos::RCP<solution_t> vectorsolution( new solution_t(env_, comm2, 1, mj));
    mj->partition(vectorsolution);

    // Transform the solution
    Teuchos::ArrayRCP<part_t> parts(myNumVertices);
    for(size_t i = 0; i < myNumVertices; i++) {
      parts[i] = vectorsolution->getPartListView()[i];
    }
    solution->setParts(parts);

  }

} // namespace Zoltan2

#endif
