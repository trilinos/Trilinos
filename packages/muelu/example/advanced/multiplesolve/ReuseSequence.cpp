// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <chrono>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_Factory.hpp>

#include <MueLu_RAPFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TransPFactory.hpp>

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
#include <MueLu_ExplicitInstantiation.hpp>
#endif

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>  // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>   // => This header defines Belos::MueLuOp
#endif

#include <MueLu_CreateXpetraPreconditioner.hpp>

#ifdef HAVE_MUELU_PAMGEN
#include "RTC_FunctionRTC.hh"
#include <MueLu_TestHelpers_Common.hpp>
#endif

// This example demonstrates how to reuse some parts of a classical SA multigrid setup between runs.
//
// In this example, we suppose that the pattern of the matrix does not change between runs so that:
// - Aggregates can be reused
// - The tentative prolongator of Smoothed-Aggregation does not change (as it derived directly from the aggregate information).
// - The pattern of coarse grid A can be reused during its computation
//
// The resulting preconditioners are identical to multigrid preconditioners built without recycling the parts described above.
// This can be verified by using the --no-recycling option.

template <class Scalar>
class Tensor {
 private:
  typedef Scalar SC;
  typedef Teuchos::ScalarTraits<SC> STS;

 public:
  Tensor()
    : useSigmaRTC_(false)
    , is3D_(true) {}

#ifdef HAVE_MUELU_PAMGEN
  Tensor(const std::string& rtcString, bool is3D = true)
    : useSigmaRTC_(true)
    , is3D_(is3D) {
    sigmaRTC_ = Teuchos::rcp(new PG_RuntimeCompiler::Function);
    std::string variableType;
    if (TYPE_EQUAL(Scalar, float) || TYPE_EQUAL(Scalar, std::complex<float>))
      variableType = "float";
    else
      variableType = "double";

    if (!sigmaRTC_->addVar(variableType, "x")) throw std::runtime_error("Error setting RTC input argument \"x\"");
    if (!sigmaRTC_->addVar(variableType, "y")) throw std::runtime_error("Error setting RTC input argument \"y\"");
    if (is3D_ &&
        !sigmaRTC_->addVar(variableType, "z")) throw std::runtime_error("Error setting RTC input argument \"z\"");
    if (!sigmaRTC_->addVar(variableType, "t")) throw std::runtime_error("Error setting RTC input argument \"t\"");
    if (!sigmaRTC_->addVar(variableType, "sigmax")) throw std::runtime_error("Error setting RTC input argument \"sigmax\"");
    if (!sigmaRTC_->addVar(variableType, "sigmay")) throw std::runtime_error("Error setting RTC input argument \"sigmay\"");
    if (is3D_ &&
        !sigmaRTC_->addVar(variableType, "sigmaz")) throw std::runtime_error("Error setting RTC input argument \"sigmaz\"");

    if (!sigmaRTC_->addBody(rtcString)) throw std::runtime_error("Error in RTC function compilation");
  }
#endif

  SC operator()(char c, SC x, SC y, SC z = 0.0) const {
#ifdef HAVE_MUELU_PAMGEN
    if (useSigmaRTC_)
      return tensorRTC(c, x, y, z);
#endif
    return tensorDefault(c, x, y, z);
  }

  void setT(double t) { t_ = t; }

  void operator=(const Tensor<Scalar>& tensor) {
    useSigmaRTC_ = tensor.useSigmaRTC_;
    is3D_        = tensor.is3D_;
    t_           = tensor.t_;
#ifdef HAVE_MUELU_PAMGEN
    sigmaRTC_ = tensor.sigmaRTC_;
#endif
  }

 private:
  SC tensorDefault(char c, SC x, SC y, SC z) const {
    // isotropic tensor
    return STS::one();
  }

#ifdef HAVE_MUELU_PAMGEN
  SC tensorRTC(char c, SC x, SC y, SC z) const {
    SC sigmax, sigmay, sigmaz;

    int cnt = 0;
    if (!sigmaRTC_->varValueFill(cnt++, x)) throw std::runtime_error("Could not fill \"x\"");
    if (!sigmaRTC_->varValueFill(cnt++, y)) throw std::runtime_error("Could not fill \"y\"");
    if (is3D_ &&
        !sigmaRTC_->varValueFill(cnt++, z)) throw std::runtime_error("Could not fill \"z\"");
    if (!sigmaRTC_->varValueFill(cnt++, t_)) throw std::runtime_error("Could not fill \"t\"");
    if (!sigmaRTC_->varAddrFill(cnt++, &sigmax)) throw std::runtime_error("Could not fill \"sigmax\"");
    if (!sigmaRTC_->varAddrFill(cnt++, &sigmay)) throw std::runtime_error("Could not fill \"sigmay\"");
    if (is3D_ &&
        !sigmaRTC_->varAddrFill(cnt++, &sigmaz)) throw std::runtime_error("Could not fill \"sigmaz\"");

    sigmaRTC_->execute();

    if (c == 'x') return sigmax;
    if (c == 'y') return sigmay;
    if (is3D_ &&
        c == 'z') return sigmaz;

    throw std::runtime_error("Wrong c");
  }
#endif

 private:
  bool useSigmaRTC_;
  bool is3D_;
  double t_;

#ifdef HAVE_MUELU_PAMGEN
  mutable Teuchos::RCP<PG_RuntimeCompiler::Function> sigmaRTC_;
#endif
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Map, class Matrix, class MultiVector>
Teuchos::RCP<Matrix> BuildMatrix(bool is3D, const Tensor<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& tensor, Teuchos::ParameterList& list,
                                 const Teuchos::RCP<const Map>& map, const Teuchos::RCP<const MultiVector>& coords) {
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal LO;
  typedef Scalar SC;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  typedef typename MultiVector::scalar_type Real;

  GO nx = list.get("nx", (GO)-1);
  GO ny = list.get("ny", (GO)-1);
  GO nz = -1;
  if (is3D) {
    // 3D
    nz = list.get("nz", (GO)-1);

    if (nx == -1 || ny == -1 || nz == -1) {
      GO n = map->getGlobalNumElements();
      nx   = (GO)Teuchos::ScalarTraits<double>::pow(n, 0.33334);
      ny   = nx;
      nz   = nx;
      TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
    }
  } else {
    // 2D
    if (nx == -1 || ny == -1) {
      GO n = map->getGlobalNumElements();
      nx   = (GO)Teuchos::ScalarTraits<double>::pow(n, 0.5);
      ny   = nx;
      TEUCHOS_TEST_FOR_EXCEPTION(nx * ny != n, std::logic_error, "You need to specify nx, ny, and nz");
    }
  }

  double one  = 1.0;
  SC stretchx = list.get("stretchx", one);
  SC stretchy = list.get("stretchy", one);
  SC stretchz = list.get("stretchz", one);

  // bool keepBCs = list.get("keepBCs", false);

  LO nnz = (is3D ? 7 : 5);

  Teuchos::RCP<Matrix> A = Galeri::Xpetra::MatrixTraits<Map, Matrix>::Build(map, nnz);

  LO numMyElements = map->getLocalNumElements();
  GO indexBase     = map->getIndexBase();

  ArrayView<const GO> myGlobalElements = map->getLocalElementList();

  std::vector<GO> inds(nnz);
  std::vector<SC> vals(nnz);

  ArrayRCP<const Real> x = coords->getData(0);
  ArrayRCP<const Real> y = coords->getData(1);
  ArrayRCP<const Real> z = (is3D ? coords->getData(2) : Teuchos::null);

  //    e
  //  b a c
  //    d
  // + f bottom and g top
  GO center, left, right, bottom, top, front, back;
  for (LO i = 0; i < numMyElements; ++i) {
    size_t n = 0;

    center = myGlobalElements[i] - indexBase;

    if (is3D) {
      // 3D
      Galeri::Xpetra::GetNeighboursCartesian3d(center, nx, ny, nz, left, right, front, back, bottom, top);

      SC w;

      w    = tensor('x', x[center], y[center], z[center]);
      SC b = -((left != -1) ? tensor('x', 0.5 * (x[center] + x[left]), 0.5 * (y[center] + y[left]), 0.5 * (z[center] + z[left])) : w) / (stretchx * stretchx);
      SC c = -((right != -1) ? tensor('x', 0.5 * (x[center] + x[right]), 0.5 * (y[center] + y[right]), 0.5 * (z[center] + z[right])) : w) / (stretchx * stretchx);
      w    = tensor('y', x[center], y[center], z[center]);
      SC d = -((front != -1) ? tensor('y', 0.5 * (x[center] + x[front]), 0.5 * (y[center] + y[front]), 0.5 * (z[center] + z[front])) : w) / (stretchy * stretchy);
      SC e = -((back != -1) ? tensor('y', 0.5 * (x[center] + x[back]), 0.5 * (y[center] + y[back]), 0.5 * (z[center] + z[back])) : w) / (stretchy * stretchy);
      w    = tensor('z', x[center], y[center], z[center]);
      SC f = -((bottom != -1) ? tensor('z', 0.5 * (x[center] + x[bottom]), 0.5 * (y[center] + y[bottom]), 0.5 * (z[center] + z[bottom])) : w) / (stretchz * stretchz);
      SC g = -((top != -1) ? tensor('z', 0.5 * (x[center] + x[top]), 0.5 * (y[center] + y[top]), 0.5 * (z[center] + z[top])) : w) / (stretchz * stretchz);

      SC a = -(b + c + d + e + f + g);

      if (left != -1) {
        inds[n]   = left;
        vals[n++] = b;
      }
      if (right != -1) {
        inds[n]   = right;
        vals[n++] = c;
      }
      if (front != -1) {
        inds[n]   = front;
        vals[n++] = d;
      }
      if (back != -1) {
        inds[n]   = back;
        vals[n++] = e;
      }
      if (bottom != -1) {
        inds[n]   = bottom;
        vals[n++] = f;
      }
      if (top != -1) {
        inds[n]   = top;
        vals[n++] = g;
      }

      if (bottom != -1 && Galeri::Xpetra::IsBoundary3d(center, nx, ny, nz)) {
        // Neumann boundary unknown (diagonal = sum of all offdiagonal)
        a = Teuchos::ScalarTraits<Scalar>::zero();
        for (size_t j = 0; j < n; j++)
          a -= vals[j];
      }
      inds[n]   = center;
      vals[n++] = a;

    } else {
      // 2D
      Galeri::Xpetra::GetNeighboursCartesian2d(center, nx, ny, left, right, front, back);

      SC w;

      w    = tensor('x', x[center], y[center]);
      SC b = -((left != -1) ? tensor('x', 0.5 * (x[center] + x[left]), 0.5 * (y[center] + y[left])) : w) / (stretchx * stretchx);
      SC c = -((right != -1) ? tensor('x', 0.5 * (x[center] + x[right]), 0.5 * (y[center] + y[right])) : w) / (stretchx * stretchx);
      w    = tensor('y', x[center], y[center]);
      SC d = -((front != -1) ? tensor('y', 0.5 * (x[center] + x[front]), 0.5 * (y[center] + y[front])) : w) / (stretchy * stretchy);
      SC e = -((back != -1) ? tensor('y', 0.5 * (x[center] + x[back]), 0.5 * (y[center] + y[back])) : w) / (stretchy * stretchy);

      SC a = -(b + c + d + e);

      if (left != -1) {
        inds[n]   = left;
        vals[n++] = b;
      }
      if (right != -1) {
        inds[n]   = right;
        vals[n++] = c;
      }
      if (front != -1) {
        inds[n]   = front;
        vals[n++] = d;
      }
      if (back != -1) {
        inds[n]   = back;
        vals[n++] = e;
      }

      if (front != -1 && Galeri::Xpetra::IsBoundary2d(center, nx, ny)) {
        // Neumann boundary unknown (diagonal = sum of all offdiagonal)
        a = Teuchos::ScalarTraits<Scalar>::zero();
        for (size_t j = 0; j < n; j++)
          a -= vals[j];
      }
      inds[n]   = center;
      vals[n++] = a;
    }

    for (size_t j = 0; j < n; j++)
      inds[j] += indexBase;

    ArrayView<GO> iv(&inds[0], n);
    ArrayView<SC> av(&vals[0], n);
    A->insertGlobalValues(myGlobalElements[i], iv, av);
  }

  A->fillComplete();

  return A;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ConstructData(bool is3D, const Tensor<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& tensor, const std::string& matrixType, Teuchos::ParameterList& galeriList,
                   Xpetra::UnderlyingLib lib, Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                   Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                   Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
                   Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>>& coordinates,
                   Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& nullspace) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef typename Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  if (is3D) {
    // 3D
    map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian3D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<real_type, LO, GO, Map, RealValuedMultiVector>("3D", map, galeriList);

  } else {
    // 2D
    map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian2D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<real_type, LO, GO, Map, RealValuedMultiVector>("2D", map, galeriList);
  }

  A = BuildMatrix<SC, LO, GO, Map, CrsMatrixWrap, RealValuedMultiVector>(is3D, tensor, galeriList, map, coordinates);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor& clp, Xpetra::UnderlyingLib& lib, int argc, char* argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  using namespace std::chrono;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

  // =========================================================================
  // Convenient definitions
  // =========================================================================
  typedef Teuchos::ScalarTraits<SC> STS;
  SC one = STS::one(), zero = STS::zero();

  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out       = *fancy;
  out.setOutputToRootOnly(0);

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  GO nx = 20, ny = 20, nz = 20;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace3D");  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                                       // manage parameters of Xpetra

  std::string xmlFileName = "reuse_seq.xml";
  clp.setOption("xml", &xmlFileName, "read parameters from a file");
  std::string solveType = "cg";
  clp.setOption("solver", &solveType, "solve type: (none | cg | standalone)");
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType tol = 1e-6;
  clp.setOption("tol", &tol, "solver convergence tolerance");
  int maxIts = 200;
  clp.setOption("its", &maxIts, "maximum number of solver iterations");
  int dim = 3;
  clp.setOption("dim", &dim, "space dimension");

  clp.recogniseAllOptions(true);
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  bool is3D = (dim == 3);

  // Retrieve matrix parameters (they may have been changed on the command line)
  // [for instance, if we changed matrix type from 2D to 3D we need to update nz]
  Teuchos::ParameterList galeriList = galeriParameters.GetParameterList();

  // =========================================================================
  // Problem construction
  // =========================================================================
  // For comments, see Driver.cpp
  out << "========================================================\n"
      << xpetraParameters << galeriParameters;
  std::string matrixType = galeriParameters.GetMatrixType();

  out << "Processor subdomains in x direction: " << galeriList.get<GO>("mx") << std::endl
      << "Processor subdomains in y direction: " << galeriList.get<GO>("my") << std::endl
      << "Processor subdomains in z direction: " << galeriList.get<GO>("mz") << std::endl
      << "========================================================" << std::endl;

  // =========================================================================
  // Setups and solves
  // =========================================================================
  std::string thickSeparator = "=============================================================";
  std::string thinSeparator  = "-------------------------------------------------------------";

  Teuchos::ParameterList paramList;
  paramList.set("verbosity", "none");
  if (xmlFileName != "")
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *comm);

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef typename Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;
  Tensor<real_type> tensor;
  if (paramList.isParameter("sigma")) {
    std::string sigmaString = paramList.get<std::string>("sigma");
    paramList.remove("sigma");
#ifdef HAVE_MUELU_PAMGEN
    out << "Switching to RTC" << std::endl;
    tensor = Tensor<real_type>(sigmaString, is3D);
#else
    (void)sigmaString;  // fix compiler warning
#endif
  }

  out << "Parameter list:" << std::endl
      << paramList << std::endl;

  // =========================================================================
  // The LOOP
  // =========================================================================
  std::vector<std::string> reuseTypes, reuseNames;
  reuseTypes.push_back("none");
  reuseNames.push_back("none");
  reuseTypes.push_back("S");
  reuseNames.push_back("smoothers");
  reuseTypes.push_back("tP");
  reuseNames.push_back("tentative P");
  reuseTypes.push_back("RP");
  reuseNames.push_back("smoothed P and R");
  reuseTypes.push_back("RAP");
  reuseNames.push_back("coarse grids");

  const size_t numSteps = 8;

  high_resolution_clock::time_point tc;
  std::vector<duration<double>> setup_time(reuseTypes.size() * numSteps);
  std::vector<duration<double>> solve_time(reuseTypes.size() * numSteps);
  std::vector<int> num_its(reuseTypes.size() * numSteps);

  for (size_t k = 0; k < reuseTypes.size(); k++) {
    out << thickSeparator << " " << reuseTypes[k] << " " << thickSeparator << std::endl;

    paramList.set("reuse: type", reuseTypes[k]);

    RCP<Matrix> A;
    RCP<const Map> map;
    RCP<RealValuedMultiVector> coordinates;
    RCP<MultiVector> nullspace;

    tensor.setT(0);

    ConstructData(is3D, tensor, matrixType, galeriList, lib, comm, A, map, coordinates, nullspace);
    A->SetMaxEigenvalueEstimate(-one);

    RCP<Vector> X = VectorFactory::Build(map);
    X->setSeed(846930886);
    X->randomize();
    RCP<Vector> B = VectorFactory::Build(map);
    A->apply(*X, *B);

    Teuchos::ParameterList userParamList = paramList.sublist("user data");
    userParamList.set<RCP<RealValuedMultiVector>>("Coordinates", coordinates);
    RCP<Hierarchy> H = MueLu::CreateXpetraPreconditioner(A, paramList);

    for (size_t t = 1; t < numSteps; t++) {
      out << thinSeparator << " Step " << t << " " << thinSeparator << std::endl;

      tensor.setT(t);

      ConstructData(is3D, tensor, matrixType, galeriList, lib, comm, A, map, coordinates, nullspace);
      A->SetMaxEigenvalueEstimate(-one);

      tc = high_resolution_clock::now();
      if (solveType == "none")
        H = MueLu::CreateXpetraPreconditioner(A, paramList);
      else
        MueLu::ReuseXpetraPreconditioner(A, H);
      setup_time[k * numSteps + t] = duration_cast<duration<double>>(high_resolution_clock::now() - tc);

      X->putScalar(zero);

      if (solveType == "none") {
        // Do nothing

      } else if (solveType == "standalone") {
        H->IsPreconditioner(false);

        tc = high_resolution_clock::now();
        H->Iterate(*B, *X, tol);
        solve_time[k * numSteps + t] = duration_cast<duration<double>>(high_resolution_clock::now() - tc);

      } else if (solveType == "cg" || solveType == "gmres") {
        H->IsPreconditioner(true);

#ifdef HAVE_MUELU_BELOS
        // Operator and Multivector type that will be used with Belos
        typedef MultiVector MV;
        typedef Belos::OperatorT<MV> OP;

        // Define Operator and Preconditioner
        Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(A));  // Turns a Xpetra::Matrix object into a Belos operator
        Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp<SC, LO, GO, NO>(H));   // Turns a MueLu::Hierarchy object into a Belos operator

        // Construct a Belos LinearProblem object
        RCP<Belos::LinearProblem<SC, MV, OP>> belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
        belosProblem->setRightPrec(belosPrec);

        bool set = belosProblem->setProblem();
        if (set == false) {
          out << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
          return EXIT_FAILURE;
        }

        // Belos parameter list
        Teuchos::ParameterList belosList;
        belosList.set("Maximum Iterations", maxIts);  // Maximum number of iterations allowed
        belosList.set("Convergence Tolerance", tol);  // Relative convergence tolerance requested
        belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
        belosList.set("Output Frequency", 1);
        belosList.set("Output Style", Belos::Brief);

        // Create an iterative solver manager
        RCP<Belos::SolverManager<SC, MV, OP>> solver;
        if (solveType == "cg") {
          solver = rcp(new Belos::PseudoBlockCGSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));
        } else if (solveType == "gmres") {
          solver = rcp(new Belos::BlockGmresSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));
        }

        // Perform solve
        Belos::ReturnType ret = Belos::Unconverged;

        tc                           = high_resolution_clock::now();
        ret                          = solver->solve();
        solve_time[k * numSteps + t] = duration_cast<duration<double>>(high_resolution_clock::now() - tc);

        // Get the number of iterations for this solve.
        out << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;
        // Check convergence
        if (ret != Belos::Converged) {
          out << std::endl
              << "ERROR:  Belos did not converge! " << std::endl;
          num_its[k * numSteps + t] = -1;

        } else {
          out << std::endl
              << "SUCCESS:  Belos converged!" << std::endl;
          num_its[k * numSteps + t] = solver->getNumIters();
        }
#endif  // ifdef HAVE_MUELU_BELOS
      } else {
        throw MueLu::Exceptions::RuntimeError("Unknown solver type: \"" + solveType + "\"");
      }

      out << "residual(A) = " << Utilities::ResidualNorm(*A, *X, *B)[0] << " [reuse \"" << reuseNames[k] << "\"]" << std::endl;
    }
  }
  for (size_t t = 1; t < numSteps; t++) {
    out << thinSeparator << std::endl;
    for (size_t k = 0; k < reuseTypes.size(); k++)
      printf("step #%d reuse \"%20s\": setup = %5.2e, solve = %5.2e [%3d], total = %5.2e\n", static_cast<int>(t), reuseNames[k].c_str(),
             setup_time[k * numSteps + t].count(), solve_time[k * numSteps + t].count(), num_its[k * numSteps + t],
             setup_time[k * numSteps + t].count() + solve_time[k * numSteps + t].count());
  }

  return EXIT_SUCCESS;
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char* argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
