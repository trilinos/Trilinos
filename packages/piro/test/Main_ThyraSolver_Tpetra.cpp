// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include <iostream>
#include <string>

#include "MockModelEval_A_Tpetra.hpp"
#include "Piro_PerformSolve.hpp"
#include "Piro_SolverFactory.hpp"
#include "Piro_StratimikosUtils.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Tpetra_Core.hpp"

#ifdef HAVE_PIRO_IFPACK2
#include "Thyra_Ifpack2PreconditionerFactory.hpp"
#endif

#ifdef HAVE_PIRO_MUELU
#include "Stratimikos_MueLuHelpers.hpp"
#endif

#include "Piro_ConfigDefs.hpp"

int main(int argc, char *argv[]) {
    int status = 0;          // 0 = pass, failures are incremented
    int overall_status = 0;  // 0 = pass, failures are incremented over multiple tests
    bool success = true;

    // Initialize MPI
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);
    int Proc = mpiSession.getRank();

    Teuchos::oblackholestream blackHole;
    std::ostream &out = (Proc == 0) ? std::cout : blackHole;

    auto appComm = Tpetra::getDefaultComm();

    using Teuchos::RCP;
    using Teuchos::rcp;

    std::string inputFile;
    bool doAll = (argc == 1);
    if (argc > 1) doAll = !strcmp(argv[1], "-v");

    Piro::SolverFactory solverFactory;

    int numTests = 1;

    for (int iTest = 0; iTest < numTests; iTest++) {
        if (doAll) {
            switch (iTest) {
                case 0:
                    inputFile = "input_Solve_NOX_Thyra.xml";
                    break;
                default:
                    std::cout << "iTest logic error " << std::endl;
                    exit(-1);
            }
        } else {
            inputFile = argv[1];
            iTest = 999;
        }

        if (Proc == 0)
            std::cout << "===================================================\n"
                      << "======  Running input file " << iTest << ": " << inputFile << "\n"
                      << "===================================================\n"
                      << std::endl;

        try {
            // Create (1) a Model Evaluator and (2) a ParameterList
            const RCP<Thyra::ModelEvaluator<double>> Model = rcp(new MockModelEval_A_Tpetra(appComm));

            RCP<Teuchos::ParameterList> piroParams =
                rcp(new Teuchos::ParameterList("Piro Parameters"));
            Teuchos::updateParametersFromXmlFile(inputFile, piroParams.ptr());

            // Use these two objects to construct a Piro solved application
            RCP<const Thyra::ResponseOnlyModelEvaluatorBase<double>> piro;
            {
                const RCP<Teuchos::ParameterList> stratParams = Piro::extractStratimikosParams(piroParams);

                Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
#ifdef HAVE_PIRO_IFPACK2
                typedef Thyra::PreconditionerFactoryBase<double> Base;
                typedef Thyra::Ifpack2PreconditionerFactory<Tpetra_CrsMatrix> Impl;
                linearSolverBuilder.setPreconditioningStrategyFactory(
                    Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
#endif
                linearSolverBuilder.setParameterList(stratParams);

                RCP<Thyra::LinearOpWithSolveFactoryBase<double>> lowsFactory =
                    createLinearSolveStrategy(linearSolverBuilder);

                RCP<Thyra::ModelEvaluator<double>> thyraModel = rcp(new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(
                    Model, lowsFactory));

                piro = solverFactory.createSolver(piroParams, thyraModel);
            }

            const Teuchos::RCP<Teuchos::ParameterList> solveParams =
                Teuchos::sublist(Teuchos::sublist(piroParams, "Analysis"), "Solve");

            Teuchos::Array<RCP<Thyra::VectorBase<double>>> responses;
            Teuchos::Array<Teuchos::Array<RCP<Thyra::MultiVectorBase<double>>>> sensitivities;

            Teuchos::Array<RCP<Thyra::MultiVectorBase<double>>> directions;
            Teuchos::Array<Teuchos::Array<RCP<Thyra::MultiVectorBase<double>>>> reducedHessian;

            const int seed = 42;

            const int np = piro->Np();

            ::Thyra::seed_randomize<double>(seed);
            directions.resize(np);
            int n_directions = 2;
            for (int l = 0; l < np; l++) {
                auto p_space = piro->getNominalValues().get_p(l)->space();
                directions[l] = Thyra::createMembers(p_space, n_directions);
                for (int i_direction = 0; i_direction < n_directions; i_direction++) {
                    ::Thyra::put_scalar(0.0, directions[l]->col(i_direction).ptr());
                    ::Thyra::set_ele(i_direction, 1.0, directions[l]->col(i_direction).ptr());
                }
            }
            Piro::PerformSolve(*piro, *solveParams, responses, sensitivities, directions, reducedHessian);

            // Extract default input parameters
            const RCP<const Thyra::VectorBase<double>> p1 = piro->getNominalValues().get_p(0);

            // Extract output arguments
            const RCP<Thyra::VectorBase<double>> g1 = responses[0];
            const RCP<Thyra::VectorBase<double>> gx = responses[1];
            const RCP<Thyra::MultiVectorBase<double>> dgdp = sensitivities[0][0];

            const RCP<Thyra::MultiVectorBase<double>> hv = reducedHessian[0][0];

            double tol = 1e-5;

            // Print out everything
            out << "Finished Model Evaluation: Printing everything {Exact in brackets}"
                << "\n-----------------------------------------------------------------"
                << std::setprecision(9) << std::endl;

            if (Teuchos::nonnull(p1)) {
                out << "\nParameters! {1,1}\n"
                    << *p1 << std::endl;
                Thyra::DetachedVectorView<double> p_view(p1->clone_v());
                double p_exact[2] = {1, 1};

                double l2_diff = std::sqrt(std::pow(p_view(0) - p_exact[0], 2) + std::pow(p_view(1) - p_exact[1], 2));
                if (l2_diff > tol) {
                    status += 100;
                    out << "\nPiro_AnalysisDrvier:  Expected parameter values are: {"
                        << p_exact[0] << ", " << p_exact[1] << "}, but the values are: {"
                        << p_view(0) << ", " << p_view(1) << "}.\n"
                        << "Difference in l2 norm: " << l2_diff << " > tol: " << tol << std::endl;
                }
            }
            if (Teuchos::nonnull(g1)) {
                out << "\nResponses! {8.0}\n"
                    << *g1 << std::endl;
                Thyra::DetachedVectorView<double> g_view(g1);
                double g_computed = g_view(0);
                double g_exact = 8.0;

                double diff = std::abs(g_exact - g_computed);
                if (diff > tol) {
                    status += 100;
                    out << "\nPiro_AnalysisDrvier:  Responce value is: {"
                        << g_exact << "}, but the computed value is: {"
                        << g_computed << "}.\n"
                        << "Absolute difference: " << diff << " > tol: " << tol << std::endl;
                }
            }
            if (Teuchos::nonnull(gx)) {
                out << "\nSolution! {1,2,3,4}\n"
                    << *gx << std::endl;
                Thyra::DetachedVectorView<double> x_view(gx);
                double x_exact[4] = {1, 2, 3, 4};

                double l2_diff = std::sqrt(std::pow(x_view(0) - x_exact[0], 2) +
                                           std::pow(x_view(1) - x_exact[1], 2) +
                                           std::pow(x_view(2) - x_exact[2], 2) +
                                           std::pow(x_view(3) - x_exact[3], 2));
                if (l2_diff > tol) {
                    status += 100;
                    out << "\nPiro_AnalysisDrvier:  Expected solution vector is : {"
                        << x_exact[0] << ", " << x_exact[1] << "}, but computed solution vector is: {"
                        << x_view(0) << ", " << x_view(1) << "}.\n"
                        << "Difference in l2 norm: " << l2_diff << " > tol: " << tol << std::endl;
                }
            }
            if (Teuchos::nonnull(dgdp)) {
                out << "\nSensitivities {6.66667, -8.0}\n"
                    << *dgdp << std::endl;
                Thyra::DetachedVectorView<double> dgdp_view(dgdp->col(0));
                double dgdp_exact[2] = {20./3, -8.};

                double l2_diff = std::sqrt(std::pow(dgdp_view(0) - dgdp_exact[0], 2) +
                                           std::pow(dgdp_view(1) - dgdp_exact[1], 2));
                if (l2_diff > tol) {
                    status += 100;
                    out << "\nPiro_AnalysisDrvier:  Expected sensitivity vector is: {"
                        << dgdp_exact[0] << ", " << dgdp_exact[1] << "}, but computed sensitivity vector is: {"
                        << dgdp_view(0) << ", " << dgdp_view(1) << "}.\n"
                        << "Difference in l2 norm: " << l2_diff << " > tol: " << tol << std::endl;
                }
            }

            if (Teuchos::nonnull(hv)) {
                double hv_exact[4] = {6., -10./3, -10./3, 4.};
                for (int i_direction = 0; i_direction < n_directions; i_direction++) {
                    out << "\n hv[" << i_direction << "]\n"
                        << *(hv->col(i_direction)) << std::endl;

                    Thyra::DetachedVectorView<double> hv_view(hv->col(i_direction));
                    double l2_diff = std::sqrt(std::pow(hv_view(0) - hv_exact[2 * i_direction], 2) +
                                               std::pow(hv_view(1) - hv_exact[2 * i_direction + 1], 2));
                    if (l2_diff > tol) {
                        status += 100;
                        out << "\nPiro_AnalysisDrvier:  Expected Hessian-vector product " << i_direction << " is: {"
                            << hv_exact[2 * i_direction] << ", " << hv_exact[2 * i_direction + 1] << "}, but computed Hessian-vector product is : {"
                            << hv_view(0) << ", " << hv_view(1) << "}.\n"
                            << "Difference in l2 norm: " << l2_diff << " > tol: " << tol << std::endl;
                    }
                }
            }

            out << "\n-----------------------------------------------------------------\n";
        }
        TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
        if (!success) status += 1000;

        overall_status += status;
    }

    if (Proc == 0) {
        if (overall_status == 0)
            std::cout << "\nTEST PASSED\n"
                      << std::endl;
        else
            std::cout << "\nTEST Failed:  " << overall_status << std::endl;
    }

    return status;
}
