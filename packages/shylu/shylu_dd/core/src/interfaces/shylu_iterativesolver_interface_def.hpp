// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLU_ITERATIVESOLVER_INTERFACE_DEF_HPP
#define SHYLU_ITERATIVESOLVER_INTERFACE_DEF_HPP

#include "ShyLU_DDCore_config.h"
#include "shylu_iterativesolver_interface_decl.hpp"

#include <Teuchos_XMLParameterListHelpers.hpp>



namespace ShyLU
{
  template <class Matrix, class Vector>
  IterativeSolverInterface<Matrix,Vector>::IterativeSolverInterface()
  {}

  template <class Matrix, class Vector>
  IterativeSolverInterface<Matrix,Vector>::IterativeSolverInterface(Matrix *inA, Teuchos::ParameterList *inpList)
  {
    A = inA;
    pList = inpList;
  }

  template <class Matrix, class Vector>
  IterativeSolverInterface<Matrix,Vector>::~IterativeSolverInterface()
  {}

  template <class Matrix, class Vector>
  int IterativeSolverInterface<Matrix,Vector>::init_matrix(Matrix *inA, Teuchos::ParameterList *inpList)
  {
    A=inA;
    pList = inpList;
    return 0;
  }

  template <class Matrix, class Vector>
  int IterativeSolverInterface<Matrix,Vector>::solveAztec(Vector *b, Vector *x)
  {
    printf("**Error** Aztec only supported by Epetra \n");
    return 1;
  }


  template <>
  int IterativeSolverInterface<Epetra_CrsMatrix,Epetra_MultiVector>::solveAztec(Epetra_MultiVector *b, Epetra_MultiVector *x)
  {

    Teuchos::ParameterList subList = pList->sublist("Aztec Input");
    //Add additional sublist parameters for aztec


    e_problem.SetOperator(A);
    e_problem.SetRHS(x);
    e_problem.SetLHS(b);

    //set Aztec paramters
    solver_aztec = new AztecOO(e_problem);
    solver_aztec->SetAztecOption(AZ_precond, AZ_none);
    solver_aztec->SetAztecOption(AZ_solver,  AZ_gmres);
    solver_aztec->SetAztecOption(AZ_max_iter, 1000);

    //solve
    solver_aztec->Iterate(1000, 1e-5);

    return 0;
  }//end solveAztec

  template <class Matrix, class Vector>
  int IterativeSolverInterface<Matrix,Vector>::solveBelos(Vector *b, Vector *x)
  {

    typedef typename Matrix::scalar_type           scalar_type;
    typedef typename Matrix::local_ordinal_type    local_ordinal_type;
    typedef typename Matrix::global_ordinal_type   global_ordinal_type;
    typedef typename Matrix::node_type             node_type;
    typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type>                                     OP;

    //Input list of belos
    //Teuchos::RCP<Teuchos::ParameterList> belosList = Teuchos::parameterList();

    //Teuchos::RCP<Teuchos::ParameterList> subList = pList->sublist("Belos Input");


    Teuchos::RCP<Teuchos::ParameterList> pListRCP = Teuchos::rcp(pList,false);
    Teuchos::RCP<Teuchos::ParameterList> belosList = Teuchos::sublist(pListRCP, "Belos Input");

      //pListRCP->sublist("Belos Input");

      //pList->sublist("Belos Input");

    std::string solver_name;
    if(belosList->isParameter("Solver"))
      {
        //solver_name = Teuchos::getParameter<string>(&belosList, "Solver");
        solver_name = belosList->get<std::string>("Solver");
      }


    //printf("start\n");
    //Add additional sublipst parameter
    //belosList.set( "Num blocks", 100);
    // belosList->set( "Block Size" , 10);
    // belosList->set( "Maximum Iterations", 1000);
    //belosList->set( "Maximum Restarts", 10);
    //belosList->set( "Convergence Tolerance" , 1e-5);

    printf("load belos parameters\n");
    Teuchos::RCP<Belos::LinearProblem<scalar_type, Vector, OP> > belos_problem;
    belos_problem = Teuchos::rcp(new Belos::LinearProblem<scalar_type, Vector, OP>());

    belos_problem->setOperator(Teuchos::rcp(A,false));
    belos_problem->setRHS(Teuchos::rcp(x,false));
    belos_problem->setLHS(Teuchos::rcp(b,false));
    belos_problem->setProblem();

    //Based on input

    Belos::SolverFactory<scalar_type, Vector, OP> factory;

    Teuchos::RCP<Belos::SolverManager<scalar_type,Vector,OP> > solver_belos;

    solver_belos = factory.create(solver_name, belosList);

    solver_belos->setProblem(belos_problem);


    // FIXME (mfh 23 Aug 2015) Not using this return value results in
    // a build warning.  You should use the return value from Belos to
    // return something other than 0, in case the solve didn't succeed
    // (didn't return Belos::Converged).  I don't know what ShyLU uses
    // for a "didn't succeed, but not catastrophically so" error code,
    // or I would put it in myself.
    //
    //Belos::ReturnType ret = solver_belos->solve();
    (void) solver_belos->solve();

    return 0;
  }// end solveBelos

   template <>
  int IterativeSolverInterface<Epetra_CrsMatrix,Epetra_MultiVector>::solveBelos(Epetra_MultiVector *b, Epetra_MultiVector *x)
  {
    // FIXME (mfh 23 Aug 2015) Maybe you should throw an exception or
    // something to signal that this routine is unimplemented, instead
    // of just returning 0.

    //typedef double                                 scalar_type;
    //typedef Teuchos::ScalarTraits<scalar_type>     scalar_type_traits;
    //typedef scalar_type_traits::magnitudeType      magnitude_type; // unused
    //typedef Epetra_Operator                        OP; // unused

    //Input list of belos
    //Teuchos::RCP<Teuchos::ParameterList> belosList = Teuchos::parameterList();

    //Teuchos::RCP<Teuchos::ParameterList> subList = pList->sublist("Belos Input");

    /*
    Teuchos::RCP<Teuchos::ParameterList> belosList = Teuchos::rcp(pList)->sublist("Belos Input");

    string solver_name;
    if(belosList->isParameter("Solver"))
      {
        solver_name = belosList->getEntry("Solver");
      }


    //printf("start\n");
    //Add additional sublipst parameter
    //belosList.set( "Num blocks", 100);
    // belosList->set( "Block Size" , 10);
    // belosList->set( "Maximum Iterations", 1000);
    //belosList->set( "Maximum Restarts", 10);
    //belosList->set( "Convergence Tolerance" , 1e-5);

    printf("load belos parameters\n");
    Teuchos::RCP<Belos::LinearProblem<scalar_type, Vector, OP> > belos_problem;
    belos_problem = Teuchos::rcp(new Belos::LinearProblem<scalar_type, Vector, OP>());

    belos_problem->setOperator(Teuchos::rcp(A,false));
    belos_problem->setRHS(Teuchos::rcp(x,false));
    belos_problem->setLHS(Teuchos::rcp(b,false));
    belos_problem->setProblem();

    //Based on input

    Belos::SolverFactory<scalar_type, Vector, OP> factory;

    Teuchos::RCP<Belos::SolverManager<scalar_type,Vector,OP> > solver_belos;

    solver_belos = factory.create(solver_name, belosList);

    solver_belos->setProblem(belos_problem);


    Belos::ReturnType ret = solver_belos->solve();

    */
    return 0;
  }// end solveBelos


  template <class Matrix, class Vector>
  int IterativeSolverInterface<Matrix,Vector>::solve(Vector *b, Vector* x)
  {
    std::string solverpackage = Teuchos::getParameter<std::string>(*pList, "Iterative Solver Package");

    std::cout << "found package " << solverpackage << std::endl;

    if(solverpackage.compare("Aztec")==0)
      {
        printf("**Warning** Aztec is not preferred \n");
        solveAztec(b,x);
      }
    else if(solverpackage.compare("Belos")==0)
      {
        solveBelos(b,x);
      }
    else
      {
        printf("**Error**Iterative Solver Package Not Provided\n");
        return 1;
      }
    return 0;
  }//end solve()


}//end namespace shylu

#endif //end ifdef shylu_interative_fie

#if defined(ShyLU_DDCore_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ShyLU_DDCore package is deprecated"
#endif
#endif

