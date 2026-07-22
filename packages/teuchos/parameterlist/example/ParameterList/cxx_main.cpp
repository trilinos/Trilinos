// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main(int argc, char* argv[])
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::Array;
  using Teuchos::tuple;
  using Teuchos::as;

  bool success = false;
  bool verbose = true;
  try {

    std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

    // Creating an empty parameter list looks like:
    ParameterList myPL;

    // Setting parameters in this list can be easily done:
    myPL.set("Max Iters", 1550, "Determines the maximum number of iterations in the solver");
    myPL.set("Tolerance", 1e-10, "The tolerance used for the convergence check");

    // For the "Solver" option, create a validator that will automatically
    // create documentation for this parameter but will also help in validation.
    RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
      solverValidator = Teuchos::rcp(
          new Teuchos::StringToIntegralParameterEntryValidator<int>(
            Teuchos::tuple<std::string>( "GMRES", "CG", "TFQMR" )
            ,"Solver"
            )
          );
    myPL.set(
        "Solver"
        ,"GMRES" // This will be validated by solverValidator right here!
        ,"The type of solver to use"
        ,solverValidator
        );

    /* The templated ``set'' method should cast the input {\it value} to the
       correct data type.  However, in the case where the compiler is not casting the input
       value to the expected data type, an explicit cast can be used with the ``set'' method:
       */
    myPL.set("Tolerance", as<float>(1e-10), "The tolerance used for the convergence check");

    /* Reference-counted pointers can also be passed through a ParameterList.
       To illustrate this we will use the Array class to create an array of 10 doubles
       representing an initial guess for a linear solver, whose memory is being managed by a
       RCP.
       */

    myPL.set<Array<double> >("Initial Guess", tuple<double>( 10, 0.0 ),
        "The initial guess as a RCP to an array object.");

    /* A hierarchy of parameter lists can be constructed using {\tt ParameterList}.  This
       means another parameter list is a valid {\it value} in any parameter list.  To create a sublist
       in a parameter list and obtain a reference to it:
       */
    ParameterList& Prec_List = myPL.sublist("Preconditioner", false,
        "Sublist that defines the preconditioner.");

    // Now this parameter list can be filled with values:
    Prec_List.set("Type", "ILU", "The tpye of preconditioner to use");
    Prec_List.set("Drop Tolerance", 1e-3,
        "The tolerance below which entries from the\n""factorization are left out of the factors.");

    // The parameter list can be queried about the existance of a parameter, sublist, or type:
    // Has a solver been chosen?
    bool solver_defined = false, prec_defined = false, dtol_double = false;
    solver_defined = myPL.isParameter("Solver");
    TEUCHOS_ASSERT_EQUALITY(solver_defined, true);
    // Has a preconditioner been chosen?
    prec_defined = myPL.isSublist("Preconditioner");
    TEUCHOS_ASSERT_EQUALITY(prec_defined, true);
    // Has a tolerance been chosen and is it a double-precision number?
    bool tol_double = false;
    tol_double = myPL.INVALID_TEMPLATE_QUALIFIER isType<double>("Tolerance");
    TEUCHOS_ASSERT_EQUALITY(tol_double, false); // It is 'float'!
    // Has a drop tolerance been chosen and is it a double-precision number?
    dtol_double = Teuchos::isParameterType<double>(Prec_List, "Drop Tolerance");
    TEUCHOS_ASSERT_EQUALITY(dtol_double, true);

    // Parameters can be retrieved from the parameter list in quite a few ways:
    // Get method that creates and sets the parameter if it doesn't exist.
    int its = -1;
    its = myPL.get("Max Iters", 1200);
    TEUCHOS_ASSERT_EQUALITY(its, 1550); // Was already ste
    // Get method that retrieves a parameter of a particular type that must exist.
    float tol = -1.0;
    tol = myPL.get<float>("Tolerance");
    TEUCHOS_ASSERT_EQUALITY(tol, as<float>(1e-10));
    // Get the "Solver" value and validate!
    std::string
      solver = solverValidator->validateString(
          Teuchos::getParameter<std::string>(myPL,"Solver")
          );

    // We can use this same syntax to get arrays out, like the initial guess.
    Array<double> init_guess = myPL.get<Array<double> >("Initial Guess");

    std::cout << "\n# Printing this parameter list using opeator<<(...) ...\n\n";
    std::cout << myPL << std::endl;

    std::cout << "\n# Printing the parameter list only showing documentation fields ...\n\n";
    myPL.print(std::cout,
        ParameterList::PrintOptions().showDoc(true).indent(2).showTypes(true));

    /* It is important to note that mispelled parameters
       (with additional space characters, capitalizations, etc.) may be ignored.
       Therefore, it is important to be aware that a given parameter has not been used.
       Unused parameters can be printed with method:
       */
    std::cout << "\n# Showing unused parameters ...\n\n";
    myPL.unused( std::cout );

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
