
#include "Trilinos_LinearProblem.h"


//=============================================================================
Trilinos_LinearProblem::Trilinos_LinearProblem(Petra_RDP_CRS_Matrix * A, 
					       Petra_RDP_MultiVector * X,
					       Petra_RDP_MultiVector * B) 
  : A_(A),
    X_(X),
    B_(B),
    ProblemSymmetric_(false),
    PDL_(unsure)
{
  
}

//=============================================================================
Trilinos_LinearProblem::Trilinos_LinearProblem(const Trilinos_LinearProblem& Problem) 
  : A_(Problem.A_),
    X_(Problem.X_),
    B_(Problem.B_),
    ProblemSymmetric_(Problem.ProblemSymmetric_),
    PDL_(Problem.PDL_)
{
}
//=============================================================================
Trilinos_LinearProblem::~Trilinos_LinearProblem(void)  
{
}

