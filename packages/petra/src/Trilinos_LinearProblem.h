#ifndef _TRILINOS_LINEARPROBLEM_H_
#define _TRILINOS_LINEARPROBLEM_H_

//! Trilinos_LinearProblem:  The Trilinos Linear Problem Class.
/*! The Trilinos_LinearProblem class is a wrapper that encapsulates the 
  general information needed for describing a linear system of equations.  
  Currently it accepts a Petra matrix, initial guess and RHS. It allows the
  user to set a ProblemDifficultyLevel (easy, moderate, hard, unsure), from
  which a set of solution parameters (solver choice etc.) will be inferred by
  Aztec_OO. (Aztec_OO is the wrapper by which a Trilinos_LinearProblem can be
  passed to Aztec's solution methods.)<p>
  Simple example (pseudo-code) usage:<p>
<code>
  Petra_RDP_CRS_Matrix A;<br>
  Petra_RDP_MultiVector X, B;<p>
  //...code to assemble/fill the matrix and vectors...<p>
  Trilinos_LinearProblem prob(&A, &X, &B);<br>
  prob.SetPDL(hard);<br>
  prob.AssertSymmetric(); //...if we know that A is symmetric...<p>
  Aztec_OO aztec(prob);<br>
  aztec.Iterate(500, 1.e-9); //...maxIters, tolerance...<p>
</code>
  Note:
<ul>
<li> Users may also choose from the full set of Aztec options and params 
directly using the Aztec_OO class.
<li> Also, there will ultimately be a "Trilinos_Strategy" class which will
  provide more sophistication in managing the interaction between
  the problem and the solver library.
</ul>
*/

#include "Petra_Petra.h"
#include "Petra_Comm.h"
#include "Petra_RDP_MultiVector.h"
#include "Petra_RDP_CRS_Matrix.h"
#ifndef DOXYGEN_SHOULD_SKIP_THIS
enum ProblemDifficultyLevel {easy, moderate, hard, unsure};
#endif
class Trilinos_LinearProblem {
    
  public:
  //!  Trilinos_LinearProblem Constructor.
  /*! Creates a Trilinos_LinearProblem instance. 
  */
  Trilinos_LinearProblem(Petra_RDP_CRS_Matrix * A, Petra_RDP_MultiVector * X,
			 Petra_RDP_MultiVector * B);

  //! Trilinos_LinearProblem Copy Constructor.
  /*! Makes copy of an existing Trilinos_LinearProblem instance.
  */
  Trilinos_LinearProblem(const Trilinos_LinearProblem& Problem);


  // Solver strategy assertions

  void AssertSymmetric(){ProblemSymmetric_ = true;};

  //! Set problem difficulty level.
  /*! This sets info that will be used by the solution object (e.g., Aztec_OO)
     that this Trilinos_LinearProblem is ultimately passed to.
  */
  void SetPDL(ProblemDifficultyLevel PDL) {PDL_ = PDL; }

  //! Trilinos_LinearProblem Destructor.
  /*! Completely deletes a Trilinos_LinearProblem object.  
  */
  virtual ~Trilinos_LinearProblem(void);

 private:

  friend class Aztec_OO;

  Petra_RDP_CRS_Matrix * A_;
  Petra_RDP_MultiVector * X_;
  Petra_RDP_MultiVector * B_;

  bool ProblemSymmetric_;
  ProblemDifficultyLevel PDL_;
    
};

#endif /* _TRILINOS_LINEARPROBLEM_H_ */
