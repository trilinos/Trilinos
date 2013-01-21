//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
                                                                                
#include "NOX_Common.H"
#include "petscsnes.h"
#include "petscda.h"

/*
   User-defined routines
*/
int FormJacobian(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
int FormFunction(SNES,Vec,Vec,void*);
int MatrixFreePreconditioner(void*,Vec,Vec);

#include "FiniteDifference.H"

// Constructor - creates the Petsc objects (maps and vectors) 
FiniteDifference::FiniteDifference(SNES* snes_, void* ctx_) :
  ctx(ctx_),
  snes(snes_),
  matStruct(SAME_NONZERO_PATTERN)
{ }

// Destructor
FiniteDifference::~FiniteDifference()
{ 
  // Nothing currently owned by this class
}


// Matrix and Residual Fills
bool FiniteDifference::evaluate(FillType f, 
			      const Vec* soln, 
			      Vec* tmp_rhs, 
			      Mat* tmp_matrix)
{
  flag = f;
  int ierr = 0;

  // Set the incoming linear objects
  if (flag == RHS_ONLY) {
    rhs = tmp_rhs;
  } 
  else if (flag == MATRIX_ONLY) {
    A = tmp_matrix;
  } 
  else if (flag == ALL) { 
    rhs = tmp_rhs;
    A = tmp_matrix;
  } 
  else {
    std::cout << "ERROR: FiniteDifference::fillMatrix() - No such flag as " 
	 << flag << std::endl;
    throw;
  }

  // Begin RHS fill
  if((flag == RHS_ONLY) || (flag == ALL)) {
    ierr = FormFunction(*snes, *soln, *rhs, ctx);CHKERRQ(ierr);
    PetscScalar minusOne = 1.0;
    ierr = VecScale( *rhs, minusOne );CHKERRQ(ierr);
  }

  // Begin Jacobian fill
  if((flag == MATRIX_ONLY) || (flag == ALL)) {
    ierr = FormJacobian(*snes, *soln, A, A, &matStruct, ctx);CHKERRQ(ierr);
  } 

  return true;
}

Mat& FiniteDifference::getJacobian()
{
  return *A;
}


/* ------------------------------------------------------------------- */
/*
   FormFunction - Performs function (residual) evaluation

   Input Parameters:
   user - user-defined application context
   X - vector

   Output Parameter:
   f - function vector
 */
int FormFunction(SNES snes,Vec x,Vec f,void *ctx)
{
  ApplicationCtx *user = (ApplicationCtx*) ctx;
  DA             da = user->da;
  PetscScalar    *xx,*ff,d;
  int            i,ierr,M,xs,xm;
  Vec            xlocal;

  PetscFunctionBegin;
  ierr = DAGetLocalVector(da,&xlocal);CHKERRQ(ierr);
  /*
     Scatter ghost points to local vector, using the 2-step process
        DAGlobalToLocalBegin(), DAGlobalToLocalEnd().
     By placing code between these two statements, computations can
     be done while messages are in transition.
  */
  ierr = DAGlobalToLocalBegin(da,x,INSERT_VALUES,xlocal);CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(da,x,INSERT_VALUES,xlocal);CHKERRQ(ierr);

  /*
     Get pointers to vector data.
       - The vector xlocal includes ghost point; the vectors x and f do
         NOT include ghost points.
       - Using DAVecGetArray() allows accessing the values using global ordering  */
  ierr = DAVecGetArray(da,xlocal,&xx);CHKERRQ(ierr);
  ierr = DAVecGetArray(da,f,&ff);CHKERRQ(ierr);

  /*
     Get local grid boundaries (for 1-dimensional DA):
       xs, xm  - starting grid index, width of local grid (no ghost points)
  */
  ierr = DAGetCorners(da,&xs,PETSC_NULL,PETSC_NULL,&xm,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DAGetInfo(da,PETSC_NULL,&M,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                   PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);

  /*
     Set function values for boundary points; define local interior grid point range:
        xsi - starting interior grid index
        xei - ending interior grid index
  */
  d = 1.0/(user->h*user->h);
  if (xs == 0) { /* left boundary */
    ff[0] = xx[0] - 1.0;
    xs++;xm--;
  }
  if (xs+xm == M) {  /* right boundary */
    ff[xs+xm-1] = d*(2.0*xx[xs+xm-2] - 2.0*xx[xs+xm-1]) -
                           12.0*xx[xs+xm-1]*xx[xs+xm-1];
    xm--;
  }

  /*
     Compute function over locally owned part of the grid (interior points only)  */
  for (i=xs; i<xs+xm; i++) {
    ff[i] = d*(xx[i-1] - 2.0*xx[i] + xx[i+1]) - 12.0*xx[i]*xx[i];
  }

  /*
     Restore vectors
  */
  ierr = DAVecRestoreArray(da,xlocal,&xx);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(da,f,&ff);CHKERRQ(ierr);
  ierr = DARestoreLocalVector(da,&xlocal);CHKERRQ(ierr);
  //PetscFunctionReturn(0);

  return 0;
}

/* ------------------------------------------------------------------- */
/*
   FormJacobian - Evaluates Jacobian matrix.

   Input Parameters:
.  snes - the SNES context
.  x - input vector
.  dummy - optional user-defined context (not used here)

   Output Parameters:
.  jac - Jacobian matrix
.  B - optionally different preconditioning matrix
.  flag - flag indicating matrix structure
*/
int FormJacobian(SNES snes,Vec x,Mat *jac,Mat *B,MatStructure*flag,void *ctx)
{
  ApplicationCtx *user = (ApplicationCtx*) ctx;
  PetscScalar    *xx,d,A[3];
  int            i,j[3],ierr,M,xs,xm;
  DA             da = user->da;

  PetscFunctionBegin;
  /*
     Get pointer to vector data
  */
  ierr = DAVecGetArray(da,x,&xx);CHKERRQ(ierr);
  ierr = DAGetCorners(da,&xs,PETSC_NULL,PETSC_NULL,&xm,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);

  /*
    Get range of locally owned matrix
  */
  ierr = DAGetInfo(da,PETSC_NULL,&M,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
                   PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);

  /*
     Determine starting and ending local indices for interior grid points.
     Set Jacobian entries for boundary points.
  */

  d = 1.0/(user->h*user->h);

  if (xs == 0) {  /* left boundary */
    i = 0; A[0] = 1.0;
    ierr = MatSetValues(*jac,1,&i,1,&i,A,INSERT_VALUES);CHKERRQ(ierr);
    xs++;xm--;
  }
  if (xs+xm == M) { /* right boundary */
    i = M-1; A[0] = 1.0;
    j[0] = M-2;     j[1] = M-1;
    A[0] = 2.0*d;   A[1] = -2.0*d - 24.0*xx[M-1];
    ierr = MatSetValues(*jac,1,&i,2,j,A,INSERT_VALUES);CHKERRQ(ierr);
    xm--;
  }

  /*
     Interior grid points
      - Note that in this case we set all elements for a particular
        row at once.
  */
  for (i=xs; i<xs+xm; i++) {
    j[0] = i - 1; j[1] = i; j[2] = i + 1;
    A[0] = A[2] = d; A[1] = -2.0*d - 24.0*xx[i];
    ierr = MatSetValues(*jac,1,&i,3,j,A,INSERT_VALUES);CHKERRQ(ierr);
  }

  /*
     Assemble matrix, using the 2-step process:
       MatAssemblyBegin(), MatAssemblyEnd().
     By placing code between these two statements, computations can be
     done while messages are in transition.

     Also, restore vector.
  */

  ierr = MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(da,x,&xx);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  return 0;
}

/* ------------------------------------------------------------------- */
/*
   MatrixFreePreconditioner - This routine demonstrates the use of a
   user-provided preconditioner.  This code implements just the null
   preconditioner, which of course is not recommended for general use.

   Input Parameters:
.  ctx - optional user-defined context, as set by PCShellSetApply()
.  x - input vector

   Output Parameter:
.  y - preconditioned vector
*/
int MatrixFreePreconditioner(void *ctx,Vec x,Vec y)
{
  int ierr;
  ierr = VecCopy(x,y);CHKERRQ(ierr);
  return 0;
}

