#include "Thyra.hpp"
#include "ThyraCG.cpp"
#include "Teuchos_Array.hpp"


using namespace Thyra;
using namespace Tuechos;

/** Set up a matrix and a right hand side to test the ThyraCG
    solver */

int main(int argc, void *argv[])
{
 
  /* Set up a the vector space.  We'll first use Epetra.  The easiest
     way to do this is to set up a VectorType and then use the
     createSpace method.  */

  VectorSpace<double> type = new EpetraVectorType();

  /* Set up the dimension and an array with the locally owned indices.
     In this case, all of the indices are locally owned. */

  int dimension = 100;
  int[dimension] indices;

  for (int i = 0; i < dimension; i++)
    {
      indices[i] = i;
    }
  
  VectorSpace<double> space = 
    type.createSpace(dimension, dimension, indices); 

  /* Set up the matrix:
    
     Since an Epetra matrix has a sophisticated structure associated
     with it, it is easiest to use a special tool in Thyra, called a
     TSFEpetraMatrixFactory, to help to construct the actual
     matrix. This allows us to specify the row structure of the
     matrix. We specify the factory with the domain and range spaces,
     which, in our case, are both the same.  */

  EpetraMatrixFactory factory(space, space);

  /* Now we use the factory to specify the row structure.  We'll do an
     easy case and create a banded matrix  */

  int bandwidth = 3; // bandwidth
  Teuchos::Array<Array<int> > colInd;
  colInd.resize(dimension);

  for (row = 0; row < dimension; row++)
    {
      colInd[row].resize(0);
      /* set up the column indices.  */
      int start = row - bandwidth / 2;
      int end = row + bandwidth / 2;
      for (int i = start; i < end; i++)
	{
	  if (i < 0 || i >= dimension) continue;
	  colInd[row].append(i);
	}
      /* set the row structure */
      int numElementsToInsert = colInd.size();
      factory.initializeNonzerosInRow(row, numElementsToInsert, colInd); 
    }

  /* Freeze the structure  */
  factory.finalize();

  /* Create the matrix */
  
  LinearOperator<double> A = factory.createMatrix();

  /* Fill the matrix.  First, we have to cast it to a Loadable
     matrix  */

  RCP<LoadableMatrix<double> > LoadA = A.matrix();
  LoadA.zero();  //initializes all the elements to zero

  /* Fill the matrix row by row.  Each row will have -1.0 off the
     diagonal element and the value of bandwidth on the diagonal.
     This ensures that the matrix is symmetric and positive definite
     and hence suitable for CG.  */

  Teuchos::Array<double> vals;  //An array to hold the values 
  for (int row = 0; row < dimension; row++)
    {
      vals.resize(0);
      int start = row - bandwidth / 2;
      int end = row + bandwidth / 2;
      for (int i = start; i < end; i++)
	{
	  if (i < 0 || i >= dimension) continue;
	  if (i == row) 
	    {
	      vals.append(double(bandwidth));
	    }
	  else
	    {
	      vals.append(-1.0);
	    }
	  addToRow(row, vals.size(), *colInd[row], *vals); 
	}
    }

  /* Now create the right hand side.  We'll do this by just assigning
     values to a vector x and multiplying A by x to obtain an
     appropriate right hand side, b */

  Vector<double> xTrue = space.createMember();
  xTrue.setToConstant(3.2);
  Vector<double> b = A * xTrue;

  /* Now, solve the system useing ThyraCG  */

  Vector<double> x = space.createMember();
  ThyraCG thyraCG();
  
  /* set tolerance and max number of iterations  */
  double tol = 1.0e-07;
  int maxit = 1000;
  thyraCG.run(A, x, b, tol, maxit); 


  /* check the answer  */
  Vector<double> diff = x - xTrue;
  cout << "The norm of the difference between the computed "
       << "and the true solution is " << diff.norm2() << endl;

  exit(0);

}
