#include "have_interface.hpp"

  int HAVE_PART=0;
  int HAVE_ITER_SOLVER=0;
  int HAVE_P_METHOD=0;
  
  #if defined HAVE_ZOLTAN2  
  #if defined HAVE_ISORROPIA
  HAVE_PART=3;
  #else
  HAVE_PART=2;
  #endif
  #elseif defined HAVE_ISORROPIA
  HAVE_PART=1;
  #endif

  #if defined HAVE_BELOS
  #if defined HAVE_AZTECOO
  HAVE_ITER_SOLVER=3;
  #else
  HAVE_ITER_SOLVER=2;
  #endif
  #elseif HAVE_AZTECOO
  HAVE_ITER_SOLVER=1;
  #endif

  #if defined HAVE_SCOTCH
  #if defined HAVE_PARMETIS
  HAVE_P_METHOD=3;
  #else
  HAVE_P_METHOD=2;
  #endif
  #elseif defined HAVE_PARMETIS
  HAVE_P_METHOD=1;
  #endif


int a = 3;
  template <class MatrixType>
  int have_interface (MatrixType A, ParameterList shyLUList)
  {
    
    
    cout << "Checking Paramter List" << endl;
    cout << a << endl;
    shyLUList.print(std::cout, 2, true, true);
    //Check if matrix is of right typeasd 
    
    //-------------------------------------------Check Partitioning Package
    string partitioningPackage = Teuchos::getParameter<string>(shyLUList, 
						       "Partitioning Package");

    cout << "Have_part_value" << endl;
    
    if(partitioningPackage.compare("Zoltan2")==0)
      {
	if((HAVE_PART!=3)||(HAVE_PART!=2))
	  { 
	    cout << "Zoltan2 has been selected, but is not supported" << endl;
	    exit(1);
	  }
      }
    else if(partitioningPackage.compare("Isorropia") ==0)
      {
	if((HAVE_PART!=3)||(HAVE_PART!=1))
	  {
	    cout << "Isorropia has been selected, but is not supported" << endl;
	    exit(1);
	  }
      }
    else
      {
	cout << "Defaulting to Zoltan2" << endl;
	/*Note: Add defaulting code */
      }

    //---------------------------------------Check Outer Solver Package
    string outerSolverPackage = Teuchos::getParameter<string>(shyLUList, 
							      "Outer Solver Library");

    if(outerSolverPackage.compare("Belos")==0)
      {
	if((HAVE_ITER_SOLVER!=3)||(HAVE_ITER_SOLVER!=2))
	  {
	    cout << "Belos has been selected, but is not supported" << endl;
	    exit(1);
	  }
      }
    else if(outerSolverPackage.compare("Aztec00")==0)
      {
	if((HAVE_ITER_SOLVER!=3)||(HAVE_ITER_SOLVER!=1))
	  {
	    cout << "AztecOO has been selected, but is not supported" << endl;
	    exit(1);
	  }
	else
	  {
	    cout << "**WARNING** AztecOO will be used for outer solver" << endl;
	    cout << "**WARNING** AztecOO is not recommended" << endl;
	  }
      }
    else
      {
	cout << "Defaulting to Belos for Outer Iterative Solver" << endl;
	/*Note: add defaulting code*/

      }

    //--------------------------- SUBLIST PARTITIONING ----------------*/
        
  }

template <class scale>
int test(scale a)
{
  return (int) a;
}
   
template int test(double);
template int have_interface(Epetra_CrsMatrix, Teuchos::ParameterList);



    
  
