#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include "basker.hpp"

using namespace std;

int main(int argc, char* argv[])
{

  typedef int Int;
  typedef double Entry;
  
  Int lnnz; //left guess
  Int unnz; //right guess

  
 /*load martrix*/
  Int annz;
  Int anrow, ancol;
  Int *Ap, *Ai;
  Entry *Ax;
  

 
  //Read in Matrix file from CSR file
  string temp;
  ifstream fp;
  fp.open(argv[1]);
  fp >> temp;
  anrow = atoi(temp.c_str());
  ancol =anrow;
  fp >> temp;
  annz = atoi(temp.c_str());
  
  cout << "Size: " << anrow << " nnz: " << annz << endl;
  Ap = new Int[anrow+1];
  Ai = new Int[annz];
  Ax = new Entry[annz];
  for(int i=0; i < anrow+1; i++)
    {
      string t;
      fp >> t;
      Ap[i] = atoi(t.c_str())-1;
    }
  for(int i=0; i < annz; i++)
    {
      string t;
      fp >> t;
      Ai[i] = atoi(t.c_str())-1;
      if(i == 0)
	{
	  cout << "First index: " << Ai[i] << endl;
	}
    }
  for(int i=0; i < annz; i++)
    {
      string t;
      fp >> t;
      Ax[i] = atof(t.c_str());
    }

  
  if(argc > 2)
    {
      lnnz = atoi(argv[2]);
       unnz = atoi(argv[3]);
    }
  else
    {
      lnnz = 2*annz;
      unnz = 2*annz; 
    }

 

  //Allocate some work space.
  Int* ws = (Int *)calloc((ancol)+(4*anrow), sizeof(Int));
  Entry* X = (Entry *)calloc(2*anrow, sizeof(Entry));
  Int* pinv = (Int *)calloc(ancol, sizeof(Int));
  Int* Lp = (Int *)calloc(ancol+1, sizeof(Int));
  Int* Li = (Int *)calloc(lnnz, sizeof(Int));
  Entry* Lx = (Entry *)calloc(lnnz, sizeof(Entry));
  Int* Up = (Int *)calloc(ancol+1, sizeof(Int));
  Int* Ui = (Int *)calloc(unnz, sizeof(Int));
  Entry* Ux = (Entry *)calloc(unnz,sizeof(Int));
  
  pinv[0] = -1;

  cout << "Done allocating space" << endl;
  Int result = basker::basker<Int, Entry>(Ap, Ai, Ax, anrow, ancol, ws, X, Lp, &Li, &Lx, Up, &Ui, &Ux, &lnnz, &unnz, pinv);

}

    
				      

