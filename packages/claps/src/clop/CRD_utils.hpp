#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "Epetra_CrsMatrix.h"

namespace CRD_utils {
  int find_index(int a[], int n, int gdof);
  void sort_and_cull(int a[], int n, int & m);
  void Epetra_datfile(const Epetra_CrsMatrix* A, char fname[]);
  void Epetra_datfile(int* A, int N, char fname[]);
  void spmat_datfile(int nrow, int rowbegp [], int colidxp [],
		     double val[], char fname[]);
  class Graph_class 
  {
  public:
    Graph_class(int N_, int A1_[], int A2_[]);
    ~Graph_class();
    void Components(int component_[], int & comp_num);
  private:
    void DFS(const int v, const int comp_num);
    int N;
    int *A1, *A2;
    int *component;
  };
}
