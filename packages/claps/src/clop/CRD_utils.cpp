#include "CRD_utils.hpp"

namespace CRD_utils {

int find_index(int a[], int n, int gdof)
{
  int *loc, loc_col;
  loc = std::lower_bound(&a[0], &a[n], gdof);
  loc_col = (int) (loc - &a[0]);
  if ((loc_col < 0) || (loc_col > n-1)) return -1;
  if (a[loc_col] != gdof) return -1;
  return loc_col;
}
  
void sort_and_cull(int a[], int n, int & m)
{
  if (n == 0) {
    m = 0;
    return;
  }
  std::sort(&a[0], &a[n]);
  int previd = a[0];
  m = 1;
  for (int i=1; i<n; i++) {
    if (a[i] != previd) {
      previd = a[i];
      a[m] = a[i];
      m++;
    }
  }
}

void Epetra_datfile(const Epetra_CrsMatrix* A, char fname[])
{
  int i, j, NumEntries, *Indices;
  double *Values;
  ofstream fout;
  sprintf(fname, "%s.dat", fname);
  fout.open(fname);
  for (i=0; i<A->NumMyRows(); i++) {
    A->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      fout << A->GRID(i) << " " << A->GCID(Indices[j]) << setw(22) << setprecision(15)
	   << Values[j] << endl;
    }
  }
  fout.close();
}

void spmat_datfile(int nrow, int rowbegp [], int colidxp [],
				 double val[], char fname[])
{
  int i, j;
  ofstream fout;
  sprintf(fname, "%s.dat", fname);
  fout.open(fname);
  for (i=0; i<nrow; i++) {
    for (j=rowbegp[i]; j<rowbegp[i+1]; j++) {
      fout << i << " " << colidxp[j] << setw(22) << setprecision(15)
	   << val[j] << endl;
    }
  }
  fout.close();
}

void Epetra_datfile(int* A, int N, char fname[])
{
  int i;
  ofstream fout;
  sprintf(fname, "%s.dat", fname);
  fout.open(fname);
  for (i=0; i<N; i++) fout << A[i] << endl;
  fout.close();
}

Graph_class::Graph_class(int N_, int A1_[], int A2_[])
{
  N = N_;
  A1 = A1_;
  A2 = A2_;
  component = 0;
  component = new int[N];
  for (int i=0; i<N; i++) component[i] = -1;
}

Graph_class::~Graph_class()
{
  delete [] component; component = 0;
}

void Graph_class::Components(int component_[], int & comp_num)
{
  comp_num = 0;
  for (int i=0; i<N; i++) {
    if (component[i] == -1) {
      DFS(i, comp_num);
      comp_num++;
    }
  }
  for (int i=0; i<N; i++) component_[i] = component[i];
}

void Graph_class::DFS(const int v, const int comp_num)
{
  component[v] = comp_num;
  for (int i=A2[v];i<A2[v+1];i++) {
    int adj_vertex = A1[i];
    if (component[adj_vertex] == -1) DFS(adj_vertex, comp_num);
  }
}

} // end of namespace CRD_utils
