#ifndef _IFP_SPARSEUTIL_H_
#define _IFP_SPARSEUTIL_H_

void shell_sort(
  const int n,
  int x[]);

void allocate_ilu(
  const int levfill,
  const int n,
  int *nzl, int *nzu,
  const int ia[], const int ja[],
  int *ial[], int *jal[],
  int *iau[], int *jau[],
  int growthl, int growthu);

int symbolic_ilu(
  const int levinc,
  const int n,
  int *nzl,
  int *nzu,
  const int ia[], const int ja[],
  int ial[], int jal[],
  int iau[], int jau[]);

#endif // _IFP_SPARSEUTIL_H_
