#define SWAP(type, r,s)  do{type t=r; r=s; s=t; } while(0)

static void siftDown( double *a, int start, int count);

void sort_double(int count, double ra[])
{
  int start, end;
 
  /* heapify */
  for (start = (count-2)/2; start >=0; start--) {
    siftDown( ra, start, count);
  }
 
  for (end=count-1; end > 0; end--) {
    SWAP(double, ra[end],ra[0]);
    siftDown(ra, 0, end);
  }
}
 
static void siftDown( double *a, int start, int end)
{
  int root = start;
 
  while ( root*2+1 < end ) {
    int child = 2*root + 1;
    if ((child + 1 < end) && (a[child] < a[child+1])) {
      child += 1;
    }
    if (a[root] < a[child]) {
      SWAP(double, a[child], a[root] );
      root = child;
    }
    else
      return;
  }
}

void      shell_sort(int n, double *arr)
{
  sort_double(n, arr);
}

