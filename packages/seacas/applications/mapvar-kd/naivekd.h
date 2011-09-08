
typedef struct TreeRec   /* Tree Record */
{
  int     Discrim; /* Stores which dimension this tree partitions data  */ 
  int     Index;   /* Stores index of point in the points array */
  struct  TreeRec   *Left, *Right; /*Pointers to sons of tree node */
} TREEREC, TreeRec, *TreePtr;
