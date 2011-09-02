#define BUCKETSIZE 100

typedef struct optkdnode {
  int bucket;  /* true if the node is a bucket node */
  int discrim; /* discriminator of node */
  real cutval; 
  struct optkdnode *loson, *hison;
  int lopt,hipt;
} optkdNode;

