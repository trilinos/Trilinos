#include "migoct_const.h"

/* static void tag_regions(pOctant *octs, int *newpids, int nocts, 
                        LB_TAG **tags, int *ntags, int **tag_pids); */
static void tag_regions(pOctant *octs, int *newpids, int nocts, 
			LB_TAG **export_tags, int *nsentags, int **tag_pids,
			LB_TAG **kept_tags, int *nkeptags);
static void malloc_new_octants(int ntags, LB_TAG *tags, int *tag_pids,
                               int *nrectags, LB_TAG **import_tags,
			       LB_TAG *kept_tags, int nkeptags);
