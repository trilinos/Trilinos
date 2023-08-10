#pragma once
#include "../internal/internal.h"
#include "../refine_map/refine_map.h"
#include "params.h"
#include "structs.h"
#include <stdio.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* assign/assign.c */
extern void assign(struct vtx_data **graph, double **yvecs, int nvtxs, int ndims, int cube_or_mesh,
                   int nsets, double *wsqrt, int *sets, int *active, int mediantype, double *goal,
                   int vwgt_max);

/* assign/assign_out.c */
extern void assign_out(int nvtxs, int *sets, int nsets, char *outname);

/* assign/mapper.c */
extern void mapper(struct vtx_data **graph, double **xvecs, int nvtxs, int *active, int *sets,
                   int ndims, int cube_or_mesh, int nsets, int mediantype, double *goal,
                   int vwgt_max);

/* assign/median.c */
extern void median_assign(struct vtx_data **graph, double *vals, int nvtxs, double *goal,
                          int using_vwgts, int *sets, double wlow, double whigh, double guess);
extern void median(struct vtx_data **graph, double *vals, int nvtxs, int *active, double *goal,
                   int using_vwgts, int *sets);

/* assign/merge_assign.c */
extern void merge_assignments(int *assignment, int *subassign, int *subsets, int subnvtxs,
                              int *loc2glob);

/* assign/rec_median.c */
extern void rec_median_1(struct vtx_data **graph, double *vals, int nvtxs, int *active,
                         int cube_or_mesh, int nsets, double *goal, int using_vwgts, int *assign,
                         int top);
extern void rec_median_k(struct vtx_data **graph, double **vals, int nvtxs, int *active, int ndims,
                         int cube_or_mesh, double *goal, int using_vwgts, int *assign);

/* assign/rotate.c */
extern void rotate2d(double **yvecs, int nmyvtxs, double theta);
extern void rotate3d(double **yvecs, int nmyvtxs, double theta, double phi, double gamma2);

/* assign/y2x.c */
extern void y2x(double **xvecs, int ndims, int nmyvtxs, double *wsqrt);
extern void x2y(double **yvecs, int ndims, int nmyvtxs, double *wsqrt);

/* bpmatch/checkbp.c */
extern void checkbp(struct vtx_data **graph, double **xvecs, int *sets, double *dists, int nvtxs,
                    int ndims);

/* bpmatch/genvals2d.c */
extern void genvals2d(double **xvecs, double *vals[4][MAXSETS], int nvtxs);

/* bpmatch/genvals3d.c */
extern void genvals3d(double **xvecs, double *vals[8][MAXSETS], int nvtxs);

/* bpmatch/inits2d.c */
extern int  findindex(int *indices, double *vals, double target, int nvals);
extern void inits2d(struct vtx_data **graph, double **xvecs, double *vals[4][MAXSETS],
                    int *indices[4][MAXSETS], int nvtxs, double *dist, int startvtx[4][MAXSETS],
                    double *size, int *sets);

/* bpmatch/inits3d.c */
extern void inits3d(struct vtx_data **graph, double **xvecs, double *vals[8][MAXSETS],
                    int *indices[8][MAXSETS], int nvtxs, double *dist, int startvtx[8][MAXSETS],
                    double *size, int *sets);

/* bpmatch/map2d.c */
extern void map2d(struct vtx_data **graph, double **xvecs, int nvtxs, int *sets, double *goal,
                  int vwgt_max);

/* bpmatch/map3d.c */
extern void map3d(struct vtx_data **graph, double **xvecs, int nvtxs, int *sets, double *goal,
                  int vwgt_max);

/* bpmatch/movevtxs.c */
extern void movevtxs(struct vtx_data **graph, int nvtxs, int nsets, double *dist,
                     int *indices[][MAXSETS], double *vals[][MAXSETS], int startvtx[][MAXSETS],
                     int *sets, double *size, double *goal, int vwgt_max);

/* bpmatch/sorts2d.c */
extern void sorts2d(double *vals[4][MAXSETS], int *indices[4][MAXSETS], int nvtxs);

/* bpmatch/sorts3d.c */
extern void sorts3d(double *vals[8][MAXSETS], int *indices[8][MAXSETS], int nvtxs);

/* coarsen/coarsen.c */
extern void coarsen(struct vtx_data **graph, int nvtxs, int nedges, int using_vwgts,
                    int using_ewgts, float *term_wgts[], int igeom, float **coords, double **yvecs,
                    int ndims, int solver_flag, int vmax, double eigtol, int nstep, int step,
                    int give_up);

/* coarsen/coarsen1.c */
extern void coarsen1(struct vtx_data **graph, int nvtxs, int nedges, struct vtx_data ***pcgraph,
                     int *pcnvtxs, int *pcnedges, int **pv2cv, int igeom, float **coords,
                     float **ccoords, int using_ewgts);

/* coarsen/countcedges.c */
extern void countcedges(struct vtx_data **graph, int nvtxs, int *start, int *seenflag, int *mflag,
                        int *v2cv, int *pcnedges);

/* coarsen/interpolate.c */
extern void ch_interpolate(double **vecs, double **cvecs, int ndims, struct vtx_data **graph,
                           int nvtxs, int *v2cv, int using_ewgts);

/* coarsen/makeccoords.c */
extern void makeccoords(struct vtx_data **graph, int cnvtxs, int *cv2v_ptrs, int *cv2v_vals,
                        int igeom, float **coords, float **ccoords);

/* coarsen/makecgraph.c */
extern void makecgraph(struct vtx_data **graph, int nvtxs, struct vtx_data ***pcgraph, int *pcnvtxs,
                       int *pcnedges, int *mflag, int *v2cv, int nmerged, int using_ewgts,
                       int igeom, float **coords, float **ccoords);

/* coarsen/makecgraph2.c */
extern void makecgraph2(struct vtx_data **graph, int nvtxs, int nedges, struct vtx_data ***pcgraph,
                        int *pcnvtxs, int *pcnedges, int *mflag, int *v2cv, int nmerged,
                        int using_ewgts, int igeom, float **coords, float **ccoords);

/* coarsen/makefgraph.c */
extern void makefgraph(struct vtx_data **graph, int nvtxs, int nedges, struct vtx_data ***pcgraph,
                       int cnvtxs, int *pcnedges, int *v2cv, int using_ewgts, int igeom,
                       float **coords, float **ccoords);

/* coarsen/makev2cv.c */
extern void makev2cv(int *mflag, int nvtxs, int *v2cv);

/* coarsen/maxmatch.c */
extern int maxmatch(struct vtx_data **graph, int nvtxs, int nedges, int *mflag, int using_ewgts,
                    int igeom, float **coords);

/* coarsen/maxmatch1.c */
extern int maxmatch1(struct vtx_data **graph, int nvtxs, int *mflag, int using_ewgts);

/* coarsen/maxmatch2.c */
extern int maxmatch2(struct vtx_data **graph, int nvtxs, int *mflag, int using_ewgts);

/* coarsen/maxmatch3.c */
extern int maxmatch3(struct vtx_data **graph, int nvtxs, int *mflag, int using_ewgts);

/* coarsen/maxmatch4.c */
extern int maxmatch4(struct vtx_data **graph, int nvtxs, int nedges, int *mflag, int using_ewgts);

/* coarsen/maxmatch5.c */
extern int maxmatch5(struct vtx_data **graph, int nvtxs, int *mflag, int igeom, float **coords);

/* coarsen/maxmatch9.c */
extern int maxmatch9(struct vtx_data **graph, int nvtxs, int *mflag, int using_ewgts);

/* connect/add_edges.c */
extern void add_edges(struct vtx_data **graph, struct edgeslist *new_edges,
                      struct ilists **old_edges, struct flists **old_ewgts, int using_ewgts);

/* connect/connect_enforce.c */
extern void connect_enforce(struct vtx_data **graph, int nvtxs, int using_ewgts, int *assignment,
                            double *goal, int nsets_tot, int *total_move, int *max_move);

/* connect/connected.c */
extern void make_connected(struct vtx_data **graph, int nvtxs, int *nedges, int *mark, int *vtxlist,
                           struct connect_data **cdata, int using_ewgts);
extern void make_unconnected(struct vtx_data **graph, int *nedges, struct connect_data **cdata,
                             int using_ewgts);
extern void print_connected(struct connect_data *cdata);
extern void free_edgeslist(struct edgeslist *edge_list);

/* connect/find_comps.c */
extern int find_comps(struct vtx_data **graph, int nvtxs, int *mark, int *vtxlist);
extern int find_edges(struct vtx_data **graph, int nvtxs, int *mark, int *vtxlist,
                      struct edgeslist **edges);

/* connect/heap.c */
extern void   heapify(struct heap *heap, int index, int nvals, int *map);
extern void   heap_build(struct heap *heap, int nvals, int *map);
extern double heap_extract_max(struct heap *heap, int nvals, int *ptag, int *map);
extern void   heap_update_val(struct heap *heap, int index, double newval, int *map);

/* eigen/Tevec.c */
extern double Tevec(double *alpha, double *beta, int j, double ritz, double *s);

/* eigen/bidir.c */
extern double bidir(double *alpha, double *beta, int j, double ritz, double *s, double hurdle);

/* eigen/bisect.c */
extern int bisect(double *alpha, double *beta, int j, double Anorm, double *workj, double *ritz,
                  int nevals_left, int nevals_right, double tol, double *ritz_sav, int max_steps);

/* eigen/checkeig.c */
extern double checkeig(double *err, struct vtx_data **A, double *y, int n, double lambda,
                       double *vwsqrt, double *work);

/* eigen/checkeig_ext.c */
extern double checkeig_ext(double *err, double *work, struct vtx_data **A, double *y, int n,
                           double extval, double *vwsqrt, double *gvec, double eigtol,
                           int warnings);

/* eigen/checkorth.c */
extern void checkorth(double **mat, int n, int dim);
extern void checkorth_float(float **mat, int n, int dim);

/* eigen/cksturmcnt.c */
extern void cksturmcnt(double *vec, int beg, int end, double x1, double x2, int *x1ck, int *x2ck,
                       int *numck);

/* eigen/eigensolve.c */
extern void eigensolve(struct vtx_data **graph, int nvtxs, int nedges, double maxdeg, int vwgt_max,
                       double *vwsqrt, int using_vwgts, int using_ewgts, float *term_wgts[],
                       int igeom, float **coords, double **yvecs, double *evals, int architecture,
                       int *assignment, double *goal, int solver_flag, int rqi_flag, int vmax,
                       int ndims, int mediantype, double eigtol);

/* eigen/get_extval.c */
extern void get_extval(double *alpha, double *beta, int j, double ritzval, double *s, double eigtol,
                       double wnorm_g, double sigma, double *extval, double *v, double *work1,
                       double *work2);

/* eigen/get_ritzvals.c */
extern int get_ritzvals(double *alpha, double *beta, int j, double Anorm, double *workj,
                        double *ritz, int d, int left_goodlim, int right_goodlim, double eigtol,
                        double bis_safety);

/* eigen/lanc_seconds.c */
extern double lanc_seconds(void);

/* eigen/lanczos_FO.c */
extern void lanczos_FO(struct vtx_data **A, int n, int d, double **y, double *lambda, double *bound,
                       double eigtol, double *vwsqrt, double maxdeg, int version);

/* eigen/lanczos_SO.c */
extern void lanczos_SO(struct vtx_data **A, int n, int d, double **y, double *lambda, double *bound,
                       double eigtol, double *vwsqrt, double maxdeg, int version, int cube_or_mesh,
                       int nsets, int *assignment, int *active, int mediantype, double *goal,
                       int vwgt_max);

/* eigen/lanczos_SO_float.c */
extern void lanczos_SO_float(struct vtx_data **A, int n, int d, double **y, double *lambda,
                             double *bound, double eigtol, double *vwsqrt, double maxdeg,
                             int version, int cube_or_mesh, int nsets, int *assignment, int *active,
                             int mediantype, double *goal, int vwgt_max);

/* eigen/lanczos_ext.c */
extern int lanczos_ext(struct vtx_data **A, int n, int d, double **y, double eigtol, double *vwsqrt,
                       double maxdeg, int version, double *gvec, double sigma);

/* eigen/lanczos_ext_float.c */
extern int lanczos_ext_float(struct vtx_data **A, int n, int d, double **y, double eigtol,
                             double *vwsqrt, double maxdeg, int version, double *gvec,
                             double sigma);

/* eigen/lanpause.c */
extern int lanpause(int j, int lastpause, int interval, double **q, int n, int *pausemode,
                    int version, double beta);
extern int lanpause_float(int j, int lastpause, int interval, float **q, int n, int *pausemode,
                          int version, double beta);

/* eigen/makeorthlnk.c */
extern struct orthlink       *makeorthlnk(void);
extern struct orthlink_float *makeorthlnk_float(void);

/* eigen/mkeigvecs.c */
extern void mkeigvecs(struct scanlink *scanlist, double *lambda, double *bound, int *index,
                      double *bj, int d, double *Sres_max, double *alpha, double *beta, int j,
                      double *s, double **y, int n, double **q);

/* eigen/mkscanlist.c */
extern struct scanlink *mkscanlist(int depth);

/* eigen/orthog1.c */
extern void orthog1(double *x, int beg, int end);
extern void orthog1_float(float *x, int beg, int end);

/* eigen/orthogonalize.c */
extern void orthogonalize(double *vec, int n, struct orthlink *orthlist);

/* eigen/orthogvec.c */
extern void orthogvec(double *vec1, int beg, int end, double *vec2);
extern void orthogvec_float(float *vec1, int beg, int end, float *vec2);

/* eigen/ql.c */
extern int ql(double d[], double e[], int n);

/* eigen/rqi.c */
extern void rqi(struct vtx_data **A, double **yvecs, int index, int n, double *r1, double *r2,
                double *v, double *w, double *x, double *y, double *work, double tol,
                double initshift, double *evalest, double *vwsqrt, struct orthlink *orthlist,
                int cube_or_mesh, int nsets, int *assignment, int *active, int mediantype,
                double *goal, int vwgt_max, int ndims);

/* eigen/rqi_ext.c */
extern void rqi_ext(void);

/* eigen/scale_diag.c */
extern void scale_diag(double *vec, int beg, int end, double *diag);
extern void scale_diag_float(float *vec, int beg, int end, float *diag);

/* eigen/scanmax.c */
extern void scanmax(double *vec, int beg, int end, struct scanlink **scanlist);

/* eigen/scanmin.c */
extern void scanmin(double *vec, int beg, int end, struct scanlink **scanlist);

/* eigen/solistout.c */
extern void solistout(struct orthlink **solist, int ngood, int j);
extern void solistout_float(struct orthlink_float **solist, int ngood, int j);

/* eigen/sorthog.c */
extern void sorthog(double *vec, int n, struct orthlink **solist, int ngood);
extern void sorthog_float(float *vec, int n, struct orthlink_float **solist, int ngood);

/* eigen/splarax.c */
extern void splarax(double *result, struct vtx_data **mat, int n, double *vec, double *vwsqrt,
                    double *work);
extern void splarax_float(float *result, struct vtx_data **mat, int n, float *vec, float *vwsqrt,
                          float *work);

/* eigen/sturmcnt.c */
extern int sturmcnt(double *alpha, double *beta, int j, double mu, double *p);

/* eigen/tri_solve.c */
extern void tri_solve(double *alpha, double *beta, int j, double lambda, double *v, double b,
                      double *d, double *e);

/* eigen/warnings.c */
void warnings(double *workn, struct vtx_data **A, double **y, int n, double *lambda, double *vwsqrt,
              double *Ares, double *bound, int *index, int d, int j, int maxj, double Sres_max,
              double eigtol, double *u, double Anorm, FILE *out_file);

/* graph/check_graph.c */
extern int check_graph(struct vtx_data **graph, int nvtxs, int nedges);
extern int is_an_edge(struct vtx_data *vertex, int v2, float *weight2);

/* graph/free_graph.c */
extern void free_graph(struct vtx_data **graph);

/* graph/graph_out.c */
extern void graph_out(struct vtx_data **graph, int nvtxs, int using_ewgts, char *tag,
                      char *file_name);

/* graph/mm_out.c */
extern void mm_out(struct vtx_data **graph, int nvtxs, int using_ewgts, char *tag, char *file_name);

/* graph/reformat.c */
extern int reformat(int *start, int *adjacency, int nvtxs, int *pnedges, int *vwgts, float *ewgts,
                    struct vtx_data ***pgraph);

/* graph/subgraph.c */
extern void make_subgraph(struct vtx_data **graph, struct vtx_data **subgraph, int subnvtxs,
                          int *psubnedges, int *assignment, int set, int *glob2loc, int *loc2glob,
                          int *degree, int using_ewgts);
extern void remake_graph(struct vtx_data **subgraph, int subnvtxs, int *loc2glob, int *degree,
                         int using_ewgts);

/* inertial/eigenvec2.c */
extern void evals2(double H[2][2], double *eval1, double *eval2);
extern void eigenvec2(double A[2][2], double eval, double evec[2], double *res);

/* inertial/eigenvec3.c */
extern void ch_evals3(double H[3][3], double *eval1, double *eval2, double *eval3);
extern void kramer3(double A[3][3], double b[3], double x[3]);
extern void ch_eigenvec3(double A[3][3], double eval, double evec[3], double *res);

/* inertial/inertial.c */
extern void inertial(struct vtx_data **graph, int nvtxs, int cube_or_mesh, int nsets, int igeom,
                     float **coords, int *sets, double *goal, int using_vwgts);

/* inertial/inertial1d.c */
extern void inertial1d(struct vtx_data **graph, int nvtxs, int cube_or_mesh, int nsets, float *x,
                       int *sets, double *goal, int using_vwgts);

/* inertial/inertial2d.c */
extern void inertial2d(struct vtx_data **graph, int nvtxs, int cube_or_mesh, int nsets, float *x,
                       float *y, int *sets, double *goal, int using_vwgts);

/* inertial/inertial3d.c */
extern void inertial3d(struct vtx_data **graph, int nvtxs, int cube_or_mesh, int nsets, float *x,
                       float *y, float *z, int *sets, double *goal, int using_vwgts);

/* inertial/make_subgeom.c */
extern void make_subgeom(int igeom, float **coords, float **subcoords, int subnvtxs, int *loc2glob);

/* input/check_input.c */
extern int check_input(struct vtx_data **graph, int nvtxs, int nedges, int igeom, float **coords,
                       char *graphname, int *assignment, double *goal, int architecture,
                       int ndims_tot, int mesh_dims[3], int global_method, int local_method,
                       int rqi_flag, int *vmax, int ndims, double eigtol);

/* input/input.c */

/* input/input_assign.c */

/* input/input_geom.c */

/* input/input_graph.c */

/* input/read_params.c */
extern void read_params(FILE *pfile);

extern double read_val(FILE *, int *);
extern int    read_int(FILE *, int *);

/* input/read_val.c */

/* input/reflect_input.c */

/* internal/check_internal.c */
extern void check_internal(struct vtx_data **graph, int nvtxs, struct bidint *int_list,
                           struct bidint *set_list, struct bidint *vtx_elems, int *total_vwgt,
                           int *assign, int nsets_tot);

/* internal/force_internal.c */
extern void force_internal(struct vtx_data **graph, int nvtxs, int using_ewgts, int *assign,
                           double *goal, int nsets_tot, int npasses_max);

/* internal/improve_internal.c */
extern int improve_internal(struct vtx_data **graph, int nvtxs, int *assign, double *goal,
                            struct bidint *int_list, struct bidint *set_list,
                            struct bidint *vtx_elems, int set1, int *locked, int *nlocked,
                            int using_ewgts, int vwgt_max, int *total_vwgt);

/* klspiff/bilistops.c */
extern void add2bilist(struct bilist *lptr, struct bilist **list);
extern void removebilist(struct bilist *lptr, struct bilist **list);
extern void movebilist(struct bilist *lptr, struct bilist **oldlist, struct bilist **newlist);

/* klspiff/buckets.c */
extern void bucketsorts(struct vtx_data **graph, int nvtxs, struct bilist ****buckets,
                        struct bilist **listspace, int **dvals, int *sets, float *term_wgts[],
                        int maxdval, int nsets, int parity, int (*hops)[MAXSETS], int *bspace,
                        int list_length, int npass, int using_ewgts);

/* klspiff/buckets1.c */
extern void bucketsort1(struct vtx_data **graph, int vtx, struct bilist ****buckets,
                        struct bilist **listspace, int **dvals, int *sets, float *term_wgts[],
                        int maxdval, int nsets, int (*hops)[MAXSETS], int using_ewgts);

/* klspiff/buckets_bi.c */
extern void bucketsorts_bi(struct vtx_data **graph, int nvtxs, struct bilist ****buckets,
                           struct bilist **listspace, int **dvals, int *sets, float *term_wgts[],
                           int maxdval, int nsets, int parity, int (*hops)[MAXSETS], int *bspace,
                           int list_length, int npass, int using_ewgts);

/* klspiff/coarsen_kl.c */
extern void coarsen_kl(struct vtx_data **graph, int nvtxs, int nedges, int using_vwgts,
                       int using_ewgts, float *term_wgts[], int igeom, float **coords, int vwgt_max,
                       int *assignment, double *goal, int architecture, int (*hops)[MAXSETS],
                       int solver_flag, int ndims, int nsets, int vmax, int mediantype,
                       int mkconnected, double eigtol, int nstep, int step, int **pbndy_list,
                       double *weights, int give_up);

/* klspiff/compress_ewgts.c */
extern void compress_ewgts(struct vtx_data **graph, int nvtxs, int nedges, double ewgt_max,
                           int using_ewgts);
extern void restore_ewgts(struct vtx_data **graph, int nvtxs);

/* klspiff/count_weights.c */
extern void count_weights(struct vtx_data **graph, int nvtxs, int *sets, int nsets, double *weights,
                          int using_vwgts);

/* klspiff/kl_init.c */
extern int kl_init(struct bilist *****bucket_ptrs, struct bilist ***listspace, int ***dvals,
                   int ***tops, int nvtxs, int nsets, int maxchange);

/* klspiff/kl_output.c */
extern void pbuckets(struct bilist ****buckets, struct bilist **listspace, int maxdeg, int nsets);
extern void p1bucket(struct bilist **bucket, struct bilist *lptr, int maxdeg);

/* klspiff/klspiff.c */
extern void klspiff(struct vtx_data **graph, int nvtxs, int *sets, int nsets, int (*hops)[MAXSETS],
                    double *goal, float *term_wgts[], int max_dev, double maxdeg, int using_ewgts,
                    int **bndy_list, double *weights);

/* klspiff/make_bndy_list.c */
extern void make_bndy_list(struct vtx_data **graph, struct bilist *movelist,
                           struct bilist ****buckets, struct bilist **listspace, int *sets,
                           int nsets, int *bspace, int **tops, int **bndy_list);

/* klspiff/make_kl_list.c */
extern int make_kl_list(struct vtx_data **graph, struct bilist *movelist, struct bilist ****buckets,
                        struct bilist **listspace, int *sets, int nsets, int *bspace, int **dvals,
                        int maxdval);

/* klspiff/nway_kl.c */
extern int nway_kl(struct vtx_data **graph, int nvtxs, struct bilist ****buckets,
                   struct bilist **listspace, int **tops, int **dvals, int *sets, int maxdval,
                   int nsets, double *goal, float *term_wgts[], int (*hops)[MAXSETS], int max_dev,
                   int using_ewgts, int **bndy_list, double *startweight);

/* klvspiff/bpm_improve.c */
extern void bpm_improve(struct vtx_data **graph, int *sets, double *goal, int max_dev,
                        int **bndy_list, double *weights, int using_vwgts);

/* klvspiff/bucketsv.c */
extern void bucketsortsv(struct vtx_data **graph, int nvtxs, struct bilist **lbuckets,
                         struct bilist **rbuckets, struct bilist *llistspace,
                         struct bilist *rlistspace, int *ldvals, int *rdvals, int *sets,
                         int maxdval, int parity, int *bspace, int list_length);

/* klvspiff/clear_dvals.c */
extern void clear_dvals(struct vtx_data **graph, int nvtxs, int *ldvals, int *rdvals, int *bspace,
                        int list_length);

/* klvspiff/coarsen_klv.c */
extern void coarsen_klv(struct vtx_data **graph, int nvtxs, int nedges, int using_vwgts,
                        int using_ewgts, float *term_wgts[], int igeom, float **coords,
                        int vwgt_max, int *assignment, double *goal, int architecture,
                        int (*hops)[MAXSETS], int solver_flag, int ndims, int nsets, int vmax,
                        int mediantype, int mkconnected, double eigtol, int nstep, int step,
                        int **pbndy_list, double *weights, int give_up);
extern void print_sep_size(int *list, struct vtx_data **graph);

/* klvspiff/countup_vtx_sep.c */
extern void countup_vtx_sep(struct vtx_data **graph, int nvtxs, int *sets);

/* klvspiff/find_bndy.c */
extern int find_bndy(struct vtx_data **graph, int nvtxs, int *assignment, int new_val,
                     int **pbndy_list);
extern int find_side_bndy(struct vtx_data **graph, int nvtxs, int *assignment, int side,
                          int new_val, int **pbndy_list);

/* klvspiff/flatten.c */
extern int  flatten(struct vtx_data **graph, int nvtxs, int nedges, struct vtx_data ***pcgraph,
                    int *pcnvtxs, int *pcnedges, int **pv2cv, int using_ewgts, int igeom,
                    float **coords, float **ccoords);
extern void find_flat(struct vtx_data **graph, int nvtxs, int *pcnvtxs, int *v2cv);
extern int  SameStructure(int node1, int node2, struct vtx_data **graph, int *scatter);

/* klvspiff/flow.c */
extern void wbpcover(int n_left, int n_right, int *pointers, int *indices, int *vweight,
                     int *psep_size, int *psep_weight, int **psep_nodes);
extern void confirm_cover(int n_left, int n_right, int *pointers, int *indices, int *flow,
                          int *vweight, int *resid, int sep_size, int *sep_nodes);
extern int  count_flow(int n_left, int n_right, int *pointers, int *flow);
extern void count_resid(int n_left, int n_right, int *resid, int *vweight, int *marked);
extern void check_resid(int n_left, int n_right, int *vweight, int *resid, int *pointers,
                        int *indices, int *flow);

/* klvspiff/klv_init.c */
extern int klv_init(struct bilist ***lbucket_ptr, struct bilist ***rbucket_ptr,
                    struct bilist **llistspace, struct bilist **rlistspace, int **ldvals,
                    int **rdvals, int nvtxs, int maxchange);

/* klvspiff/klvspiff.c */
extern void klvspiff(struct vtx_data **graph, int nvtxs, int *sets, double *goal, int max_dev,
                     int **bndy_list, double *weights);

/* klvspiff/make_bpgraph.c */
extern void make_bpgraph(struct vtx_data **graph, int *sets, int *bndy_list, int sep_size,
                         int set_match, int **ppointers, int **pindices, int **pvweight,
                         int **ploc2glob, int *pnleft, int *pnright, int using_vwgts);
extern void check_bpgraph(int n_left, int n_right, int *pointers, int *indices);
extern void print_bpgraph(int nleft, int nright, int *pointers, int *indices, int *vwgts);

/* klvspiff/make_sep_list.c */
extern int make_sep_list(int *bspace, int list_length, int *sets);

/* klvspiff/matching.c */
extern void bpcover(int n_left, int n_right, int *pointers, int *indices, int *sep_size,
                    int *sep_nodes);
extern void confirm_match(int n_left, int n_right, int *pointers, int *indices, int *matching,
                          int sep_size, int *sep_nodes);
extern int  match_size(int *matching, int nleft);

/* klvspiff/nway_klv.c */
extern int nway_klv(struct vtx_data **graph, int nvtxs, struct bilist **lbuckets,
                    struct bilist **rbuckets, struct bilist *llistspace, struct bilist *rlistspace,
                    int *ldvals, int *rdvals, int *sets, int maxdval, double *goal, int max_dev,
                    int **bndy_list, double *weightsum);

/* main/interface.c */
extern int INTER_FACE(int nvtxs, int *start, int *adjacency, int *vwgts, float *ewgts, float *x,
                      float *y, float *z, char *outassignname, char *outfilename, int *assignment,
                      int architecture, int ndims_tot, int mesh_dims[3], double *goal,
                      int global_method, int local_method, int rqi_flag, int vmax, int ndims,
                      double eigtol, long seed);

/* main/main.c */
extern int main(void);

/* main/user_params.c */

/* misc/count.c */
extern void count(struct vtx_data **graph, int nvtxs, int *sets, int nsets, int (*hops)[MAXSETS],
                  int dump, int using_ewgts);

/* misc/countup.c */
void countup(struct vtx_data **graph,        /* graph data structure */
             int               nvtxs,        /* number of vtxs in graph */
             int              *assignment,   /* set number of each vtx (length nvtxs+1) */
             int               ndims,        /* number of cuts at each level */
             int               architecture, /* what's the target parallel machine? */
             int               ndims_tot,    /* total number of hypercube dimensions */
             int               mesh_dims[3], /* extent of mesh in each dimension */
             int               print_lev,    /* level of output */
             FILE             *outfile,      /* output file if not NULL */
             int               using_ewgts   /* are edge weights being used? */
);

/* misc/countup_cube.c */
void countup_cube(struct vtx_data **graph,      /* graph data structure */
                  int               nvtxs,      /* number of vtxs in graph */
                  int              *assignment, /* set number of each vtx (length nvtxs+1) */
                  int               ndims,      /* number of cuts at each level */
                  int               ndims_tot,  /* total number of divisions of graph */
                  int               print_lev,  /* level of output */
                  FILE             *outfile,    /* output file if not NULL */
                  int               using_ewgts /* are edge weights being used? */
);

/* misc/countup_mesh.c */
void countup_mesh(struct vtx_data **graph,        /* graph data structure */
                  int               nvtxs,        /* number of vtxs in graph */
                  int              *assignment,   /* set number of each vtx (length nvtxs+1) */
                  int               mesh_dims[3], /* extent of mesh in each dimension */
                  int               print_lev,    /* level of output */
                  FILE             *outfile,      /* output file if not NULL */
                  int               using_ewgts   /* are edge weights being used? */
);

/* misc/define_subcubes.c */
extern int define_subcubes(int nsets_real, int ndims_tot, int ndims, struct set_info *set,
                           struct set_info *set_info, int *subsets, int inert, int *pstriping,
                           int hop_mtx_special[MAXSETS][MAXSETS]);

/* misc/define_submeshes.c */
extern int define_submeshes(int nsets, int cube_or_mesh, int *mesh_dims, struct set_info *set,
                            struct set_info *set_info, int *subsets, int inert, int *striping,
                            int *dir, int hop_mtx_special[MAXSETS][MAXSETS]);

/* misc/divide_procs.c */
extern int divide_procs(int architecture, int ndims, int ndims_tot, struct set_info *info_set,
                        struct set_info *divide_set, int *subsets, int inert, int *pndims_real,
                        int *pnsets_real, int *pstriping, int *cut_dirs, int *mesh_dims,
                        int hops_special[][MAXSETS]);

/* misc/find_maxdeg.c */
extern double find_maxdeg(struct vtx_data **graph, int nvtxs, int using_ewgts, float *pmax_ewgt);

/* misc/make_maps.c */
extern int  make_maps(int *setlists, int *list_ptrs, int set, int *glob2loc, int *loc2glob);
extern void make_maps2(int *assignment, int nvtxs, int set, int *glob2loc, int *loc2glob);

/* misc/make_setlists.c */
extern void make_setlists(int *setlists, int *list_ptrs, int nsets, int *subsets, int *subassign,
                          int *loc2glob, int subnvtxs, int first);

/* misc/make_subgoal.c */
extern void make_subgoal(double *goal, double *subgoal, int nsets, int cube_or_mesh, int nsets_tot,
                         int mesh_dims[3], int set, double sub_vwgt_sum);

/* misc/make_term_props.c */
extern void make_term_props(struct vtx_data **graph, int sub_nvtxs, int *loc2glob, int *assignment,
                            int architecture, int ndims_tot, int ndims, struct set_info *set_info,
                            int setnum, int nsets, int set_max, int *subsets, float *term_wgts[],
                            int using_ewgts);

/* misc/merge_goals.c */
extern void merge_goals(double *goal, double *merged_goal, struct set_info *set_info, int *subsets,
                        int nsets, int ndims_tot, int cube_or_mesh, int mesh_dims[3],
                        double vwgt_sum);

/* misc/perturb.c */
extern void perturb_init(int n);
extern void perturb_clear(void);
extern void perturb(double *result, double *vec);
extern void perturb_float(float *result, float *vec);

/* misc/sequence.c */
extern void sequence(struct vtx_data **graph, int nvtxs, int nedges, int using_ewgts,
                     double *vwsqrt, int solver_flag, int rqi_flag, int vmax, double eigtol);

/* misc/simple_part.c */
extern void simple_part(struct vtx_data **graph, int nvtxs, int *sets, int nsets, int simple_type,
                        double *goal);

/* misc/time_kernels.c */
extern void time_kernels(struct vtx_data **A, int n, double *vwsqrt);

/* misc/timing.c */
extern void clear_timing(void);
extern void time_out(FILE *outfile);

/* optimize/determinant.c */
extern double determinant(double M[3][3], int ndims);

/* optimize/func2d.c */
extern double func2d(double coeffs[5], double theta);
extern double grad2d(double coeffs[5], double theta);
extern double hess2d(double coeffs[5]);

/* optimize/func3d.c */
extern double func3d(double coeffs[15], double theta, double phi, double gamma2);
extern void   grad3d(double coeffs[15], double grad[3], double theta, double phi, double gamma2);
extern void   hess3d(double coeffs[15], double hess[3][3]);
extern double constraint(double *coeffs2);
extern void   gradcon(double *coeffs2, double grad[3]);
extern void   hesscon(double *coeffs2, double hess[3][3]);

/* optimize/opt2d.c */
extern double opt2d(struct vtx_data **graph, double **yvecs, int nvtxs, int nmyvtxs);

/* optimize/opt3d.c */
extern void opt3d(struct vtx_data **graph, double **yvecs, int nvtxs, int nmyvtxs, double *vwsqrt,
                  double *ptheta, double *pphi, double *pgamma, int using_vwgts);

/* refine_map/compute_cube_edata.c */
extern double compute_cube_edata(struct refine_edata *edata, struct refine_vdata *vdata,
                                 int nsets_tot, struct vtx_data **comm_graph, int *node2vtx);

/* refine_map/compute_cube_vdata.c */
extern void compute_cube_vdata(struct refine_vdata *vdata, struct vtx_data **comm_graph, int vtx,
                               int mask, int *vtx2node);

/* refine_map/compute_mesh_edata.c */
extern double compute_mesh_edata(struct refine_edata *edata, struct refine_vdata *vdata,
                                 int mesh_dims[3], struct vtx_data **comm_graph, int *node2vtx);

/* refine_map/compute_mesh_vdata.c */
extern void compute_mesh_vdata(struct refine_vdata *vdata, struct vtx_data **comm_graph, int vtx,
                               int *vtx2node, int mesh_dims[3], int dim);

/* refine_map/find_edge_cube.c */
extern struct refine_edata *find_edge_cube(int node, int dim, struct refine_edata *edata,
                                           int nsets_tot);

/* refine_map/find_edge_mesh.c */
extern struct refine_edata *find_edge_mesh(int vertex, int dim, struct refine_edata *edata,
                                           int *mesh_dims, int *vtx2node);

/* refine_map/init_cube_edata.c */
extern void init_cube_edata(struct refine_edata *edata, int node1, int dim, int mask);

/* refine_map/init_mesh_edata.c */
extern void init_mesh_edata(struct refine_edata *edata, int mesh_dims[3]);

/* refine_map/make_comm_graph.c */
extern int make_comm_graph(struct vtx_data ***pcomm_graph, struct vtx_data **graph, int nvtxs,
                           int using_ewgts, int *assign, int nsets_tot);

/* refine_map/refine_cube.c */
extern int refine_cube(struct vtx_data **comm_graph, int ndims_tot, double maxdesire, int *vtx2node,
                       int *node2vtx);

/* refine_map/refine_map.c */
extern void refine_map(struct vtx_data **graph, int nvtxs, int using_ewgts, int *assign,
                       int cube_or_mesh, int ndims_tot, int mesh_dims[3]);

/* refine_map/refine_mesh.c */
extern int refine_mesh(struct vtx_data **comm_graph, int cube_or_mesh, int mesh_dims[3],
                       double maxdesire, int *vtx2node, int *node2vtx);

/* refine_map/update_cube_edata.c */
extern void update_cube_edata(int vertex, int dim, struct refine_edata *edata,
                              struct refine_vdata *vdata, struct vtx_data **comm_graph,
                              int *node2vtx, int *vtx2node, int nsets_tot, double *best_desire,
                              int imax, struct refine_edata **desire_ptr);

/* refine_map/update_cube_vdata.c */
extern void update_cube_vdata(int old_side, int mask, int neighbor_node, double ewgt,
                              struct refine_vdata *vdata);

/* refine_map/update_mesh_edata.c */
extern void update_mesh_edata(int vertex, int dim, struct refine_edata *edata,
                              struct refine_vdata *vdata, struct vtx_data **comm_graph,
                              int mesh_dims[3], int *node2vtx, int *vtx2node, double *best_desire,
                              int imax, struct refine_edata **desire_ptr);

/* refine_map/update_mesh_vdata.c */
extern void update_mesh_vdata(int old_loc, int new_loc, int dim, double ewgt,
                              struct refine_vdata *vdata, int mesh_dims[3], int neighbor,
                              int *vtx2node);

/* refine_part/kl_refine.c */
extern int kl_refine(struct vtx_data **graph, struct vtx_data **subgraph, struct bilist *set_list,
                     struct bilist *vtx_elems, int *new_assign, int set1, int set2, int *glob2loc,
                     int *loc2glob, int *sub_assign, int *old_sub_assign, int *degrees,
                     int using_ewgts, int (*hops)[MAXSETS], double *goal, int *sizes,
                     float *term_wgts[], int architecture, int mesh_dims[3]);

/* refine_part/make_maps_ref.c */
extern void make_maps_ref(struct vtx_data **graph, struct bilist *set_list,
                          struct bilist *vtx_elems, int *assignment, int *sub_assign, int set1,
                          int set2, int *glob2loc, int *loc2glob, int *psub_nvtxs, int *pvwgt_max,
                          int *pvwgt_sum1, int *pvwgt_sum2);

/* refine_part/make_terms_ref.c */
extern void make_terms_ref(struct vtx_data **graph, int using_ewgts, int subnvtxs, int *loc2glob,
                           int set0, int set1, int *assignment, int architecture, int mesh_dims[3],
                           float *term_wgts[]);

/* refine_part/refine_part.c */
extern int refine_part(struct vtx_data **graph, int nvtxs, int using_ewgts, int *assign,
                       int architecture, int ndims_tot, int mesh_dims[3], double *goal);

/* submain/balance.c */
extern void balance(struct vtx_data **graph, int nvtxs, int nedges, int using_vwgts,
                    int using_ewgts, double *vwsqrt, int igeom, float **coords, int *assignment,
                    double *goal, int architecture, int ndims_tot, int *mesh_dims,
                    int global_method, int local_method, int rqi_flag, int vmax, int ndims,
                    double eigtol, int (*hops)[MAXSETS]);

/* submain/divide.c */
extern void chaco_divide(struct vtx_data **graph, int nvtxs, int nedges, int using_vwgts,
                         int using_ewgts, double *vwsqrt, int igeom, float **coords,
                         int *assignment, double *goal, int architecture, float *term_wgts[],
                         int global_method, int local_method, int rqi_flag, int vmax, int ndims,
                         double eigtol, int (*hop_mtx)[MAXSETS], int nsets, int striping);

/* submain/submain.c */
extern int submain(struct vtx_data **graph, int nvtxs, int nedges, int using_vwgts, int using_ewgts,
                   int igeom, float **coords, char *outassignname, char *outfilename,
                   int *assignment, double *goal, int architecture, int ndims_tot, int mesh_dims[3],
                   int global_method, int local_method, int rqi_flag, int vmax, int ndims,
                   double eigtol, long seed);

/* symmlq/aprod.c */
extern int aprod(long *lnvtxs, double *x, double *y, double *dA, double *vwsqrt, double *work,
                 double *dorthlist);

/* symmlq/msolve.c */
extern int msolve(long int *nvtxs, double *x, double *y);

/* symmlq/symmlq.c */
extern int symmlq(long int *n, double *b, double *r1, double *r2, double *v, double *w, double *x,
                  double *y, double *work, long int *checka, long int *goodb, long int *precon,
                  double *shift, long int *nout, long int *itnlim, double *rtol, long int *istop,
                  long int *itn, double *anorm, double *acond, double *rnorm, double *ynorm,
                  double *a, double *vwsqrt, double *orthlist, double *macheps, double *normxlim,
                  long int *itnmin);

/* symmlq/symmlqblas.c */
extern int chdaxpy(long int *n, double *da, double *dx, long int *incx, double *dy, long int *incy);
extern int chdcopy(long int *n, double *dx, long int *incx, double *dy, long int *incy);
extern double ch_ddot(long int *n, double *dx, long int *incx, double *dy, long int *incy);
extern double chdnrm2(long int *n, double *dx, long int *incx);

/* tinvit/tinvit.c */
extern int tinvit(long int *nm, long int *n, double *d, double *e, double *e2, long int *m,
                  double *w, long int *ind, double *z, long int *ierr, double *rv1, double *rv2,
                  double *rv3, double *rv4, double *rv6);

/* util/affirm.c */
extern int affirm(char *prompt);

/* util/array_alloc_2D.c */
extern void *array_alloc_2D_ret(size_t dim1, size_t dim2, size_t size);

/* util/bail.c */
extern void bail(char *msg, int status);

/* util/bit_reverse.c */
extern int bit_reverse(int val, int nbits);

/* util/checkpnt.c */
extern void checkpnt(char *tag);

/* util/cpvec.c */
extern void cpvec(double *copy, int beg, int end, double *vec);
extern void float_to_double(double *copy, int beg, int end, float *vec);
extern void double_to_float(float *copy, int beg, int end, double *vec);

/* util/dot.c */
extern double dot(double *vec1, int beg, int end, double *vec2);
extern double dot_float(float *vec1, int beg, int end, float *vec2);

/* util/doubleout.c */
extern void doubleout(double number, int mode);
extern void doubleout_file(FILE *outfile, double number, int mode);

/* util/gray.c */
extern int gray(int i);

/* util/input_int.c */
extern int input_int(void);

/* util/machine_params.c */
extern void machine_params(double *double_epsilon, double *double_max);

/* util/makevwsqrt.c */
extern void makevwsqrt(double *vwsqrt, struct vtx_data **graph, int nvtxs);
extern void make_subvector(double *vec, double *subvec, int subnvtxs, int *loc2glob);

/* util/mergesort.c */
extern void ch_mergesort(double *vals, int nvals, int *indices, int *space);

/* util/mkvec.c */
extern double *mkvec(int nl, int nh);
extern double *mkvec_ret(int nl, int nh);
extern void    frvec(double *v, int nl);
extern float  *mkvec_float(int nl, int nh);
extern float  *mkvec_ret_float(int nl, int nh);
extern void    frvec_float(float *v, int nl);

/* util/norm.c */
extern double ch_norm(double *vec, int beg, int end);
extern double norm_float(float *vec, int beg, int end);

/* util/normalize.c */
extern double ch_normalize(double *vec, int beg, int end);
extern double sign_normalize(double *vec, int beg, int end, int k);
extern double normalize_float(float *vec, int beg, int end);

/* util/random.c */
extern long   init_rand_port(long seed);
extern long   get_init_rand_port(void);
extern long   genr_rand_port(long init_rand);
extern long   rand_port(void);
extern double rand_rect_port(void);
extern int    main(void);

/* util/randomize.c */
extern void   randomize(int *array, int n);
extern double drandom(void);
extern void   setrandom(long int seed);

/* util/scadd.c */
extern void scadd(double *vec1, int beg, int end, double fac, double *vec2);
extern void scadd_float(float *vec1, int beg, int end, float fac, float *vec2);
extern void scadd_mixed(double *vec1, int beg, int end, double fac, float *vec2);

/* util/seconds.c */
extern double seconds(void);

/* util/setvec.c */
extern void setvec(double *vec, int beg, int end, double setval);
extern void setvec_float(float *vec, int beg, int end, float setval);

/* util/shell_sort.c */
extern void sort_double(int count, double ra[]);
extern void shell_sort(int n, double *arr);

/* util/smalloc.c */
extern void *smalloc(size_t n);
extern void *smalloc_ret(size_t n);
extern void *srealloc(void *ptr, size_t n);
extern void *srealloc_ret(void *ptr, size_t n);
extern void  sfree(void *ptr);
extern void  smalloc_stats(void);
extern int   smalloc_num(void);

/* util/strout.c */
extern void strout(char *msg);

/* util/tri_prod.c */
extern double tri_prod(double *v1, double *v2, double *v3, double *wsqrt, int n);

/* util/true_or_false.c */
extern char *true_or_false(int flag);

/* util/update.c */
extern void update(double *vec1, int beg, int end, double *vec2, double fac, double *vec3);
extern void update_float(float *vec1, int beg, int end, float *vec2, float fac, float *vec3);

/* util/vecout.c */
extern void vecout(double *vec, int beg, int end, char *tag, char *file_name);
extern void vecnorm(double *vec, int beg, int end);

/* util/vecran.c */
extern void vecran(double *vec, int beg, int end);
extern void vecran_float(float *vec, int beg, int end);

/* util/vecscale.c */
extern void vecscale(double *vec1, int beg, int end, double alpha, double *vec2);
extern void vecscale_float(float *vec1, int beg, int end, float alpha, float *vec2);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
