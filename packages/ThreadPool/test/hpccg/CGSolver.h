
struct cgsolve_data {
  int     nRow ; 
  int   * A_pc ; 
  int   * A_ia ; 
  float * A_a ; 
  int     max_iter ; 
  int     print_iter ; 
  float   tolerance ; 

  int     np ; 
  int     ip ; 
  int   * recv_pc ; 
  int   * send_pc ; 
  int   * send_id ; 
}; 

void cgsolve_set_lhs( const struct cgsolve_data * data ,
                      const double * const x ,
                            double * const b );

void cgsolve( const struct cgsolve_data * data ,
              const double * const b ,
                    double * const x ,
                    int    * const iter_count ,
                    double * const norm_resid ,
                    double * const dt_mxv ,
                    double * const dt_axpby ,
                    double * const dt_dot );

