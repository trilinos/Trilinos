extern void   set_method(double method_number);
extern void   *my_malloc(int size);
extern int    in_box(double pt[3], double lower[3], double upper[3]);
extern void   bounds_to_origin_size(double min[3], double max[3],
				    double origin[3], double size[3]);
extern void   bounds_to_origin(double min[3], double max[3], 
			       double origin[3]);
extern void   child_bounds_wrapper(pOctant oct, int input, 
				   double cmin[3], double cmax[3]);
extern void   child_bounds(double pmin[3], double pmax[3], double porigin[3],
			   int cnum, double cmin[3], double cmax[3]);

extern int    compare(unsigned *x, unsigned *y);
extern int    hilbert_bounds(double min[3], double max[3], int cnum);
extern int    hilbert2d_bounds(double min[3], double max[3], int cnum);

extern int    child_which_wrapper(pOctant oct, double point[3]);
extern int    child_which(double origin[3], double point[3]);
extern double dist_point_box(double point[3], double min[3], double max[3]);
extern int    convert_to_gray(int input);
extern int    convert_from_gray(int input);

extern int    convert_to_hilbert(int n, int o);
extern int    convert_from_hilbert(int n, int o);
extern int    child_orientation(int o, int cnum);
extern int    change_to_hilbert2d(double min[3], double max[3], 
				  double origin[3], int cnum);
extern int    change_to_hilbert(double min[3], double max[3], double origin[3],
				int cnum);
