extern void     *my_malloc(int size);
extern int      in_box(double pt[3], double lower[3], double upper[3]);
extern void     bounds_to_origin_size(double min[3], double max[3],
                                      double origin[3], double size[3]);
extern void     bounds_to_origin(double min[3], double max[3], 
				 double origin[3]);
extern void     child_bounds(double pmin[3], double pmax[3], double porigin[3],
			     int cnum, double cmin[3], double cmax[3]);
extern int      child_which(double origin[3], double point[3]);
extern double   dist_point_box(double point[3], double min[3], double max[3]);

