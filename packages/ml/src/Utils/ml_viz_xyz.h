#ifndef ML_VIZ_XYZ_H
#define ML_VIZ_XYZ_H

#include "ml_include.h"
#include "ml_viz_stats.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

int ML_Aggregate_VisualizeXYZ( ML_Aggregate_Viz_Stats info,
			      char base_filename[],
			      ML_Comm *comm,
			      double * vector);

int ML_PlotXYZ(int Npoints, double* x, double* y, double* z,
	       char base_filename[],
	       USR_COMM comm, double * vector);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* #ifndef ML_VIZ_XYZ_H */

