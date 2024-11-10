/********************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/********************************************************************* */
/*          visualization routines                                     */
/********************************************************************* */

#ifndef __MLVIZOPENDX__
#define __MLVIZOPENDX__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_viz_stats.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

 extern int ML_Aggregate_VisualizeWithOpenDX( ML_Aggregate_Viz_Stats info,
					      char base_filename[],
					      ML_Comm * comm);
 extern int ML_Aggregate_VisualizeXYZ( ML_Aggregate_Viz_Stats info,
				      char base_filename[],
				      ML_Comm * comm, double * vector);


#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* #ifndef __MLVIZOPENDX__ */
