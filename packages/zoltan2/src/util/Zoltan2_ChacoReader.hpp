// @HEADER
// ***********************************************************************
//      Copyright message goes here.
// ***********************************************************************

/*! \file ChacoReader.hpp
 *  \brief declarations for chaco file reader
 */

namespace Zoltan2 {

int chaco_input_graph( FILE *fin, char *inname, int **start,
int **adjacency, int  *nvtxs, int  *vwgt_dim, float **vweights,
int  *ewgt_dim, float **eweights );

int chaco_input_geom( FILE *fingeom, char *geomname, int nvtxs,
int  *igeom, float **x, float **y, float **z);

} // namespace Zoltan2
