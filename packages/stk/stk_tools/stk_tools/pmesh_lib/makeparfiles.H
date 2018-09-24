#ifndef SALINAS_MAKEPARFILE_H
#define SALINAS_MAKEPARFILE_H

namespace stk_tools {

void MakeParFile(
		 const int& my_proc_id,
		 const int& num_procs,
		 const int& ncuts_x,
		 const int& ncuts_y,
		 const int& ncuts_z,
		 const int& nelem_per_edge,
		 const double& lenx,
		 const double& leny,
		 const double& lenzi,
		 const int& numraid,
		 const char* rootdir,
		 const char* subdir );

}

#endif

