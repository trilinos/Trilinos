#ifdef PETRA_MPI
#ifndef _GSCOMM_COMM_H_
#define _GSCOMM_COMM_H_

//
// This plan is constructed for gather/scatter operations
// OO version derived from Zoltan's LB_comm_create
//

#include "Petra_Petra.h"
#include "GSComm_Plan.h"

class GSComm_Comm {
    
  public:

    GSComm_Comm();

    ~GSComm_Comm();

    bool Do( GSComm_Plan & comm_plan,
	     const int & tag,
	     char * export_objs,
	     const int & obj_size,
	     char * import_objs );

    bool Do_Posts( GSComm_Plan & comm_plan,
	               const int & tag,
	               char * export_objs,
	               const int & obj_size,
	               char * import_objs );

    bool Do_Waits( GSComm_Plan & comm_plan,
	           const int & tag,
	           char * export_objs,
	           const int & obj_size,
	           char * import_objs );

    bool DoReverse( GSComm_Plan & comm_plan,
	            const int & tag,
	            char * export_objs,
	            const int & obj_size,
	            char * import_objs );

    bool DoReverse_Posts( GSComm_Plan & comm_plan,
	                      const int & tag,
	                      char * export_objs,
	                      const int & obj_size,
	                      char * import_objs );

    bool DoReverse_Waits( GSComm_Plan & comm_plan,
	                  const int & tag,
	                  char * export_objs,
	                  const int & obj_size,
	                  char * import_objs );

  private:

    char * recv_array_;
    char * send_array_;

    GSComm_Plan * comm_plan_reverse_;

};

// GSComm_Comm constructor
inline GSComm_Comm::GSComm_Comm()
	: recv_array_(0),
	  send_array_(0),
	  comm_plan_reverse_(0)
{
}

// GSComm_Comm destructor
inline GSComm_Comm::~GSComm_Comm()
{
  if( send_array_ != 0 ) delete [] send_array_;

  if( comm_plan_reverse_ != 0 ) delete comm_plan_reverse_;
}

// GSComm_Comm Do method
inline bool GSComm_Comm::Do( GSComm_Plan & comm_plan,
		            const int & tag,
		            char * export_objs,
		            const int & obj_size,
		            char * import_objs )
{
  bool comm_flag = true;

  Do_Posts( comm_plan, tag, export_objs, obj_size, import_objs );
  Do_Waits( comm_plan, tag, export_objs, obj_size, import_objs );

  return comm_flag;
}

// GSComm_Comm DoReverse method
inline bool GSComm_Comm::DoReverse( GSComm_Plan & comm_plan,
		                   const int & tag,
		                   char * export_objs,
		                   const int & obj_size,
		                   char * import_objs )
{
  bool comm_flag = true;

  DoReverse_Posts( comm_plan, tag, export_objs, obj_size, import_objs );
  DoReverse_Waits( comm_plan, tag, export_objs, obj_size, import_objs );

  return comm_flag;
}

#endif /* _GSCOMM_COMM_H_ */
#endif /* PETRA_MPI */
