// $Id$

#ifndef pamgen_mesh_specification_LT_H
#define pamgen_mesh_specification_LT_H

#include "pamgen_element_dictionary.h"
#include <string>
#include <sstream>
#include <ios>
using namespace PAMGEN_NEVADA;
namespace ms_lt{


/*****************************************************************************/
class Mesh_Specification
/*****************************************************************************/
// This abstract class represents a mesh specification.
//
// A note on numbering.  Most file formats for specifying meshes were
// developed in the FORTRAN world, where the default convention is for array 
// indices to begin at 1 rather than 0.  In addition, when a mesh is
// decomposed for parallel computations, its topological entities are
// assigned a global number.  As a result, there is considerable potential
// for confusion about numbering.  Unless otherwise specified, local indexing 
// starts at 0.  When indices returned by a function are based on the 
// FORTRAN convention (indices starting at 1), this will be noted.  Global 
// numbering will always be explicitly noted when not obvious.
{
  public:

  enum MSPPA {
    ELMT_NODE_LINKAGE,
    NODE_SET_NODES,
    SIDE_SET_ELEMENTS,
    SIDE_SET_FACES,
    SIDE_SET_NODES,
    SIDE_SET_NODE_COUNTER,
    COMM_NODE_IDS,
    COMM_NODE_PROC_IDS,
    COMM_ELEM_IDS,
    COMM_ELEM_PROC_IDS,
    COMM_SIDE_IDS,
    NUM_MSPPA
  };

  enum MSPPDA{
    ATTRIBUTES,
    NODE_SET_DF,
    SIDE_SET_DF,
    NUM_MSPPDA
  };

  enum MSPSA{
    INFO_STRINGS,
    COORDINATE_NAMES,
    ELEMENT_TYPES,
    NUM_MSPSA
  };

  enum MSPA {
    ELEM_ORDER_MAP=0,
    BLOCK_ID,
    ELEMENTS_IN_BLOCK,
    NODES_PER_ELEMENT,
    ELEMENT_ATTRIBUTES,
    NODE_SET_ID,
    NUM_NODES_IN_NODE_SET,
    NUM_DF_IN_NODE_SET,
    SIDE_SET_ID,
    NUM_ELEMENTS_IN_SIDE_SET,
    NUM_NODES_IN_SIDE_SET,
    NUM_DF_IN_SIDE_SET,
    GLOBAL_ELEMENT_NUMBERS,
    GLOBAL_NODE_NUMBERS,
    NBR_PROC_LIST,
    ELEM_BLK_IDS_GLOBAL,
    ELEM_BLK_CNTS_GLOBAL,
    NS_IDS_GLOBAL,
    NS_CNTS_GLOBAL,
    NS_DF_CNTS_GLOBAL,
    SS_IDS_GLOBAL,
    SS_CNTS_GLOBAL,
    SS_DF_CNTS_GLOBAL,
    INTERNAL_ELEMENTS,
    BORDER_ELEMENTS,
    INTERNAL_NODES,
    BORDER_NODES,
    EXTERNAL_NODES,
    NODE_CMAP_NODE_CNTS,
    NODE_CMAP_IDS,
    ELEM_CMAP_ELEM_CNTS,
    ELEM_CMAP_IDS,
    NUM_MSPA};

  enum MSIA {
    DIM=0,
    PROC_ID,
    NUM_QA_RECORDS,
    NUM_INFO_RECORDS,
    NUM_TOTAL_PROC,
    NUM_PROC_IN_FILE,
    NUM_NODES,
    NUM_ELEMENTS,
    NUM_EDGES,
    NUM_FACES,
    NUM_BLOCKS,
    NUM_NODE_SETS,
    NUM_SIDE_SET_NODES,
    NUM_SIDE_SETS,
    NUM_NODES_GLOBAL,
    NUM_ELEMS_GLOBAL,
    NUM_ELM_BLKS_GLOBAL,
    NUM_NODE_SETS_GLOBAL,
    NUM_SIDE_SETS_GLOBAL,
    NUM_INTERNAL_NODES,
    NUM_BORDER_NODES,
    NUM_EXTERNAL_NODES,
    NUM_INTERNAL_ELEMS,
    NUM_BORDER_ELEMS,
    NUM_NODE_COMM_MAPS,
    NUM_ELEM_COMM_MAPS,
    NUM_NBR_PROCS,
    NUM_MSIA};

  long long getMSI(MSIA ind){return msia[ind];}
  void setMSI(MSIA ind,long long the_int){msia[ind] = the_int;}

        long long * getMSP(MSPA ind)       {return mspa[ind];}
  const long long  * getMSP(MSPA ind) const {return mspa[ind];}

        std::string * getMSPSA(MSPSA ind)       {return mspsa[ind];}
  const std::string * getMSPSA(MSPSA ind) const {return mspsa[ind];}

        long long * const * getMSPP(MSPPA ind)       {return msppa[ind];}
  const long long  * const * getMSPP(MSPPA ind) const {return msppa[ind];}

        double * const * getMSPPD(MSPPDA ind)       {return msppda[ind];}
  const double * const * getMSPPD(MSPPDA ind) const {return msppda[ind];}

  static Mesh_Specification * static_storage;

  Mesh_Specification();
  Mesh_Specification( long long pid){
    Zero_Set();
    msia[PROC_ID] = pid;}
  virtual ~Mesh_Specification();
  
  std::string getErrorString()  {return error_stream.str();}
  std::string getWarningString(){return warning_stream.str();}
  

  // Access functions for local data

    const  std::string& Title()    const {return title;}

    // Nodal coordinates
    
    const double* Coord() const {return coord;}
          double* Coord()       {return coord;}
      // Stored by node, then by coordinate component.  Thus, Coord()[n] returns 
      // the X coordinate of the nth node, while Coord()[n+Number_Of_Nodes()]
      // returns the Y coordinate of the nth node if Dimensionality()>1.

	  virtual const std::string& File_Type() {return file_type;}
    
    // Information records

    typedef std::string QA_Record[4];
    const QA_Record *QA_Records() const {return qa_strings;}
          QA_Record *QA_Records()       {return qa_strings;}
      // The QA records give information on every code that has "touched"
      // the data in the mesh specification.  Each QA record consists
      // of four parts enumerated below:
    enum { 
           QA_CODE_NAME = 0,
           QA_CODE_DESCRIPTOR = 1,
           QA_ANALYSIS_DATE = 2,
           QA_ANALYSIS_TIME = 3
    };


    bool Are_Warnings_Suppressed()                  const;

    void Parallel_Data_Size(long long,long long,long long,long long,long long,long long,long long);
    void Allocate_Locational_Data();
    void Allocate_LoadBal_Data();
    void Allocate_Global_Data();
    void Global_Data_Size( long long, long long, long long, long long, long long, long long, long long);
    void Allocate_Parallel_Data();
    void Free_Parallel_Data();
    void Free_Locational_Data();
    void Free_Global_Data();

// Definition functions

    void Specify_Global_Information(const std::string &title,
                                    long long dimensionality,
                                    long long number_of_nodes,
                                    long long number_of_elements,
                                    long long number_of_element_blocks,
                                    long long number_of_node_sets,
                                    long long number_of_side_sets,
                                    long long number_of_qa_records,
                                    long long number_of_info_records);
      // Specifies the overall dimensions of the mesh, and allocates
      // storage for next level of information. 
    
    void Specify_Block_Information(long long index,
                                   long long block_id,
                                   long long number_of_block_elements,
                                   long long number_of_nodes_per_element,
                                   long long number_of_element_attributes,
                                   Element_Type block_element_type);
      // Specifies the dimensions of a particular element block, and
      // allocates storage for the next level of information for the block.
    

    
    void Specify_Node_Set_Information(long long index,
                                      long long node_set_id,
                                      long long number_of_nodes_in_node_set,
                                      long long number_of_df_in_node_set);
      // Specifies the dimensions of a particular node set, and allocates
      // storage for the next level of information for the node set.
    
    void Specify_Side_Set_Information(long long index,
                                      long long side_set_id,
                                      long long number_of_faces_in_side_set,
                                      long long number_of_nodes_in_side_set,
                                      long long number_of_df_in_side_set);
      // Specifies the dimensions of a particular side set, and allocates
      // storage for the next level of information for the side set.
    
    void Resize_Info_Store(long long number_of_info_records);
      // Increase the number of information records.

    void Suppress_Warnings(long long);


    void Free_NonTransient_Storage();
      // Free all storage except that which contains data required for
      // time step dumps.  After a call to Free_NonTransient_Storage(),
      // Block_ID() and Number_Of_Block_Elements() will return pointers
      // to meaningful data, but all other functions returning a pointer
      // will return a null pointer.



  protected:

    std::string title;
    std::stringstream error_stream;
    std::stringstream warning_stream;
 
    double*     coord;
    
    Element_Type* block_element_type;
    
    QA_Record*  qa_strings;

    bool         suppress_warnings;

    //nem data
    std::string file_type;

    //Arrays for storing ints pointers,
    // automatically sized by the enums
    long long msia[NUM_MSIA];
    long long * mspa[NUM_MSPA];
    long long * * msppa[NUM_MSPPA];
    double * * msppda[NUM_MSPPDA];
    std::string * mspsa[NUM_MSPSA];
  private:

    Mesh_Specification(const Mesh_Specification &);
    Mesh_Specification &operator=(const Mesh_Specification &);
    
    void Zero_Set();
    void Free();
};

}
#endif
