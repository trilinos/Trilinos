#ifndef STK_DEBUGGING_MACROS_HPP
#define STK_DEBUGGING_MACROS_HPP

//----------------------------------------------------------------------

// Use macro below to activate counters that track calls to mesh-modification routines
//#define STK_MESH_MODIFICATION_COUNTERS

#ifdef STK_MESH_MODIFICATION_COUNTERS
#define INCREMENT_MODIFICATION_COUNTER(METHOD_TYPE, MOD_TYPE) {++m_modification_counters[METHOD_TYPE][MOD_TYPE];}
#define INCREMENT_ENTITY_MODIFICATION_COUNTER(METHOD_TYPE, RANK,MOD_TYPE) {++m_entity_modification_counters[METHOD_TYPE][RANK][MOD_TYPE];}
#else
#define INCREMENT_MODIFICATION_COUNTER(METHOD_TYPE, MOD_TYPE) {}
#define INCREMENT_ENTITY_MODIFICATION_COUNTER(METHOD_TYPE, RANK,MOD_TYPE) {}
#endif

#endif
