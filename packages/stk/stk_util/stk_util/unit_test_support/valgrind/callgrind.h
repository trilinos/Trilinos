#ifndef stk_perf_callgrind_h
#define stk_perf_callgrind_h

// A dummy file used so that performance unit tests will compile even if the
// current system does not have valgrind

#define CALLGRIND_START_INSTRUMENTATION ThrowRequireMsg(false, "Not able to include correct callgrind headers")
#define CALLGRIND_STOP_INSTRUMENTATION ThrowRequireMsg(false, "Not able to include correct callgrind headers")
#define CALLGRIND_TOGGLE_COLLECT ThrowRequireMsg(false, "Not able to include correct callgrind headers")

#define __VALGRIND_MAJOR__    0
#define __VALGRIND_MINOR__    0

#endif
