/*
 * io_timer.h
 *
 *  Created on: Apr 20, 2009
 *      Author: thkorde
 */

#ifndef IO_TIMER_H_
#define IO_TIMER_H_

#if defined(USE_TIMERS)
#define Start_Timer(timer) { timer = MPI_Wtime(); }
#define Stop_Timer(name, timer)  { timer = MPI_Wtime() - timer; log_debug(LOG_ALL, "%s Time = %10.8f", name, timer); }
#else
#define Start_Timer(timer)  {}
#define Stop_Timer(name, timer)   {}
#endif

#if defined(USE_TIMERS) && defined(USE_TIMER_BARRIERS)
#define Timer_Barrier(comm) { MPI_Barrier(comm); }
#else
#define Timer_Barrier(comm) {}
#endif

#endif /* IO_TIMER_H_ */
