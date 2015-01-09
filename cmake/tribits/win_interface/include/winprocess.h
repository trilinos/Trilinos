#ifdef _MSC_VER
# define NOMINMAX
# include <Winsock2.h>
# include <process.h>
# define getpid _getpid
inline void sleep(int sec)
{
  Sleep(sec * 1000);
}
#pragma comment(lib, "Ws2_32.lib") 
#endif
