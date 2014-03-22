#ifdef DEBUG
  #define PRINTF(...) printf(__VA_ARGS__)
#else
  #define PRINTF(...)
#endif

#ifndef DEBUG
  #define DEBUG 0
#endif
