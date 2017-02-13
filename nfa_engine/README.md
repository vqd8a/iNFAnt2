Note:
  + use negative numbers for states, which are also accepting states, in NFA transition table
  + multi-byte fectch input string into register before processing instead of reading bytes one by one
  + use 32 bits (unsigned int) to represent state id -> each thread reads only one transition instead of two
  + store NFA transition table in texture memory is better than in global memory