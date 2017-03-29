// Dear emacs, this is -*- c++ -*-                                                             
#ifndef GLUONBB_ERRORCHECK_H
#define GLUONBB_ERRORCHECK_H

// Error checking macro                                                                        
#define CHECK( NAME, ARG )                           \
  do {                                               \
  const bool result = ARG;                         \
  const char* name = NAME;                         \
  if(!result) {                                    \
  ::Error(name, "Failed to execute: \"%s\"",     \
	  #ARG );                                \
  return 1;                                      \
  }                                                \
  } while( false )
#endif 

