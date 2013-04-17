/*
 * tserror.h --
 *     Provide thread-safe strerror calls
 *
 * Copyright (c) 2012, SRI International.  All Rights Reserved.
 */

#ifndef tserror_h
#define tserror_h

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NO_TLS
char *ts_strerror(int errnum);
#else
#define ts_strerror strerror
#endif

void tserror_freeThread();

#ifdef __cplusplus
}
#endif

#endif /* tserror_h */

