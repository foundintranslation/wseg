/*
 * Counts.h --
 *	Utility functions for counts
 *
 * Copyright (c) 2006 SRI International.  All Rights Reserved.
 *
 * @(#)$Header: /home/srilm/devel/lm/src/RCS/Counts.h,v 1.4 2006/09/05 01:03:47 stolcke Exp $
 *
 */
#pragma once
#ifndef _Counts_h_
#define _Counts_h_

#ifdef PRE_ISO_CXX
# include <iostream.h>
#else
# include <iostream>
using namespace std;
#endif
#include <stdio.h>

#include "Boolean.h"
#include "XCount.h"
#include "File.h"

#ifdef USE_LONGLONG_COUNTS
typedef unsigned long long Count;	/* a count of something */
#else
typedef unsigned long Count;		/* a count of something */
#endif
typedef double FloatCount;		/* a fractional count */

/*
 * Type-dependent count <--> string conversions
 */
extern char ctsBuffer[100];

inline const char *
countToString(unsigned count)
{
    sprintf(ctsBuffer, "%u", count);
    return ctsBuffer;
}

inline const char *
countToString(int count)
{
    sprintf(ctsBuffer, "%d", count);
    return ctsBuffer;
}

inline const char *
countToString(long count)
{
    sprintf(ctsBuffer, "%ld", count);
    return ctsBuffer;
}

inline const char *
countToString(unsigned long count)
{
    sprintf(ctsBuffer, "%lu", count);
    return ctsBuffer;
}

inline const char *
countToString(unsigned long long count)
{
    sprintf(ctsBuffer, "%llu", count);
    return ctsBuffer;
}

inline const char *
countToString(XCount count)
{
    return countToString((XCountValue)count);
}

template <class CountT>
inline const char *
countToString(CountT count)
{
    sprintf(ctsBuffer, "%lg", (double)count);
    return ctsBuffer;
}

inline Boolean
stringToCount(const char *str, unsigned int &count)
{
    /*
     * scanf("%u") doesn't check for a positive sign, so we have to ourselves.
     */
    return (*str != '-' && sscanf(str, "%u", &count) == 1);
}

inline Boolean
stringToCount(const char *str, int &count)
{
    return (sscanf(str, "%d", &count) == 1);
}

inline Boolean
stringToCount(const char *str, unsigned short &count)
{
    /*
     * scanf("%u") doesn't check for a positive sign, so we have to ourselves.
     */
    return (*str != '-' && sscanf(str, "%hu", &count) == 1);
}

inline Boolean
stringToCount(const char *str, unsigned long &count)
{
    /*
     * scanf("%lu") doesn't check for a positive sign, so we have to ourselves.
     */
    return (*str != '-' && sscanf(str, "%lu", &count) == 1);
}

inline Boolean
stringToCount(const char *str, unsigned long long &count)
{
    /*
     * scanf("%lu") doesn't check for a positive sign, so we have to ourselves.
     */
    return (*str != '-' && sscanf(str, "%llu", &count) == 1);
}

inline Boolean
stringToCount(const char *str, long &count)
{
    return (sscanf(str, "%ld", &count) == 1);
}

inline Boolean
stringToCount(const char *str, XCount &count)
{
    XCountValue x;
    if (stringToCount(str, x)) {
    	count = x;
	return true;
    } else {
    	return false;
    }
}

template <class CountT>
static inline Boolean
stringToCount(const char *str, CountT &count)
{
    double x;
    if (sscanf(str, "%lf", &x) == 1) {
	count = x;
	return true;
    } else {
	return false;
    }
}

/*
 * Binary count I/O
 * 	Functions return 0 on failure,  number of bytes read/written otherwise
 */

unsigned writeBinaryCount(FILE *fp, unsigned long long count,
						    unsigned minBytes = 0);
unsigned writeBinaryCount(FILE *fp, float count);
unsigned writeBinaryCount(FILE *fp, double count);

inline unsigned
writeBinaryCount(FILE *fp, unsigned long count) {
    return writeBinaryCount(fp, (unsigned long long)count);
}

inline unsigned
writeBinaryCount(FILE *fp, unsigned count)
{
    return writeBinaryCount(fp, (unsigned long long)count);
}

inline unsigned
writeBinaryCount(FILE *fp, unsigned short count)
{
    return writeBinaryCount(fp, (unsigned long long)count);
}

inline unsigned
writeBinaryCount(FILE *fp, XCount count)
{
    return writeBinaryCount(fp, (unsigned long long)count);
}

//unsigned readBinaryCount(FILE *fp, unsigned long long &count);
//unsigned readBinaryCount(FILE *fp, float &count);
//unsigned readBinaryCount(FILE *fp, double &count);

//inline unsigned
//readBinaryCount(FILE *fp, unsigned long &count)
//{
//    unsigned long long lcount;
//    unsigned result = readBinaryCount(fp, lcount);
//    if (result > 0) {
//	count = lcount;
//    }
//    return result;
//}

//inline unsigned
//readBinaryCount(FILE *fp, unsigned &count)
//{
//    unsigned long long lcount;
//    unsigned result = readBinaryCount(fp, lcount);
//    if (result > 0) {
//	count = lcount;
//    }
//    return result;
//}

//inline unsigned
//readBinaryCount(FILE *fp, unsigned short &count)
//{
//    unsigned long long lcount;
//    unsigned result = readBinaryCount(fp, lcount);
//    if (result > 0) {
//	count = lcount;
//    }
//    return result;
//}

//inline unsigned
//readBinaryCount(FILE *fp, XCount &count)
//{
//    unsigned long long lcount;
//    unsigned result = readBinaryCount(fp, lcount);
//    if (result > 0) {
//	count = lcount;
//    }
//    return result;
//}

#endif /* _Counts_h_ */
