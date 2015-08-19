/*************************************************************************
    > File Name: msg.h
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Fri 31 Jul 2015 03:18:29 PM CDT
 ************************************************************************/

#ifndef SIMPLEMPL_MSG_H
#define SIMPLEMPL_MSG_H

#include <cstdarg>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include "namespace.h"

SIMPLEMPL_BEGIN_NAMESPACE

/// message type for print functions 
enum MessageType {
	kNONE = 0, 
	kINFO = 1, 
	kWARN = 2, 
	kERROR = 3, 
	kDEBUG = 4, 
    kASSERT = 5
};

/// print to screen (stdout)
int mplPrint(MessageType m, const char* format, ...);
/// print to stream 
int mplPrintStream(MessageType m, FILE* stream, const char* format, ...);
/// core function to print formatted data from variable argument list 
int mplVPrintStream(MessageType m, FILE* stream, const char* format, va_list args);
/// format to a buffer 
int mplSPrint(MessageType m, char* buf, const char* format, ...);
/// core function to format a buffer 
int mplVSPrint(MessageType m, char* buf, const char* format, va_list args);
/// format prefix 
int mplSPrintPrefix(MessageType m, char* buf);

/// assertion 
void mplPrintAssertMsg(const char* expr, const char* fileName, unsigned lineNum, const char* funcName, const char* format, ...);
void mplPrintAssertMsg(const char* expr, const char* fileName, unsigned lineNum, const char* funcName);

#define mplAssertMsg(condition, args...) do {\
    if (!(condition)) \
    {\
        ::SIMPLEMPL_NAMESPACE::mplPrintAssertMsg(#condition, __FILE__, __LINE__, __PRETTY_FUNCTION__, args); \
        abort(); \
    }\
} while (false)
#define mplAssert(condition) do {\
    if (!(condition)) \
    {\
        ::SIMPLEMPL_NAMESPACE::mplPrintAssertMsg(#condition, __FILE__, __LINE__, __PRETTY_FUNCTION__); \
        abort(); \
    }\
} while (false)

/// static assertion 
template <bool>
struct mplStaticAssert;
template <>
struct mplStaticAssert<true> 
{
    mplStaticAssert(const char* = NULL) {}
};


SIMPLEMPL_END_NAMESPACE

#endif
