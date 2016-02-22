/*************************************************************************
    > File Name: Msg.cpp
    > Author: Yibo Lin
    > Mail: yibolin@utexas.edu
    > Created Time: Fri 31 Jul 2015 03:20:14 PM CDT
 ************************************************************************/

#include "Msg.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>

SIMPLEMPL_BEGIN_NAMESPACE

int mplPrint(MessageType m, const char* format, ...)
{
	va_list args;
	va_start(args, format);
	int ret = mplVPrintStream(m, stdout, format, args);
	va_end(args);

	return ret;
}

int mplPrintStream(MessageType m, FILE* stream, const char* format, ...)
{
	va_list args;
	va_start(args, format);
	int ret = mplVPrintStream(m, stream, format, args);
	va_end(args);

	return ret;
}

int mplVPrintStream(MessageType m, FILE* stream, const char* format, va_list args)
{
	// print prefix 
    char prefix[8];
    mplSPrintPrefix(m, prefix);
    // merge prefix and format 
    char formatBuf[256];
    sprintf(formatBuf, "%s%s", prefix, format);

	// print message 
    // only print once to ensure multi-thread safe 
    int ret = vfprintf(stream, formatBuf, args);
	
	return ret;
}

int mplSPrint(MessageType m, char* buf, const char* format, ...)
{
	va_list args;
	va_start(args, format);
	int ret = mplVSPrint(m, buf, format, args);
	va_end(args);

	return ret;
}

int mplVSPrint(MessageType m, char* buf, const char* format, va_list args)
{
	// print prefix 
    char prefix[8];
    mplSPrintPrefix(m, prefix);
	sprintf(buf, "%s", prefix);

	// print message 
	int ret = vsprintf(buf+strlen(prefix), format, args);
	
	return ret;
}

int mplSPrintPrefix(MessageType m, char* prefix)
{
	switch (m)
	{
		case kNONE:
            return sprintf(prefix, "%c", '\0');
		case kINFO:
			return sprintf(prefix, "(I) ");
		case kWARN:
            return sprintf(prefix, "(W) ");
		case kERROR:
            return sprintf(prefix, "(E) ");
		case kDEBUG:
            return sprintf(prefix, "(D) ");
        case kASSERT:
            return sprintf(prefix, "(A) ");
		default:
			mplAssertMsg(0, "unknown message type");
	}
    return 0;
}

void mplPrintAssertMsg(const char* expr, const char* fileName, unsigned lineNum, const char* funcName, const char* format, ...)
{
    // construct message 
    char buf[1024];
    va_list args;
	va_start(args, format);
    vsprintf(buf, format, args);
    va_end(args);

    // print message 
    mplPrintStream(kASSERT, stderr, "%s:%u: %s: Assertion `%s' failed: %s\n", fileName, lineNum, funcName, expr, buf);
}

void mplPrintAssertMsg(const char* expr, const char* fileName, unsigned lineNum, const char* funcName)
{
    // print message
    mplPrintStream(kASSERT, stderr, "%s:%u: %s: Assertion `%s' failed\n", fileName, lineNum, funcName, expr);
}

SIMPLEMPL_END_NAMESPACE
