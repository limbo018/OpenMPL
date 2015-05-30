#==========================================================================
#                    Configuration under darwin environment
# ==========================================================================

# detect compiler 
ifneq ($(shell which clang++),)
	CXX = clang++
	AR = ar
else
	CXX = g++
	AR = ar
endif
endif

CXXFLAGS_BASIC = -ferror-limit=1 -W -Wall -Wextra -Wreturn-type -ansi -m64 -Wno-deprecated
CXXFLAGS_DEBUG = -g -DDEBUG $(CXXFLAGS_BASIC) 
CXXFLAGS_RELEASE = -O2 -fopenmp $(CXXFLAGS_BASIC) 

ARFLAGS = rvs

