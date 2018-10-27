#==========================================================================
#                    Configuration under unix environment
# ==========================================================================

ifeq ("x","y")
# detect compiler 
ifneq ($(shell which g++48),)
	CXX = g++48
	AR = ar
else
ifneq ($(shell which g++47),)
	CXX = g++47
	AR = ar
else 
	CXX = g++
	AR = ar
endif
endif
endif

CXX = g++ -std=gnu++11
CXXFLAGS_BASIC =  -fmax-errors=1 -W -Wall -Wextra -Wreturn-type -ansi -m64 -Wno-deprecated -Wno-unused-local-typedefs
CXXFLAGS_DEBUG = -g -rdynamic -DQDEBUG $(CXXFLAGS_BASIC) 
CXXFLAGS_RELEASE = -O3 -fopenmp $(CXXFLAGS_BASIC) 
ARFLAGS = rvs

# gcc linker provides fine link control, while clang does not
STATIC_LINK_FLAG = -Wl,-Bstatic
DYNAMIC_LINK_FLAG = -Wl,-Bdynamic
