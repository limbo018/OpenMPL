#==========================================================================
#                    Configuration under darwin environment
# ==========================================================================

# detect compiler 
ifneq ($(shell which clang++),)
	CXX = clang-omp++
	AR = ar
else
	CXX = g++
	AR = ar
endif

CXXFLAGS_BASIC = -ferror-limit=1 -W -Wall -Wextra -Wreturn-type -m64 -Wno-deprecated -Wno-unused-parameter -stdlib=libstdc++
CXXFLAGS_DEBUG = -g -DDEBUG $(CXXFLAGS_BASIC) 
CXXFLAGS_RELEASE = -O3 -fopenmp $(CXXFLAGS_BASIC) 

ARFLAGS = rvs

# dummy flags for linker
STATIC_LINK_FLAG = 
DYNAMIC_LINK_FLAG = 
