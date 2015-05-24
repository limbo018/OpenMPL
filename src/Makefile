#==========================================================================
#                         C++ Compiler
# ==========================================================================
# search for suitable g++
# c++11 is required
CXX49 = $(shell which g++-4.9)
CXX48 = $(shell which g++48)
CXX47 = $(shell which g++47)

# search for g++-4.9, g++48, then g++47, then g++
ifneq ($(CXX49),)
	CXX = g++-4.9
	# g++-4.9 support colored error information
	CXXFLAGS += -fdiagnostics-color=always #-static # no static option under mac 
else 
	ifneq ($(CXX48),)
		CXX = g++48
		CXXFLAGS += -fmax-errors=1 #-static
	else
		ifneq ($(CXX47),)
			CXX = g++47
			CXXFLAGS += -fmax-errors=1 #-static
		else 
			CXX = g++
			# no fmax-errors for gcc-4.4 
			# the g++ in Oracle is not built staticly 
		endif
	endif
endif

#==========================================================================
#                         Special Library
# ==========================================================================

INCLUDE = \
		  -I $(BOOST_DIR)/include \
		  -I $(LIMBO_DIR) \
		  -I $(GUROBI_HOME)/include

# only gurobi60 is dynamically linked 
# it may have problem in linking boost library (.a or .so) without -rpath
LIBS = \
	   -L $(BOOST_DIR)/lib -lboost_system -lboost_timer -Wl,-rpath=$(BOOST_DIR)/lib \
	   -L $(LIMBO_DIR)/lib -lgdsparser \
	   -L $(GUROBI_HOME)/lib -lgurobi_c++ \
	   -L $(GUROBI_HOME)/lib -lgurobi60


# ==========================================================================
#                         Standard Setting
# ==========================================================================
# compilation flags
# STATWINDOW_POLYGON: if defined, store pointer of polygons in each window; 
#                     otherwise store pointer of decomposed rectangles from these polygons.
#                     The latter can provide around 1 second faster in window route area calculation than the former in small benchmark.
# ENABLE_ILP

CXXFLAGS += -W -Wall -Wextra -Wreturn-type -ansi -g #-fopenmp # -funroll-all-loops -m64 -fPIC -Wno-deprecated
CXXFLAGS += -fmax-errors=1 -m64 -Wno-deprecated -DDEBUG

MKDIR = if [ ! -d $(@D) ]; then mkdir -p $(@D); fi
BPATH = ../bin/
OPATH = ../obj/
SRCS	= $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=$(OPATH)%.o)

EXE	 = $(BPATH)SimpleMPL

.SUFFIXES: .cc .cpp .o .x

all: $(EXE) tags

$(EXE): $(OBJS)
	@$(MKDIR)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS) 

$(OPATH)%.o:%.cpp
	@$(MKDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

tags: $(SRCS)
	ctags -R *

clean:
	rm -f *.o *~ $(EXE) $(OPATH)*
	
