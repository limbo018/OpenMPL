#==========================================================================
#                         Compilation Flags
# ==========================================================================

# recommend make command: make -j CSDP=1 
# be careful that no additional space is allowed
# default DBG is off
DBG = 1
# default GPROF is off, used to enable runtime profiling
GPROF = 0
# default GUROBI is off 
GUROBI = 1
# default LEMONCBC is off 
LEMONCBC = 0
# default CSDP is off 
CSDP = 0

UNAME_S = $(shell uname -s)
ifeq ($(UNAME_S), Linux)
	-include ../Unix.mk
endif 
ifeq ($(UNAME_S), Darwin)
	-include ../Darwin.mk
endif 

ifeq ($(DBG), 1)
	CXXFLAGS = $(CXXFLAGS_DEBUG) #-D_DGOUT -DDEBUG_SDPCOLORING -DDEBUG_ILPCOLORING  -DDEBUG_LPCOLORING 
else
	CXXFLAGS = $(CXXFLAGS_RELEASE) #-DCOMPONENTS -DDEBUG_NONINTEGERS -DDEBUG_NOANCHOR
endif

ifeq ($(GPROF), 1)
	CXXFLAGS += -pg
endif

CXXFLAGS += -g -DGUROBI=$(GUROBI) -DLEMONCBC=$(LEMONCBC) -DCSDP=$(CSDP)

ifeq ($(UNAME_S), Linux)
ifeq ($(GUROBI), 0)
ifeq ($(LEMONCBC), 0)
	CXXFLAGS += -static
endif
endif
endif

#==========================================================================
#                         Special Library
# ==========================================================================

INCLUDE = \
		  -I $(BOOST_DIR)/include \
		  -I $(LIMBO_DIR)/include

# only gurobi60 is dynamically linked 
# static libraries have .a as suffix in both unix and darwin 
# but dynamic libraries have different suffixes 
LIBS = $(STATIC_LINK_FLAG) \
	   -L $(LIMBO_DIR)/lib -lgdsparser -lprogramoptions \
	   -L $(BOOST_DIR)/lib -lboost_system -lboost_timer -lboost_chrono -lboost_iostreams \
	   -L $(ZLIB_DIR) -lz

# boost chrono links to to librt.a under Linux
ifeq ($(UNAME_S), Linux)
#LIBS += $(STATIC_LINK_FLAG) -lrt 
LIBS += -L $(BOOST_DIR)/lib -lboost_system
endif 

# GUROBI 
ifeq ($(GUROBI), 1)
INCLUDE += -I $(GUROBI_HOME)/include
LIBS += -L $(GUROBI_HOME)/lib \
		$(STATIC_LINK_FLAG) -lgurobi_c++ \
		$(DYNAMIC_LINK_FLAG) -lgurobi81
endif 
# Cbc lemon api (Cbc must be installed with static flag and lemon must be installed with Cbc as MIP solver)
ifeq ($(LEMONCBC), 1)
INCLUDE += -I $(CBC_DIR)/include \
		   -I $(LEMON_DIR)/include
LIBS += -L $(CBC_DIR)/lib $(DYNAMIC_LINK_FLAG) -lCbc -lCbcSolver -lCgl -lClp -lCoinUtils -lOsiCbc -lOsiClp -lOsi \
		-L $(LEMON_DIR)/lib $(STATIC_LINK_FLAG) -lemon \
		$(DYNAMIC_LINK_FLAG) -lz -lbz2 \
		-Wl,-rpath,$(CBC_DIR)/lib:$(LEMON_DIR)/lib
endif
# Csdp api (OpenBLAS and Fortran library are necessary)
# it looks like pthead must be linked dynamically when GUROBI is also on 
ifeq ($(CSDP), 1)
LIBS += -L $(LIMBO_DIR)/lib $(STATIC_LINK_FLAG) -lsdp \
	   -L $(LIMBO_DIR)/lib $(STATIC_LINK_FLAG) -lopenblas-st -lm 
# CSDP links to Fortran under Linux 
ifeq ($(UNAME_S), Linux)
LIBS += $(DYNAMIC_LINK_FLAG) -lrt -lgfortran 
endif
ifeq ($(GUROBI), 1) 
LIBS += $(DYNAMIC_LINK_FLAG) -lpthread 
else
ifeq ($(LEMONCBC), 1)
LIBS += $(DYNAMIC_LINK_FLAG) -lpthread 
else
LIBS += $(STATIC_LINK_FLAG) -lpthread
endif
endif
endif

# BOOST_REG_INTERSECTION
ifeq ($(BOOST_REG_INTERSECTION), 1)
	CXXFLAGS += -DBOOST_REG_INTERSECTION
endif

# ==========================================================================
#                         Standard Setting
# ==========================================================================

PROGRAM_ROOT_DIR = ../..
MKDIR = if [ ! -d $(@D) ]; then mkdir -p $(@D); fi
EXEDIR = $(PROGRAM_ROOT_DIR)/bin
OBJDIR = $(PROGRAM_ROOT_DIR)/obj
LIBDIR = $(PROGRAM_ROOT_DIR)/lib
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=$(OBJDIR)/%.o)
OBJS_NOMAIN = $(filter-out $(OBJDIR)/main.o, $(OBJS)) # filter out main.o 
DEPS = $(OBJS:%.o=%.d) 	# one dependency file for each source

LIB = $(LIBDIR)/OpenMPL.a
EXE	= $(EXEDIR)/OpenMPL

.SUFFIXES: .cc .cpp .o .x

all: $(EXE) $(LIB) tags

$(EXE): $(OBJDIR)/main.o $(LIB)
	@$(MKDIR)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIB) $(LIBS) 

$(LIB): $(OBJS)
	@$(MKDIR)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS_NOMAIN)

$(OBJDIR)/%.o:%.cpp
	@$(MKDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

# rule to generate a dep file by using the C preprocessor
# (see man cpp for details on the -MM and -MT options)
$(OBJDIR)/%.d:%.cpp
	$(CXX) $(CXXFLAGS) $< -MM -MT $(@:%.d=%.o) >$@ $(INCLUDE)

tags: $(SRCS)
	ctags -R *


-include $(DEPS)

.PHONY: clean
clean: cleandep
	rm -f $(LIB) $(EXE) $(OBJS)
	
.PHONY: cleandep
cleandep:
	rm -f $(DEPS)