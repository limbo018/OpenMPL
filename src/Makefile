#==========================================================================
#                         Compilation Flags
# ==========================================================================

# default DBG is off
DBG = 0

UNAME_S = $(shell uname -s)
ifeq ($(UNAME_S), Linux)
	-include Unix.mk
endif 
ifeq ($(UNAME_S), Darwin)
	-include Darwin.mk
endif 

ifeq ($(DBG), 1)
	CXXFLAGS = $(CXXFLAGS_DEBUG)
else
	CXXFLAGS = $(CXXFLAGS_RELEASE)
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

MKDIR = if [ ! -d $(@D) ]; then mkdir -p $(@D); fi
EXEDIR = ../bin
OBJDIR = ../obj
LIBDIR = ../lib
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=$(OBJDIR)/%.o)
OBJS_NOMAIN = $(filter-out $(OBJDIR)/main.o, $(OBJS)) # filter out main.o 
DEPS = $(OBJS:%.o=%.d) 	# one dependency file for each source

LIB = $(LIBDIR)/SimpleMPL.a
EXE	= $(EXEDIR)/SimpleMPL

.SUFFIXES: .cc .cpp .o .x

all: $(EXE) $(LIB)

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
	$(CXX) $(CXXFLAGS) $< -MM -MT $(@:%.d=%.o) >$@

-include $(DEPS)

.PHONY: clean
clean: cleandep
	rm -f $(LIB) $(EXE) $(OBJS)
	
.PHONY: cleandep
cleandep:
	rm -f $(DEPS)
