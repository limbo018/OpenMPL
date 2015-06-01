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
# static libraries have .a as suffix in both unix and darwin 
# but dynamic libraries have different suffixes 
LIBS = \
	   $(BOOST_DIR)/lib/libboost_system.a \
	   $(BOOST_DIR)/lib/libboost_timer.a \
	   $(BOOST_DIR)/lib/libboost_chrono.a \
	   $(LIMBO_DIR)/lib/libgdsparser.a \
	   $(GUROBI_HOME)/lib/libgurobi_c++.a \
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
