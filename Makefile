## Basic options + ROOT
CXXFLAGS = -Wall -g -O2 $(shell root-config --cflags)
#CXXFLAGS = -Wall -fPIC $(shell root-config --cflags)
LDFLAGS  = $(shell root-config --libs) -lMinuit2
#LDFLAGS  = $(shell root-config --libs) -lRooFit -lMinuit2

## LHAPDF
#CXXFLAGS += $(shell /data2/analysis/kenichi/local/bin/lhapdf-config --cflags)
#LDFLAGS  += $(shell /data2/analysis/kenichi/local/bin/lhapdf-config --libs) -Wl,-rpath,$(shell /data2/analysis/kenichi/local/bin/lhapdf-config --libdir)

## BOOST
#CXXFLAGS += -I$(BOOST)

## MySQL
# CXXFLAGS += $(shell mysql_config --include)
# LDFLAGS  += $(shell mysql_config --libs)

## Fortran
# LDFLAGS += -lgfortran
# F77=gfortran
# #F77=f77
# #F77_OPTS= -fno-automatic -fno-second-underscore -falign-commons -fd-lines-as-comments
# F77_OPTS= -fno-automatic -fno-second-underscore -fd-lines-as-comments

ROOT_DIC = obj/root_dict

## Files to be compiled
##    A *.cc file is compiled into an executable 
##    if the corresponding header file (*.h) doesn't exist.
##    Otherwise a *.cc file is compiled and linked to executables.
HDRS = $(filter-out LinkDef.h, $(wildcard *.h))
SRCS = $(wildcard *.cc)

OBJS_ALL = $(patsubst %.cc, obj/%.o, $(SRCS))
OBJS_SUB = $(patsubst  %.h, obj/%.o, $(HDRS))
OBJS_EXE = $(filter-out $(OBJS_SUB), $(OBJS_ALL))
EXES = $(patsubst obj/%.o, %, $(OBJS_EXE))
DEPS = $(patsubst %.cc, obj/%.d, $(SRCS))

# SRCS_F = $(wildcard *.f)
# OBJS_F = $(patsubst %.f, obj/%.o, $(SRCS_F))

DEPFILE = depend.mk

.PHONY : all
all : $(EXES)
$(EXES) : % : obj/%.o $(OBJS_SUB) $(ROOT_DIC).o
#$(EXES) : % : obj/%.o $(OBJS_SUB) $(OBJS_F)
	@echo "Create  $@"
	@$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

# obj/%.o : %.f 
# 	@echo "Compile $<"
# 	@$(F77) $(F77_OPTS) -c -o $@ $< 

obj/%.o : %.cc
	@echo "Compile $<"
	@$(CXX) -c $(CXXFLAGS) -MMD -MP -o $@ $<

$(ROOT_DIC).o : Event.h LinkDef.h
	@echo "Rootcint"
	@rootcint -f $(ROOT_DIC).cc -c $(CXXFLAGS) $^
	@$(CXX) -c $(CXXFLAGS) -I. -o $@ $(ROOT_DIC).cc

.PHONY : clean clean-obj clean-exe
clean : clean-obj clean-exe
clean-obj : 
	rm -f $(OBJS_ALL) $(DEPS) $(OBJS_F) $(ROOT_DIC).*
clean-exe :
	rm -f $(EXES)

-include $(DEPS)
