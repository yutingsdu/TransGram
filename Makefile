#-- for now these MUST point to the included "samtools-0.x.x" and "gclib" sub-directories
BAM  := ./samtools-0.1.18
GDIR := ./gclib
#--

INCDIRS := -I. -I${GDIR} -I${BAM}

CXX   := $(if $(CXX),$(CXX),g++)

BASEFLAGS := -Wall -Wextra ${INCDIRS} -fsigned-char -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -std=c++11 -fno-strict-aliasing -fno-exceptions -fno-rtti
#for gcc 8+ add: -Wno-class-memaccess
GCCVER5 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 5)
ifeq "$(GCCVER5)" "1"
 BASEFLAGS += -Wno-implicit-fallthrough
endif

GCCVER8 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 8)
ifeq "$(GCCVER8)" "1"
  BASEFLAGS += -Wno-class-memaccess
endif

LINKER  := $(if $(LINKER),$(LINKER),g++)

LDFLAGS := $(if $(LDFLAGS),$(LDFLAGS),-g)

LDFLAGS += -L${BAM}

LIBS    := -lbam -lz

ifneq (,$(filter %nothreads %prof %profile, $(MAKECMDGOALS)))
 NOTHREADS=1 #yuting New LR
endif

#detect MinGW (Windows environment)
ifneq (,$(findstring mingw,$(shell ${CXX} -dumpmachine)))
 WINDOWS=1
endif

# Misc. system commands
#ifdef WINDOWS ##<-- use MSYS
# RM = del /Q
#else
RM = rm -f
#endif

# File endings
ifdef WINDOWS
 EXE = .exe
else
 EXE =
endif

# Non-windows systems need pthread
ifndef WINDOWS
 ifndef NOTHREADS
   LIBS := -pthread ${LIBS}
   BASEFLAGS += -pthread
 endif
endif

ifdef NOTHREADS
  BASEFLAGS += -DNOTHREADS
endif

DMACH := $(shell ${CXX} -dumpmachine)

ifneq (,$(filter %release %static, $(MAKECMDGOALS)))
  # -- release build
  RELEASE_BUILD=1
  CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O3)
  CXXFLAGS += -DNDEBUG $(BASEFLAGS)
else
  ifneq (,$(filter %memcheck %memdebug %tsan %tcheck %thrcheck, $(MAKECMDGOALS)))
     #use sanitizer in gcc 4.9+
     GCCVER49 := $(shell expr `${CXX} -dumpversion | cut -f1,2 -d.` \>= 4.9)
     ifeq "$(GCCVER49)" "0"
       $(error gcc version 4.9 or greater is required for this build target)
     endif
     CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O0)
     SANLIBS :=
     ifneq (,$(filter %tsan %tcheck %thrcheck, $(MAKECMDGOALS)))
        # thread sanitizer only (incompatible with address sanitizer)
        CXXFLAGS += -fno-omit-frame-pointer -fsanitize=thread -fsanitize=undefined $(BASEFLAGS)
        SANLIBS := -ltsan
     else
        # address sanitizer
        CXXFLAGS += -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address $(BASEFLAGS)
        SANLIBS := -lasan
     endif
     ifeq "$(GCCVER5)" "1"
       CXXFLAGS += -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=vptr
       CXXFLAGS += -fsanitize=float-cast-overflow -fsanitize=object-size
       #CXXFLAGS += -fcheck-pointer-bounds -mmpx
     endif
     CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG -fno-common -fstack-protector
     LIBS := ${SANLIBS} -lubsan -ldl ${LIBS}
  else
     ifneq (,$(filter %prof %profile, $(MAKECMDGOALS)))
     ## profiling build
       CXXFLAGS := -DNDEBUG $(BASEFLAGS) -g -pg
       LDFLAGS += -g -pg
     else
        #just plain debug build
        DEBUG_BUILD=1
        CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O0)
        ifneq (, $(findstring darwin, $(DMACH)))
           CXXFLAGS += -gdwarf-3
        endif
        CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG $(BASEFLAGS)
     endif
  endif
endif

ifdef RELEASE_BUILD
 ifneq (,$(findstring static, $(MAKECMDGOALS)))
    STATIC_CLIB=1
 endif
endif

ifdef STATIC_CLIB
 LDFLAGS += -static-libgcc -static-libstdc++
endif

ifdef DEBUG_BUILD
  #$(warning Building DEBUG version of transgram-graph.. )
  DBG_WARN=@echo
  DBG_WARN+='WARNING: built DEBUG version [much slower], use "make clean release" for a faster, optimized version of the program.'
endif

OBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o ${GDIR}/GBam.o \
 ${GDIR}/gdna.o ${GDIR}/codons.o ${GDIR}/GFastaIndex.o ${GDIR}/GFaSeqGet.o ${GDIR}/gff.o 

ifneq (,$(filter %memtrace %memusage %memuse, $(MAKECMDGOALS)))
    CXXFLAGS += -DGMEMTRACE
    OBJS += ${GDIR}/proc_mem.o
endif

ifndef NOTHREADS
 OBJS += ${GDIR}/GThreads.o 
endif

%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

OBJS3 := $(OBJS)
OBJS += build_graph.o tablemaker.o tmerge.o common.o

OBJS1 = simplify-graph.o
OBJS2 = trans2path.o

OBJS3 += rlink.o tablemakerE.o tmergeE.o

all release static debug: transgram-graph${EXE} transgram-path-search${EXE} transgram-trans2pathsinfo${EXE} transgram-filter${EXE} transgram-expression${EXE} 
memcheck memdebug tsan tcheck thrcheck: transgram-graph${EXE} transgram-path-search${EXE} transgram-trans2pathsinfo${EXE} transgram-filter${EXE} transgram-expression${EXE}
memuse memusage memtrace: transgram-graph${EXE} transgram-path-search${EXE} transgram-trans2pathsinfo${EXE} transgram-filter${EXE} transgram-expression${EXE}
prof profile: transgram-graph${EXE} transgram-path-search${EXE} transgram-trans2pathsinfo${EXE} transgram-filter${EXE} transgram-expression${EXE}
nothreads: transgram-graph${EXE} transgram-path-search${EXE} transgram-trans2pathsinfo${EXE} transgram-filter${EXE} transgram-expression${EXE}

transgram-path-search.o : simplify-graph.h
simplify-graph.o : simplify-graph.h

transgram-graph.o : $(GDIR)/GBitVec.h $(GDIR)/GHashMap.hh $(GDIR)/GBam.h
build_graph.o : build_graph.h tablemaker.h $(GDIR)/GBam.h $(GDIR)/GBitVec.h
tmerge.o : build_graph.h tmerge.h 
tablemaker.o : tablemaker.h build_graph.h
common.o : common.h

transgram-expression.o : $(GDIR)/GBitVec.h $(GDIR)/GHashMap.hh $(GDIR)/GBam.h
rlink.o : rlink.h tablemakerE.h $(GDIR)/GBam.h $(GDIR)/GBitVec.h
tmergeE.o : rlink.h tmergeE.h
tablemakerE.o : tablemakerE.h rlink.h

${BAM}/libbam.a: 
	cd ${BAM} && make lib
transgram-graph${EXE}: ${BAM}/libbam.a $(OBJS) transgram-graph.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}

transgram-expression${EXE}: ${BAM}/libbam.a $(OBJS3) transgram-expression.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}

transgram-path-search${EXE}: $(OBJS1) transgram-path-search.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}

transgram-trans2pathsinfo${EXE}: $(OBJS2) transgram-trans2pathsinfo.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}

transgram-filter${EXE}: transgram-filter.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}
#test demo tests: transgram-graph${EXE}
#	@./run_tests.sh

.PHONY : clean cleanall cleanAll allclean

COPTS  = -ggdb -O2 -fopenmp -std=gnu++0x -fstack-protector-all
CFLAGS =
COMPILE = $(CXX) $(INCS) $(CFLAGS) $(COPTS) ${LDFLAGS}
VPATH = .

#simplify-graph.o : simplify-graph.h
#transgram-path-search : transgram-path-search.o simplify-graph.o
#	gcc transgram-path-search.cpp -c -o simplify-graph.o

#%.o: %.cpp %.h
#	$(COMPILE) $< -c -o $@

# list executable file names
#
#EXECS1 = transgram-path-search
#
# compile and link
#
default:
	@echo
	@echo " to build:"
	@echo "    make all"
	@echo
	@echo " to clean:"
	@echo "    make clean"
	@echo "    make realclean"
	@echo

#all: directories $(EXECS1)

#$(EXECS1): $(OBJS1)
#	$(foreach EX, $(EXECS1), $(COMPILE) $(EX).cpp  -c -o $(DIR_OBJ)/$(EX).o;)
#	$(foreach EX, $(EXECS1), $(COMPILE) $(OBJS1) $(DIR_OBJ)/$(EX).o -o $(DIR_BIN)/$(EX) $(LIBST);)


#$(DIR_OBJ)/%.o: %.cpp %.h
#	$(COMPILE) $< -c -o $@


# target for removing all object files

#	echo $(PATH)
clean:
	${RM} transgram-graph${EXE} transgram-path-search${EXE} transgram-graph.o*  transgram-path-search.o $(OBJS) $(OBJS1) $(OBJS2)
	${RM} transgram-trans2pathsinfo${EXE} transgram-trans2pathsinfo.o 
	${RM} transgram-filter${EXE} transgram-filter.o
	${RM} transgram-expression${EXE} transgram-expression.o
	${RM} core.* rlink.o tablemakerE.o tmergeE.o trans2path.o
allclean cleanAll cleanall:
	cd ${BAM} && make clean
	${RM} transgram-graph${EXE} transgram-path-search${EXE} transgram-graph.o* transgram-path-search.o $(OBJS) $(OBJS1)
	${RM} transgram-trans2pathsinfo${EXE} transgram-trans2pathsinfo.o
	${RM} transgram-filter${EXE} transgram-filter.o
	${RM} transgram-expression${EXE} transgram-expression.o
	${RM} core.* rlink.o tablemakerE.o tmergeE.o trans2path.o
