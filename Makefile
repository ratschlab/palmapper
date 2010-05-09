SrcDir = ./src
OutDir = ./out
ObjDir = $(OutDir)/o
ExeDir = $(OutDir)/a

SVNVERSION = $(shell svnversion)

CC = g++
#CFLAGS = -Wall -g # debug
#CFLAGS = -Wno-unused-parameter -Wformat -Wformat-security -Wimplicit -Wparentheses -Wshadow -O9 -fexpensive-optimizations -frerun-cse-after-loop -fcse-follow-jumps -finline-functions -fschedule-insns2 -fthread-jumps -fforce-addr -fstrength-reduce -funroll-loops -march=native -mtune=native -pthread # linux amd64 optimized
CFLAGS = -O9 -Wall -g -Wno-unused-parameter -Wformat -Wformat-security -Wimplicit -Wparentheses -Wshadow # generic
GMFLAGS = -DGM
INCLUDE =  -Ishogun/ -Idyn_prog/ -Isrc
LDFLAGS = 

SHOGUN_OBJ = $(ObjDir)/pm/genomemapper/shogun/init.o \
	$(ObjDir)/pm/genomemapper/shogun/Mathematics.o \
	$(ObjDir)/pm/genomemapper/shogun/io.o \
	$(ObjDir)/pm/genomemapper/shogun/Parallel.o \
	$(ObjDir)/pm/genomemapper/shogun/Version.o \
	$(ObjDir)/pm/genomemapper/shogun/SGObject.o \
	$(ObjDir)/pm/genomemapper/shogun/ShogunException.o \
	$(ObjDir)/pm/genomemapper/shogun/Signal.o

DYNPROG_OBJ = $(ObjDir)/pm/genomemapper/dyn_prog/Mathmatics_dp.o \
	$(ObjDir)/pm/genomemapper/dyn_prog/io_dp.o \
	$(ObjDir)/pm/genomemapper/dyn_prog/qpalma_dp.o \
	$(ObjDir)/pm/genomemapper/dyn_prog/debug_tools.o \
	$(ObjDir)/pm/genomemapper/dyn_prog/penalty_info_dp.o \
	$(ObjDir)/pm/genomemapper/dyn_prog/result_align.o \
	$(ObjDir)/pm/genomemapper/dyn_prog/fill_matrix.o

GM_OBJ = $(ObjDir)/gm/genomemapper/GenomeMaps.o \
	$(ObjDir)/gm/genomemapper/QPalma.o \
	$(ObjDir)/gm/genomemapper/align.o \
	$(ObjDir)/gm/genomemapper/TopAlignments.o \
	$(ObjDir)/gm/genomemapper/IntervalQuery.o \
	$(ObjDir)/gm/genomemapper/genomemapper.o \
	$(ObjDir)/gm/genomemapper/init.o \
	$(ObjDir)/gm/genomemapper/print.o \
	$(ObjDir)/gm/genomemapper/Chromosome.o \
	$(ObjDir)/gm/genomemapper/Config.o \
	$(ObjDir)/gm/genomemapper/Genome.o \
	$(ObjDir)/gm/genomemapper/Hits.o \
	$(ObjDir)/gm/genomemapper/Read.o \
	$(ObjDir)/gm/genomemapper/Statistics.o \
	$(ObjDir)/gm/genomemapper/Util.o \
	$(SHOGUN_OBJ) $(DYNPROG_OBJ)

PM_OBJ = $(ObjDir)/pm/genomemapper/GenomeMaps.o \
	$(ObjDir)/pm/genomemapper/QPalma.o \
	$(ObjDir)/pm/genomemapper/align.o \
	$(ObjDir)/pm/genomemapper/TopAlignments.o \
	$(ObjDir)/pm/genomemapper/IntervalQuery.o \
	$(ObjDir)/pm/genomemapper/genomemapper.o \
	$(ObjDir)/pm/genomemapper/init.o \
	$(ObjDir)/pm/genomemapper/print.o \
	$(ObjDir)/pm/genomemapper/Chromosome.o \
	$(ObjDir)/pm/genomemapper/Config.o \
	$(ObjDir)/pm/genomemapper/Genome.o \
	$(ObjDir)/pm/genomemapper/Hits.o \
	$(ObjDir)/pm/genomemapper/Read.o \
	$(ObjDir)/pm/genomemapper/Statistics.o \
	$(ObjDir)/pm/genomemapper/Util.o \
	$(SHOGUN_OBJ) $(DYNPROG_OBJ)

PMIDX_OBJ = $(ObjDir)/pm/pmindex/init.o \
	$(ObjDir)/pm/pmindex/printindex.o \
	$(ObjDir)/pm/pmindex/usage.o \
	$(ObjDir)/pm/pmindex/write.o \
	$(ObjDir)/pm/pmindex/load.o \
	$(ObjDir)/pm/pmindex/index.o \
	$(ObjDir)/pm/pmindex/alloc.o \
	$(ObjDir)/pm/pmindex/pmindex.o

GMIDX_OBJ = $(ObjDir)/gm/pmindex/init.o \
	$(ObjDir)/gm/pmindex/printindex.o \
	$(ObjDir)/gm/pmindex/usage.o \
	$(ObjDir)/gm/pmindex/write.o \
	$(ObjDir)/gm/pmindex/load.o \
	$(ObjDir)/gm/pmindex/index.o \
	$(ObjDir)/gm/pmindex/alloc.o \
	$(ObjDir)/gm/pmindex/pmindex.o

CurrentDir := $(shell pwd)

all: palmapper pmindex genomemapper gmindex

palmapper: $(PM_OBJ) src/genomemapper/*.h 
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o palmapper $(PM_OBJ) -lpthread -lz -lm

pmindex:  $(PMIDX_OBJ) src/pmindex/*.h src/pmindex/pmindex_symbols.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o pmindex $(PMIDX_OBJ) 

genomemapper: $(GM_OBJ) src/genomemapper/*.h
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o genomemapper $(GM_OBJ) -lpthread -lz -lm

gmindex: $(GMIDX_OBJ) src/pmindex/*.h src/pmindex/pmindex_symbols.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o gmindex $(GMIDX_OBJ)

clean:
	rm -rf $(OutDir) palmapper pmindex genomemapper gmindex

test:
	(cd testcase; make test)

release_pm:
	make clean; mkdir -p ../release; cd ..; rsync -av $(CurrentDir) release; cd release; rm -rf */.settings */.cproject */.project */.svn */*/.svn */*/*/.svn; tar czvf ../release.$(SVNVERSION).tar.gz .; cd ..; rm -rf release

#todo also delete shogun...
release_gm:
	make clean; mkdir -p ../release; cd ..; rsync -av $(CurrentDir) release; cd release; rm -rf */.settings */.cproject */.project */.svn */*/.svn */*/*/.svn; tar czvf ../genomemapper.$(SVNVERSION).tar.gz .; cd ..; rm -rf release

# generic rule for compiling c++
$(ObjDir)/pm/%.o : $(SrcDir)/%.cpp
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

$(ObjDir)/pm/%.o : $(SrcDir)/%.cpp $(SrcDir)/%.h
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

$(ObjDir)/pm/%.o : $(SrcDir)/%.c
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

$(ObjDir)/gm/%.o : $(SrcDir)/%.cpp
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) $(GMFLAGS) $(INCLUDE) -o $@ -c $<

$(ObjDir)/gm/%.o : $(SrcDir)/%.cpp $(SrcDir)/%.h
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) $(GMFLAGS) $(INCLUDE) -o $@ -c $<

$(ObjDir)/gm/%.o : $(SrcDir)/%.c
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) $(GMFLAGS) $(INCLUDE) -o $@ -c $<

