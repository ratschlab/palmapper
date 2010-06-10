SrcDir = ./src
OutDir = ./out
ObjDir = $(OutDir)/o
ExeDir = $(OutDir)/a

SVNVERSION = $(shell svnversion)

CC = g++
#CFLAGS = -Wall -g # debug
#CFLAGS = -Wno-unused-parameter -Wformat -Wformat-security -Wimplicit -Wparentheses -Wshadow -O9 -fexpensive-optimizations -frerun-cse-after-loop -fcse-follow-jumps -finline-functions -fschedule-insns2 -fthread-jumps -fforce-addr -fstrength-reduce -funroll-loops -march=native -mtune=native -pthread # linux amd64 optimized
CFLAGS = -O9 -Wall -g -Wno-unused-parameter -Wformat -Wformat-security -Wimplicit -Wparentheses -Wshadow # generic
#CFLAGS = -O9 -Wall -g -pg -Wno-unused-parameter -Wformat -Wformat-security -Wimplicit -Wparentheses -Wshadow # generic
GMFLAGS = -DGM
INCLUDE =  -Ishogun/ -Idyn_prog/ -Isrc
LDFLAGS = 

SHOGUN_OBJ = $(ObjDir)/pm/palmapper/shogun/init.o \
	$(ObjDir)/pm/palmapper/shogun/Mathematics.o \
	$(ObjDir)/pm/palmapper/shogun/io.o \
	$(ObjDir)/pm/palmapper/shogun/Parallel.o \
	$(ObjDir)/pm/palmapper/shogun/Version.o \
	$(ObjDir)/pm/palmapper/shogun/SGObject.o \
	$(ObjDir)/pm/palmapper/shogun/ShogunException.o \
	$(ObjDir)/pm/palmapper/shogun/Signal.o

DYNPROG_OBJ = $(ObjDir)/pm/palmapper/dyn_prog/Mathmatics_dp.o \
	$(ObjDir)/pm/palmapper/dyn_prog/io_dp.o \
	$(ObjDir)/pm/palmapper/dyn_prog/qpalma_dp.o \
	$(ObjDir)/pm/palmapper/dyn_prog/debug_tools.o \
	$(ObjDir)/pm/palmapper/dyn_prog/penalty_info_dp.o \
	$(ObjDir)/pm/palmapper/dyn_prog/result_align.o \
	$(ObjDir)/pm/palmapper/dyn_prog/fill_matrix.o
	
LANG_OBJ = $(ObjDir)/pm/lang/Thread.o

GM_OBJ = $(ObjDir)/gm/palmapper/GenomeMaps.o \
	$(ObjDir)/gm/palmapper/QPalma.o \
	$(ObjDir)/gm/palmapper/align.o \
	$(ObjDir)/gm/palmapper/TopAlignments.o \
	$(ObjDir)/gm/palmapper/IntervalQuery.o \
	$(ObjDir)/gm/palmapper/palmapper.o \
	$(ObjDir)/gm/palmapper/init.o \
	$(ObjDir)/gm/palmapper/print.o \
	$(ObjDir)/gm/palmapper/Chromosome.o \
	$(ObjDir)/gm/palmapper/Config.o \
	$(ObjDir)/gm/palmapper/Genome.o \
	$(ObjDir)/gm/palmapper/Hits.o \
	$(ObjDir)/gm/palmapper/Mapper.o \
	$(ObjDir)/pm/palmapper/QueryFile.o \
	$(ObjDir)/gm/palmapper/Read.o \
	$(ObjDir)/gm/palmapper/Statistics.o \
	$(ObjDir)/gm/palmapper/Util.o \
	$(SHOGUN_OBJ) $(DYNPROG_OBJ) $(LANG_OBJ)

PM_OBJ = $(ObjDir)/pm/palmapper/GenomeMaps.o \
	$(ObjDir)/pm/palmapper/QPalma.o \
	$(ObjDir)/pm/palmapper/align.o \
	$(ObjDir)/pm/palmapper/TopAlignments.o \
	$(ObjDir)/pm/palmapper/IntervalQuery.o \
	$(ObjDir)/pm/palmapper/palmapper.o \
	$(ObjDir)/pm/palmapper/init.o \
	$(ObjDir)/pm/palmapper/print.o \
	$(ObjDir)/pm/palmapper/Chromosome.o \
	$(ObjDir)/pm/palmapper/Config.o \
	$(ObjDir)/pm/palmapper/Genome.o \
	$(ObjDir)/pm/palmapper/Hits.o \
	$(ObjDir)/pm/palmapper/Mapper.o \
	$(ObjDir)/pm/palmapper/QueryFile.o \
	$(ObjDir)/pm/palmapper/Read.o \
	$(ObjDir)/pm/palmapper/Statistics.o \
	$(ObjDir)/pm/palmapper/Util.o \
	$(SHOGUN_OBJ) $(DYNPROG_OBJ) $(LANG_OBJ)

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

palmapper: $(PM_OBJ) src/palmapper/*.h 
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o palmapper $(PM_OBJ) -lpthread -lz -lm

pmindex:  $(PMIDX_OBJ) src/pmindex/*.h src/pmindex/pmindex_symbols.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o pmindex $(PMIDX_OBJ) 

genomemapper: $(GM_OBJ) src/palmapper/*.h
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o genomemapper $(GM_OBJ) -lpthread -lz -lm

gmindex: $(GMIDX_OBJ) src/pmindex/*.h src/pmindex/pmindex_symbols.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o gmindex $(GMIDX_OBJ)

clean:
	rm -rf $(OutDir) palmapper pmindex genomemapper gmindex

test:
	(cd testcase; make test)

release_pm:
	make clean; mkdir -p ../release; cd ..; rsync -av $(CurrentDir) release; cd release; rm -rf */.settings */.cproject */.project */.svn */*/.svn */*/*/.svn */*/*/*/.svn; tar czvf ../release.$(SVNVERSION).tar.gz .; cd ..; rm -rf release

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

