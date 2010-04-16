SrcDir = ./src
OutDir = ./out
ObjDir = $(OutDir)/o
ExeDir = $OutDir/a

CC = g++
#CFLAGS = -Wall -g # debug
CFLAGS = -Wno-unused-parameter -Wformat -Wformat-security -Wimplicit -Wparentheses -Wshadow -O9 -fexpensive-optimizations -frerun-cse-after-loop -fcse-follow-jumps -finline-functions -fschedule-insns2 -fthread-jumps -fforce-addr -fstrength-reduce -funroll-loops -march=native -mtune=native -pthread # linux amd64 optimized
#CFLAGS = -O9 -Wall -g -Wno-unused-parameter -Wformat -Wformat-security -Wimplicit -Wparentheses -Wshadow # generic
INCLUDE =  -Ishogun/ -Idyn_prog/ -Isrc
LDFLAGS = 

SHOGUN_OBJ = $(ObjDir)/genomemapper/shogun/init.o \
	$(ObjDir)/genomemapper/shogun/Mathematics.o \
	$(ObjDir)/genomemapper/shogun/io.o \
	$(ObjDir)/genomemapper/shogun/Parallel.o \
	$(ObjDir)/genomemapper/shogun/Version.o \
	$(ObjDir)/genomemapper/shogun/SGObject.o \
	$(ObjDir)/genomemapper/shogun/ShogunException.o \
	$(ObjDir)/genomemapper/shogun/Signal.o

DYNPROG_OBJ = $(ObjDir)/genomemapper/dyn_prog/Mathmatics_dp.o \
	$(ObjDir)/genomemapper/dyn_prog/io_dp.o \
	$(ObjDir)/genomemapper/dyn_prog/qpalma_dp.o \
	$(ObjDir)/genomemapper/dyn_prog/debug_tools.o \
	$(ObjDir)/genomemapper/dyn_prog/penalty_info_dp.o \
	$(ObjDir)/genomemapper/dyn_prog/result_align.o \
	$(ObjDir)/genomemapper/dyn_prog/fill_matrix.o

GM_OBJ = $(ObjDir)/genomemapper/GenomeMaps.o \
	$(ObjDir)/genomemapper/QPalma.o \
	$(ObjDir)/genomemapper/align.o \
	$(ObjDir)/genomemapper/TopAlignments.o \
	$(ObjDir)/genomemapper/IntervalQuery.o \
	$(ObjDir)/genomemapper/genomemapper.o \
	$(ObjDir)/genomemapper/init.o \
	$(ObjDir)/genomemapper/print.o \
	$(ObjDir)/genomemapper/Chromosome.o \
	$(ObjDir)/genomemapper/Config.o \
	$(ObjDir)/genomemapper/Genome.o \
	$(ObjDir)/genomemapper/Hits.o \
	$(ObjDir)/genomemapper/Read.o \
	$(ObjDir)/genomemapper/Statistics.o \
	$(ObjDir)/genomemapper/Util.o \
	$(SHOGUN_OBJ) $(DYNPROG_OBJ)

IDX_OBJ = $(ObjDir)/pmindex/init.o \
	$(ObjDir)/pmindex/printindex.o \
	$(ObjDir)/pmindex/usage.o \
	$(ObjDir)/pmindex/write.o \
	$(ObjDir)/pmindex/load.o \
	$(ObjDir)/pmindex/index.o \
	$(ObjDir)/pmindex/alloc.o \
	$(ObjDir)/pmindex/pmindex.o

CurrentDir := $(shell pwd)

all: palmapper pmindex

palmapper: $(GM_OBJ) src/genomemapper/*.h 
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o palmapper $(GM_OBJ) -lpthread -lz -lm

pmindex:  $(IDX_OBJ) src/pmindex/*.h src/pmindex/pmindex_symbols.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o pmindex $(IDX_OBJ) 

clean:
	rm -rf $(OutDir) palmapper pmindex

test:
	(cd testcase; make test)

release:
	make clean; mkdir -p ../release; cd ..; rsync -av $(CurrentDir) release; cd release; rm -rf */.settings */.cproject */.project */.svn */*/.svn */*/*/.svn; tar czvf ../release.tar.gz .; cd ..; rm -rf release

# generic rule for compiling c++
$(ObjDir)/%.o : $(SrcDir)/%.cpp
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

$(ObjDir)/%.o : $(SrcDir)/%.cpp $(SrcDir)/%.h
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

$(ObjDir)/%.o : $(SrcDir)/%.c
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<
