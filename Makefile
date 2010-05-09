SrcDir = ./src
OutDir = ./out
ObjDir = $(OutDir)/o
ExeDir = $OutDir/a

CC = g++
#CFLAGS = -Wall -g # debug
#CFLAGS = -Wno-unused-parameter -Wformat -Wformat-security -Wimplicit -Wparentheses -Wshadow -O9 -fexpensive-optimizations -frerun-cse-after-loop -fcse-follow-jumps -finline-functions -fschedule-insns2 -fthread-jumps -fforce-addr -fstrength-reduce -funroll-loops -march=native -mtune=native -pthread # linux amd64 optimized
CFLAGS = -O9 -Wall -g -Wno-unused-parameter -Wformat -Wformat-security -Wimplicit -Wparentheses -Wshadow # generic
INCLUDE =  -Ishogun/ -Idyn_prog/ -Isrc
LDFLAGS = 

SHOGUN_OBJ = $(ObjDir)/palmapper/shogun/init.o \
	$(ObjDir)/palmapper/shogun/Mathematics.o \
	$(ObjDir)/palmapper/shogun/io.o \
	$(ObjDir)/palmapper/shogun/Parallel.o \
	$(ObjDir)/palmapper/shogun/Version.o \
	$(ObjDir)/palmapper/shogun/SGObject.o \
	$(ObjDir)/palmapper/shogun/ShogunException.o \
	$(ObjDir)/palmapper/shogun/Signal.o

DYNPROG_OBJ = $(ObjDir)/palmapper/dyn_prog/Mathmatics_dp.o \
	$(ObjDir)/palmapper/dyn_prog/io_dp.o \
	$(ObjDir)/palmapper/dyn_prog/qpalma_dp.o \
	$(ObjDir)/palmapper/dyn_prog/debug_tools.o \
	$(ObjDir)/palmapper/dyn_prog/penalty_info_dp.o \
	$(ObjDir)/palmapper/dyn_prog/result_align.o \
	$(ObjDir)/palmapper/dyn_prog/fill_matrix.o

GM_OBJ = $(ObjDir)/palmapper/GenomeMaps.o \
	$(ObjDir)/palmapper/QPalma.o \
	$(ObjDir)/palmapper/align.o \
	$(ObjDir)/palmapper/TopAlignments.o \
	$(ObjDir)/palmapper/IntervalQuery.o \
	$(ObjDir)/palmapper/palmapper.o \
	$(ObjDir)/palmapper/init.o \
	$(ObjDir)/palmapper/print.o \
	$(ObjDir)/palmapper/Chromosome.o \
	$(ObjDir)/palmapper/Config.o \
	$(ObjDir)/palmapper/Genome.o \
	$(ObjDir)/palmapper/Hits.o \
	$(ObjDir)/palmapper/Read.o \
	$(ObjDir)/palmapper/Statistics.o \
	$(ObjDir)/palmapper/Util.o \
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

palmapper: $(GM_OBJ) src/palmapper/*.h 
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
