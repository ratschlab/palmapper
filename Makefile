SrcDir = ./src
OutDir = ./out
ObjDir = $(OutDir)/o
ExeDir = $OutDir/a

CC = g++
# CFLAGS = -Wall -g 
#CFLAGS = -Wno-sign-compare -Wno-unused-parameter -Wformat -Wformat-security -Wimplicit -Wparentheses -Wshadow -O9 -fexpensive-optimizations -frerun-cse-after-loop -fcse-follow-jumps -finline-functions -fschedule-insns2 -fthread-jumps -fforce-addr -fstrength-reduce -funroll-loops -march=native -mtune=native -pthread
CFLAGS = -Wall -g -Wno-unused-parameter -Wformat -Wformat-security -Wimplicit -Wparentheses -Wshadow 
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

GM_OBJ = $(ObjDir)/genomemapper/report_maps.o \
	$(ObjDir)/genomemapper/qpalma_alignment.o \
	$(ObjDir)/genomemapper/qpalma_filter.o \
	$(ObjDir)/genomemapper/qpalma_init.o \
	$(ObjDir)/genomemapper/align.o \
	$(ObjDir)/genomemapper/alloc.o \
	$(ObjDir)/genomemapper/TopAlignments.o \
	$(ObjDir)/genomemapper/IntervalQuery.o \
	$(ObjDir)/genomemapper/genomemapper.o \
	$(ObjDir)/genomemapper/index.o \
	$(ObjDir)/genomemapper/init.o \
	$(ObjDir)/genomemapper/print.o \
	$(ObjDir)/genomemapper/Chromosome.o \
	$(ObjDir)/genomemapper/Config.o \
	$(ObjDir)/genomemapper/Genome.o \
	$(ObjDir)/genomemapper/Hit.o \
	$(ObjDir)/genomemapper/Read.o \
	$(ObjDir)/genomemapper/Statistics.o \
	$(ObjDir)/genomemapper/Util.o \
	$(SHOGUN_OBJ) $(DYNPROG_OBJ)

#	$(ObjDir)/genomemapper/Aligner.o \
#	$(ObjDir)/genomemapper/Config.o \
#	$(ObjDir)/genomemapper/Read.o \

IDX_OBJ = $(ObjDir)/mkindex/init.o \
	$(ObjDir)/mkindex/printindex.o \
	$(ObjDir)/mkindex/usage.o \
	$(ObjDir)/mkindex/write.o \
	$(ObjDir)/mkindex/load.o \
	$(ObjDir)/mkindex/index.o \
	$(ObjDir)/mkindex/alloc.o \
	$(ObjDir)/mkindex/mkindex.o

genomemapper: $(GM_OBJ) src/genomemapper/*.h src/genomemapper/genomemapper_symbols.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o genomemapper $(GM_OBJ) -lpthread -lz -lm

gmindex:  $(IDX_OBJ) src/mkindex/*.h src/mkindex/mkindex_symbols.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o gmindex $(IDX_OBJ) 

clean:
	rm -rf $(OutDir) genomemapper gmindex

all: genomemapper gmindex

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
	
#src/genomemapper/%.o: src/genomemapper/%.c
#	$(CC) $(CFLAGS) $(INCLUDE) -c $? -o $@
#
#src/genomemapper/%.o: src/genomemapper/%.cpp
#	$(CC) $(CFLAGS) $(INCLUDE) -c $? -o $@
#
#src/mkindex/%.o: src/mkindex/%.c
#	$(CC) $(CFLAGS) $(INCLUDE) -c $? -o $@


