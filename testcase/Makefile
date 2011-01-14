dataDir = data
resultDir ?= $(dataDir)
ext = 
binDir = ..
#nrThreads ?= $(shell grep processor /proc/cpuinfo | wc -l)
nrThreads=1

Unmapped = split_1m.000.unmapped$(ext)
MappedSam = split_1m.000.sam$(ext)
MappedBed = split_1m.000.bedx$(ext)
Coverage = split_1m.000.wiggle$(ext)
#Spliced = split_1m.000.spliced$(ext)
ResultFilesSam = $(resultDir)/$(MappedSam) $(resultDir)/$(Unmapped) $(resultDir)/$(Coverage)
ResultFilesBed = $(resultDir)/$(MappedBed) $(resultDir)/$(Unmapped) $(resultDir)/$(Coverage)

result = result.md5
resultUnspliced = resultUnspliced.md5
resultRtrim = resultRtrim.md5
resultPolytrim = resultPolytrim.md5
resultNoss = resultNoss.md5
resultBed = resultBed.md5

DnaFile = $(dataDir)/c_elegans.WS209.dna.fa
SpliceAccFile = $(dataDir)/C_elegans_SpliceSitePred_WS209/acc_pred.bspf
SpliceDonFile = $(dataDir)/C_elegans_SpliceSitePred_WS209/don_pred.bspf
# Since make is unable to handle tools creating multiple products, we use only one representative
SpliceFiles = $(SpliceAccFile) # $(SpliceDonFile)
# Since make is unable to handle tools creating multiple products, we use only one representative
IndexFiles := $(DnaFile).cid # $(DnaFile).mfd $(DnaFile)/mrc $(DnaFile).mta
BwaIndexFiles := $(DnaFile).bwt # $(DnaFile).mfd $(DnaFile)/mrc $(DnaFile).mta

CurrentDir := $(shell pwd)

inFiles = -i $(dataDir)/c_elegans.WS209.dna.fa -q $(dataDir)/split_1m.000
accDon = -acc $(SpliceAccFile)/contig_%i%c -don $(SpliceDonFile)/contig_%i%c 
outFiles =  -u $(resultDir)/$(Unmapped) #-H $(resultDir)/$(Spliced)
outSam= -f sam -o $(resultDir)/$(MappedSam)
outBed= -f bedx -o $(resultDir)/$(MappedBed)
outCoverage =  -report-coverage-wig $(resultDir)/$(Coverage)
qpalma =  -L 35 -K 12 -C 55 -S -qpalma $(dataDir)/parameters.qpalma 
other = -M 4 -G 1 -E 4 -z 1 -l 20 -seed-hit-cancel-threshold 1000
polytrim = -polytrim 40
rtrim = -rtrim 10
report =
threads = -threads $(nrThreads)

allFlags=$(inFiles) $(outFiles) $(outCoverage) $(other) $(threads) $(report)
MapCmd = $(binDir)/palmapper $(allFlags) #-no-ss-pred -min-spliced-segment-len 20  -rlim 10000  -don-consensus ?? -acc-consensus ??

test: map

$(dataDir)/split_1m.000: 
	bzip2 -d $(dataDir)/split_1m.000.bz2

$(DnaFile): $(SpliceFiles)

$(IndexFiles): $(DnaFile) # pmindex
	mkdir -p $(dataDir)
	bash -c "time $(binDir)/pmindex -i $(DnaFile) -v"

$(BwaIndexFiles): $(DnaFile) # pmindex
	mkdir -p $(dataDir)
	bash -c "time $(binDir)/bwa index $(DnaFile) "

$(SpliceFiles):
	mkdir -p $(dataDir)
	wget -c ftp://ftp.tuebingen.mpg.de/pub/fml/raetsch-lab/predictions/splice/C_elegans/C_elegans_SpliceSitePred_WS209.tar.bz2
	(cd $(dataDir) && tar xvf ../C_elegans_SpliceSitePred_WS209.tar.bz2 && rm ../C_elegans_SpliceSitePred_WS209.tar.bz2) 
	bzip2 -d $(dataDir)/C_elegans_SpliceSitePred_WS209/genome/C_elegans_WS209_genome.fasta.bz2
	mv $(dataDir)/C_elegans_SpliceSitePred_WS209/genome/C_elegans_WS209_genome.fasta $(DnaFile)

valgrind: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "valgrind $(MapCmd)"

cachegrind: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "valgrind --tool=cachegrind $(MapCmd)"

result.md5: $(ResultFilesSam)
	mkdir -p $(dataDir)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesSam))) > $(result)


resultUnspliced.md5: $(ResultFilesSam)
	mkdir -p $(dataDir)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesSam))) > $(resultUnspliced)


resultRtrim.md5: $(ResultFilesSam)
	mkdir -p $(dataDir)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesSam))) > $(resultRtrim) 


resultPolytrim.md5: $(ResultFilesSam)
	mkdir -p $(dataDir)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesSam))) >  $(resultPolytrim) 


resultNoss.md5: $(ResultFilesSam)
	mkdir -p $(dataDir)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesSam))) >  $(resultNoss) 

resultBed.md5: $(ResultFilesBed)
	mkdir -p $(dataDir)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesBed))) >  $(resultBed) 

clean:
	rm -rf $(dataDir)/C_elegans_SpliceSitePred_WS209/ $(dataDir)/c_elegans.WS209.dna.fa* $(dataDir)/split_1m.000.*

map: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma)"
#	valgrind $(MapCmd)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesSam))) 
	(cd $(resultDir) ; md5sum -c $(CurrentDir)/$(result))

unsplicedmap:$(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "time $(MapCmd)  $(outSam)"
#	valgrind $(MapCmd)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesSam))) 
	(cd $(resultDir) ; md5sum -c $(CurrentDir)/$(resultUnspliced))

rtrimmap: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) $(rtrim)"
#	valgrind $(MapCmd)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesSam))) 
	(cd $(resultDir) ; md5sum -c $(CurrentDir)/$(resultRtrim))

polytrimmap: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) $(polytrim)"
#	valgrind $(MapCmd)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesSam))) 
	(cd $(resultDir) ; md5sum -c $(CurrentDir)/$(resultPolytrim))

nossmap: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "time $(MapCmd) $(outSam) $(qpalma) -no-ss-pred"
#	valgrind $(MapCmd)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesSam))) 
	(cd $(resultDir) ; md5sum -c $(CurrentDir)/$(resultNoss))

bedmap: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "time $(MapCmd)  $(accDon) $(outBed) $(qpalma)"
#	valgrind $(MapCmd)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesBed))) 
	(cd $(resultDir) ; md5sum -c $(CurrentDir)/$(resultBed))


bwamap: $(binDir)/palmapper $(DnaFile) $(BwaIndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma)"
#	valgrind $(MapCmd)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesSam))) 
	(cd $(resultDir) ; md5sum -c $(CurrentDir)/$(result))

bwaunsplicedmap: $(binDir)/palmapper $(DnaFile) $(BwaIndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "time $(MapCmd)  $(outSam) -bwa"
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesSam))) 
	(cd $(resultDir) ; md5sum -c $(CurrentDir)/$(resultUnspliced))

bwartrimmap: $(binDir)/palmapper $(DnaFile) $(BwaIndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) $(rtrim) -bwa"
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesSam))) 
	(cd $(resultDir) ; md5sum -c $(CurrentDir)/$(resultRtrim))

bwapolytrimmap: $(binDir)/palmapper $(DnaFile) $(BwaIndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) $(polytrim) -bwa"
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesSam))) 
	(cd $(resultDir) ; md5sum -c $(CurrentDir)/$(resultPolytrim))

bwanossmap: $(binDir)/palmapper $(DnaFile) $(BwaIndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "time $(MapCmd) $(outSam) $(qpalma) -no-ss-pred -bwa"
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesSam))) 
	(cd $(resultDir) ; md5sum -c $(CurrentDir)/$(resultNoss))

bwabedmap: $(binDir)/palmapper $(DnaFile) $(BwaIndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "time $(MapCmd)  $(accDon) $(outBed) $(qpalma) -bwa"
	(cd $(resultDir) ; md5sum $(notdir $(ResultFilesBed))) 
	(cd $(resultDir) ; md5sum -c $(CurrentDir)/$(resultBed))