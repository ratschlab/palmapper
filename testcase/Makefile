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


#Checksum file for each testcase
result = result.md5
resultUnspliced =resultUnspliced.md5
resultPolytrim = resultPolytrim.md5
resultRtrim = resultRtrim.md5
resultNoss = resultNoss.md5
resultBed = resultBed.md5

#Output directory for each testcase
resultDirTest = $(resultDir)/map/
resultDirUnspliced = $(resultDir)/mapUnspliced/
resultDirPolytrim = $(resultDir)/mapPolytrim/
resultDirRtrim = $(resultDir)/mapRtrim/
resultDirNoss = $(resultDir)/mapNoSpliceSitePred/
resultDirBed = $(resultDir)/mapBed/

#Output file names for each testcase
ResultFilesTest = $(resultDirTest)/$(MappedSam) $(resultDirTest)/$(Unmapped) $(resultDirTest)/$(Coverage)
ResultFilesUnspliced = $(resultDirUnspliced)/$(MappedSam) $(resultDirUnspliced)/$(Unmapped) $(resultDirUnspliced)/$(Coverage)
ResultFilesPolytrim = $(resultDirPolytrim)/$(MappedSam) $(resultDirPolytrim)/$(Unmapped) $(resultDirPolytrim)/$(Coverage)
ResultFilesRtrim = $(resultDirRtrim)/$(MappedSam) $(resultDirRtrim)/$(Unmapped) $(resultDirRtrim)/$(Coverage)
ResultFilesNoss = $(resultDirNoss)/$(MappedSam) $(resultDirNoss)/$(Unmapped) $(resultDirNoss)/$(Coverage)
ResultFilesBed = $(resultDirBed)/$(MappedBed) $(resultDirBed)/$(Unmapped) $(resultDirBed)/$(Coverage)


DnaFile = $(dataDir)/c_elegans.WS209.dna.fa
SpliceAccFile = $(dataDir)/C_elegans_SpliceSitePred_WS209/acc_pred.bspf
SpliceDonFile = $(dataDir)/C_elegans_SpliceSitePred_WS209/don_pred.bspf
# Since make is unable to handle tools creating multiple products, we use only one representative
SpliceFiles = $(SpliceAccFile) # $(SpliceDonFile)
# Since make is unable to handle tools creating multiple products, we use only one representative
IndexFiles := $(DnaFile).cid # $(DnaFile).mfd $(DnaFile)/mrc $(DnaFile).mta
BwaIndexFiles := $(DnaFile).bwt # $(DnaFile).mfd $(DnaFile)/mrc $(DnaFile).mta

CurrentDir := $(shell pwd)

#PALMapper parameters
inFiles = -i $(dataDir)/c_elegans.WS209.dna.fa -q $(dataDir)/split_1m.000
accDon = -acc $(SpliceAccFile)/contig_%i%c -don $(SpliceDonFile)/contig_%i%c 
outSam= -f sam 
outBed= -f bedx
qpalma =  -L 35 -K 12 -C 55 -S -qpalma $(dataDir)/parameters.qpalma 
other = -M 4 -G 1 -E 4 -z 1 -l 20 -seed-hit-cancel-threshold 1000
polytrim = -polytrim 40
rtrim = -rtrim 25
report =
threads = -threads $(nrThreads)

#Parameter base for all testcases
allFlags=$(inFiles) $(other) $(threads) $(report) -rlim 10000
MapCmd = $(binDir)/palmapper $(allFlags) 

#outFiles =  -u $(resultDir)/$(Unmapped) #-H $(resultDir)/$(Spliced)
#-o $(resultDir)/$(MappedSam)
#outCoverage =  -report-coverage-wig $(resultDir)/$(Coverage)

test: map

mapall: map unsplicedmap polytrimmap rtrimmap nossmap bedmap
bwamapall: bwamap bwaunsplicedmap bwapolytrimmap bwartrimmap bwanossmap bwabedmap

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

clean:
	rm -rf $(dataDir)/C_elegans_SpliceSitePred_WS209/ $(dataDir)/c_elegans.WS209.dna.fa* $(dataDir)/$(resultDirTest) $(dataDir)/$(resultDirUnspliced) $(dataDir)/$(resultDirRtrim) $(dataDir)/$(resultDirPolytrim) $(dataDir)/$(resultDirNoss) $(dataDir)/$(resultDirBed)

allmd5: result.md5 resultUnspliced.md5 resultPolytrim.md5 resultRtrim.md5 resultNoss.md5 resultBed.md5

#Checksum file creation rules
result.md5: $(ResultFilesTest)
	mkdir -p $(resultDirTest)
	(cd $(resultDirTest) ; md5sum $(notdir $(ResultFilesTest))) > $(result)

resultUnspliced.md5: $(ResultFilesUnspliced)
	mkdir -p $(resultDirUnspliced)
	(cd $(resultDirUnspliced) ; md5sum $(notdir $(ResultFilesUnspliced))) > $(resultUnspliced)

resultRtrim.md5: $(ResultFilesRtrim)
	mkdir -p $(resultDirRtrim)
	(cd $(resultDirRtrim) ; md5sum $(notdir $(ResultFilesRtrim))) > $(resultRtrim) 

resultPolytrim.md5: $(ResultFilesPolytrim)
	mkdir -p $(resultDirPolytrim)
	(cd $(resultDirPolytrim) ; md5sum $(notdir $(ResultFilesPolytrim))) >  $(resultPolytrim) 

resultNoss.md5: $(ResultFilesNoss)
	mkdir -p $(resultDirNoss)
	(cd $(resultDirNoss) ; md5sum $(notdir $(ResultFilesNoss))) >  $(resultNoss) 

resultBed.md5: $(ResultFilesBed)
	mkdir -p $(resultDirBed)
	(cd $(resultDirBed) ; md5sum $(notdir $(ResultFilesBed))) >  $(resultBed) 


#Testcases with array-based index
map: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(resultDirTest)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) -u $(resultDirTest)/$(Unmapped) -o $(resultDirTest)/$(MappedSam) -report-coverage-wig $(resultDirTest)/$(Coverage)"
	(cd $(resultDirTest) ; md5sum $(notdir $(ResultFilesTest))) 
	(cd $(resultDirTest) ; md5sum -c $(CurrentDir)/$(result)) || true

unsplicedmap: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(resultDirUnspliced)
	bash -c "time $(MapCmd) $(outSam) -u $(resultDirUnspliced)/$(Unmapped) -o $(resultDirUnspliced)/$(MappedSam) -report-coverage-wig $(resultDirUnspliced)/$(Coverage)"
	(cd $(resultDirUnspliced) ; md5sum $(notdir $(ResultFilesUnspliced))) 
	(cd $(resultDirUnspliced) ; md5sum -c $(CurrentDir)/$(resultUnspliced)) || true

rtrimmap: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(resultDirRtrim)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) $(rtrim) -u $(resultDirRtrim)/$(Unmapped) -o $(resultDirRtrim)/$(MappedSam) -report-coverage-wig $(resultDirRtrim)/$(Coverage)"
	(cd $(resultDirRtrim) ; md5sum $(notdir $(ResultFilesRtrim))) 
	(cd $(resultDirRtrim) ; md5sum -c $(CurrentDir)/$(resultRtrim)) || true

polytrimmap: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(resultDirPolytrim)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) $(polytrim) -u $(resultDirPolytrim)/$(Unmapped) -o $(resultDirPolytrim)/$(MappedSam) -report-coverage-wig $(resultDirPolytrim)/$(Coverage)"
	(cd $(resultDirPolytrim) ; md5sum $(notdir $(ResultFilesPolytrim))) 
	(cd $(resultDirPolytrim) ; md5sum -c $(CurrentDir)/$(resultPolytrim)) || true

nossmap: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(resultDirNoss)
	bash -c "time $(MapCmd) -no-ss-pred $(outSam) $(qpalma) -u $(resultDirNoss)/$(Unmapped) -o $(resultDirNoss)/$(MappedSam) -report-coverage-wig $(resultDirNoss)/$(Coverage)"
	(cd $(resultDirNoss) ; md5sum $(notdir $(ResultFilesNoss))) 
	(cd $(resultDirNoss) ; md5sum -c $(CurrentDir)/$(resultNoss)) || true

bedmap: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(resultDirBed)
	bash -c "time $(MapCmd) $(accDon)  $(outBed) $(qpalma) -u $(resultDirBed)/$(Unmapped) -o $(resultDirBed)/$(MappedBed) -report-coverage-wig $(resultDirBed)/$(Coverage)"
	(cd $(resultDirBed) ; md5sum $(notdir $(ResultFilesBed))) 
	(cd $(resultDirBed) ; md5sum -c $(CurrentDir)/$(resultBed)) || true

#Testcases with bwa-based index
bwamap: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(resultDirTest)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) -u $(resultDirTest)/$(Unmapped) -o $(resultDirTest)/$(MappedSam) -report-coverage-wig $(resultDirTest)/$(Coverage) -bwa"
	(cd $(resultDirTest) ; md5sum $(notdir $(ResultFilesTest))) 
	(cd $(resultDirTest) ; md5sum -c $(CurrentDir)/$(result)) || true

bwaunsplicedmap: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(resultDirUnspliced)
	bash -c "time $(MapCmd) $(outSam) -u $(resultDirUnspliced)/$(Unmapped) -o $(resultDirUnspliced)/$(MappedSam) -report-coverage-wig $(resultDirUnspliced)/$(Coverage) -bwa"
	(cd $(resultDirUnspliced) ; md5sum $(notdir $(ResultFilesUnspliced))) 
	(cd $(resultDirUnspliced) ; md5sum -c $(CurrentDir)/$(resultUnspliced)) || true

bwartrimmap: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(resultDirRtrim)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) $(rtrim) -u $(resultDirRtrim)/$(Unmapped) -o $(resultDirRtrim)/$(MappedSam) -report-coverage-wig $(resultDirRtrim)/$(Coverage) -bwa"
	(cd $(resultDirRtrim) ; md5sum $(notdir $(ResultFilesRtrim))) 
	(cd $(resultDirRtrim) ; md5sum -c $(CurrentDir)/$(resultRtrim)) || true

bwapolytrimmap: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(resultDirPolytrim)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) $(polytrim) -u $(resultDirPolytrim)/$(Unmapped) -o $(resultDirPolytrim)/$(MappedSam) -report-coverage-wig $(resultDirPolytrim)/$(Coverage) -bwa"
	(cd $(resultDirPolytrim) ; md5sum $(notdir $(ResultFilesPolytrim))) 
	(cd $(resultDirPolytrim) ; md5sum -c $(CurrentDir)/$(resultPolytrim)) || true

bwanossmap: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(resultDirNoss)
	bash -c "time $(MapCmd) -no-ss-pred $(outSam) $(qpalma) -u $(resultDirNoss)/$(Unmapped) -o $(resultDirNoss)/$(MappedSam) -report-coverage-wig $(resultDirNoss)/$(Coverage) -bwa"
	(cd $(resultDirNoss) ; md5sum $(notdir $(ResultFilesNoss))) 
	(cd $(resultDirNoss) ; md5sum -c $(CurrentDir)/$(resultNoss)) || true

bwabedmap: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(resultDirBed)
	bash -c "time $(MapCmd) $(accDon)  $(outBed) $(qpalma) -u $(resultDirBed)/$(Unmapped) -o $(resultDirBed)/$(MappedBed) -report-coverage-wig $(resultDirBed)/$(Coverage) -bwa"
	(cd $(resultDirBed) ; md5sum $(notdir $(ResultFilesBed))) 
	(cd $(resultDirBed) ; md5sum -c $(CurrentDir)/$(resultBed)) || true


