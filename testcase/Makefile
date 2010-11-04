dataDir = data
resultDir ?= $(dataDir)
ext = 
binDir = ..
#nrThreads ?= $(shell grep processor /proc/cpuinfo | wc -l)
nrThreads=1

Unmapped = split_1m.000.unmapped$(ext)
Mapped = split_1m.000.sam$(ext)
Coverage = split_1m.000.wiggle$(ext)
#Spliced = split_1m.000.spliced$(ext)
ResultFiles = $(resultDir)/$(Mapped) $(resultDir)/$(Unmapped) $(resultDir)/$(Coverage)

DnaFile = $(dataDir)/c_elegans.WS209.dna.fa
SpliceAccFile = $(dataDir)/C_elegans_SpliceSitePred_WS209/acc_pred.bspf
SpliceDonFile = $(dataDir)/C_elegans_SpliceSitePred_WS209/don_pred.bspf
# Since make is unable to handle tools creating multiple products, we use only one representative
SpliceFiles = $(SpliceAccFile) # $(SpliceDonFile)
# Since make is unable to handle tools creating multiple products, we use only one representative
IndexFiles := $(DnaFile).cid # $(DnaFile).mfd $(DnaFile)/mrc $(DnaFile).mta

CurrentDir := $(shell pwd)

inFiles = -i $(dataDir)/c_elegans.WS209.dna.fa -q $(dataDir)/split_1m.000
accDon = -acc $(SpliceAccFile)/contig_%i%c -don $(SpliceDonFile)/contig_%i%c 
outFiles = -o $(resultDir)/$(Mapped) -u $(resultDir)/$(Unmapped) #-H $(resultDir)/$(Spliced)

output = -f sam -report-coverage-wig $(resultDir)/$(Coverage)
qpalma = -qpalma $(dataDir)/parameters.qpalma 
report = 
other = -M 4 -G 1 -E 4 -l 20 -L 35 -K 12 -C 55 -z 1 -S -seed-hit-cancel-threshold 1000
threads = -threads $(nrThreads)

allFlags=$(inFiles) $(accDon) $(outFiles) $(output) $(qpalma) $(report) $(other) $(threads)

MapCmd = $(binDir)/palmapper $(allFlags) #-no-ss-pred -min-spliced-segment-len 20  -rlim 10000  -don-consensus ?? -acc-consensus ??

test: map

$(dataDir)/split_1m.000: 
	bzip2 -d $(dataDir)/split_1m.000.bz2

$(DnaFile): $(SpliceFiles)

$(IndexFiles): $(DnaFile) # pmindex
	mkdir -p $(dataDir)
	bash -c "time $(binDir)/pmindex -i $(DnaFile) -v"

$(SpliceFiles):
	mkdir -p $(dataDir)
	wget -c ftp://ftp.tuebingen.mpg.de/pub/fml/raetsch-lab/predictions/splice/C_elegans/C_elegans_SpliceSitePred_WS209.tar.bz2
	(cd $(dataDir) && tar xvf ../C_elegans_SpliceSitePred_WS209.tar.bz2 && rm ../C_elegans_SpliceSitePred_WS209.tar.bz2) 
	bzip2 -d $(dataDir)/C_elegans_SpliceSitePred_WS209/genome/C_elegans_WS209_genome.fasta.bz2
	mv $(dataDir)/C_elegans_SpliceSitePred_WS209/genome/C_elegans_WS209_genome.fasta $(DnaFile)

map: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "time $(MapCmd)"
#	valgrind $(MapCmd)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFiles))) 
	(cd $(resultDir) ; md5sum -c $(CurrentDir)/result.md5)

result.md5: $(ResultFiles)
	mkdir -p $(dataDir)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFiles))) > result.md5 

clean:
	rm -rf $(dataDir)/C_elegans_SpliceSitePred_WS209/ $(dataDir)/c_elegans.WS209.dna.fa* $(dataDir)/split_1m.000.*
