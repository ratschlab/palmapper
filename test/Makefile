dataDir = data
resultDir ?= $(dataDir)
ext =
binDir = ..

Unmapped = split_1m.000.unmapped$(ext)
Mapped = split_1m.000.mapped$(ext)
Spliced = split_1m.000.spliced$(ext)
ResultFiles = $(resultDir)/$(Unmapped) $(resultDir)/$(Spliced) $(resultDir)/$(Unmapped)

DnaFile = $(dataDir)/c_elegans.WS200.dna.fa
SpliceAccFile = $(dataDir)/C_elegans_SpliceSitePred_WS200/acc_pred.bspf
SpliceDonFile = $(dataDir)/C_elegans_SpliceSitePred_WS200/don_pred.bspf
SpliceFiles = $(SpliceAccFile) $(SpliceDonFile)
# Since make is unable to handle tools creating multiple products, we use only one representative
IndexFiles := $(DnaFile).cid # $(DnaFile).mfd $(DnaFile)/mrc $(DnaFile).mta

CurrentDir := $(shell pwd)

inFiles = -i $(dataDir)/c_elegans.WS200.dna.fa -q $(dataDir)/split_1m.000
accDon = -acc $(SpliceAccFile)/contig_%i%c -don $(SpliceDonFile)/contig_%i%c 
outFiles = -o $(resultDir)/$(Mapped) -u $(resultDir)/$(Unmapped) -H $(resultDir)/$(Spliced)

filter = -filter-max-mismatches 3 -filter-max-gaps 0 -filter-splice-min-edit 2 -filter-splice-region 5 -f bedx
qpalma = -qpalma $(dataDir)/parameters.qpalma -qpalma-use-map-max-len 2000
report = -report-map-read -report-spliced-read  -report-map-region -report-splice-sites 0.9
other = -M 6 -G 2 -E 6 -l 18 -L 35 -K 12 -C 55 -I 25000 -NI 2 -SA 5 -CT 10 -z 10 -S -seed-hit-cancel-threshold 10000

allFlags=$(inFiles) $(accDon) $(outFiles) $(filter) $(qpalma) $(report) $(other)

MapCmd = $(binDir)/genomemapper $(allFlags)

dnaFile : $(DnaFile)

$(dataDir)/split_1m.000: 
	bzip2 -d $(dataDir)/split_1m.000.bz2

$(DnaFile): $(SpliceFiles)
	bzip2 -d $(dataDir)/C_elegans_SpliceSitePred_WS200/genome/C_elegans_WS200_genome.fasta.bz2
	mv $(dataDir)/C_elegans_SpliceSitePred_WS200/genome/C_elegans_WS200_genome.fasta.bz2 $(dataDir)/c_elegans.WS200.dna.fa

$(IndexFiles): $(DnaFile) # gmindex
	mkdir -p $(dataDir)
	bash -c "time $(binDir)/gmindex -i $(DnaFile) -v"

$(SpliceFiles):
	mkdir -p $(dataDir)
	wget ftp://ftp.tuebingen.mpg.de/pub/fml/raetsch-lab/predictions/splice/C_elegans/C_elegans_SpliceSitePred_WS200.tar.bz2
	(cd $(dataDir) && tar xvf ../C_elegans_SpliceSitePred_WS200.tar.bz2 && rm ../C_elegans_SpliceSitePred_WS200.tar.bz2) 

map: $(binDir)/genomemapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "time $(MapCmd)"
	(cd $(resultDir) ; md5sum -c $(CurrentDir)/result.md5)

result.md5: $(ResultFiles)
	mkdir -p $(dataDir)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFiles))) > result.md5 
