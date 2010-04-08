dataDir = ../data
resultDir ?= $(dataDir)
ext =
binDir = ..

Unmapped = split_1m.000.unmapped$(ext)
Mapped = split_1m.000.mapped$(ext)
Spliced = split_1m.000.spliced$(ext)
ResultFiles = $(resultDir)/$(Unmapped) $(resultDir)/$(Spliced) $(resultDir)/$(Unmapped)

DnaFile = $(dataDir)/c_elegans.WS200.dna.fa
# Since make is unable to handle tools creating multiple products, we use only one representative
IndexFiles := $(DnaFile).cid # $(DnaFile).mfd $(DnaFile)/mrc $(DnaFile).mta

CurrentDir := $(shell pwd)

inFiles = -i $(dataDir)/c_elegans.WS200.dna.fa -q $(dataDir)/split_1m.000
accDon = -acc $(dataDir)/acc_pred.bspf/pred/contig_%i%c -don $(dataDir)/don_pred.bspf/pred/contig_%i%c 
outFiles = -o $(resultDir)/$(Mapped) -u $(resultDir)/$(Unmapped) -H $(resultDir)/$(Spliced)

filter = -filter-max-mismatches 3 -filter-max-gaps 0 -filter-splice-min-edit 2 -filter-splice-region 5 -f bedx
qpalma = -qpalma $(dataDir)/parameters_final.qpalma -qpalma-use-map-max-len 2000
report = -report-map-read -report-spliced-read  -report-map-region -report-splice-sites 0.9
other = -M 6 -G 2 -E 6 -l 18 -L 35 -K 12 -C 55 -I 25000 -NI 2 -SA 5 -CT 10 -z 10 -S -seed-hit-cancel-threshold 10000

allFlags=$(inFiles) $(accDon) $(outFiles) $(filter) $(qpalma) $(report) $(other)

MapCmd = $(binDir)/genomemapper $(allFlags)

dnaFile : $(DnaFile)

$(DnaFile) :
	wget -O $(DnaFile).gz ftp://ftp.wormbase.org/pub/wormbase/genomes/c_elegans/sequences/dna/c_elegans.WS200.dna.fa.gz
	gunzip $(DnaFile).gz

$(IndexFiles): $(DnaFile) # gmindex
	bash -c "time $(binDir)/gmindex -i $(DnaFile) -v"

map: $(binDir)/genomemapper $(DnaFile) $(IndexFiles)
	bash -c "time $(MapCmd)"
	(cd $(resultDir) ; md5sum -c $(CurrentDir)/result.md5)
	

result.md5: $(ResultFiles)
	(cd $(resultDir) ; md5sum $(notdir $(ResultFiles))) > result.md5 