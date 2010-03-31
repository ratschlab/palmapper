
There are some useful defines that are useful to save memory:
It's tricky, but necessary to get it to run efficiently on 32GB RAM
(total saving 7GB/4 + 2*3GB*3/4: 6.25GB)

in interval_query.cpp (7GB):
----------------------------

#define CONVERT_SCORE_MAPS // save 25% 

in genome.h (3GB):
------------------

#define USE_CHR_BIN 

and

#define USE_CHR_BIN_CLASS CDNAArray2 // saves 50%
or
#define USE_CHR_BIN_CLASS CDNAArray4 // saves 75%, but non-ACGT letters are converted to A


in report_maps.h (3GB):
-----------------------

#define CHR_MAP_DNAARRAY

and

#define CHR_MAP_DNAARRAY_CLASS CDNAArray2 // saves 50%, no mapped_region reporting
or
#define CHR_MAP_DNAARRAY_CLASS CDNAArray4 // save 75%, only spliced_reads_best and mapped_reads_best
#define CHR_MAP_DNAARRAY_2BIT // this has to be switched on for compatibility with saved maps


without optimizations (index depth 13; node148):
37.8GB vmem
real    4m25.004s
user    3m37.286s
sys     1m3.848s
# 336 (100) unspliced, 54 (30) spliced alignments, 410 unmapped (spliced 13.8%, unmapped 51.2%)

with all optimizations (index depth 13; node148):
31.8GB vmem
real    6m8.100s
user    5m16.932s
sys     1m6.396s
# 418 (128) unspliced, 65 (36) spliced alignments, 505 unmapped (spliced 13.5%, unmapped 51.1%)
# [capture_hits] timing: 3.8900, 33.5600, 99.6400 (42.3200 ss access; 4779nt and 3.0 threads per alignment)

with all optimizations (index depth 12; node182):
31.1GB vmem
real    4m56.098s
user    4m16.616s
sys     0m48.299s
# 343 (106) unspliced, 57 (27) spliced alignments, 437 unmapped (spliced 14.2%, unmapped 52.2%)


some other tests:
-----------------

1. 
with -seed-hit-cancel-threshold 10000 (was -seed-hit-cancel-threshold 1000)

with all optimizations (index depth 12; node182):
30.9GB vmem
real    5m55.345s
user    4m32.549s
sys     0m49.803s
# 409 (124) unspliced, 65 (37) spliced alignments, 466 unmapped (spliced 13.7%, unmapped 49.6%)
# [capture_hits] timing: 3.3800, 29.3800, 61.8900 (34.3400 ss access; 4508nt and 3.2 threads per alignment)

2. 
with -report-ro ../GM12878_2x75.map.us (10.5%)
with -seed-hit-cancel-threshold 10000

with all optimizations (index depth 13; node148):
31.8GB vmem
real    6m7.464s
user    6m39.561s
sys     1m49.483s
# 416 (127) unspliced, 90 (48) spliced alignments, 463 unmapped (spliced 17.8%, unmapped 47.8%)
# [capture_hits] timing: 4.8000, 227.0200, 1128.2300 (47.9800 ss access; 25963nt and 3.1 threads per alignment)


with all optimizations (index depth 12; node182):
real    4m54.695s
user    5m9.199s
sys     1m8.444s
# 394 (119) unspliced, 85 (45) spliced alignments, 437 unmapped (spliced 17.7%, unmapped 47.7%)
# [capture_hits] timing: 3.9700, 168.3700, 581.5800 (30.9800 ss access; 25582nt and 3.1 threads per alignment)


3. 
with -qpalma-use-map-max-len 5000 (was -qpalma-use-map-max-len 10000)

with all optimizations (index depth 12; node182):
real    5m37.613s
user    5m29.557s
sys     0m52.559s
# 411 (126) unspliced, 85 (44) spliced alignments, 453 unmapped (spliced 17.1%, unmapped 47.7%)
# [capture_hits] timing: 3.7800, 91.0900, 249.3900 (34.6800 ss access; 16414nt and 3.1 threads per alignment)

4. 
with -rlim 5000 (was -rlim 1000)

with all optimizations (index depth 13; node148):
32.2GB vmem
real    15m22.346s
user    26m49.621s
sys     6m28.120s
# 2092 (699) unspliced, 464 (250) spliced alignments, 2372 unmapped (spliced 18.2%, unmapped 48.1%)
# [capture_hits] timing: 28.2800, 1590.0000, 9258.4395 (64.8600 ss access; 26051nt and 3.5 threads per alignment)

15.3 min - 6.1min (exp 2) = 0.14s real time/read

5. 
with -rlim 5000 (was -rlim 1000)
with qpalma filter stats (was off)

with all optimizations (index depth 13; node148):
[filter] reason 0:      32 spliced      148 unspliced   17.78%
[filter] reason 1:      2 spliced       600 unspliced   0.33%
[filter] reason 2:      96 spliced      211 unspliced   31.27%
[filter] reason 3:      22 spliced      78 unspliced    22.00%


Decision:
* use all optimizations
* use index of depth 12
* use -qpalma-use-map-max-len 10000
* align all reads (not just the unmapped from unspliced alignments)
* align about 10% of all experiments in a first step, then do later as many as possible

