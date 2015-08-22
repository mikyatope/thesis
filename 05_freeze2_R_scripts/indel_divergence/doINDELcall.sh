
time find chr4/ -type f | parallel --eta --gnu -j1 "python callFixedINDELs.py {} >> ALL.fixed_indels.csv"

time find chr2L/ -type f | parallel --eta --gnu -j1 "python callFixedINDELs.py {} >> ALL.fixed_indels.csv"

time find chr2R/ -type f | parallel --eta --gnu -j1 "python callFixedINDELs.py {} >> ALL.fixed_indels.csv"

time find chr3L/ -type f | parallel --eta --gnu -j1 "python callFixedINDELs.py {} >> ALL.fixed_indels.csv"

time find chr3R/ -type f | parallel --eta --gnu -j1 "python callFixedINDELs.py {} >> ALL.fixed_indels.csv"

time find chrX/ -type f | parallel --eta --gnu -j1 "python callFixedINDELs.py {} >> ALL.fixed_indels.csv"
