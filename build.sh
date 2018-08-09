## Visit https://github.com/simongog/sdsl-lite to get a copy of SDSL
inc_sdsl=<Path to the include folder of SDSL>
lib_sdsl=<Path to the lib folder of SDSL>

## Visit https://github.com/jlabeit/parallel-range-lite to get a copy of
## Parallel Range Lite
inc_sa=<Path to the include folder of Parallel Range Lite>
lib_sa=<Path to the lib folder of Parallel Range Lite>

## Integer sorting of the Problem Based Benchmar Suite (PBBS -
## http://www.cs.cmu.edu/~pbbs/) 
g++ -O2 -c external/sortPBBS/blockRadixSort.cpp -o external/sortPBBS/blockRadixSort.o -fcilkplus

echo "[BWT] Compiling space efficient BWT"
g++ -std=c++14 -O3 main.cpp external/sortPBBS/blockRadixSort.o -o bwt_se \
-I$inc_sdsl -L$lib_sdsl -I$inc_sa -L$lib_sa -lsdsl -llibprange -ldl -fcilkplus \
-lcilkrts 
