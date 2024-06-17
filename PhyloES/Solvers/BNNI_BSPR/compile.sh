

# cd src
# rm fast_me_
echo "$PWD"

g++  -c -g test_main.cc fastme.cc SPR.cc bme.cc bNNI.cc graph.cc heap.cc traverse.cc initialiser.cc utils.cc -fopenmp
g++  -o  fast_me   bme.o bNNI.o graph.o heap.o SPR.o traverse.o utils.o initialiser.o fastme.o test_main.o -lm -fopenmp -lpthread

rm *o


# echo "$PWD"



# cd ..
# ./src/fast_me_ -i mat.mat -m b -n -s -f 17
# ./fastme -i mat.mat -m b -n -s -u init_topology.nwk -f 17

