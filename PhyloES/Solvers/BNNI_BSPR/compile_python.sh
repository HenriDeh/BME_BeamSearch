


g++ -c -Ofast -fopenmp -fPIC bridge.cc fastme.cc SPR.cc bme.cc bNNI.cc graph.cc heap.cc traverse.cc initialiser.cc utils.cc
g++ -shared -Wl,-soname,bridge.so -o bridge.so bridge.o bme.o bNNI.o graph.o heap.o SPR.o traverse.o utils.o initialiser.o fastme.o -fopenmp -lm -lpthread



# g++ -c -fopenmp -fPIC bridge.cc 
# g++ -shared -Wl,-soname,bridge.so -o bridge.so. bridge.o 

rm *.o

