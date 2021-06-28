
g++ -std=c++11 -shared -fPIC src/LSODA.cpp src/wrapper.cpp -o liblsoda.so -O3 -Wall