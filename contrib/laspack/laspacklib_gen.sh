gcc -g -c *.c -lm
gcc -shared -fPIC -g -o liblaspack-dbg.so  *.c -lm
rm *.o
gcc -g -c *.c -lm
gcc -shared -fPIC  -o liblaspack-opt.so  *.c -lm
rm *.o