GCC = gcc-12

.PHONY: test clean all

all: a.so b.so

a.so: a.c cgeneric.h
	$(GCC) -Wall -fpic -O -c -o a.o a.c -lm
	$(GCC) -shared -o a.so a.o -lm

b.so: b.c cgeneric.h smat-operations.c smat-operations.h spde-precision.c spde-precision.h
	$(GCC) -Wall -fpic -O -c -o b.o b.c -lm
	$(GCC) -shared -o b.so b.o -lm

test: test1.c test2.c smat-operations.o spde-precision.o test-funcs.o cgeneric.h
	$(GCC) -o test1.o test1.c smat-operations.o spde-precision.o test-funcs.o -lm
	$(GCC) -o test2.o test2.c smat-operations.o spde-precision.o test-funcs.o -lm

%.o: %.c %.h
	$(GCC) -Wall -fpic -O -c -o $@ $< -lm

clean:
	rm *.o
	rm *.so
