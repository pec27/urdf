# -fPIC process independent code (for shared libraries)
# -std=c99 complex numbers

LIB = build/liburdf.so
DEBUG_LIB = build/liburdf_debug.so
CC = gcc -Wall -fPIC -shared -std=c99 -O2
HDR =
SRC = src/urdfs.c

all: build ${LIB} ${DEBUG_LIB}

test: all
	python -m urdf.test

build:
	mkdir -p build

${LIB}: ${SRC} ${HDR} Makefile
	${CC} -DNDEBUG -o ${LIB} ${SRC}

${DEBUG_LIB}: ${SRC} ${HDR} Makefile
	${CC} -o ${DEBUG_LIB} ${SRC}
clean:
	rm -rf build

