CC=g++
GCCOPTS = -Wall
OPTS = $(GCCOPTS)

all: opt dbg

opt: main.cpp solver.o
	$(CC) main.cpp solver.o -O3 -o df.exe

dbg: main.cpp solver_dbg.o
	$(CC) main.cpp solver_dbg.o -g -o df_dbg.exe

solver.o: solver.cpp
	$(CC) -c solver.cpp -O3 $(GCCOPTS) -o solver.o

solver_dbg.o: solver.cpp
	$(CC) -c solver.cpp -g $(GCCOPTS) -o solver_dbg.o

asm:
	$(CC) main.cpp -S -fverbose-asm -o main.s
	$(CC) main.cpp -S -fverbose-asm -O3 -o main_opt.s