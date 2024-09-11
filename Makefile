CC=g++

all: opt dbg

opt:
	$(CC) main.cpp -O3 -o df.exe

dbg:
	$(CC) main.cpp -g -o df_dbg.exe

asm:
	$(CC) main.cpp -S -fverbose-asm -o main.s
	$(CC) main.cpp -S -fverbose-asm -O3 -o main_opt.s