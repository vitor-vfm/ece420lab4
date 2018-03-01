N_PROC=4

.PHONY: datatrim serialtester main

all: datatrim serialtester main

datatrim:
	gcc datatrim.c -o datatrim

serialtester:
	gcc serialtester.c Lab4_IO.c -o serialtester -lm

main:
	mpicc -std=c11 -g -Wall main.c Lab4_IO.c -o main -lm 

run: main
	mpirun -np $(N_PROC) ./main

run_generate: main datatrim
	./datatrim
	mpirun -np $(N_PROC) ./main
