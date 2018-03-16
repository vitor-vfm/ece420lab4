N_PROC=4

.PHONY: datatrim serialtester main main_2

all: datatrim serialtester main main_2

datatrim:
	gcc datatrim.c -o datatrim

serialtester:
	gcc serialtester.c Lab4_IO.c -o serialtester -lm

main:
	mpicc -std=c11 -g -Wall main.c Lab4_IO.c -o main -lm

main_2:
	mpicc -std=c11 -g -Wall main_2.c Lab4_IO.c -o main_2 -lm

run: main
	mpirun -np $(N_PROC) ./main

run_2: main_2
	mpirun -np $(N_PROC) ./main_2

run_generate: main datatrim
	./datatrim
	mpirun -np $(N_PROC) ./main
