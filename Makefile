CC=qcc
CFLAGS=-O2 -std=c99 -Wall -D_XOPEN_SOURCE=700
MPIGEN=-source -D_MPI=1
MPICC=mpicc
MPIFLAGS=-O2 -std=c99 -Wall -D_MPI=1
# 2024-04-03: OSMesa not needed? LINKFLAGS=-lfb_osmesa -lOSMesa -lm
LINKFLAGS=-lm

.PHONY: source
SHELL:=/bin/bash
LMOD=source scripts/load_modules.sh

SRC_DIR:=src
BIN_DIR:=bin
OBJ_DIR:=obj
OUT_DIR:=out
SLURM_DIR:=slurm-out
INC_DIR:=$(SRC_DIR)/include

toml: $(INC_DIR)/toml.h $(INC_DIR)/toml.c
	gcc $(CFLAGS) -I$(INC_DIR) -c $(INC_DIR)/$@.c -o $(OBJ_DIR)/$@.o

drop: $(SRC_DIR)/drop.c $(OBJ_DIR)/toml.o
	$(LMOD); wait && \
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $(OBJ_DIR)/$@.o
	$(CC) $(CFLAGS) -I$(INC_DIR) -autolink $(OBJ_DIR)/toml.o $(OBJ_DIR)/$@.o \
		-o $(BIN_DIR)/$@ $(LINKFLAGS)
	@echo "Complete"
	@echo "Run: ./bin/$@"

.PHONY: run
run:
	source scripts/choose_job.sh

.PHONY: clean
clean:
	rm $(SLURM_DIR)/*
	rm -r $(OUT_DIR)/*

.PHONY: deepclean
deepclean:
	rm $(OBJ_DIR)/*
	rm $(BIN_DIR)/*
	rm $(SLURM_DIR)/*
	rm -r $(OUT_DIR)/*

.PHONY: all
all: | toml drop

