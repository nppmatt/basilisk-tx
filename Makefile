CC=qcc
CFLAGS=-O2 -std=c99 -Wall -D_XOPEN_SOURCE=700
MPICC=mpicc
MPIFLAGS=-O2 -std=c99 -Wall -D_MPI=1 -D_XOPEN_SOURCE=700
# 2024-04-03: OSMesa not needed? LINKFLAGS=-lfb_osmesa -lOSMesa -lm
LINKFLAGS=-lm

.PHONY: source
SHELL:=/usr/bin/env bash
LMOD=source scripts/load_modules.sh
MPIGEN=source scripts/mpigen.sh

SRC_DIR:=src
BIN_DIR:=bin
OBJ_DIR:=obj
OUT_DIR:=out
SLURM_DIR:=slurm-out
INC_DIR:=$(SRC_DIR)/include

toml: $(INC_DIR)/toml.h $(INC_DIR)/toml.c
	gcc $(CFLAGS) -I$(INC_DIR) -c $(INC_DIR)/$@.c -o $(OBJ_DIR)/$@.o

drop: $(SRC_DIR)/drop.c $(OBJ_DIR)/toml.o
	rm $(SRC_DIR)/_$@.c
	$(MPIGEN) $@.c
	$(LMOD); wait && \
	$(MPICC) $(MPIFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/_$@.c -o $(OBJ_DIR)/$@.o
	$(MPICC) $(MPIFLAGS) -I$(INC_DIR) $(OBJ_DIR)/toml.o $(OBJ_DIR)/$@.o \
		-o $(BIN_DIR)/$@ $(LINKFLAGS)
	@echo "Complete"

.PHONY: run
run:
	source scripts/choose_job.sh

.PHONY: clean
clean:
	rm -f *.mp4
	rm -f $(SLURM_DIR)/*
	rm -rf $(OUT_DIR)/*

.PHONY: deepclean
deepclean:
	rm -f *.mp4
	rm -f $(SLURM_DIR)/*
	rm -rf $(OUT_DIR)/*
	rm -f $(OBJ_DIR)/*
	rm -f $(BIN_DIR)/*

