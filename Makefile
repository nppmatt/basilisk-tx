CC = qcc
CFLAGS = -O2 -std=c99 -Wall -D_XOPEN_SOURCE=700
LINKFLAGS = -lfb_osmesa -lOSMesa -lm

SRC_DIR := src
OUT_DIR := out

fall_test: $(SRC_DIR)/fall_test.c
	$(CC) $(CFLAGS) -autolink $< -o $(OUT_DIR)/$@ $(LINKFLAGS)

all: fall_test

.PHONY: clean
clean:
	rm -r $(OUT_DIR)/*
