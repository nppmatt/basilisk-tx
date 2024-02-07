CC=qcc
CFLAGS=-O2 -std=c99 -Wall -D_XOPEN_SOURCE=700
LINKFLAGS=-lfb_osmesa -lOSMesa -lm

.PHONY: source
SHELL:=/bin/bash
LMOD=source load_modules.sh


SRC_DIR:=src
OUT_DIR:=out

fall_test: $(SRC_DIR)/fall_test.c
	$(LMOD); wait && \
	$(CC) $(CFLAGS) -autolink $< -o $(OUT_DIR)/$@ $(LINKFLAGS)
	@echo "Finished!"

fall_test2: $(SRC_DIR)/fall_test2.c
	$(LMOD); wait && \
	$(CC) $(CFLAGS) -autolink $< -o $(OUT_DIR)/$@ $(LINKFLAGS)
	@echo "Finished!"

partial_wetting: $(SRC_DIR)/partial_wetting.c
	$(LMOD); wait && \
	$(CC) $(CFLAGS) -autolink $< -o $(OUT_DIR)/$@ $(LINKFLAGS)
	@echo "Finished!"

partial_wetting_fast: $(SRC_DIR)/partial_wetting_fast.c
	$(LMOD); wait && \
	$(CC) $(CFLAGS) -autolink $< -o $(OUT_DIR)/$@ $(LINKFLAGS)
	@echo "Finished!"

pw_fine: $(SRC_DIR)/pw_fine.c
	$(LMOD); wait && \
	$(CC) $(CFLAGS) -autolink $< -o $(OUT_DIR)/$@ $(LINKFLAGS)
	@echo "Finished!"

pw_fast_fine: $(SRC_DIR)/pw_fast_fine.c
	$(LMOD); wait && \
	$(CC) $(CFLAGS) -autolink $< -o $(OUT_DIR)/$@ $(LINKFLAGS)
	@echo "Finished!"

all: fall_test fall_test2 partial_wetting partial_wetting_fast pw_fine pw_fast_fine

.PHONY: clean
clean:
	rm $(OUT_DIR)/*
