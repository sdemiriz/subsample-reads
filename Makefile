.PHONY: run-example run-example-prg example example-map example-prg example-map-prg activate install

PYTHON := python
PIP := $(PYTHON) -m pip
SR := python -m subsample_reads
VENV_ACTIVATE := source venv/bin/activate
EXAMPLE_DIR := examples

install: requirements.txt
	@echo "Creating virtual environment and installing dependencies..."
	python3 -m venv venv
	$(VENV_ACTIVATE) && $(PIP) install -r requirements.txt

example: install
	@echo "Generating example BAM file..."
	$(VENV_ACTIVATE) && $(PYTHON) $(EXAMPLE_DIR)/make_example.py

example-map: install
	@echo "Generating example BAM file used for mapping..."
	$(VENV_ACTIVATE) && $(PYTHON) $(EXAMPLE_DIR)/make_example_map.py

example-prg: actiinstallvate
	@echo "Generating example BAM file used for HLA-LA PRG..."
	$(VENV_ACTIVATE) && $(PYTHON) $(EXAMPLE_DIR)/make_example_prg.py

example-map-prg: install
	@echo "Generating example BAM file used for mapping with HLA-LA PRG..."
	$(VENV_ACTIVATE) && $(PYTHON) $(EXAMPLE_DIR)/make_example_map_prg.py


run-example: IN_BAM = $(EXAMPLE_DIR)/example.bam
run-example: MAP_BAM = $(EXAMPLE_DIR)/example-map.bam
run-example: SUB_BAM = $(EXAMPLE_DIR)/example-out.bam
run-example: BED = $(EXAMPLE_DIR)/example-map.bed
run-example: PLT = $(EXAMPLE_DIR)/example-out.png
run-example: example example-map
	
	@echo "[STEP 1] Mapping example..."
	$(VENV_ACTIVATE) && $(SR) map --in-bam $(MAP_BAM) --contig chr1 --start 100 --end 1100 --interval-count 10 --bed-dir $(EXAMPLE_DIR)

	@echo "[STEP 2] Sampling example..."
	$(VENV_ACTIVATE) && $(SR) sample --in-bam $(IN_BAM) --bed $(BED) --out-bam $(SUB_BAM)

	@echo "[STEP 3] Plotting example..."
	$(VENV_ACTIVATE) && $(SR) plot --in-bam $(IN_BAM) --map-bam $(MAP_BAM) --sub-bam $(SUB_BAM) --bed $(BED) --out-plt $(PLT)

run-example-prg: IN_BAM = $(EXAMPLE_DIR)/example-prg.bam
run-example-prg: MAP_BAM = $(EXAMPLE_DIR)/example-prg-map.bam
run-example-prg: SUB_BAM = $(EXAMPLE_DIR)/example-prg-out.bam
run-example-prg: BED = $(EXAMPLE_DIR)/example-prg-map.bed
run-example-prg: PLT = $(EXAMPLE_DIR)/example-prg-out.png
run-example-prg: example-prg example-map-prg

	@echo "[STEP 1] Mapping PRG example..."
	$(VENV_ACTIVATE) && $(SR) map --in-bam $(MAP_BAM) --contig chr6 --start 29941260 --end 29942260 --interval-count 10 --bed-dir $(EXAMPLE_DIR)

	@echo "[STEP 2] Sampling PRG example..."
	$(VENV_ACTIVATE) && $(SR) sample --prg GRCh38 --in-bam $(IN_BAM) --bed $(BED) --out-bam $(SUB_BAM)

	@echo "[STEP 3] Plotting PRG example..."
	$(VENV_ACTIVATE) && $(SR) plot --in-bam $(IN_BAM) --map-bam $(MAP_BAM) --sub-bam $(SUB_BAM) --bed $(BED) --out-plt $(PLT)

algo-demo: IN_BAM = $(EXAMPLE_DIR)/algo-demo.bam
algo-demo: BED = $(EXAMPLE_DIR)/algo-demo.bed
algo-demo: SUB_BAM = $(EXAMPLE_DIR)/algo-demo-out.bam
algo-demo:
	@echo "Running algorithm demo..."
	$(VENV_ACTIVATE) && $(PYTHON) $(EXAMPLE_DIR)/algo_demo.py
	$(VENV_ACTIVATE) && $(SR) sample --in-bam $(IN_BAM) --bed $(BED) --out-bam $(SUB_BAM)
	
	@echo "The following step will fail without samtools being installed"
	@echo "Reads shown in the Application Notes figure (r2, r4, r5, r7, r10, r15, r14) should be present: "
	samtools view $(SUB_BAM) | cut -f 1 | tr '\n' ' '

all: run-example run-example-prg algo-demo

clean:
	@echo "Cleaning up created example files..."
	rm -rf examples/*.bam
	rm -rf examples/*.bam.bai
	rm -rf examples/*.bed
	rm -rf examples/*.png
