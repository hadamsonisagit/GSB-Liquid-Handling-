# opentrons-gsb-oligo-cloning

Opentrons OT-2 automation scripts for high-throughput oligo pool cloning using the [Salis Lab Genetic Systems Builder (GSB)](https://www.denovodna.com/build_genetic_systems) pipeline. Automates Steps 4 (PCR setup + thermocycling) and Step 6 (Golden Gate assembly setup + thermocycling) across a full 96-well plate.

> **Note:** These scripts were developed for protocol documentation and portfolio purposes. Logic has been validated against the Opentrons Python API v2.15 but has not yet been tested on physical hardware. Use with appropriate caution and verify volumes in simulation before running.

---

## Background

The GSB pipeline converts a synthesized oligonucleotide pool into sequence-verified plasmid assemblies in 96-well format. Each well receives a unique pair of PCR primers that amplify a specific DNA fragment from the shared pool. Those fragments are then assembled into plasmids via Golden Gate (Type IIS restriction + ligation). The result is a plate of unique constructs built in parallel.

These scripts automate the two most liquid-handling-intensive steps:

- **Step 4** — PCR plate setup and thermocycling
- **Step 6** — Golden Gate assembly plate setup and thermocycling

Steps 1–3 (design, purchasing, oligo pool reconstitution) and Step 5 (Zymo ZR-96 DNA purification) are performed manually between the two scripts.

---

## Hardware Requirements

| Component | Details |
|---|---|
| Robot | Opentrons OT-2 |
| Left pipette | P300 Single-Channel GEN2 |
| Right pipette | P20 Multi-Channel GEN2 |
| Module | Thermocycler Module GEN2 (occupies Slots 7 + 10) |
| Tip racks | 300 µL (×1) and 20 µL (×3) |
| Tube rack | NEST 15-well tube rack for 15 mL conicals |
| Plates | 96-well PCR plates (NEST or equivalent) |

---

## Deck Layout

Both scripts share the same deck configuration:

```
┌───────────┬───────────┬───────────┐
│  THERMO   │  20µL     │  20µL     │
│  CYCLER   │  Tips #2  │  Tips #3  │
│ (7 + 10)  ├───────────┼───────────┤
│           │  20µL     │   empty   │
│           │  Tips #1  │           │
├───────────┼───────────┼───────────┤
│  Slot 4   │  Slot 5   │  300µL    │
│ (varies)  │ (varies)  │  Tips     │
├───────────┼───────────┼───────────┤
│ Tube Rack │   empty   │   empty   │
│ (mastermix│           │           │
│  tube A1) │           │           │
└───────────┴───────────┴───────────┘
        (front of robot)

Step 4: Slot 4 = Forward Primer Plate | Slot 5 = Reverse Primer Plate
Step 6: Slot 4 = PCR Product Eluate Plate | Slot 5 = empty
```

---

## Script Overview

### `gsb_step4_pcr_setup.py` — PCR Setup + Thermocycler

Sets up a 96-well PCR plate and runs the thermocycler program automatically.

**Reaction composition** (default, 2.5 µM primers):

| Reagent | Volume |
|---|---|
| NEBNext Ultra II Q5 Master Mix (2X) | 12.5 µL |
| Oligopool Template (0.5 ng/µL) | 1.0 µL |
| Nuclease-free water | 1.5 µL |
| Forward primer (well-specific) | 5.0 µL |
| Reverse primer (well-specific) | 5.0 µL |
| **Total** | **25.0 µL** |

**Pipetting strategy:**
- P300 single-channel dispenses mastermix from a single 15 mL conical tube into each well individually. Single-tube preparation ensures every well receives identically composed mastermix with minimal inter-well variability.
- P20 multi-channel adds forward then reverse primers column by column with new tips every column.
- Liquid level tracking adjusts aspiration height in the conical tube as volume depletes across 96 dispenses.

**Thermocycler program:**

| Step | Temperature | Time |
|---|---|---|
| Initial denaturation | 98°C | 2 min |
| Denaturation | 98°C | 20 sec |
| Annealing | 61°C | 30 sec |
| Extension | 72°C | 6 sec |
| Cycles | — | ×26 |
| Final extension | 72°C | 60 sec |
| Hold | 4°C | ∞ |

---

### `gsb_step6_assembly_setup.py` — Golden Gate Assembly + Thermocycler

Sets up a 96-well Golden Gate assembly plate and runs the digestion/ligation thermocycler program automatically.

**Reaction composition** (default settings):

| Reagent | Volume |
|---|---|
| T4 DNA Ligase Buffer (10X) | 2.5 µL |
| T4 DNA Ligase (2000 U/µL) | 0.5 µL |
| BsaI-HFv2 (20 U/µL) | 1.5 µL |
| Vector fragment (calculated) | variable |
| Nuclease-free water | variable |
| PCR product (Step 5 eluate) | 10.0 µL |
| **Total** | **25.0 µL** |

Vector volume is calculated to achieve equimolar ratio with PCR fragments:

```
vector_vol = (pcr_conc × pcr_vol × vector_bp) / (vector_conc × fragment_bp × n_fragments)
```

Capped at a maximum of 75 ng total vector per well.

**Pipetting strategy:**
- P300 single-channel dispenses assembly mastermix from a single 15 mL conical tube well-by-well. Liquid level tracking active throughout.
- P20 multi-channel transfers well-specific PCR products column by column with new tips every column to prevent cross-contamination between unique assemblies.

**Thermocycler program:**

| Step | Temperature | Time |
|---|---|---|
| Digestion | 37°C | 5 min |
| Ligation | 16°C | 5 min |
| Cycles | — | ×60 |
| Hold | 4°C | ∞ |

---

## Configuration

All user-adjustable parameters are grouped at the top of each script under `USER CONFIGURATION`. Key parameters:

**Step 4:**

```python
NUM_COLUMNS          = 12     # 1–12, set lower for partial plate runs
PRIMER_CONC_UM       = 2.5    # Primer concentration in µM
PRIMER_START_COLUMN  = 1      # Primer plate column that maps to destination column 1
INSPECT_PAUSE        = False  # Set True to pause after mastermix for visual QC
DRY_RUN              = False  # Set True to rehearse without moving liquid
```

**Step 6:**

```python
NUM_COLUMNS          = 12
PCR_PRODUCT_VOL      = 10.0   # µL of eluate per well (8–16 µL range)
PCR_CONC_NG_UL       = 10.0   # Measured by Nanodrop or plate reader
NUM_FRAGMENTS        = 5      # Average fragments per assembly
VECTOR_CONC_NG_UL    = 30.0   # Vector stock concentration
VECTOR_LENGTH_BP     = 1500   # Vector length in base pairs
ELUATE_START_COLUMN  = 1      # Eluate plate column that maps to assembly column 1
ENZYME               = 'BsaI-HFv2'
ENZYME_VOL           = 1.5    # Use 3.0 µL for BsmBI, Esp3I, or SapI
DRY_RUN              = False
```

---

## Reagents

### Step 4

| Reagent | Catalog # | Notes |
|---|---|---|
| NEBNext Ultra II Q5 Master Mix (2X) | NEB M0544X | |
| Nuclease-free water | NEB B1500 | |
| IDT custom primers in TE buffer | — | 96-well plate format, 2.5 µM |
| Synthesized oligopool | Twist / IDT / Genscript | Reconstituted to 0.5 ng/µL |

### Step 6

| Reagent | Catalog # | Notes |
|---|---|---|
| T4 DNA Ligase Buffer (10X) | NEB B0202S | |
| T4 DNA Ligase (2000 U/µL) | NEB M0202M | |
| BsaI-HFv2 (20 U/µL) | NEB R3733 | Default enzyme |
| BbsI-HF (20 U/µL) | NEB R3539L | Alternative |
| BsmBI-v2 (10 U/µL) | NEB R0739L | Use 3.0 µL |
| Esp3I (10 U/µL) | NEB R0734L | Use 3.0 µL |
| SapI (10 U/µL) | NEB R0569L | Use 3.0 µL |
| Nuclease-free water | NEB B1500 | |

---

## Usage

### 1. Prepare mastermix

**Step 4 mastermix** — combine in a single 15 mL conical tube, vortex thoroughly just before loading:

The script prints exact volumes to the run log at startup based on your configuration. Example for a full 96-well plate at default settings:

```
NEBNext Q5 Master Mix (2X) : 1330.6 µL
Oligopool Template         :  106.5 µL
Nuclease-free water        :  159.7 µL
TOTAL                      : 1596.8 µL
```

**Step 6 mastermix** — combine on ice, mix gently by flicking (do not vortex enzymes):

```
T4 DNA Ligase Buffer (10X) : 264.0 µL
T4 DNA Ligase (2000 U/µL)  :  52.8 µL
BsaI-HFv2                  : 158.4 µL
Vector Fragment             : variable
Nuclease-free water        : variable
```

### 2. Load the deck

Follow the pre-run checklist printed by the script at startup. The robot will pause and wait for confirmation before any liquid movement begins.

### 3. Run

Upload the script to the Opentrons app, calibrate labware, and press Run. The robot will:
- Print a full run summary to the log
- Pause for pre-run checklist confirmation
- Execute all liquid transfers
- Run the thermocycler program automatically
- Open the lid and pause when complete

### 4. Between steps

After Step 4, proceed to **Step 5 (Zymo ZR-96 DNA purification)** manually before running the Step 6 script. The purification step is centrifuge-based and not automated in this pipeline.

---

## Labware Notes

Scripts use NEST 96-well PCR plates (`nest_96_wellplate_100ul_pcr_full_skirt`) as the default labware definition. If using a different plate brand, substitute the appropriate Opentrons labware definition at the lines marked:

```python
# Swap for your labware def if needed
```

Custom labware definitions can be created at [Opentrons Labware Creator](https://labware.opentrons.com/create/) and loaded alongside the script in the Opentrons app.

---

## Features

- **Two-pipette strategy** — P300 single-channel for mastermix (minimal inter-well variability), P20 multi-channel for primers and PCR product transfer (speed)
- **Liquid level tracking** — aspiration height in the conical tube adjusts dynamically as volume depletes, preventing air aspiration across 96 dispenses
- **Automatic tip management** — shared tips for common reagents, new tips per column for well-specific transfers
- **Integrated thermocycler** — full PCR and Golden Gate cycling programs run automatically without user intervention
- **Column offset support** — `PRIMER_START_COLUMN` and `ELUATE_START_COLUMN` parameters allow non-standard plate layouts
- **Dry run mode** — rehearse all robot movements without transferring liquid
- **Pre-run checklists** — robot pauses at start for user confirmation
- **Reagent calculator** — exact load volumes printed to run log with 10% overage
- **Optional inspect pause** — pause after mastermix dispense for visual QC before primers are added (Step 4)

---

## Acknowledgments

Protocol based on the Salis Lab Genetic Systems Builder pipeline:

> Salis, H.A., Adamson, H., & Vakoun, A. *Oligopool-to-Plasmids Protocol v1.1.* Salis Lab, Penn State University. Available at: https://www.denovodna.com/build_genetic_systems

Scripts developed independently for automation and portfolio purposes. Not affiliated with the Salis Lab or De Novo DNA.

---

## License

MIT License — free to use, modify, and distribute with attribution.

---

## Author

[Your Name]
[Your institution / LinkedIn URL]
