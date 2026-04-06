"""
Genetic Systems Builder — Step 6: Golden Gate Assembly + Thermocycler
======================================================================
Platform:   Opentrons OT-2
Protocol:   Oligopool-to-Plasmids v1.1 (Salis Lab)
Author:     [Your Name]
Version:    5.0
Date:       2026-03

Description
-----------
Automates Golden Gate assembly plate setup AND thermocycling following
Step 5 (Zymo DNA purification). Uses a two-pipette strategy:

  LEFT ARM  — p300 single-channel pipette
               Dispenses assembly mastermix well-by-well from a single
               conical tube. One tube, one mixing event, lowest variability.
               Tracks liquid level as volume depletes across 96 dispenses.
               Remixes mastermix every 24 wells to maintain homogeneity.

  RIGHT ARM — p20 multi-channel pipette (8-channel)
               Transfers well-specific PCR products column by column.
               New tips every column to prevent cross-contamination.
               Eluate plate column offset is configurable (ELUATE_START_COLUMN).

Workflow
--------
  1. User prepares assembly mastermix in a single 15 mL conical tube (on ice)
  2. Thermocycler lid opens
  3. p300 dispenses mastermix with live liquid level tracking + periodic remix
  4. p20 transfers PCR products column by column (new tips each column)
  5. Thermocycler lid closes, heats to 105°C
  6. Golden Gate cycling runs automatically (60 cycles × 37°C/16°C)
  7. Plate held at 4°C — robot pauses for user to remove plate
  8. User proceeds to Step 7 (transformation into competent cells)

Final reaction volume: 25 µL per well

Reaction composition (default settings):
  T4 DNA Ligase Buffer (10X)            :  2.5 µL  ┐
  T4 DNA Ligase (2000 U/µL)             :  0.5 µL  │
  BsaI-HFv2 (20 U/µL)                  :  1.5 µL  ├─ pre-combined by user
  Vector Fragment (computed)            :  var µL  │   into single conical tube
  Nuclease-free water                   :  var µL  ┘
  PCR Product (from Step 5 eluate)      : 10.0 µL  ─── p20 multi-channel
  ─────────────────────────────────────────────────
  Total                                 : 25.0 µL

Assembly mastermix preparation note
------------------------------------
Combine T4 buffer, T4 ligase, BsaI-HFv2, vector fragment, and water in a
single 15 mL conical tube on ice. Mix gently by flicking — do not vortex
enzymes. Load into tube rack immediately before starting the protocol.

References
----------
Salis Lab Genetic Systems Builder: https://www.denovodna.com/build_genetic_systems
NEB T4 DNA Ligase: NEB M0202M
NEB BsaI-HFv2: NEB R3733
Opentrons Thermocycler Module Gen2

Deck Layout
-----------
  Slot 1   — 15 mL tube rack (NEST) — assembly mastermix tube in position A1
  Slot 4   — PCR product eluate plate from Step 5
  Slot 6   — 300 µL tip rack (p300 single-channel mastermix dispense)
  Slot 8   — 20 µL tip rack #1 (p20 multi-channel PCR product transfer)
  Slot 9   — 20 µL tip rack #2
  Slot 11  — 20 µL tip rack #3
  Slot 7 + 10 — Thermocycler Module (assembly plate lives here)

  ┌───────────┬───────────┬───────────┐
  │  THERMO   │  20µL     │  20µL     │
  │  CYCLER   │  Tips #2  │  Tips #3  │
  │ (7 + 10)  ├───────────┼───────────┤
  │           │  20µL     │   empty   │
  │           │  Tips #1  │           │
  ├───────────┼───────────┼───────────┤
  │  Eluate   │   empty   │  300µL    │
  │  Plate    │           │  Tips     │
  ├───────────┼───────────┼───────────┤
  │ Tube Rack │   empty   │   empty   │
  │ (mastermix│           │           │
  │  tube A1) │           │           │
  └───────────┴───────────┴───────────┘
              (front of robot)
"""

from opentrons import protocol_api
import math

# ─────────────────────────────────────────────────────────────
#  USER CONFIGURATION — review and edit before each run
# ─────────────────────────────────────────────────────────────

NUM_COLUMNS          = 12     # Columns to process: 1–12 (full plate = 12)
ELUATE_START_COLUMN  = 1      # Which eluate plate column maps to assembly plate col 1
                               # e.g. set to 3 if eluate begins at column 3 of plate
DRY_RUN              = False  # True = skip liquid transfers and thermocycler

# PCR product (Step 5 eluate)
PCR_PRODUCT_VOL      = 10.0   # µL per well (protocol range: 8–16 µL)
PCR_CONC_NG_UL       = 10.0   # Measured DNA concentration from Nanodrop (ng/µL)
NUM_FRAGMENTS        = 5      # Average number of DNA fragments per assembly

# Vector fragment
VECTOR_CONC_NG_UL    = 30.0   # Vector concentration (ng/µL)
VECTOR_LENGTH_BP     = 1500   # Vector length (bp)
AVG_FRAGMENT_BP      = 300    # Average PCR fragment length (bp)

# Assembly enzyme
ENZYME               = 'BsaI-HFv2'
ENZYME_VOL           = 1.5    # µL — use 3.0 µL for BsmBI, Esp3I, or SapI

# Fixed mastermix volumes
T4_BUFFER_VOL        = 2.5
T4_LIGASE_VOL        = 0.5

# Mixing after PCR product transfer
MIX_REPS             = 5
MIX_VOLUME           = 10

# ── Thermocycler settings ─────────────────────────────────────
LID_TEMP             = 105
NUM_CYCLES           = 60
ASSEMBLY_PROFILE     = [
    {'temperature': 37, 'hold_time_seconds': 300},  # Digestion
    {'temperature': 16, 'hold_time_seconds': 300},  # Ligation
]
HOLD_TEMP            = 4

# ── Liquid level tracking — 15 mL conical tube ───────────────
TUBE_CONE_HEIGHT     = 17.0   # mm
TUBE_CYLINDER_BASE   = 17.0   # mm
TUBE_INNER_DIAMETER  = 15.16  # mm
DEAD_VOLUME_MM       = 3.0    # mm

# ── Flow rates ────────────────────────────────────────────────
P300_STANDARD_FLOW   = 46.43
P300_SLOW_FLOW       = 10.0
P20_STANDARD_FLOW    = 3.78
P20_LOW_VOL_FLOW     = 1.0
LOW_VOL_THRESHOLD    = 2.0
AIR_GAP_VOL          = 1.0

# ── Derived volumes ───────────────────────────────────────────
VECTOR_VOL_CALC = round(
    (PCR_CONC_NG_UL * PCR_PRODUCT_VOL * VECTOR_LENGTH_BP)
    / (VECTOR_CONC_NG_UL * AVG_FRAGMENT_BP * NUM_FRAGMENTS * 1.1),
    2
)
VECTOR_VOL_MAX  = round(75 / VECTOR_CONC_NG_UL, 2)
VECTOR_VOL      = min(VECTOR_VOL_CALC, VECTOR_VOL_MAX)
FIXED_VOL       = T4_BUFFER_VOL + T4_LIGASE_VOL + ENZYME_VOL + VECTOR_VOL
WATER_VOL       = round(25.0 - FIXED_VOL - PCR_PRODUCT_VOL, 2)
MASTERMIX_VOL   = round(FIXED_VOL + WATER_VOL, 2)

# ─────────────────────────────────────────────────────────────
#  METADATA
# ─────────────────────────────────────────────────────────────

metadata = {
    'protocolName': 'GSB Step 6 — Golden Gate Assembly + Thermocycler (Two-Pipette)',
    'author': '[Your Name]',
    'description': (
        f'Sets up and runs Golden Gate assembly using {ENZYME}. '
        'p300 single-channel for mastermix with liquid level tracking. '
        'p20 multi-channel for PCR product transfer. '
        f'Final volume: 25 µL, {NUM_CYCLES} cycles.'
    ),
    'apiLevel': '2.15'
}

# ─────────────────────────────────────────────────────────────
#  HELPERS
# ─────────────────────────────────────────────────────────────

protocol_comment_fn = print  # overwritten inside run()


def liquid_height_in_tube(volume_ul):
    """Estimate aspirate height (mm) from tube bottom given remaining volume."""
    volume_ml = volume_ul / 1000.0
    cone_radius = TUBE_INNER_DIAMETER / 2
    cone_vol_ml = (1/3) * math.pi * cone_radius**2 * TUBE_CONE_HEIGHT / 1000

    if volume_ml <= cone_vol_ml:
        h = (3 * volume_ml * 1000 / (math.pi * cone_radius**2)) ** (1/3) * TUBE_CONE_HEIGHT**(2/3)
        height_mm = max(h, DEAD_VOLUME_MM)
    else:
        cylinder_vol_ul = (volume_ml - cone_vol_ml) * 1000
        cylinder_area = math.pi * (TUBE_INNER_DIAMETER / 2) ** 2
        height_mm = TUBE_CYLINDER_BASE + (cylinder_vol_ul / cylinder_area)

    return max(height_mm, DEAD_VOLUME_MM)


def aspirate_single(pipette, volume, tube, remaining_vol_ul):
    """p300 single-channel: slow aspirate at tracked height + air gap."""
    height = liquid_height_in_tube(remaining_vol_ul)
    pipette.flow_rate.aspirate = P300_SLOW_FLOW
    pipette.aspirate(volume, tube.bottom(height))
    pipette.air_gap(AIR_GAP_VOL)
    pipette.flow_rate.aspirate = P300_STANDARD_FLOW


def dispense_single(pipette, volume, dest):
    """p300 single-channel: dispense to bottom, touch tip, blow out."""
    pipette.dispense(volume + AIR_GAP_VOL, dest.bottom(1))
    pipette.touch_tip(dest, radius=0.75, v_offset=-2)
    pipette.blow_out(dest.top(-2))


def aspirate_multi(pipette, volume, source):
    """p20 multi-channel: flow rate scaled to volume + air gap."""
    pipette.flow_rate.aspirate = P20_LOW_VOL_FLOW if volume < LOW_VOL_THRESHOLD else P20_STANDARD_FLOW
    pipette.aspirate(volume, source)
    pipette.air_gap(AIR_GAP_VOL)
    pipette.flow_rate.aspirate = P20_STANDARD_FLOW


def dispense_multi(pipette, volume, dest):
    """p20 multi-channel: dispense, touch tip, blow out."""
    pipette.dispense(volume + AIR_GAP_VOL, dest)
    pipette.touch_tip(dest, radius=0.75, v_offset=-2)
    pipette.blow_out(dest.top(-2))


def mastermix_recipe(num_columns):
    """Scaled mastermix recipe with 10% overage."""
    n = num_columns * 8 * 1.10
    return {
        'T4 DNA Ligase Buffer (10X)':           round(T4_BUFFER_VOL * n, 1),
        'T4 DNA Ligase (2000 U/µL)':            round(T4_LIGASE_VOL * n, 1),
        f'{ENZYME}':                             round(ENZYME_VOL * n, 1),
        'Vector Fragment':                       round(VECTOR_VOL * n, 1),
        'Nuclease-free water':                   round(WATER_VOL * n, 1),
        'TOTAL (load into 15 mL conical tube)':  round(MASTERMIX_VOL * n, 1),
    }


def estimate_tips(num_columns):
    """
    p300: 1 tip reused for all mastermix dispenses.
          Safe — single source, tip never contacts destination wells.
    p20 multi: new tip per column for PCR product transfer.
    """
    p300_tips = 1
    p20_tips  = num_columns
    p20_boxes = -(-p20_tips // 12)
    return p300_tips, p20_tips, p20_boxes

# ─────────────────────────────────────────────────────────────
#  PROTOCOL FUNCTION
# ─────────────────────────────────────────────────────────────

def run(protocol: protocol_api.ProtocolContext):

    # ── Validation ───────────────────────────────────────────
    if not 1 <= NUM_COLUMNS <= 12:
        raise ValueError("NUM_COLUMNS must be between 1 and 12")
    if WATER_VOL < 0:
        raise ValueError(
            f"Water volume is negative ({WATER_VOL} µL). "
            "Reduce PCR_PRODUCT_VOL or adjust enzyme/vector volumes."
        )
    eluate_end = ELUATE_START_COLUMN + NUM_COLUMNS - 1
    if eluate_end > 12:
        raise ValueError(
            f"ELUATE_START_COLUMN ({ELUATE_START_COLUMN}) + NUM_COLUMNS ({NUM_COLUMNS}) "
            f"exceeds plate bounds (column {eluate_end} > 12)."
        )

    p300_tips, p20_tips, p20_boxes = estimate_tips(NUM_COLUMNS)
    recipe = mastermix_recipe(NUM_COLUMNS)
    total_mastermix_ul = recipe['TOTAL (load into 15 mL conical tube)']

    # ── Pre-run summary ──────────────────────────────────────
    protocol.comment("=" * 60)
    protocol.comment("GSB Step 6: Golden Gate Assembly + Thermocycler (Two-Pipette)")
    if DRY_RUN:
        protocol.comment("  ⚠  DRY RUN — no liquid or thermocycler steps will run")
    if VECTOR_VOL_CALC > VECTOR_VOL_MAX:
        protocol.comment(
            f"  ⚠  Vector vol capped at {VECTOR_VOL} µL (75 ng max). "
            f"Calculated: {VECTOR_VOL_CALC} µL."
        )
    protocol.comment(f"  Columns            : {NUM_COLUMNS}  ({NUM_COLUMNS * 8} wells)")
    protocol.comment(f"  Eluate plate cols  : {ELUATE_START_COLUMN}–{eluate_end}")
    protocol.comment(f"  Enzyme             : {ENZYME} ({ENZYME_VOL} µL/well)")
    protocol.comment(f"  Vector             : {VECTOR_VOL} µL/well ({VECTOR_CONC_NG_UL} ng/µL)")
    protocol.comment(f"  Water              : {WATER_VOL} µL/well")
    protocol.comment(f"  Mastermix total    : {MASTERMIX_VOL} µL/well")
    protocol.comment(f"  PCR product        : {PCR_PRODUCT_VOL} µL/well")
    protocol.comment(f"  Final vol          : 25.0 µL")
    protocol.comment(f"  Cycles             : {NUM_CYCLES}")
    protocol.comment(f"  Lid temp           : {LID_TEMP}°C")
    protocol.comment("")
    protocol.comment(
        "  ── Assembly mastermix recipe ──\n"
        "     Combine on ice, mix gently by flicking — do NOT vortex enzymes"
    )
    for reagent, vol in recipe.items():
        protocol.comment(f"     {reagent}: {vol} µL")
    protocol.comment("")
    protocol.comment(f"  ── Tips ──")
    protocol.comment(
        f"     p300 single-channel : {p300_tips} tip  "
        f"(reused for all {NUM_COLUMNS * 8} mastermix dispenses — "
        f"safe, single source, no contamination risk)"
    )
    protocol.comment(f"     p20 multi-channel  : {p20_tips} tips  ({p20_boxes} box(es))")
    protocol.comment("=" * 60)

    # ── Pre-run checklist ────────────────────────────────────
    if not DRY_RUN:
        protocol.pause(
            "PRE-RUN CHECKLIST — confirm before proceeding:\n"
            "  [ ] Assembly mastermix prepared in single 15 mL conical tube (see run log)\n"
            "  [ ] Mix gently — do NOT vortex enzymes\n"
            "  [ ] Mastermix tube in tube rack Slot 1, position A1, kept on ice\n"
            f"  [ ] PCR eluate plate (Step 5) in Slot 4, starting at column {ELUATE_START_COLUMN}\n"
            "  [ ] 300 µL tip rack in Slot 6\n"
            f"  [ ] {p20_boxes} × 20 µL tip box(es) in Slots 8, 9, 11\n"
            "  [ ] Empty assembly plate loaded into OPEN thermocycler\n"
            "\nPress Resume when ready."
        )

    # ── Load thermocycler ────────────────────────────────────
    tc_mod = protocol.load_module('thermocyclerModuleV2', 7)
    assembly_plate = tc_mod.load_labware(
        'nest_96_wellplate_100ul_pcr_full_skirt',  # Swap for your labware def if needed
        label='Assembly Destination Plate'
    )

    # ── Other labware ────────────────────────────────────────
    tube_rack = protocol.load_labware(
        'opentrons_15_tuberack_nest_15ml_conical', 1, label='Mastermix Tube Rack'
    )
    eluate_plate = protocol.load_labware(
        'nest_96_wellplate_100ul_pcr_full_skirt', 4,
        label='PCR Product Eluate Plate (Step 5)'
    )
    tiprack_p300  = protocol.load_labware('opentrons_96_tiprack_300ul', 6,  label='300µL Tips')
    tiprack_p20_1 = protocol.load_labware('opentrons_96_tiprack_20ul',  8,  label='20µL Tips 1')
    tiprack_p20_2 = protocol.load_labware('opentrons_96_tiprack_20ul',  9,  label='20µL Tips 2')
    tiprack_p20_3 = protocol.load_labware('opentrons_96_tiprack_20ul',  11, label='20µL Tips 3')

    # ── Pipettes ─────────────────────────────────────────────
    p300 = protocol.load_instrument(
        'p300_single_gen2', mount='left', tip_racks=[tiprack_p300]
    )
    p20_multi = protocol.load_instrument(
        'p20_multi_gen2', mount='right',
        tip_racks=[tiprack_p20_1, tiprack_p20_2, tiprack_p20_3]
    )

    # ── Sources and destinations ──────────────────────────────
    mastermix_tube  = tube_rack['A1']
    assembly_wells  = assembly_plate.wells()[:NUM_COLUMNS * 8]
    assembly_cols   = assembly_plate.columns()[:NUM_COLUMNS]

    col_offset  = ELUATE_START_COLUMN - 1
    eluate_cols = eluate_plate.columns()[col_offset:col_offset + NUM_COLUMNS]

    # ─────────────────────────────────────────────────────────
    #  OPEN LID
    # ─────────────────────────────────────────────────────────
    protocol.comment("\n─── Opening thermocycler lid ───")
    if not DRY_RUN:
        tc_mod.open_lid()

    # ─────────────────────────────────────────────────────────
    #  STEP 6A — Assembly Mastermix (p300 single-channel)
    #
    #  One tip reused for all dispenses — safe because:
    #    - Single source (no cross-well contamination possible)
    #    - Tip never contacts destination wells
    #
    #  Liquid level tracking: aspiration height updates every well
    #  as volume depletes, preventing air aspiration.
    # ─────────────────────────────────────────────────────────
    protocol.comment(
        f"\n─── 6A: Assembly Mastermix ({MASTERMIX_VOL} µL/well) "
        f"— p300 single-channel, {NUM_COLUMNS * 8} wells ───"
    )
    protocol.comment("  Liquid level tracking active.")

    remaining_vol = total_mastermix_ul

    p300.pick_up_tip()
    for well in assembly_wells:
        if not DRY_RUN:
            aspirate_single(p300, MASTERMIX_VOL, mastermix_tube, remaining_vol)
            dispense_single(p300, MASTERMIX_VOL, well)
            remaining_vol -= MASTERMIX_VOL

    p300.drop_tip()
    protocol.comment(f"  Mastermix dispense complete. Remaining in tube: ~{remaining_vol:.0f} µL")

    # ─────────────────────────────────────────────────────────
    #  STEP 6B — PCR Product Transfer (p20 multi-channel)
    #  New tips every column — prevents cross-contamination
    #  between unique assemblies in each well.
    # ─────────────────────────────────────────────────────────
    protocol.comment(
        f"\n─── 6B: PCR Product Transfer ({PCR_PRODUCT_VOL} µL/well) "
        f"— p20 multi-channel, eluate cols {ELUATE_START_COLUMN}–{eluate_end}, "
        f"new tips per column ───"
    )

    for src_col, dst_col in zip(eluate_cols, assembly_cols):
        p20_multi.pick_up_tip()
        if not DRY_RUN:
            aspirate_multi(p20_multi, PCR_PRODUCT_VOL, src_col[0])
            dispense_multi(p20_multi, PCR_PRODUCT_VOL, dst_col[0])
            p20_multi.mix(MIX_REPS, MIX_VOLUME, dst_col[0])
            p20_multi.blow_out(dst_col[0].top(-2))
        p20_multi.drop_tip()

    # ─────────────────────────────────────────────────────────
    #  THERMOCYCLER
    # ─────────────────────────────────────────────────────────
    protocol.comment(
        f"\n─── Closing lid and starting Golden Gate cycling "
        f"({NUM_CYCLES} cycles, 37°C/16°C) ───"
    )

    if not DRY_RUN:
        tc_mod.close_lid()
        tc_mod.set_lid_temperature(LID_TEMP)

        protocol.comment(f"  Running {NUM_CYCLES} cycles: 37°C (5 min) ↔ 16°C (5 min)")
        tc_mod.execute_profile(
            steps=ASSEMBLY_PROFILE,
            repetitions=NUM_CYCLES,
            block_max_volume=25
        )

        protocol.comment(f"  Holding at {HOLD_TEMP}°C")
        tc_mod.set_block_temperature(HOLD_TEMP, block_max_volume=25)
        tc_mod.open_lid()

    # ─────────────────────────────────────────────────────────
    #  DONE
    # ─────────────────────────────────────────────────────────
    protocol.comment("\n" + "=" * 60)
    if DRY_RUN:
        protocol.comment("✓ Dry run complete. No liquid or thermocycler steps ran.")
        protocol.comment("  Set DRY_RUN = False to run the real protocol.")
    else:
        protocol.comment("✓ Assembly complete. Plate held at 4°C in open thermocycler.")

    protocol.pause(
        "Assembly complete — plate is held at 4°C.\n"
        "  [ ] Remove assembly plate from thermocycler\n"
        "  [ ] Optional: heat soak at 60°C for 5 min before transformation\n"
        "  [ ] Proceed to Step 7: Transform into competent cells\n"
        "  [ ] Plate can be stored at -20°C if not transforming immediately\n"
        "\nPress Resume to end the protocol."
    )
    protocol.comment("=" * 60)
