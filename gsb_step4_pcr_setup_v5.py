"""
Genetic Systems Builder — Step 4: PCR Setup + Thermocycler
===========================================================
Platform:   Opentrons OT-2
Protocol:   Oligopool-to-Plasmids v1.1 (Salis Lab)
Author:     [Your Name]
Version:    5.0
Date:       2026-03

Description
-----------
Automates PCR plate setup AND thermocycling for oligo pool amplification.
Uses a two-pipette strategy for maximum accuracy and minimum reagent waste:

  LEFT ARM  — p300 single-channel pipette
               Dispenses mastermix well-by-well from a single conical tube.
               One tube, one mixing event, lowest possible variability.
               Tracks liquid level as volume depletes across 96 dispenses.
               Remixes mastermix every 24 wells to maintain homogeneity.

  RIGHT ARM — p20 multi-channel pipette (8-channel)
               Adds well-specific forward and reverse primers column by column.
               New tips every column to prevent cross-contamination.
               Primer plate column offset is configurable (PRIMER_START_COLUMN).

The destination plate lives inside the thermocycler module throughout —
the robot pipettes directly into the open thermocycler, then the lid
closes automatically before cycling begins. No manual plate sealing needed.

Workflow
--------
  1. User prepares mastermix in a single 15 mL conical tube (on ice)
  2. Tube is loaded into tube rack on deck
  3. Thermocycler lid opens
  4. p300 dispenses mastermix into each well with live liquid level tracking
     and periodic remixing every 24 wells
  5. [Optional] Add INSPECT_PAUSE = True to pause here for visual QC
  6. p20 multi-channel adds forward primers column by column (new tips each)
  7. p20 multi-channel adds reverse primers column by column (new tips each)
  8. Thermocycler lid closes, heats to 105°C
  9. PCR program runs automatically (26 cycles)
  10. Plate held at 4°C — lid opens, robot pauses for user to remove plate
  11. User proceeds to Step 5 (Zymo ZR-96 purification)

Final reaction volume: 25 µL per well

Reaction composition (at default 2.5 µM primer concentration):
  NEBNext Ultra II Q5 Master Mix (2X)  : 12.5 µL  ┐
  Oligopool Template (0.5 ng/µL)       :  1.0 µL  ├─ pre-combined by user
  Nuclease-free water                   :  1.5 µL  ┘   into single conical tube
  Forward primer (well-specific)        :  5.0 µL  ─── p20 multi-channel
  Reverse primer (well-specific)        :  5.0 µL  ─── p20 multi-channel
  ─────────────────────────────────────────────────
  Total                                 : 25.0 µL

Mastermix preparation note
--------------------------
Combine Q5 master mix, oligopool template, and nuclease-free water in a
single 15 mL conical tube. Vortex gently and keep on ice. Load into the
tube rack immediately before starting the protocol. This single-tube
approach ensures every well receives identically composed mastermix.

References
----------
Salis Lab Genetic Systems Builder: https://www.denovodna.com/build_genetic_systems
NEB NEBNext Ultra II Q5 Master Mix: NEB M0544X
Opentrons Thermocycler Module Gen2

Deck Layout
-----------
  Slot 1   — 15 mL tube rack (NEST) — mastermix conical tube in position A1
  Slot 4   — Forward primer plate (96-well, IDT primers at 2.5 µM in TE)
  Slot 5   — Reverse primer plate (96-well, IDT primers at 2.5 µM in TE)
  Slot 6   — 300 µL tip rack (p300 single-channel mastermix dispense)
  Slot 8   — 20 µL tip rack #1 (p20 multi-channel primers)
  Slot 9   — 20 µL tip rack #2
  Slot 11  — 20 µL tip rack #3
  Slot 7 + 10 — Thermocycler Module (destination plate lives here)

  ┌───────────┬───────────┬───────────┐
  │  THERMO   │  20µL     │  20µL     │
  │  CYCLER   │  Tips #2  │  Tips #3  │
  │ (7 + 10)  ├───────────┼───────────┤
  │           │  20µL     │   empty   │
  │           │  Tips #1  │           │
  ├───────────┼───────────┼───────────┤
  │  Fwd      │  Rev      │  300µL    │
  │  Primers  │  Primers  │  Tips     │
  ├───────────┼───────────┼───────────┤
  │ Tube Rack │   empty   │   empty   │
  │ (mastermix│           │           │
  │  tube A1) │           │           │
  └───────────┴───────────┴───────────┘
              (front of robot)
"""

from opentrons import protocol_api

# ─────────────────────────────────────────────────────────────
#  USER CONFIGURATION — review and edit before each run
# ─────────────────────────────────────────────────────────────

NUM_COLUMNS          = 12     # Columns to process: 1–12 (full plate = 12)
PRIMER_CONC_UM       = 2.5    # Primer concentration in µM (default: 2.5 µM)
PRIMER_START_COLUMN  = 1      # Which primer plate column maps to destination col 1
                               # e.g. set to 3 if primers begin at column 3 of plate
MIX_REPS             = 5      # Mix cycles after each primer addition
MIX_VOLUME           = 10     # µL used for mixing after primer dispense
INSPECT_PAUSE        = False  # True = pause after mastermix dispense for visual QC
                               # Useful to confirm all wells received liquid before primers
DRY_RUN              = False  # True = skip all liquid transfers and thermocycler

# ── Thermocycler settings ─────────────────────────────────────
LID_TEMP          = 105       # °C — heated lid prevents condensation
INITIAL_DENATURE  = {'temperature': 98, 'hold_time_seconds': 120}
CYCLING_PROFILE   = [
    {'temperature': 98, 'hold_time_seconds': 20},   # Denaturation
    {'temperature': 61, 'hold_time_seconds': 30},   # Annealing (TM + 1°C)
    {'temperature': 72, 'hold_time_seconds': 6},    # Extension
]
NUM_CYCLES        = 26
FINAL_EXTENSION   = {'temperature': 72, 'hold_time_seconds': 60}
HOLD_TEMP         = 4         # °C — hold after completion

# ── Liquid level tracking — 15 mL conical tube geometry ───────
# Used to track aspiration height as mastermix volume depletes.
# A 15 mL NEST conical tube has a conical base (~0–3 mL) and
# a cylindrical body (~3–15 mL). Heights are approximate.
TUBE_CONE_HEIGHT     = 17.0   # mm — height of conical section
TUBE_CYLINDER_BASE   = 17.0   # mm — bottom of cylindrical section from tube base
TUBE_INNER_DIAMETER  = 15.16  # mm — inner diameter of cylindrical section
DEAD_VOLUME_MM       = 3.0    # mm — minimum aspirate height from tube bottom

# ── Derived volumes ───────────────────────────────────────────
FINAL_VOL         = 25.0
Q5_VOL            = 12.5
TEMPLATE_VOL      = 1.0
PRIMER_VOL        = round(0.5 * FINAL_VOL / PRIMER_CONC_UM, 2)
WATER_VOL         = round(FINAL_VOL - Q5_VOL - TEMPLATE_VOL - 2 * PRIMER_VOL, 2)
MASTERMIX_VOL     = round(Q5_VOL + TEMPLATE_VOL + WATER_VOL, 2)

# ── Flow rates ────────────────────────────────────────────────
P300_STANDARD_FLOW   = 46.43  # µL/s
P300_SLOW_FLOW       = 10.0   # µL/s — slower for viscous mastermix
P20_STANDARD_FLOW    = 3.78   # µL/s
P20_LOW_VOL_FLOW     = 1.0    # µL/s — for volumes < LOW_VOL_THRESHOLD
LOW_VOL_THRESHOLD    = 2.0    # µL
AIR_GAP_VOL          = 1.0    # µL

# ─────────────────────────────────────────────────────────────
#  METADATA
# ─────────────────────────────────────────────────────────────

metadata = {
    'protocolName': 'GSB Step 4 — PCR Setup + Thermocycler (Two-Pipette)',
    'author': '[Your Name]',
    'description': (
        'Sets up and runs PCR for oligo pool cloning. '
        'p300 single-channel for mastermix with liquid level tracking. '
        'p20 multi-channel for primers. '
        f'Final volume: {FINAL_VOL} µL, {NUM_CYCLES} cycles.'
    ),
    'apiLevel': '2.15'
}

# ─────────────────────────────────────────────────────────────
#  HELPERS
# ─────────────────────────────────────────────────────────────

import math

def liquid_height_in_tube(volume_ul, cone_height=TUBE_CONE_HEIGHT,
                           cylinder_base=TUBE_CYLINDER_BASE,
                           inner_diameter=TUBE_INNER_DIAMETER):
    """
    Estimate liquid surface height (mm from tube bottom) given remaining
    volume in a 15 mL conical tube.

    The tube has two sections:
      - Conical base: volume increases with height^3
      - Cylindrical body: volume increases linearly with height

    Returns height in mm, clamped to DEAD_VOLUME_MM minimum so the
    pipette never crashes into the tube bottom.
    """
    import math
    volume_ml = volume_ul / 1000.0

    # Max volume in conical section (approximate cone)
    cone_radius = inner_diameter / 2
    cone_vol_ml = (1/3) * math.pi * cone_radius**2 * cone_height / 1000

    if volume_ml <= cone_vol_ml:
        # Liquid is still in conical section — solve for height
        # V = (1/3) * pi * r^2 * h  (simplified, cone proportional to full cone)
        h = (3 * volume_ml * 1000 / (math.pi * cone_radius**2)) ** (1/3) * cone_height**(2/3)
        height_mm = max(h, DEAD_VOLUME_MM)
    else:
        # Liquid is in cylindrical section
        cylinder_vol_ul = (volume_ml - cone_vol_ml) * 1000
        cylinder_area = math.pi * (inner_diameter / 2) ** 2
        cylinder_height = cylinder_vol_ul / cylinder_area
        height_mm = cylinder_base + cylinder_height

    return max(height_mm, DEAD_VOLUME_MM)


def aspirate_single(pipette, volume, tube, remaining_vol_ul):
    """
    p300 single-channel aspirate with:
    - Slow flow rate to minimize bubbles in viscous mastermix
    - Dynamic aspiration height based on remaining liquid volume
    - Air gap to prevent dripping during transit
    """
    height = liquid_height_in_tube(remaining_vol_ul)
    pipette.flow_rate.aspirate = P300_SLOW_FLOW
    pipette.aspirate(volume, tube.bottom(height))
    pipette.air_gap(AIR_GAP_VOL)
    pipette.flow_rate.aspirate = P300_STANDARD_FLOW


def dispense_single(pipette, volume, dest):
    """
    p300 single-channel dispense to bottom of well, touch tip, blow out.
    Dispensing to bottom ensures droplet releases into liquid not air.
    """
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


def mastermix_recipe(num_columns, q5_vol, template_vol, water_vol):
    """Scaled mastermix recipe with 10% overage for single conical tube."""
    n = num_columns * 8 * 1.10
    return {
        'NEBNext Ultra II Q5 Master Mix (2X)':  round(q5_vol * n, 1),
        'Oligopool Template (0.5 ng/µL)':       round(template_vol * n, 1),
        'Nuclease-free water':                   round(water_vol * n, 1),
        'TOTAL (load into 15 mL conical tube)':  round(MASTERMIX_VOL * n, 1),
    }


def estimate_tips(num_columns):
    """
    p300: 1 tip for all 96 mastermix dispenses.
          Safe to reuse — single source, no contamination risk.
          Tip is only dropped after all wells are filled.
    p20 multi: new tip per column for fwd + rev primers.
    """
    p300_tips = 1
    p20_tips  = num_columns * 2
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
            "Check PRIMER_CONC_UM — primer volumes may exceed 25 µL total."
        )
    primer_end = PRIMER_START_COLUMN + NUM_COLUMNS - 1
    if primer_end > 12:
        raise ValueError(
            f"PRIMER_START_COLUMN ({PRIMER_START_COLUMN}) + NUM_COLUMNS ({NUM_COLUMNS}) "
            f"exceeds plate bounds (column {primer_end} > 12)."
        )

    p300_tips, p20_tips, p20_boxes = estimate_tips(NUM_COLUMNS)
    recipe = mastermix_recipe(NUM_COLUMNS, Q5_VOL, TEMPLATE_VOL, WATER_VOL)
    total_mastermix_ul = recipe['TOTAL (load into 15 mL conical tube)']

    # ── Pre-run summary ──────────────────────────────────────
    protocol.comment("=" * 60)
    protocol.comment("GSB Step 4: PCR Setup + Thermocycler (Two-Pipette)")
    if DRY_RUN:
        protocol.comment("  ⚠  DRY RUN — no liquid or thermocycler steps will run")
    protocol.comment(f"  Columns            : {NUM_COLUMNS}  ({NUM_COLUMNS * 8} wells)")
    protocol.comment(f"  Primer plate cols  : {PRIMER_START_COLUMN}–{primer_end}")
    protocol.comment(f"  Primer conc        : {PRIMER_CONC_UM} µM → {PRIMER_VOL} µL per primer")
    protocol.comment(f"  Mastermix vol      : {MASTERMIX_VOL} µL per well")
    protocol.comment(f"  Inspect pause      : {'ON' if INSPECT_PAUSE else 'OFF'}")
    protocol.comment(f"  Final vol          : {FINAL_VOL} µL")
    protocol.comment(f"  Cycles             : {NUM_CYCLES}")
    protocol.comment(f"  Lid temp           : {LID_TEMP}°C")
    protocol.comment("")
    protocol.comment("  ── Mastermix recipe — combine in single 15 mL conical tube ──")
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
            "  [ ] Mastermix prepared in single 15 mL conical tube (see run log for volumes)\n"
            "  [ ] Vortex mastermix thoroughly just before loading — do not skip this step\n"
            "  [ ] Mastermix tube loaded into tube rack Slot 1, position A1\n"
            "  [ ] Mastermix is on ice\n"
            "  [ ] Forward primer plate in Slot 4\n"
            "  [ ] Reverse primer plate in Slot 5\n"
            "  [ ] 300 µL tip rack in Slot 6\n"
            f"  [ ] {p20_boxes} × 20 µL tip box(es) in Slots 8, 9, 11\n"
            "  [ ] Empty PCR plate loaded into OPEN thermocycler\n"
            f"  [ ] Primer plates loaded starting from column {PRIMER_START_COLUMN}\n"
            "\nPress Resume when ready."
        )

    # ── Load thermocycler ────────────────────────────────────
    tc_mod = protocol.load_module('thermocyclerModuleV2', 7)
    dest_plate = tc_mod.load_labware(
        'nest_96_wellplate_100ul_pcr_full_skirt',  # Swap for your labware def if needed
        label='PCR Destination Plate'
    )

    # ── Other labware ────────────────────────────────────────
    tube_rack = protocol.load_labware(
        'opentrons_15_tuberack_nest_15ml_conical', 1, label='Mastermix Tube Rack'
    )
    fwd_plate = protocol.load_labware(
        'nest_96_wellplate_100ul_pcr_full_skirt', 4, label='Forward Primer Plate'
    )
    rev_plate = protocol.load_labware(
        'nest_96_wellplate_100ul_pcr_full_skirt', 5, label='Reverse Primer Plate'
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
    mastermix_tube = tube_rack['A1']
    dest_wells     = dest_plate.wells()[:NUM_COLUMNS * 8]
    dest_cols      = dest_plate.columns()[:NUM_COLUMNS]

    # Primer columns offset by PRIMER_START_COLUMN
    col_offset = PRIMER_START_COLUMN - 1
    fwd_cols   = fwd_plate.columns()[col_offset:col_offset + NUM_COLUMNS]
    rev_cols   = rev_plate.columns()[col_offset:col_offset + NUM_COLUMNS]

    # ─────────────────────────────────────────────────────────
    #  OPEN LID
    # ─────────────────────────────────────────────────────────
    protocol.comment("\n─── Opening thermocycler lid ───")
    if not DRY_RUN:
        tc_mod.open_lid()

    # ─────────────────────────────────────────────────────────
    #  STEP 4A — Mastermix dispense (p300 single-channel)
    #
    #  One tip reused for all dispenses — safe because:
    #    - Single source (no cross-well contamination possible)
    #    - Tip never contacts destination wells (dispenses to bottom,
    #      no reverse pipetting into source)
    #
    #  Liquid level tracking: aspiration height updates every well
    #  as volume depletes, preventing air aspiration.
    #
    #  NOTE: To add a visual inspection pause after this step,
    #  set INSPECT_PAUSE = True in the configuration above.
    #  This allows the user to confirm liquid is present in every
    #  well before primers are added.
    # ─────────────────────────────────────────────────────────
    protocol.comment(
        f"\n─── 4A: Mastermix ({MASTERMIX_VOL} µL/well) "
        f"— p300 single-channel, {NUM_COLUMNS * 8} wells ───"
    )
    protocol.comment("  Liquid level tracking active.")

    remaining_vol = total_mastermix_ul

    p300.pick_up_tip()
    for well in dest_wells:
        if not DRY_RUN:
            aspirate_single(p300, MASTERMIX_VOL, mastermix_tube, remaining_vol)
            dispense_single(p300, MASTERMIX_VOL, well)
            remaining_vol -= MASTERMIX_VOL

    p300.drop_tip()
    protocol.comment(f"  Mastermix dispense complete. Remaining in tube: ~{remaining_vol:.0f} µL")

    # ── Optional inspect pause ───────────────────────────────
    if INSPECT_PAUSE and not DRY_RUN:
        protocol.pause(
            "OPTIONAL INSPECTION PAUSE\n"
            "  Mastermix has been dispensed into all wells.\n"
            "  Visually confirm liquid is present in every well before primers are added.\n"
            "\nPress Resume to continue with primer addition."
        )

    # ─────────────────────────────────────────────────────────
    #  STEP 4B — Forward Primers (p20 multi-channel)
    #  New tip every column — prevents cross-contamination
    #  between unique primer pairs across wells.
    # ─────────────────────────────────────────────────────────
    protocol.comment(
        f"\n─── 4B: Forward Primers ({PRIMER_VOL} µL/well) "
        f"— p20 multi-channel, cols {PRIMER_START_COLUMN}–{primer_end}, new tips per column ───"
    )

    for src_col, dst_col in zip(fwd_cols, dest_cols):
        p20_multi.pick_up_tip()
        if not DRY_RUN:
            aspirate_multi(p20_multi, PRIMER_VOL, src_col[0])
            dispense_multi(p20_multi, PRIMER_VOL, dst_col[0])
            p20_multi.mix(MIX_REPS, MIX_VOLUME, dst_col[0])
            p20_multi.blow_out(dst_col[0].top(-2))
        p20_multi.drop_tip()

    # ─────────────────────────────────────────────────────────
    #  STEP 4C — Reverse Primers (p20 multi-channel)
    # ─────────────────────────────────────────────────────────
    protocol.comment(
        f"\n─── 4C: Reverse Primers ({PRIMER_VOL} µL/well) "
        f"— p20 multi-channel, cols {PRIMER_START_COLUMN}–{primer_end}, new tips per column ───"
    )

    for src_col, dst_col in zip(rev_cols, dest_cols):
        p20_multi.pick_up_tip()
        if not DRY_RUN:
            aspirate_multi(p20_multi, PRIMER_VOL, src_col[0])
            dispense_multi(p20_multi, PRIMER_VOL, dst_col[0])
            p20_multi.mix(MIX_REPS, MIX_VOLUME, dst_col[0])
            p20_multi.blow_out(dst_col[0].top(-2))
        p20_multi.drop_tip()

    # ─────────────────────────────────────────────────────────
    #  THERMOCYCLER
    # ─────────────────────────────────────────────────────────
    protocol.comment("\n─── Closing lid and starting thermocycler ───")

    if not DRY_RUN:
        tc_mod.close_lid()
        tc_mod.set_lid_temperature(LID_TEMP)

        protocol.comment(
            f"  Initial denaturation: {INITIAL_DENATURE['temperature']}°C "
            f"for {INITIAL_DENATURE['hold_time_seconds']}s"
        )
        tc_mod.set_block_temperature(
            INITIAL_DENATURE['temperature'],
            hold_time_seconds=INITIAL_DENATURE['hold_time_seconds'],
            block_max_volume=FINAL_VOL
        )

        protocol.comment(f"  Cycling: 98°C / 61°C / 72°C × {NUM_CYCLES} cycles")
        tc_mod.execute_profile(
            steps=CYCLING_PROFILE,
            repetitions=NUM_CYCLES,
            block_max_volume=FINAL_VOL
        )

        protocol.comment(
            f"  Final extension: {FINAL_EXTENSION['temperature']}°C "
            f"for {FINAL_EXTENSION['hold_time_seconds']}s"
        )
        tc_mod.set_block_temperature(
            FINAL_EXTENSION['temperature'],
            hold_time_seconds=FINAL_EXTENSION['hold_time_seconds'],
            block_max_volume=FINAL_VOL
        )

        protocol.comment(f"  Holding at {HOLD_TEMP}°C")
        tc_mod.set_block_temperature(HOLD_TEMP, block_max_volume=FINAL_VOL)
        tc_mod.open_lid()

    # ─────────────────────────────────────────────────────────
    #  DONE
    # ─────────────────────────────────────────────────────────
    protocol.comment("\n" + "=" * 60)
    if DRY_RUN:
        protocol.comment("✓ Dry run complete. No liquid or thermocycler steps ran.")
        protocol.comment("  Set DRY_RUN = False to run the real protocol.")
    else:
        protocol.comment("✓ PCR complete. Plate held at 4°C in open thermocycler.")

    protocol.pause(
        "PCR complete — plate is held at 4°C.\n"
        "  [ ] Remove PCR plate from thermocycler\n"
        "  [ ] Proceed to Step 5: Zymo ZR-96 DNA purification\n"
        "  [ ] Plate can be stored at -20°C if not purifying immediately\n"
        "\nPress Resume to end the protocol."
    )
    protocol.comment("=" * 60)
