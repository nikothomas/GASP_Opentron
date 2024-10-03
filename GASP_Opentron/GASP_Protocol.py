# Import the required libraries from the Opentrons API
import csv
import os
from collections import namedtuple
import math

import opentrons
from opentrons import protocol_api

# Metadata
metadata = {
    'protocolName': 'GASP Protocol',
    'author': 'Nikolas Yanek-Chrones <research@icarai.io>',
    'description': 'Automated setup for GASP PCR reactions using Opentron OT-2 with variable DNA volumes.',
    'apiLevel': '2.20'  # Use the latest API level supported by your OT-2
}

# Constants
# Well Positions
ROWS = ['A', 'B', 'C', 'D']
COLUMNS = ['1', '2', '3', '4', '5', '6']
SAMPLE_POSITIONS = [f"{row}{col}" for row in ROWS for col in COLUMNS]
NUCLEASE_FREE_LOC_CODE = 'A1'
ONETAQ_LOC_CODE = 'A2'
REVERSE_PRIMER_LOC_CODE = 'B1'
FORWARD_PRIMER_LOC_CODE = 'B2'
# Calculation Constants
TOTAL_REACTION_VOLUME = 25.0  # µL
DNA_SAMPLE_INITIAL_VOLUME = 20.0  # µL
DNA_MASS_TO_VOLUME_CONVERSION_FACTOR = 10.0  # µL·ng
FORWARD_PRIMER_VOLUME_PER_SAMPLE = 1.0  # µL
REVERSE_PRIMER_VOLUME_PER_SAMPLE = 1.0  # µL
ONETAQ_VOLUME_PER_SAMPLE = 12.5  # µL
MIX_REPETITIONS = 3
MIX_VOLUME = 10.0  # µL

# String Constants
ONETAQ_MASTER_MIX = 'onetaq_master_mix'
NUCLEASE_FREE_WATER = 'nuclease_free_water'
FORWARD_PRIMER = 'forward_primer'
REVERSE_PRIMER = 'reverse_primer'

# Dictionary Keys
INITIAL_VOLUME_DICT_KEY = 'initial_volume'
USED_VOLUME_DICT_KEY = 'used_volume'
ASPIRATION_VOLUME_DICT_KEY = 'aspiration_volume'
WATER_VOLUME_DICT_KEY = 'water_volume'

# Parameter Defaults
DEFAULT_DNA_MASS = 5.0
required_vol_tuple = namedtuple('required_vols', ['forward_primer', 'reverse_primer', 'onetaq_master_mix', 'nuclease_free_water'])

# Helper function to get the approximate volumes we will require based on the number of samples
def get_volumes_needed(num_samples: int, total_water_volume: float) -> required_vol_tuple:
    # Calculate required volumes with a small excess (+1 sample) for pipetting errors
    required_forward_primer_vol = (num_samples + 1) * FORWARD_PRIMER_VOLUME_PER_SAMPLE
    required_reverse_primer_vol = (num_samples + 1) * REVERSE_PRIMER_VOLUME_PER_SAMPLE
    required_onetaq_vol = (num_samples + 1) * ONETAQ_VOLUME_PER_SAMPLE
    # Round up total_water_volume to the nearest hundred microliters
    required_nuclease_free_water_vol = math.ceil(total_water_volume / 100.0) * 100
    return required_vol_tuple(required_forward_primer_vol,
                              required_reverse_primer_vol,
                              required_onetaq_vol,
                              required_nuclease_free_water_vol)

def add_parameters(parameters):
    for pos in SAMPLE_POSITIONS[1:8]:
        parameters.add_float(
            variable_name=f"mass_{pos}",
            display_name=f"{pos} Sample Mass",
            description="The Qubit quantified mass of the sample, 0 means no sample present.",
            default=DEFAULT_DNA_MASS,
            minimum=0,
            maximum=20,
            unit="ng"
        )
    for pos in SAMPLE_POSITIONS[8:]:
        parameters.add_float(
            variable_name=f"mass_{pos}",
            display_name=f"{pos} Sample Mass",
            description="The Qubit quantified mass of the sample, 0 means no sample present.",
            default=0.0,
            minimum=0,
            maximum=20,
            unit="ng"
        )

def run(protocol: protocol_api.ProtocolContext):
    # Load tip racks
    tip_rack_20ul = protocol.load_labware(load_name='opentrons_96_tiprack_20ul', location=4)

    # Define liquids
    nuclease_free_water = protocol.define_liquid(
        name="Nuclease Free Water", description="", display_color="#052FFF"
    )
    reverse_primer = protocol.define_liquid(
        name="Reverse Primer", description="", display_color="#FF05FF"
    )
    forward_primer = protocol.define_liquid(
        name="Forward Primer", description="", display_color="#E1FF05"
    )
    onetaq_master_mix = protocol.define_liquid(
        name="OneTaq Master Mix", description="", display_color="#CCAB68"
    )

    # Define Template DNA liquids
    template_dna_liquids = {}
    for pos in SAMPLE_POSITIONS:
        template_dna_liquids[pos] = protocol.define_liquid(
            name=f"Sample {pos} liquid", description="", display_color="#163D20",
        )

    # Load pipettes
    p20_single = protocol.load_instrument(
        instrument_name='p20_single_gen2',
        mount='right',
        tip_racks=[tip_rack_20ul]
    )
    p20_multi = protocol.load_instrument(
        instrument_name='p20_multi_gen2',
        mount='left',
        tip_racks=[tip_rack_20ul]
    )

    # Load labware
    pcr_plate = protocol.load_labware(
        load_name='eppendorflobind_96_wellplate_150ul',
        location=1,
        label='PCR Plate'
    )
    tube_rack_15_50ml = protocol.load_labware(
        load_name='opentrons_10_tuberack_nest_4x50ml_6x15ml_conical',
        location=2,
        label='15/50 mL Tube Rack'
    )
    tube_rack_2ml = protocol.load_labware(
        load_name='opentrons_24_tuberack_nest_1.5ml_snapcap',
        location=5,
        label='1.5 mL Tube Rack'
    )

    # Initialize tracking variables
    num_samples = 0
    pcr_well_compositions = {}
    total_master_mix_used = 0.0
    total_forward_primer_used = 0.0
    total_reverse_primer_used = 0.0
    total_water_used = 0.0
    total_dna_used_per_sample = []
    dna_info = {}  # key: tube_location, value: volume info

    # Build sample_masses dictionary
    sample_masses = {}
    for pos in SAMPLE_POSITIONS:
        mass = getattr(protocol.params, f"mass_{pos}")
        sample_masses[pos] = mass

    # Process samples and calculate volumes
    active_samples = []
    total_water_volume_needed = 0.0
    for pos in SAMPLE_POSITIONS:
        mass = sample_masses[pos]
        if mass > 0.01:
            num_samples += 1
            initial_volume = DNA_SAMPLE_INITIAL_VOLUME
            aspiration_volume = DNA_MASS_TO_VOLUME_CONVERSION_FACTOR / mass
            water_volume = TOTAL_REACTION_VOLUME - (
                    ONETAQ_VOLUME_PER_SAMPLE + FORWARD_PRIMER_VOLUME_PER_SAMPLE + REVERSE_PRIMER_VOLUME_PER_SAMPLE + aspiration_volume)
            total_water_volume_needed += water_volume
            active_samples.append(pos)
        else:
            initial_volume = 0.0
            aspiration_volume = 0.0
            water_volume = 0.0
        dna_info[pos] = {
            INITIAL_VOLUME_DICT_KEY: initial_volume,
            USED_VOLUME_DICT_KEY: 0.0,
            ASPIRATION_VOLUME_DICT_KEY: aspiration_volume,
            WATER_VOLUME_DICT_KEY: water_volume
        }

    # Calculate required reagent volumes using the helper function
    required_volumes = get_volumes_needed(num_samples, total_water_volume_needed)

    protocol.pause(f"\nPlease load the following to the tube rack:\n"
                   f"{required_volumes.nuclease_free_water} µL {nuclease_free_water.name} into the tube rack at position {NUCLEASE_FREE_LOC_CODE}\n"
                   f"{required_volumes.onetaq_master_mix} µL {onetaq_master_mix.name} into the tube rack at position {ONETAQ_LOC_CODE}\n"
                   f"{required_volumes.forward_primer} µL {forward_primer.name} into the tube rack at position {FORWARD_PRIMER_LOC_CODE}\n"
                   f"{required_volumes.reverse_primer} µL {reverse_primer.name} into the tube rack at position {REVERSE_PRIMER_LOC_CODE}\n")

    # Assign reagents to specific wells and load calculated volumes
    nuclease_free_water_wells = tube_rack_15_50ml.wells_by_name()[NUCLEASE_FREE_LOC_CODE]
    nuclease_free_water_wells.load_liquid(liquid=nuclease_free_water, volume=required_volumes.nuclease_free_water)

    onetaq_master_mix_wells = tube_rack_15_50ml.wells_by_name()[ONETAQ_LOC_CODE]
    onetaq_master_mix_wells.load_liquid(liquid=onetaq_master_mix, volume=required_volumes.onetaq_master_mix)

    reverse_primer_wells = tube_rack_15_50ml.wells_by_name()[REVERSE_PRIMER_LOC_CODE]
    reverse_primer_wells.load_liquid(liquid=reverse_primer, volume=required_volumes.reverse_primer)

    forward_primer_wells = tube_rack_15_50ml.wells_by_name()[FORWARD_PRIMER_LOC_CODE]
    forward_primer_wells.load_liquid(liquid=forward_primer, volume=required_volumes.forward_primer)

    # Load Template DNA liquids into wells
    template_dna_wells = {}
    for pos in SAMPLE_POSITIONS:
        template_dna_wells[pos] = tube_rack_2ml.wells_by_name()[pos]
        volume = dna_info[pos][INITIAL_VOLUME_DICT_KEY]
        if volume > 0:
            template_dna_wells[pos].load_liquid(liquid=template_dna_liquids[pos], volume=volume)

    # Define destination wells on the PCR plate
    all_destination_wells = [f"{row}{col}" for col in ['10', '11', '12'] for row in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']]
    destination_wells = all_destination_wells[:len(active_samples)]

    # Step 1: Add OneTaq 2X Master Mix to each PCR well
    protocol.comment("\n---------------------------------------\n"
                     "Adding OneTaq 2X Master Mix to PCR wells\n"
                     "----------------------------------------\n")
    p20_single.pick_up_tip()
    for well_name in destination_wells:
        dest_well = pcr_plate.wells_by_name()[well_name]
        p20_single.aspirate(ONETAQ_VOLUME_PER_SAMPLE, onetaq_master_mix_wells.bottom(z=1))
        p20_single.dispense(ONETAQ_VOLUME_PER_SAMPLE, dest_well.bottom(z=1))
        p20_single.touch_tip(dest_well)
        total_master_mix_used += ONETAQ_VOLUME_PER_SAMPLE
        pcr_well_compositions[well_name] = {
            'well_name': well_name,
            'master_mix': ONETAQ_VOLUME_PER_SAMPLE,
            'forward_primer': 0.0,
            'reverse_primer': 0.0,
            'dna_volume': 0.0,
            'water_volume': 0.0
        }
    p20_single.drop_tip()

    # Step 2: Add Forward Primer to each PCR well
    protocol.comment("\n---------------------------------\n"
                     "Adding Forward Primer to PCR wells\n"
                     "----------------------------------\n")
    p20_single.pick_up_tip()
    for well_name in destination_wells:
        dest_well = pcr_plate.wells_by_name()[well_name]
        p20_single.aspirate(FORWARD_PRIMER_VOLUME_PER_SAMPLE, forward_primer_wells.bottom(z=1))
        p20_single.dispense(FORWARD_PRIMER_VOLUME_PER_SAMPLE, dest_well.bottom(z=1))
        p20_single.touch_tip(dest_well)
        total_forward_primer_used += FORWARD_PRIMER_VOLUME_PER_SAMPLE
        pcr_well_compositions[well_name]['forward_primer'] = FORWARD_PRIMER_VOLUME_PER_SAMPLE
    p20_single.drop_tip()

    # Step 3: Add Reverse Primer to each PCR well
    protocol.comment("\n---------------------------------\n"
                     "Adding Reverse Primer to PCR wells\n"
                     "----------------------------------\n")
    p20_single.pick_up_tip()
    for well_name in destination_wells:
        dest_well = pcr_plate.wells_by_name()[well_name]
        p20_single.aspirate(REVERSE_PRIMER_VOLUME_PER_SAMPLE, reverse_primer_wells.bottom(z=1))
        p20_single.dispense(REVERSE_PRIMER_VOLUME_PER_SAMPLE, dest_well.bottom(z=1))
        p20_single.touch_tip(dest_well)
        total_reverse_primer_used += REVERSE_PRIMER_VOLUME_PER_SAMPLE
        pcr_well_compositions[well_name]['reverse_primer'] = REVERSE_PRIMER_VOLUME_PER_SAMPLE
    p20_single.drop_tip()

    # Step 4: Transfer Variable Volumes of Template DNA
    protocol.comment("\n-------------------------------------\n"
                     "Transferring Template DNA to PCR wells\n"
                     "--------------------------------------\n")
    for sample_pos, dest_well_name in zip(active_samples, destination_wells):
        dest_well = pcr_plate.wells_by_name()[dest_well_name]
        source_well = template_dna_wells[sample_pos]
        vol_dna = dna_info[sample_pos][ASPIRATION_VOLUME_DICT_KEY]
        if vol_dna > 0:
            p20_single.pick_up_tip()
            p20_single.aspirate(vol_dna, source_well.bottom(z=1))
            p20_single.dispense(vol_dna, dest_well.bottom(z=1))
            p20_single.touch_tip(dest_well)
            p20_single.drop_tip()
            total_dna_used_per_sample.append(vol_dna)
            pcr_well_compositions[dest_well_name]['dna_volume'] = vol_dna
            dna_info[sample_pos][USED_VOLUME_DICT_KEY] += vol_dna

    # Step 5: Add Nuclease-Free Water to each PCR well
    protocol.comment("\n--------------------------------------------\n"
                     "Transferring Nuclease Free Water to PCR wells\n"
                     "---------------------------------------------\n")
    for sample_pos, dest_well_name in zip(active_samples, destination_wells):
        dest_well = pcr_plate.wells_by_name()[dest_well_name]
        vol_water = dna_info[sample_pos][WATER_VOLUME_DICT_KEY]
        if vol_water > 0:
            p20_single.pick_up_tip()
            p20_single.aspirate(vol_water, nuclease_free_water_wells.bottom(z=1))
            p20_single.dispense(vol_water, dest_well.bottom(z=1))
            p20_single.touch_tip(dest_well)
            p20_single.drop_tip()
            total_water_used += vol_water
            pcr_well_compositions[dest_well_name]['water_volume'] = vol_water

    # Step 6: Mix the contents of each well using the multichannel pipette
    columns_to_mix = set(well_name[1:] for well_name in destination_wells)
    for column in columns_to_mix:
        dest_wells = pcr_plate.columns_by_name()[column]
        p20_multi.pick_up_tip()
        p20_multi.mix(MIX_REPETITIONS, MIX_VOLUME, dest_wells[0].bottom(z=1))
        p20_multi.touch_tip()
        p20_multi.drop_tip()

    # Verbose Output at the End of the Protocol
    protocol.comment("PCR Well Compositions:")
    for well_name in sorted(pcr_well_compositions.keys()):
        comp = pcr_well_compositions[well_name]
        total_volume = (
                comp['master_mix']
                + comp['forward_primer']
                + comp['reverse_primer']
                + comp['dna_volume']
                + comp['water_volume']
        )
        protocol.comment(f"Well {well_name}:")
        protocol.comment(f"  Master Mix: {comp['master_mix']} µL")
        protocol.comment(f"  Forward Primer: {comp['forward_primer']} µL")
        protocol.comment(f"  Reverse Primer: {comp['reverse_primer']} µL")
        protocol.comment(f"  Template DNA: {comp['dna_volume']} µL")
        protocol.comment(f"  Nuclease-Free Water: {comp['water_volume']} µL")
        protocol.comment(f"  Total Volume: {total_volume} µL")

    # Final volumes for all tubes with liquid
    protocol.comment("Final volumes for all tubes with liquid:")
    onetaq_master_mix_final = required_volumes.onetaq_master_mix - total_master_mix_used
    protocol.comment(f"OneTaq Master Mix: Final Volume = {onetaq_master_mix_final} µL")

    forward_primer_final = required_volumes.forward_primer - total_forward_primer_used
    protocol.comment(f"Forward Primer: Final Volume = {forward_primer_final} µL")

    reverse_primer_final = required_volumes.reverse_primer - total_reverse_primer_used
    protocol.comment(f"Reverse Primer: Final Volume = {reverse_primer_final} µL")

    nuclease_free_water_final = required_volumes.nuclease_free_water - total_water_used
    protocol.comment(f"Nuclease-Free Water: Final Volume = {nuclease_free_water_final} µL")

    for tube_location, volumes in dna_info.items():
        initial_volume = volumes[INITIAL_VOLUME_DICT_KEY]
        used_volume = volumes[USED_VOLUME_DICT_KEY]
        final_volume = initial_volume - used_volume
        if initial_volume > 0:
            protocol.comment(f"Template DNA Tube {tube_location}: Final Volume = {final_volume} µL")

    # End of Protocol
