# Import the required libraries from the Opentrons API
import math
from collections import namedtuple
from opentrons import protocol_api

# Metadata
metadata = {
    'protocolName': 'GASP Protocol',
    'author': 'Nikolas Yanek-Chrones <research@icarai.io>',
    'description': 'Automated setup for GASP PCR reactions using Opentron OT-2 with variable DNA volumes.',
    'apiLevel': '2.20'
}

# Constants

# Labware constants
TIP_RACK_SLOT = 4
PCR_PLATE_SLOT = 1
TUBE_RACK_SLOT = 5
TIP_RACK_TYPE = 'opentrons_96_tiprack_20ul'
PCR_PLATE_TYPE = 'eppendorflobind_96_wellplate_150ul'
TUBE_RACK_TYPE = 'opentrons_24_tuberack_nest_1.5ml_snapcap'

# Pipette constants
P20_SINGLE_TYPE = 'p20_single_gen2'
P20_MULTI_TYPE = 'p20_multi_gen2'
P20_SINGLE_MOUNT = 'right'
P20_MULTI_MOUNT = 'left'

# Row and Column definitions
ROWS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
COLUMNS = ['1', '2', '3', '4', '5', '6']
SAMPLE_POSITIONS = [f"{row}{col}" for row in ROWS[:4] for col in COLUMNS]  # Adjusted to 'A'-'D'

# Reagent positions
NUCLEASE_FREE_LOC_CODE = 'A1'
ONETAQ_LOC_CODE = 'A2'
REVERSE_PRIMER_LOC_CODE = 'B1'
FORWARD_PRIMER_LOC_CODE = 'B2'

# Reagent Volumes
TOTAL_REACTION_VOLUME = 25.0  # µL
DNA_MASS_TO_VOLUME_CONVERSION_FACTOR = 10.0  # µL·ng
FORWARD_PRIMER_VOLUME = 1.0  # µL
REVERSE_PRIMER_VOLUME = 1.0  # µL
ONETAQ_VOLUME = 12.5  # µL
MIX_REPETITIONS = 3
MIX_VOLUME = 10.0  # µL
DEFAULT_DNA_MASS = 5.0

# Thresholds and initial volumes
MASS_THRESHOLD = 0.01  # ng
INITIAL_DNA_VOLUME = 20.0  # µL

# Destination plate layout
DESTINATION_COLUMNS = ['10', '11', '12']
DESTINATION_ROWS = ROWS  # Using all rows from 'A' to 'H'

# Mixing parameters
MIXING_Z_POSITION = 1  # mm above the bottom

# Liquid colors
NUCLEASE_FREE_WATER_COLOR = '#052FFF'
REVERSE_PRIMER_COLOR = '#FF05FF'
FORWARD_PRIMER_COLOR = '#E1FF05'
ONETAQ_MASTER_MIX_COLOR = '#CCAB68'
TEMPLATE_DNA_COLOR = '#163D20'

required_vol_tuple = namedtuple('required_vols', ['forward_primer', 'reverse_primer', 'onetaq_master_mix', 'nuclease_free_water'])

# Helper function to calculate required volumes with a small excess for pipetting errors
def get_volumes_needed(num_samples: int, total_water_volume: float) -> required_vol_tuple:
    return required_vol_tuple(
        forward_primer=(num_samples + 1) * FORWARD_PRIMER_VOLUME,
        reverse_primer=(num_samples + 1) * REVERSE_PRIMER_VOLUME,
        onetaq_master_mix=(num_samples + 1) * ONETAQ_VOLUME,
        nuclease_free_water=math.ceil(total_water_volume / 100.0) * 100
    )

def run(protocol: protocol_api.ProtocolContext):
    # Load labware
    tip_rack_20ul = protocol.load_labware(TIP_RACK_TYPE, TIP_RACK_SLOT)
    pcr_plate = protocol.load_labware(PCR_PLATE_TYPE, PCR_PLATE_SLOT, 'PCR Plate')
    tube_rack_2ml = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_SLOT, '1.5 mL Tube Rack')

    # Load pipettes
    p20_single = protocol.load_instrument(P20_SINGLE_TYPE, P20_SINGLE_MOUNT, tip_racks=[tip_rack_20ul])
    p20_multi = protocol.load_instrument(P20_MULTI_TYPE, P20_MULTI_MOUNT, tip_racks=[tip_rack_20ul])

    # Define reagents
    nuclease_free_water = protocol.define_liquid("Nuclease Free Water", "", NUCLEASE_FREE_WATER_COLOR)
    reverse_primer = protocol.define_liquid("Reverse Primer", "", REVERSE_PRIMER_COLOR)
    forward_primer = protocol.define_liquid("Forward Primer", "", FORWARD_PRIMER_COLOR)
    onetaq_master_mix = protocol.define_liquid("OneTaq Master Mix", "", ONETAQ_MASTER_MIX_COLOR)
    template_dna_liquids = {pos: protocol.define_liquid(f"Sample {pos} liquid", "", TEMPLATE_DNA_COLOR) for pos in SAMPLE_POSITIONS}

    # Initialize tracking variables
    sample_masses = {pos: getattr(protocol.params, f"mass_{pos}", DEFAULT_DNA_MASS) for pos in SAMPLE_POSITIONS}
    dna_info = {}
    active_samples = []
    total_water_volume_needed = 0.0

    # Process samples and calculate volumes
    for pos, mass in sample_masses.items():
        if mass > MASS_THRESHOLD:
            aspiration_volume = DNA_MASS_TO_VOLUME_CONVERSION_FACTOR / mass
            water_volume = TOTAL_REACTION_VOLUME - (ONETAQ_VOLUME + FORWARD_PRIMER_VOLUME + REVERSE_PRIMER_VOLUME + aspiration_volume)
            total_water_volume_needed += water_volume
            active_samples.append(pos)
            initial_volume = INITIAL_DNA_VOLUME
        else:
            aspiration_volume = water_volume = initial_volume = 0.0
        dna_info[pos] = {
            'initial_volume': initial_volume,
            'aspiration_volume': aspiration_volume,
            'water_volume': water_volume
        }

    # Calculate required reagent volumes
    required_volumes = get_volumes_needed(len(active_samples), total_water_volume_needed)
    protocol.pause(f"Please load the following to the tube rack:\n"
                   f"[{required_volumes.nuclease_free_water} µL {nuclease_free_water.name} to {NUCLEASE_FREE_LOC_CODE}]\n"
                   f"[{required_volumes.onetaq_master_mix} µL {onetaq_master_mix.name} to {ONETAQ_LOC_CODE}]\n"
                   f"[{required_volumes.forward_primer} µL {forward_primer.name} to {FORWARD_PRIMER_LOC_CODE}]\n"
                   f"[{required_volumes.reverse_primer} µL {reverse_primer.name} to {REVERSE_PRIMER_LOC_CODE}]")

    # Load reagents into wells
    reagents = {
        'nuclease_free_water': (tube_rack_2ml.wells_by_name()[NUCLEASE_FREE_LOC_CODE], nuclease_free_water, required_volumes.nuclease_free_water),
        'onetaq_master_mix': (tube_rack_2ml.wells_by_name()[ONETAQ_LOC_CODE], onetaq_master_mix, required_volumes.onetaq_master_mix),
        'reverse_primer': (tube_rack_2ml.wells_by_name()[REVERSE_PRIMER_LOC_CODE], reverse_primer, required_volumes.reverse_primer),
        'forward_primer': (tube_rack_2ml.wells_by_name()[FORWARD_PRIMER_LOC_CODE], forward_primer, required_volumes.forward_primer)
    }

    for name, (well, liquid, volume) in reagents.items():
        well.load_liquid(liquid=liquid, volume=volume)

    # Load template DNA liquids into wells
    for pos in active_samples:
        volume = dna_info[pos]['initial_volume']
        if volume > 0:
            tube_rack_2ml.wells_by_name()[pos].load_liquid(liquid=template_dna_liquids[pos], volume=volume)

    # Define destination wells on the PCR plate
    destination_wells = [f"{row}{col}" for col in DESTINATION_COLUMNS for row in DESTINATION_ROWS][:len(active_samples)]
    dest_well_objs = [pcr_plate.wells_by_name()[well] for well in destination_wells]

    # Step 1: Distribute reagents
    for reagent in reagents:
        if reagent == "nuclease_free_water":
            for sample_pos, dest_well in zip(active_samples, dest_well_objs):
                vol_water = dna_info[sample_pos]['water_volume']
                if vol_water > 0:
                    p20_single.transfer(volume=vol_water, source=reagents['nuclease_free_water'][0], dest=dest_well)
        else:
            p20_single.distribute(volume=reagents[reagent][2], source=reagents[reagent][0], dest=dest_well_objs)

    # Step 2: Transfer Variable Volumes of Template DNA
    for sample_pos, dest_well in zip(active_samples, dest_well_objs):
        source_well = tube_rack_2ml.wells_by_name()[sample_pos]
        vol_dna = dna_info[sample_pos]['aspiration_volume']
        if vol_dna > 0:
            p20_single.transfer(volume=vol_dna, source=source_well, dest=dest_well)

    # Step 3: Mix the contents of each well using the multichannel pipette
    columns_to_mix = set(well[1:] for well in destination_wells)
    for column in columns_to_mix:
        dest_wells = pcr_plate.columns_by_name()[column]
        p20_multi.pick_up_tip()
        p20_multi.mix(MIX_REPETITIONS, MIX_VOLUME, dest_wells[0].bottom(z=MIXING_Z_POSITION))
        p20_multi.drop_tip()
