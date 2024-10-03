# Import the required libraries from the Opentrons API
import csv
import os
from collections import namedtuple

from opentrons import protocol_api

# Metadata
metadata = {
    'protocolName': 'GASP Protocol with Variable DNA Volumes from CSV',
    'author': 'Nikolas Yanek-Chrones <research@icarai.io>',
    'description': 'Automated setup for GASP PCR reactions using Opentron OT-2 with variable DNA volumes from CSV input',
    'apiLevel': '2.20'  # Use the latest API level supported by your OT-2
}

# Locations within wells/plates
A1_LOC_CODE = 'A1'
A2_LOC_CODE = 'A2'
A3_LOC_CODE = 'A3'
A4_LOC_CODE = 'A4'
A5_LOC_CODE = 'A5'
A6_LOC_CODE = 'A6'
B1_LOC_CODE = 'B1'
B2_LOC_CODE = 'B2'
B3_LOC_CODE = 'B3'
B4_LOC_CODE = 'B4'
B5_LOC_CODE = 'B5'
B6_LOC_CODE = 'B6'
C1_LOC_CODE = 'C1'
C2_LOC_CODE = 'C2'
C3_LOC_CODE = 'C3'
C4_LOC_CODE = 'C4'
C5_LOC_CODE = 'C5'
C6_LOC_CODE = 'C6'
D1_LOC_CODE = 'D1'
D2_LOC_CODE = 'D2'
D3_LOC_CODE = 'D3'
D4_LOC_CODE = 'D4'
D5_LOC_CODE = 'D5'
D6_LOC_CODE = 'D6'
NUCLEASE_FREE_LOC_CODE = 'A1'
ONETAQ_LOC_CODE = 'A2'
FORWARD_PRIMER_LOC_CODE = 'B1'
REVERSE_PRIMER_LOC_CODE = 'B2'

# Calculation Constants
FORWARD_PRIMER_VOLUME_PER_SAMPLE = 1  # 1 µL of forward primer per sample for 25 µL reaction
REVERSE_PRIMER_VOLUME_PER_SAMPLE = 1  # 1 µL of reverse primer per sample for 25 µL reaction
REQUIRED_ONETAQ_VOLUME_PER_SAMPLE = 12.5  # 12.5 µL of OneTaq per sample for 25 µL reaction

# String Constants
ONETAQ_MASTER_MIX = 'onetaq_master_mix'
NUCLEASE_FREE_WATER = 'nuclease_free_water'
FORWARD_PRIMER = 'forward_primer'
REVERSE_PRIMER = 'reverse_primer'

# Dictionary Keys
INITIAL_VOLUME_DICT_KEY = 'initial_volume'
USED_VOLUME_DICT_KEY = 'initial_volume'
ASPIRATION_VOLUME_DICT_KEY = 'aspiration_volume'
WATER_VOLUME_DICT_KEY = 'water_volume'

# Parameter Defaults
DEFAULT_DNA_MASS = 5.0
required_vol_tuple = namedtuple('required_vols', ['forward_primer', 'reverse_primer', 'onetaq_master_mix', 'nuclease_free_water'])

# Helper function to get the approximate volumes we will require based on the number of samples
def get_volumes_needed(num_samples: int) -> required_vol_tuple:
    # Calculate required volumes with a small excess (+1 sample) for pipetting errors
    required_forward_primer_vol = (num_samples + 1) * FORWARD_PRIMER_VOLUME_PER_SAMPLE
    required_reverse_primer_vol = (num_samples + 1) * REVERSE_PRIMER_VOLUME_PER_SAMPLE
    required_onetaq_vol = (num_samples + 1) * REQUIRED_ONETAQ_VOLUME_PER_SAMPLE
    required_nuclease_free_water_vol = 500
    return required_vol_tuple(required_forward_primer_vol,
                              required_reverse_primer_vol,
                              required_onetaq_vol,
                              required_nuclease_free_water_vol)

def add_parameters(parameters):
    # Samples A1-A6
    parameters.add_float(
        variable_name="mass_A1",
        display_name="A1 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_A2",
        display_name="A2 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_A3",
        display_name="A3 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_A4",
        display_name="A4 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_A5",
        display_name="A5 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_A6",
        display_name="A6 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    # Samples B1-B6
    parameters.add_float(
        variable_name="mass_B1",
        display_name="B1 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_B2",
        display_name="B2 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_B3",
        display_name="B3 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_B4",
        display_name="B4 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_B5",
        display_name="B5 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_B6",
        display_name="B6 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    # Samples C1-C6
    parameters.add_float(
        variable_name="mass_C1",
        display_name="C1 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_C2",
        display_name="C2 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_C3",
        display_name="C3 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_C4",
        display_name="C4 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_C5",
        display_name="C5 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_C6",
        display_name="C6 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    # Samples D1-D6
    parameters.add_float(
        variable_name="mass_D1",
        display_name="D1 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_D2",
        display_name="D2 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_D3",
        display_name="D3 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_D4",
        display_name="D4 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_D5",
        display_name="D5 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

    parameters.add_float(
        variable_name="mass_D6",
        display_name="D6 Sample Mass",
        description="The Qubit quantified mass of the sample, 0 means no sample present.",
        default=DEFAULT_DNA_MASS,
        minimum=0,
        maximum=20,
        unit="ng"
    )

def run(protocol: protocol_api.ProtocolContext):
    # Load tip racks
    tip_rack_20ul = protocol.load_labware(load_name='opentrons_96_tiprack_20ul', location=4)

    # Define liquids
    template_dna = protocol.define_liquid(
        name="Template DNA", description="", display_color="#FF3105"
    )
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
    Template_DNA_liquid_A1 = protocol.define_liquid(
        name="Sample A1 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_A2 = protocol.define_liquid(
        name="Sample A2 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_A3 = protocol.define_liquid(
        name="Sample A3 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_A4 = protocol.define_liquid(
        name="Sample A4 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_A5 = protocol.define_liquid(
        name="Sample A5 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_A6 = protocol.define_liquid(
        name="Sample A6 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_B1 = protocol.define_liquid(
        name="Sample B1 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_B2 = protocol.define_liquid(
        name="Sample B2 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_B3 = protocol.define_liquid(
        name="Sample B3 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_B4 = protocol.define_liquid(
        name="Sample B4 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_B5 = protocol.define_liquid(
        name="Sample B5 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_B6 = protocol.define_liquid(
        name="Sample B6 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_C1 = protocol.define_liquid(
        name="Sample C1 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_C2 = protocol.define_liquid(
        name="Sample C2 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_C3 = protocol.define_liquid(
        name="Sample C3 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_C4 = protocol.define_liquid(
        name="Sample C4 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_C5 = protocol.define_liquid(
        name="Sample C5 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_C6 = protocol.define_liquid(
        name="Sample C6 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_D1 = protocol.define_liquid(
        name="Sample D1 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_D2 = protocol.define_liquid(
        name="Sample D2 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_D3 = protocol.define_liquid(
        name="Sample D3 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_D4 = protocol.define_liquid(
        name="Sample D4 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_D5 = protocol.define_liquid(
        name="Sample D5 liquid", description="", display_color="#163D20",
    )
    Template_DNA_liquid_D6 = protocol.define_liquid(
        name="Sample D6 liquid", description="", display_color="#163D20",
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
        load_name='opentrons_96_wellplate_200ul_pcr_full_skirt',
        location=1,
        label='PCR Plate'
    )
    tube_rack_15_50ml = protocol.load_labware(
        load_name='opentrons_10_tuberack_nest_4x50ml_6x15ml_conical',
        location=2,
        label='15/50 mL Tube Rack'
    )
    tube_rack_2ml = protocol.load_labware(
        load_name='opentrons_24_tuberack_nest_2ml_screwcap',
        location=5,
        label='2 mL Tube Rack'
    )

    # Initialize tracking variables
    num_samples = 0
    pcr_well_compositions = {}
    total_master_mix_used = 0.0
    total_forward_primer_used = 0.0
    total_reverse_primer_used = 0.0
    total_water_used = 0.0
    total_dna_used_per_sample = []
    # For DNA tubes, create a dict to track volumes used and final volumes
    dna_info: dict[str:[dict[str:float, str:float, str:float]]] = {}  # key: tube_location, value: {'initial_volume':20.0, 'used_volume':0.0, 'aspiration_volume: 1.0'}

    sample_masses = [protocol.params.mass_A1, protocol.params.mass_A2, protocol.params.mass_A3, protocol.params.mass_A4, protocol.params.mass_A5, protocol.params.mass_A6,
                      protocol.params.mass_B1, protocol.params.mass_B2, protocol.params.mass_B3, protocol.params.mass_B4, protocol.params.mass_B5, protocol.params.mass_B6,
                      protocol.params.mass_C1, protocol.params.mass_C2, protocol.params.mass_C3, protocol.params.mass_C4, protocol.params.mass_C5, protocol.params.mass_C6,
                      protocol.params.mass_D1, protocol.params.mass_D2, protocol.params.mass_D3, protocol.params.mass_D4, protocol.params.mass_D5, protocol.params.mass_D6,]

    sample_tubes = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6',
                        'B1', 'B2', 'B3', 'B4', 'B5', 'B6',
                        'C1', 'C2', 'C3', 'C4', 'C5', 'C6',
                        'D1', 'D2', 'D3', 'D4', 'D5', 'D6']

    for i, mass in enumerate(sample_masses):
            aspiration_volume = 0.0
            water_volume = 0.0
            initial_volume = 0.0
            # Calculate water volume
            if mass > 0.01:
                num_samples += 1
                initial_volume = 20.0
                aspiration_volume = 10 / mass
                water_volume = 25.0 - (12.5 + 1.0 + 1.0 + aspiration_volume)
            dna_info[sample_tubes[i]] = {INITIAL_VOLUME_DICT_KEY: initial_volume,
                                                     USED_VOLUME_DICT_KEY: 0.0, 
                                                     ASPIRATION_VOLUME_DICT_KEY: aspiration_volume,
                                                     WATER_VOLUME_DICT_KEY: water_volume}

    # Calculate required reagent volumes using the helper function
    required_volumes = get_volumes_needed(num_samples)

    protocol.pause(f"\nPlease load the following to the tube rack:\n"
                   f"{required_volumes.nuclease_free_water}ul {nuclease_free_water.name} into the tube rack at position {NUCLEASE_FREE_LOC_CODE}\n"
                   f"{required_volumes.onetaq_master_mix}ul {onetaq_master_mix.name} into the tube rack at position {ONETAQ_LOC_CODE}\n"
                   f"{required_volumes.forward_primer}ul {forward_primer.name} into the tube rack at position {FORWARD_PRIMER_LOC_CODE}\n"
                   f"{required_volumes.reverse_primer}ul {reverse_primer.name} into the tube rack at position {REVERSE_PRIMER_LOC_CODE}\n")

    # Assign reagents to specific wells and load calculated volumes
    nuclease_free_water_wells = tube_rack_15_50ml.wells_by_name()['A1']
    nuclease_free_water_wells.load_liquid(liquid=nuclease_free_water, volume=required_volumes.nuclease_free_water)  # Adjust if necessary

    onetaq_master_mix_wells = tube_rack_15_50ml.wells_by_name()['A2']
    onetaq_master_mix_wells.load_liquid(liquid=onetaq_master_mix, volume=required_volumes.onetaq_master_mix)

    reverse_primer_wells = tube_rack_15_50ml.wells_by_name()['B1']
    reverse_primer_wells.load_liquid(liquid=reverse_primer, volume=required_volumes.reverse_primer)

    forward_primer_wells = tube_rack_15_50ml.wells_by_name()['B2']
    forward_primer_wells.load_liquid(liquid=forward_primer, volume=required_volumes.forward_primer)

    Template_DNA_well_A1 = tube_rack_2ml.wells_by_name()['A1']
    Template_DNA_well_A2 = tube_rack_2ml.wells_by_name()['A2']
    Template_DNA_well_A3 = tube_rack_2ml.wells_by_name()['A3']
    Template_DNA_well_A4 = tube_rack_2ml.wells_by_name()['A4']
    Template_DNA_well_A5 = tube_rack_2ml.wells_by_name()['A5']
    Template_DNA_well_A6 = tube_rack_2ml.wells_by_name()['A6']
    Template_DNA_well_B1 = tube_rack_2ml.wells_by_name()['B1']
    Template_DNA_well_B2 = tube_rack_2ml.wells_by_name()['B2']
    Template_DNA_well_B3 = tube_rack_2ml.wells_by_name()['B3']
    Template_DNA_well_B4 = tube_rack_2ml.wells_by_name()['B4']
    Template_DNA_well_B5 = tube_rack_2ml.wells_by_name()['B5']
    Template_DNA_well_B6 = tube_rack_2ml.wells_by_name()['B6']
    Template_DNA_well_C1 = tube_rack_2ml.wells_by_name()['C1']
    Template_DNA_well_C2 = tube_rack_2ml.wells_by_name()['C2']
    Template_DNA_well_C3 = tube_rack_2ml.wells_by_name()['C3']
    Template_DNA_well_C4 = tube_rack_2ml.wells_by_name()['C4']
    Template_DNA_well_C5 = tube_rack_2ml.wells_by_name()['C5']
    Template_DNA_well_C6 = tube_rack_2ml.wells_by_name()['C6']
    Template_DNA_well_D1 = tube_rack_2ml.wells_by_name()['D1']
    Template_DNA_well_D2 = tube_rack_2ml.wells_by_name()['D2']
    Template_DNA_well_D3 = tube_rack_2ml.wells_by_name()['D3']
    Template_DNA_well_D4 = tube_rack_2ml.wells_by_name()['D4']
    Template_DNA_well_D5 = tube_rack_2ml.wells_by_name()['D5']
    Template_DNA_well_D6 = tube_rack_2ml.wells_by_name()['D6']
    Template_DNA_well_A1.load_liquid(liquid=Template_DNA_liquid_A1, volume=dna_info[A1_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_A2.load_liquid(liquid=Template_DNA_liquid_A2, volume=dna_info[A2_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_A3.load_liquid(liquid=Template_DNA_liquid_A3, volume=dna_info[A3_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_A4.load_liquid(liquid=Template_DNA_liquid_A4, volume=dna_info[A4_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_A5.load_liquid(liquid=Template_DNA_liquid_A5, volume=dna_info[A5_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_A6.load_liquid(liquid=Template_DNA_liquid_A6, volume=dna_info[A6_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_B1.load_liquid(liquid=Template_DNA_liquid_B1, volume=dna_info[B1_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_B2.load_liquid(liquid=Template_DNA_liquid_B2, volume=dna_info[B2_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_B3.load_liquid(liquid=Template_DNA_liquid_B3, volume=dna_info[B3_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_B4.load_liquid(liquid=Template_DNA_liquid_B4, volume=dna_info[B4_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_B5.load_liquid(liquid=Template_DNA_liquid_B5, volume=dna_info[B5_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_B6.load_liquid(liquid=Template_DNA_liquid_B6, volume=dna_info[B6_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_C1.load_liquid(liquid=Template_DNA_liquid_C1, volume=dna_info[C1_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_C2.load_liquid(liquid=Template_DNA_liquid_C2, volume=dna_info[C2_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_C3.load_liquid(liquid=Template_DNA_liquid_C3, volume=dna_info[C3_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_C4.load_liquid(liquid=Template_DNA_liquid_C4, volume=dna_info[C4_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_C5.load_liquid(liquid=Template_DNA_liquid_C5, volume=dna_info[C5_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_C6.load_liquid(liquid=Template_DNA_liquid_C6, volume=dna_info[C6_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_D1.load_liquid(liquid=Template_DNA_liquid_D1, volume=dna_info[D1_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_D2.load_liquid(liquid=Template_DNA_liquid_D2, volume=dna_info[D2_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_D3.load_liquid(liquid=Template_DNA_liquid_D3, volume=dna_info[D3_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_D4.load_liquid(liquid=Template_DNA_liquid_D4, volume=dna_info[D4_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_D5.load_liquid(liquid=Template_DNA_liquid_D5, volume=dna_info[D5_LOC_CODE][INITIAL_VOLUME_DICT_KEY])
    Template_DNA_well_D6.load_liquid(liquid=Template_DNA_liquid_D6, volume=dna_info[D6_LOC_CODE][INITIAL_VOLUME_DICT_KEY])

    # Define destination wells on the PCR plate (from A10 to H12)
    all_destination_wells = [
        'A10', 'B10', 'C10', 'D10', 'E10', 'F10', 'G10', 'H10',
        'A11', 'B11', 'C11', 'D11', 'E11', 'F11', 'G11', 'H11',
        'A12', 'B12', 'C12', 'D12', 'E12', 'F12', 'G12', 'H12'
    ]

    # Re-adjust destination wells if any samples were skipped
    destination_wells = all_destination_wells[:num_samples]

    # Step 1: Add OneTaq 2X Master Mix (12.5 µL) to each PCR well
    p20_single.pick_up_tip(tip_rack_20ul.wells_by_name()['A1'])
    for well_name in destination_wells:
        dest_well = pcr_plate.wells_by_name()[well_name]
        p20_single.aspirate(12.5, onetaq_master_mix_wells.bottom(z=1))
        protocol.comment(f"Aspirated 12.5 µL of OneTaq Master Mix from {onetaq_master_mix_wells}.")
        p20_single.dispense(12.5, dest_well.bottom(z=1))
        protocol.comment(f"Dispensed 12.5 µL of OneTaq Master Mix into well {well_name}.")
        p20_single.touch_tip(dest_well)
        total_master_mix_used += 12.5
        # Initialize composition dict for this well
        well_comp = {
            'well_name': well_name,
            'master_mix': 12.5,
            'forward_primer': 0.0,
            'reverse_primer': 0.0,
            'dna_volume': 0.0,
            'water_volume': 0.0
        }
        pcr_well_compositions[well_name] = well_comp
    p20_single.drop_tip()

    # Step 2: Add Forward Primer (1 µL) to each PCR well
    p20_single.pick_up_tip(tip_rack_20ul.wells_by_name()['B1'])
    for well_name in destination_wells:
        dest_well = pcr_plate.wells_by_name()[well_name]
        p20_single.aspirate(1.0, forward_primer_wells.bottom(z=1))
        protocol.comment(f"Aspirated 1.0 µL of Forward Primer from {forward_primer_wells}.")
        p20_single.dispense(1.0, dest_well.bottom(z=1))
        protocol.comment(f"Dispensed 1.0 µL of Forward Primer into well {well_name}.")
        p20_single.touch_tip(dest_well)
        total_forward_primer_used += 1.0
        pcr_well_compositions[well_name]['forward_primer'] = 1.0
    p20_single.drop_tip()

    # Step 3: Add Reverse Primer (1 µL) to each PCR well
    p20_single.pick_up_tip(tip_rack_20ul.wells_by_name()['C1'])
    for well_name in destination_wells:
        dest_well = pcr_plate.wells_by_name()[well_name]
        p20_single.aspirate(1.0, reverse_primer_wells.bottom(z=1))
        protocol.comment(f"Aspirated 1.0 µL of Reverse Primer from {reverse_primer_wells}.")
        p20_single.dispense(1.0, dest_well.bottom(z=1))
        protocol.comment(f"Dispensed 1.0 µL of Reverse Primer into well {well_name}.")
        p20_single.touch_tip(dest_well)
        total_reverse_primer_used += 1.0
        pcr_well_compositions[well_name]['reverse_primer'] = 1.0
    p20_single.drop_tip()

    # Step 4: Transfer Variable Volumes of Template DNA from tube rack to PCR plate
    for i, (sample_tube_name, dest_well_name) in enumerate(zip(sample_tubes, destination_wells)):
        dest_well = pcr_plate.wells_by_name()[dest_well_name]
        source_well = tube_rack_2ml.wells_by_name()[sample_tube_name]
        vol_dna = dna_info[sample_tube_name][ASPIRATION_VOLUME_DICT_KEY]
        p20_single.pick_up_tip()
        p20_single.aspirate(vol_dna, source_well.bottom(z=1))
        protocol.comment(f"Aspirated {vol_dna} µL of Template DNA from {source_well}.")
        p20_single.dispense(vol_dna, dest_well.bottom(z=1))
        protocol.comment(f"Dispensed {vol_dna} µL of Template DNA into well {dest_well_name}.")
        p20_single.touch_tip(dest_well)
        p20_single.drop_tip()
        # Update volumes
        total_dna_used_per_sample.append(vol_dna)
        pcr_well_compositions[dest_well_name]['dna_volume'] = vol_dna
        dna_info[sample_tube_name][USED_VOLUME_DICT_KEY] += vol_dna

    # Step 5: Add Nuclease-Free Water to each PCR well to bring volume to 25 µL
    for i, (sample_tube_name, dest_well_name) in enumerate(zip(sample_tubes, destination_wells)):
        dest_well = pcr_plate.wells_by_name()[dest_well_name]
        vol_water = dna_info[sample_tube_name][WATER_VOLUME_DICT_KEY]
        if vol_water > 0:
            p20_single.pick_up_tip()
            p20_single.aspirate(vol_water, nuclease_free_water_wells.bottom(z=1))
            protocol.comment(f"Aspirated {vol_water} µL of Nuclease-Free Water from {nuclease_free_water_wells}.")
            p20_single.dispense(vol_water, dest_well.bottom(z=1))
            protocol.comment(f"Dispensed {vol_water} µL of Nuclease-Free Water into well {dest_well_name}.")
            p20_single.touch_tip(dest_well)
            p20_single.drop_tip()
            total_water_used += vol_water
            pcr_well_compositions[dest_well_name]['water_volume'] = vol_water

    # Step 6: Mix the contents of each well using the multichannel pipette
    # Group destination wells by column for multichannel pipetting
    columns_to_mix = set()
    for well_name in destination_wells:
        column = well_name[1:]  # Extract column number from well name
        columns_to_mix.add(column)

    for column in columns_to_mix:
        dest_wells = pcr_plate.columns_by_name()[column]
        p20_multi.pick_up_tip()
        p20_multi.mix(3, 10, dest_wells[0].bottom(z=1))  # Mix 3 times with 10 µL
        protocol.comment(f"Mixed contents of column {column} using multichannel pipette.")
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
    onetaq_master_mix_initial = required_volumes.onetaq_master_mix
    onetaq_master_mix_final = onetaq_master_mix_initial - total_master_mix_used
    protocol.comment(f"OneTaq Master Mix: Final Volume = {onetaq_master_mix_final} µL")

    forward_primer_initial = required_volumes.forward_primer
    forward_primer_final = forward_primer_initial - total_forward_primer_used
    protocol.comment(f"Forward Primer: Final Volume = {forward_primer_final} µL")

    reverse_primer_initial = required_volumes.reverse_primer
    reverse_primer_final = reverse_primer_initial - total_reverse_primer_used
    protocol.comment(f"Reverse Primer: Final Volume = {reverse_primer_final} µL")

    nuclease_free_water_final = required_volumes.nuclease_free_water - total_water_used
    protocol.comment(f"Nuclease-Free Water: Final Volume = {nuclease_free_water_final} µL")

    for tube_location, volumes in dna_info.items():
        initial_volume = volumes[INITIAL_VOLUME_DICT_KEY]
        used_volume = volumes[USED_VOLUME_DICT_KEY]
        final_volume = initial_volume - used_volume
        protocol.comment(f"Template DNA Tube {tube_location}: Final Volume = {final_volume} µL")

    # End of Protocol
