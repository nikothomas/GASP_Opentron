# Import the required libraries from the Opentrons API
import csv
import os
from opentrons import protocol_api

# Metadata
metadata = {
    'protocolName': 'GASP Protocol with Variable DNA Volumes from CSV',
    'author': 'Nikolas Yanek-Chrones <research@icarai.io>',
    'description': 'Automated setup for GASP PCR reactions using Opentron OT-2 with variable DNA volumes from CSV input',
    'apiLevel': '2.20'  # Use the latest API level supported by your OT-2
}

TESTING = 1  # Set to 0 when running on the actual robot
CSV_FILE_PATH = 'Opentron_input.csv'  # Update this path to a testing CSV file
FORWARD_PRIMER_VOLUME_PER_SAMPLE = 1  # 1 µL of forward primer per sample for 25 µL reaction
REVERSE_PRIMER_VOLUME_PER_SAMPLE = 1  # 1 µL of reverse primer per sample for 25 µL reaction
REQUIRED_ONETAQ_VOLUME_PER_SAMPLE = 12.5  # 12.5 µL of OneTaq per sample for 25 µL reaction

def get_volumes_needed(num_samples: int):
    # Calculate required volumes with a small excess (+1 sample) for pipetting errors
    required_forward_primer_vol = (num_samples + 1) * FORWARD_PRIMER_VOLUME_PER_SAMPLE
    required_reverse_primer_vol = (num_samples + 1) * REVERSE_PRIMER_VOLUME_PER_SAMPLE
    required_onetaq_vol = (num_samples + 1) * REQUIRED_ONETAQ_VOLUME_PER_SAMPLE
    return required_forward_primer_vol, required_reverse_primer_vol, required_onetaq_vol

def add_parameters(parameters):
    parameters.add_csv_file(
        variable_name="dna_masses",
        display_name="DNA Masses and Locations",
        description=(
            "Table with two columns:"
            " source slot, dna mass,"
        )
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
        name="OneTaq Master Mix", description="", display_color="#E1FF05"
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
    temp_module = protocol.load_module(
        module_name='temperature module gen2',
        location=10
    )

    # Initialize tracking variables
    pcr_well_compositions = {}
    total_master_mix_used = 0.0
    total_forward_primer_used = 0.0
    total_reverse_primer_used = 0.0
    total_water_used = 0.0
    total_dna_used_per_sample = []
    # For DNA tubes, create a dict to track volumes used and final volumes
    dna_tube_volumes = {}  # key: tube_location, value: {'initial_volume':20.0, 'used_volume':0.0}

    # Read the CSV data
    if TESTING == 1:
        if not os.path.exists(CSV_FILE_PATH):
            protocol.comment(f"CSV file not found at {CSV_FILE_PATH}.")
            return
        else:
            with open(CSV_FILE_PATH, 'r') as csvfile:
                csv_reader = csv.reader(csvfile)
                csv_data = list(csv_reader)
    else:
        csv_data = protocol.params.dna_masses.parse_as_csv()

    # Parse the CSV data
    sample_data = []
    for row in csv_data:
        if len(row) >= 2:
            tube_location = row[0].strip()
            try:
                dna_mass = float(row[1].strip())
                sample_data.append({'tube_location': tube_location, 'dna_mass': dna_mass})
            except ValueError:
                protocol.comment(f"Invalid DNA mass '{row[1]}' at tube location '{tube_location}'. Skipping.")
        else:
            protocol.comment(f"Incomplete row '{row}'. Skipping.")

    # Proceed only if there are samples
    if not sample_data:
        protocol.comment("No valid sample data provided in CSV.")
        return

    # Determine the number of samples
    num_samples = len(sample_data)

    # Calculate required reagent volumes using the helper function
    (
        required_forward_primer_vol,
        required_reverse_primer_vol,
        required_onetaq_vol
    ) = get_volumes_needed(num_samples)

    # Assign reagents to specific wells and load calculated volumes
    nuclease_free_water_wells = tube_rack_15_50ml.wells_by_name()['A1']
    # For nuclease-free water, estimate required volume
    total_water_needed = 0  # Will be updated after calculating water volumes
    nuclease_free_water_wells.load_liquid(liquid=nuclease_free_water, volume=1000)  # Adjust if necessary

    onetaq_master_mix_wells = tube_rack_15_50ml.wells_by_name()['A2']
    onetaq_master_mix_wells.load_liquid(liquid=onetaq_master_mix, volume=required_onetaq_vol)

    reverse_primer_wells = tube_rack_15_50ml.wells_by_name()['B1']
    reverse_primer_wells.load_liquid(liquid=reverse_primer, volume=required_reverse_primer_vol)

    forward_primer_wells = tube_rack_15_50ml.wells_by_name()['B2']
    forward_primer_wells.load_liquid(liquid=forward_primer, volume=required_forward_primer_vol)

    # Template DNA wells are in the 24-well plate (tube_rack_2ml)

    # Map tube locations to wells in tube_rack_2ml
    tube_rack_positions = [well.well_name for well in tube_rack_2ml.wells()]
    source_wells = []
    valid_samples = []

    for sample in sample_data:
        tube_location = sample['tube_location']
        if tube_location in tube_rack_positions:
            source_well = tube_rack_2ml.wells_by_name()[tube_location]
            source_wells.append(source_well)
            valid_samples.append(sample)
            # Assign the template_dna liquid to the source_well
            source_well.load_liquid(liquid=template_dna, volume=20)
            # Initialize DNA tube volumes tracking
            dna_tube_volumes[tube_location] = {'initial_volume': 20.0, 'used_volume': 0.0}
        else:
            protocol.comment(f"Tube location '{tube_location}' not found in tube rack. Skipping sample.")

    # Update sample_data to only include valid samples
    sample_data = valid_samples
    num_samples = len(sample_data)  # Recalculate in case some samples were skipped

    # Recalculate required reagent volumes if the number of valid samples has changed
    (
        required_forward_primer_vol,
        required_reverse_primer_vol,
        required_onetaq_vol
    ) = get_volumes_needed(num_samples)

    # Adjust loaded volumes if necessary
    onetaq_master_mix_wells.load_liquid(liquid=onetaq_master_mix, volume=required_onetaq_vol)
    reverse_primer_wells.load_liquid(liquid=reverse_primer, volume=required_reverse_primer_vol)
    forward_primer_wells.load_liquid(liquid=forward_primer, volume=required_forward_primer_vol)

    # Define destination wells on the PCR plate (from A10 to H12)
    all_destination_wells = [
        'A10', 'B10', 'C10', 'D10', 'E10', 'F10', 'G10', 'H10',
        'A11', 'B11', 'C11', 'D11', 'E11', 'F11', 'G11', 'H11',
        'A12', 'B12', 'C12', 'D12', 'E12', 'F12', 'G12', 'H12'
    ]

    # Re-adjust destination wells if any samples were skipped
    destination_wells = all_destination_wells[:num_samples]

    # Calculate DNA volumes and water volumes
    dna_volumes = []
    water_volumes = []
    for sample in sample_data:
        mass = sample['dna_mass']
        vol_dna = mass * 0.1  # in µL
        # Ensure DNA volume is within pipette limits
        if vol_dna < 1.0:
            vol_dna = 1.0
        elif vol_dna > 20.0:
            vol_dna = 20.0
            protocol.comment(f"DNA volume for mass {mass} exceeds 20 µL. Set to 20 µL.")
        dna_volumes.append(vol_dna)
        # Calculate water volume
        vol_water = 25.0 - (12.5 + 1.0 + 1.0 + vol_dna)
        if vol_water < 0:
            vol_water = 0.0
            protocol.comment(f"Total volume exceeds 25 µL for mass {mass}. No water added.")
        water_volumes.append(vol_water)
        total_water_needed += vol_water  # Update total water needed

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
    for i, (source_well, dest_well_name) in enumerate(zip(source_wells, destination_wells)):
        dest_well = pcr_plate.wells_by_name()[dest_well_name]
        vol_dna = dna_volumes[i]
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
        # Update dna_tube_volumes
        tube_location = sample_data[i]['tube_location']
        dna_tube_volumes[tube_location]['used_volume'] += vol_dna

    # Step 5: Add Nuclease-Free Water to each PCR well to bring volume to 25 µL
    for i, well_name in enumerate(destination_wells):
        dest_well = pcr_plate.wells_by_name()[well_name]
        vol_water = water_volumes[i]
        if vol_water > 0:
            p20_single.pick_up_tip()
            p20_single.aspirate(vol_water, nuclease_free_water_wells.bottom(z=1))
            protocol.comment(f"Aspirated {vol_water} µL of Nuclease-Free Water from {nuclease_free_water_wells}.")
            p20_single.dispense(vol_water, dest_well.bottom(z=1))
            protocol.comment(f"Dispensed {vol_water} µL of Nuclease-Free Water into well {well_name}.")
            p20_single.touch_tip(dest_well)
            p20_single.drop_tip()
            total_water_used += vol_water
            pcr_well_compositions[well_name]['water_volume'] = vol_water

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
    onetaq_master_mix_initial = required_onetaq_vol
    onetaq_master_mix_final = onetaq_master_mix_initial - total_master_mix_used
    protocol.comment(f"OneTaq Master Mix: Final Volume = {onetaq_master_mix_final} µL")

    forward_primer_initial = required_forward_primer_vol
    forward_primer_final = forward_primer_initial - total_forward_primer_used
    protocol.comment(f"Forward Primer: Final Volume = {forward_primer_final} µL")

    reverse_primer_initial = required_reverse_primer_vol
    reverse_primer_final = reverse_primer_initial - total_reverse_primer_used
    protocol.comment(f"Reverse Primer: Final Volume = {reverse_primer_final} µL")

    nuclease_free_water_final = 1000 - total_water_used
    protocol.comment(f"Nuclease-Free Water: Final Volume = {nuclease_free_water_final} µL")

    for tube_location, volumes in dna_tube_volumes.items():
        initial_volume = volumes['initial_volume']
        used_volume = volumes['used_volume']
        final_volume = initial_volume - used_volume
        protocol.comment(f"Template DNA Tube {tube_location}: Final Volume = {final_volume} µL")

    # End of Protocol
