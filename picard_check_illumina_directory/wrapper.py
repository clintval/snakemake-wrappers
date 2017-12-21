"""Snakemake wrapper for picard CheckIlluminaDirectory."""

__author__ = "clintval"
__copyright__ = "Copyright 2018, Clint Valentine"
__email__ = "valentine.clint@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")
read_structure - snakemake.params.get('read_structure')
data_types = snakemake.params.get('data_types', 'null')
lanes_to_check = snakemake.params.get('lanes_to_check', 'null')
tile_numbers = snakemake.params.get('tile_numbers', 'null')
fake_files = snakemake.params.get('fake_files', 'false')
link_locs = snakemake.params.get('link_locs', 'false')

command = (
    "picard CheckIlluminaDirectory"
    " {extra}"
    " BASECALLS_DIR={snakemake.input.basecalls_dir}"
    " READ_STRUCTURE={read_structure}"
    " LINK_LOCS={link_locs}"
    " FAKE_FILES={fake_files}")

if data_types != 'null':
    assert isinstance(data_types, (list, tuple)), '`data_types` must be a sequence'
    data_type_params = [" DATA_TYPES=" + str(data_type) for data_type in data_types]
    command += data_type_params
if lanes_to_check != 'null':
    assert isinstance(lanes_to_check, (list, tuple)), '`lanes_to_check` must be a sequence'
    lanes_to_check_params = [" LANES=" + str(lanes) for lanes in lanes_to_check]
    command += lanes_to_check_params
if tile_numbers != 'null':
    assert isinstance(tile_numbers, (list, tuple)), '`tile_numbers` must be a sequence'
    tile_numbers_params = [" TILE_NUMBERS=" + str(tile_number) for tile_number in tile_numbers]
    command += tile_numbers_params
command += " {log}"

shell(command)
