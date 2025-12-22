import argparse
import os
import sys
import time
import logging
import colorlog
from datetime import datetime

VERSION="0.0.0"

## Set up logger
logger = logging.getLogger('mofongo')

## format for cli
colored_format = (
    "%(green)s%(asctime)s - %(name)s - "
    "%(log_color)s%(levelname)-8s%(reset)s "
    "%(white)s%(message)s"
)

color_formatter = colorlog.ColoredFormatter(
    colored_format,
    datefmt='%H:%M:%S',
    log_colors={
        'DEBUG':    'bold_cyan', 
        'INFO':     'bold_blue', 
        'WARNING':  'bold_yellow', 
        'ERROR':    'bold_red', 
        'CRITICAL': 'bold_red,bg_white',
    }
)

## different format for text file log
normal_formater = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

def exit_mofongo_help():
    logger.error('Please use --help to see options')
    sys.exit()

def parse_args():
    """
    Main argument parsing function and monfongo runtime
    """
    parser = argparse.ArgumentParser(description=f"mofongo v{VERSION} -- Differentiate Between Engineered and Natural Plasmids Using a Linear Pangenome")
    
    other_group = parser.add_argument_group('Other Options', 'Other settings when running mofongo')
    other_group.add_argument('--debug', action='store_true', default=False, help='Print out more debug information')
    other_group.add_argument('-t', '--threads', type=int, default=3, help='Number of threads to use')
    
    input_group = parser.add_argument_group('Data Input Options', 'Options for inputs into mofongo')
    input_group.add_argument('-r', nargs='+', metavar='FILE', help="Path to one or more single-end *.fastq(.gz) files")
    input_group.add_argument('-1', '--first-paired', nargs='+', metavar='FILE', help="Path to first read of paired-end *.fastq(.gz) files (must be with -2)")
    input_group.add_argument('-2', '--second-paired', nargs='+', metavar='FILE', help="Path to second read or paired-end *.fastq(.gz) files (must be with -1)")
    input_group.add_argument('-c', '--contigs', nargs='+', metavar='FILE', help='Path to contigs or genomes *.fasta(.gz)')
    
    database_group = parser.add_argument_group('Database options', 'Options for pangenome used in mofongo')
    database_group.add_argument('-p', '--pangenome', metavar='FILE', required=True, help='Path to fasta pangenome file')
    
    output_group = parser.add_argument_group('Output options', "Output options for mofongo")
    output_group.add_argument('-o', '--output', default='mofongo_output', type=str, help='Working directory for mofongo (default=mofongo_output)')
    
    algorithm_group = parser.add_argument_group('Algorithm options', "Algorithm options for mofongo")
    
    args = parser.parse_args()
    out = args.output
    
    r = args.r 
    fwd_reads = args.first_paired
    rev_reads = args.second_paired
    contigs = args.contigs
    
    ## make sure sequence data of some form exists
    if r == None and fwd_reads == None and rev_reads == None and contigs == None:
        logger.error("No sequence data provided, must provide either reads or contigs")
        exit_mofongo_help()
    
    if (fwd_reads != None and rev_reads != None) and len(fwd_reads) != len(rev_reads):
        logger.error("Number of forward reads (-1 argument) is not equal to number of reverse reads (-2 argument)")
        exit_mofongo_help()
        
        
    ## setting up logger 
    ## file log output
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"mofongo_{timestamp}.log"
    fh = logging.FileHandler(log_filename, mode='w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(normal_formater)

    ## cli log
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(color_formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)
    logger.setLevel(logging.DEBUG)
    
    if os.path.exists(out):
        logger.error(f'Output directory `{out}` already exists, please either delete it or use another location')
        exit_mofongo_help()
    else:
        os.mkdir(out)
    
    debug = args.debug
    if debug:
        ch.setLevel(logging.DEBUG)
    
    logger.info(f'Running mofongo v{VERSION}')    
    # run_mofongo(args)
    logger.info('')
    logger.info(f'mofongo v{VERSION} completed successfully.')
    
    
if __name__=="__main__":
    parse_args()