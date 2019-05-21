# FastaPatternFinder

## Scan and find tri peptide patterns (PXX/PXP/XXP) from Fasta protein sequences. FastPatternFinder is written in python3 and requires python3 to run.

### prerequisites

     Modules required
        re
        Pandas
        BioPython
        psutil
        argparse
        
### Usage

#### python FastaPatternFinder.py -i <input_filename> -o <output_filename>
Defining output file is an option if the output file name is not provided the output will be written to output.csv in the working directory.

#### python FastaPatternFinder.py -help # for help
