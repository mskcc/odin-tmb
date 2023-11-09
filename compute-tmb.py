#!/usr/bin/env python3
import json
import argparse
import sys
import pandas as pd

"""
Script to calculated tmb, intended to take analysis.maf
 usage:
     compute-tmb inputMaf outputFile tumorId assay normalType


"""

# def get_number_of_on_target_mutations(maf_file,tumorId,assayDb):




def main():
    """
    Script to calculate TMB tumor mutational burden value
    """
    parser = argparse.ArgumentParser(description='Script to calculate TMB tumor mutational burden value')
    parser.add_argument("--maf_file",         required=True, dest = 'maf_file',         help="Input maf file")
    parser.add_argument("--output_filename",  required=False,dest = 'output_filename',  help="output_filename")
    parser.add_argument("--tumorId",          required=True, dest = 'tumorId',          help="tumorId")
    parser.add_argument("--assay",            required=True, dest = 'assay',            help="assay should be one of these IMPACT341, IMPACT410, IMPACT468, IMPACT505")
    # parser.add_argument("--normalType",       required=True, dest = 'normalType',       help="normalType")
    parser.add_argument("--assayDb_file",       required=False, dest = 'assayDb_file', default = 'data/tmbAssayDb.json', help="assayDb_file")
    
    args = parser.parse_args()

    try:
        with open(args.assayDb_file, 'r') as f:
            assayDb = json.load(f)
    except:
        print('Error: This assayDb_file file can not be found:',args.assayDb_file)
        sys.exit(1)

    if not args.assay in assayDb:
        print('Assay should be one of these: IMPACT341, IMPACT410, IMPACT468, IMPACT505')
        sys.exit(1)

    maf_df = pd.read_csv(args.maf_file, sep='\t',comment='#',usecols=['Tumor_Sample_Barcode', 'Hugo_Symbol'])

    filtered_maf_df=maf_df[maf_df['Hugo_Symbol'].isin(assayDb[args.assay]['genes'])] #only filtering out genes that are not included in the panel
    
    number_of_mutataions=maf_df.shape[0]
    number_of_on_target_mutataions=filtered_maf_df.shape[0]
    number_of_discarded=number_of_mutataions-number_of_on_target_mutataions
    tmb_val=round(number_of_on_target_mutataions/assayDb[args.assay]['genomicSize']*1000000,2)
    
    print('Total number of mutations for this sample:',number_of_mutataions)
    print('Number of on-target mutations:',number_of_on_target_mutataions)
    print('Number of off-target mutations:',number_of_discarded)
    print('TMB values is calculated by Number of on-target mutations/Assay length*1,000,000')
    print('TMB value:',tmb_val)

    if args.output_filename is not None:
        output_message = \
        'Total number of mutations for this sample: '+str(number_of_mutataions)+'\n'+ \
        'Number of on-target mutations: '+str(number_of_on_target_mutataions)+'\n'+ \
        'Number of off-target mutations: '+str(number_of_discarded)+'\n'+ \
        'TMB values is calculated by Number of on-target mutations/Assay length*1,000,000\n'+ \
        'TMB value: '+str(tmb_val)+'\n'
        
        open(args.output_filename,'w').write(output_message)


if __name__ == '__main__':
    main()