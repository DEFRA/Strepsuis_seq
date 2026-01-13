#! /usr/bin/env python3

"""
Combine outputs from all four srst2 runs to generate a single results file in csv format,
together with kraken outputs and read stats. 
"""


import pandas as pd
import argparse
import re
from datetime import date


def combineData(recNTable, MLSTTable, serotypeTable, virulenceTable, verifyCSV, krakenID, statsCSV, dataDir):

    # Get info for logging
    date_out = date.today().strftime('%d%b%y')

    # read recN data
    recN_df = pd.read_table(recNTable, sep='\t', names=list(range(2)), skiprows=1)
    recN_df.rename({0: 'Sample'}, axis=1, inplace=True)
    recN_df.rename({1: 'recN-Pos'}, axis=1, inplace=True)
    recN_df['Sample'] = recN_df['Sample'].astype(object)
    recN_df['recN-Pos'].fillna('noRecN', inplace=True)

    # read MLST data
    MLST_df = pd.read_table(MLSTTable, sep='\t')
    MLST_df[['Sample', 'ST', 'aroA', 'cpn60', 'dpr', 'gki', 'mutS', 'recA', 'thrA']] = \
        MLST_df[['Sample', 'ST', 'aroA', 'cpn60', 'dpr', 'gki', 'mutS', 'recA', 'thrA']].astype(object)

    # read serotype data
    serotype_df = pd.read_table(serotypeTable, sep='\t')
    serotype_df['Sample'] = serotype_df['Sample'].astype(object)
    serotype_df.rename({'ST': 'Serotype'}, axis=1, inplace=True)
    serotype_df['Serotype'] = serotype_df['Serotype'].astype(str)
    serotype_df['Serotype'] = serotype_df['Serotype'].str.replace(r'\D', '', regex=True)
    serotype_df['mismatches'] = serotype_df['mismatches'].astype(str)
    
    # Check the mismatch data and remove serotype assignment if snps, indels or holes above thresholds
    _st = []
    
    # regex for extracting nuber of snps, indels and holes from each sample
    for index, row in serotype_df.iterrows():
        mismatch_data = row['mismatches']
        diffs = re.split(r'/', mismatch_data)
        matches = re.search(r"(\d+snp)?(\d+indel)?(\d+holes)?", diffs[-1]).groups()
        if matches[0]:
            snp = int(re.split(r'snp?', matches[0])[0])
        else: snp = 0
        if matches[1]:
            indel = int(re.split(r'indel?', matches[1])[0])
        else: indel = 0
        if matches[2]:
            hole = int(re.split(r'holes?', matches[2])[0])
        else: hole = 0

    # Set anomoly thresholds based on depth
        depth_data = row['depth']
        if depth_data >= 20:
            max_snp = 20
            max_indel = 5
            max_hole = 10
        else:
            max_snp = 30
            max_indel = 10
            max_hole = 20

    # flag serotype assignment if anomolies exceed thresholds    
        if snp > max_snp: row['Serotype'] += ' (LowMatch)'
        if indel > max_indel: row['Serotype'] += ' (LowMatch)'
        if hole > max_hole: row['Serotype'] += ' (LowMatch)'
        r = row.to_dict()
        _st.append(r)

    parsedST_df = pd.DataFrame(_st)

    # read serotype verification
    verify_df = pd.read_csv(verifyCSV)
    verify_df[['Sample', 'Pos483', 'RefF', 'RefR', 'AltF', 'AltR']] = \
        verify_df[['Sample', 'Pos483', 'RefF', 'RefR', 'AltF', 'AltR']].astype(object)

    # update serotypes based on verification
    verified_serotype_df = pd.merge(parsedST_df, verify_df, on='Sample', how='left')
    # Quote from https://rdcu.be/d56Ee:
    # ## "The analysis revealed that all serotype 2 and all serotype 14 strains had a G nucleotide at position 483 of
    # ## the cpsK gene, while all serotype 1 and all serotype 1/2 strains contained either a C or T at that nucleotide
    # ## position."
    # if prelim serotype is 1 and base 483 is G, then Serotype is 14
    # if prelim serotype is 1 and base 483 is C or T, then Serotype is 1
    # if prelim serotype is 2 and base 483 is G, then Serotype is 2
    # if prelim serotype is 2 and base 483 is C or T, then Serotype is 1/2

    verified_serotype_df.loc[(verified_serotype_df['Serotype'] == '1') & (verified_serotype_df['Pos483'] == 'G'),
                             'Serotype'] = '14'
    verified_serotype_df.loc[(verified_serotype_df['Serotype'] == '1') & (verified_serotype_df['Pos483'] == 'C'),
                             'Serotype'] = '1'
    verified_serotype_df.loc[(verified_serotype_df['Serotype'] == '1') & (verified_serotype_df['Pos483'] == 'T'),
                             'Serotype'] = '1'
    verified_serotype_df.loc[(verified_serotype_df['Serotype'] == '2') & (verified_serotype_df['Pos483'] == 'G'),
                             'Serotype'] = '2'
    verified_serotype_df.loc[(verified_serotype_df['Serotype'] == '2') & (verified_serotype_df['Pos483'] == 'C'),
                             'Serotype'] = "'1/2"
    verified_serotype_df.loc[(verified_serotype_df['Serotype'] == '2') & (verified_serotype_df['Pos483'] == 'T'),
                             'Serotype'] = "'1/2"
    verified_serotype_df.loc[(verified_serotype_df['Serotype'] != '1') & (verified_serotype_df['Serotype'] != '2'),
                             'Serotype'] = verified_serotype_df['Serotype']

    # read virulence data
    virulence_df = pd.read_table(virulenceTable, sep='\t', names=list(range(4)), skiprows=1)
    virulence_df.rename({0: 'Sample'}, axis=1, inplace=True)
    virulence_df.fillna('-', inplace=True)
    virulence_df['epf'] = virulence_df[1].str.contains('epf') | (virulence_df[2].str.contains('epf')) | \
        (virulence_df[3].str.contains('epf')) | (virulence_df[4].str.contains('epf')) | \
        (virulence_df[5].str.contains('epf'))
    virulence_df['mrp'] = (virulence_df[1].str.contains('mrp')) | (virulence_df[2].str.contains('mrp')) | \
        (virulence_df[3].str.contains('mrp')) | (virulence_df[4].str.contains('mrp')) | \
        (virulence_df[5].str.contains('mrp'))
    virulence_df['sly'] = (virulence_df[1].str.contains('sly')) | (virulence_df[2].str.contains('sly')) | \
        (virulence_df[3].str.contains('sly')) | (virulence_df[4].str.contains('sly')) | \
        (virulence_df[5].str.contains('sly'))
    virulence_df['virA'] = (virulence_df[1].str.contains('virA')) | (virulence_df[2].str.contains('virA')) | \
        (virulence_df[3].str.contains('virA')) | (virulence_df[4].str.contains('virA')) \
        (virulence_df[5].str.contains('virA'))
    virulence_df['GI3'] = (virulence_df[1].str.contains('GI3')) | (virulence_df[2].str.contains('GI3')) | \
        (virulence_df[3].str.contains('GI3')) | (virulence_df[4].str.contains('GI3')) | \
        (virulence_df[5].str.contains('GI3'))
    virulence_df = virulence_df[['Sample', 'epf', 'mrp', 'sly', 'virA', 'GI3']]
    virulence_df = virulence_df.astype(str)
    virulence_df.replace('False', '-', inplace=True)

    # read krakenID data
    krakenID_df = pd.read_csv(krakenID)
    krakenID_df['Sample'] = krakenID_df['Sample'].astype(object)

    # read statsCSV
    readstats_df = pd.read_csv(statsCSV)
    readstats_df['Sample'] = readstats_df['Sample'].astype(object)

    # Merge dataframes
    sero_mlst_df = pd.merge(verified_serotype_df, MLST_df, on='Sample', how='outer')
    # Remove serotype assigment if coverage is anomolous compared to MLST loci
    sero_mlst_df.loc[(sero_mlst_df['depth_x'].div(sero_mlst_df['depth_y']) <= 0.25),
                     'Serotype'] += '?'
    sero_mlst_df.loc[(sero_mlst_df['depth_x'].div(sero_mlst_df['depth_y']) <= 0.12),
                     'Serotype'] = 'NoMatch'

    id_df = pd.merge(pd.merge(readstats_df, krakenID_df, on='Sample', how='outer'), 
                     recN_df, on='Sample', how='outer')

    finalout_df = pd.merge(pd.merge(id_df, sero_mlst_df, on='Sample', how='outer'),
                           virulence_df, on='Sample', how='outer')
    finalout_df.set_index('Sample', inplace=True)

    # Add additional infomation, missing values, flags for lower quality data and sort
    finalout_df['Serotype'].fillna('Undetermined', inplace=True)
    finalout_df.loc[(finalout_df['krakenID'] == 'Streptococcus suis') & (finalout_df['Flag'] == 'LowData') 
                    & (finalout_df['Serotype'] == 'Undetermined'), 'Result'] = 'RepeatRequired'
    finalout_df.loc[(finalout_df['Flag'] != 'LowData') & (finalout_df['krakenID'] != 'Streptococcus suis'),
                    'Result'] = 'Not Strep. suis'
    finalout_df.loc[(finalout_df['krakenID'] == 'Streptococcus suis') & (finalout_df['Abundance'] < 80 ),
                    'Result'] = 'Other Strep. spp.'
    finalout_df.loc[(finalout_df['krakenID'] == 'Streptococcus suis') & (finalout_df['recN-Pos'] == 'noRecN' ),
                    'Result'] = 'Possible other Strep. spp.'
    finalout_df.loc[(finalout_df['krakenID'] == 'Streptococcus suis') & (finalout_df['recN-Pos'] == 'recN-Pos'),
                    'Result'] = 'Strep. suis'
    finalout_df.loc[(finalout_df['krakenID'] == 'Streptococcus suis') & (finalout_df['recN-Pos'] == 'recN-Pos*'),
                    'Result'] = 'Strep. suis'
    finalout_df.loc[(finalout_df['krakenID'] == 'Streptococcus suis') & (finalout_df['recN-Pos'] == 'recN-Pos?'),
                    'Result'] = 'Strep. suis'
    finalout_df.loc[(finalout_df['krakenID'] == 'Streptococcus suis') & (finalout_df['recN-Pos'] == 'recN-Pos*?'),
                    'Result'] = 'Strep. suis'
    
    # TODO - ensure other combinations are covered..... 

    finalout_df['Result'].replace('nan', '', inplace=True)
    rescol = finalout_df.pop('Result')
    finalout_df.insert(0, 'Result', rescol)
    serocol = finalout_df.pop('Serotype')
    finalout_df.insert(1, 'Serotype', serocol)
    finalout_df.sort_index(inplace=True)

    # Write to csv
    finalout_df.to_csv("{}_Ssuis-FinalOut_{}.csv".format(dataDir, date_out), float_format='%.2f')

    # Append log info
    # with open("{}_FinalOut_{}.csv".format(date_out), "a") as outFile:
    #    outFile.close


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('recNTable', help='............')
    parser.add_argument('MLSTTable', help='............')
    parser.add_argument('serotypeTable', help='...............')
    parser.add_argument('virulenceTable', help='................')
    parser.add_argument('verifyCSV', help='................')
    parser.add_argument('krakenID', help='..............')
    parser.add_argument('statsCSV', help='..............')
    parser.add_argument('dataDir', help='.............')
    # parser.add_argument('commitId', help='Nextflow capture of git commit')

    args = parser.parse_args()

    combineData(**vars(args))
