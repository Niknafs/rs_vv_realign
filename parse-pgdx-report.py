#! /usr/bin/env python

# Usage: python parse-pgdx-report.py -i my-report.xlsx -o my-report-summary.tsv

import sys
import os
import xlrd
from argparse import ArgumentParser as argument_parser

#-----------------------------------------------------------------------# 
class InputParser(object):
    def __init__(self):
        parser = argument_parser()
        parser.add_argument('-i','--input', help = 'input path')
        parser.add_argument('-o','--output', help = 'output path')
        args = parser.parse_args()
        self.__dict__.update(args.__dict__)
#-----------------------------------------------------------------------# 
def encode(x):
    try:
        return str(x)
    except:
        try:
            return x.encode('utf-8')
        except:
            print >>sys.stderr, 'Unable to convert to string. Operation aborted.'
            print >>sys.stderr, x
            sys.exit()

            
#-----------------------------------------------------------------------# 
def read_report_data(path):
    book = xlrd.open_workbook(path)
    overview_sheet = book.sheet_by_index(0)
    summary_sheet = book.sheet_by_index(1)

    data = []
    # read in PGDX and lab ID of the sample
    id_string = map(encode, overview_sheet.col_values(0))[4]
    startInd = id_string.index('PGDX')
    id_string = id_string[startInd:]
    if '-' in id_string:
        pgdx_id, lab_id = id_string.split('-')
    else:
        pgdx_id, lab_id = id_string.split('_')
    data = data + [('Case ID', id_string), \
                       ('PGDX ID',pgdx_id), \
                       ('Lab ID', lab_id)]

    first_col =  map(encode,overview_sheet.col_values(0))
    last_index = first_col.index('Limitations of Approach')
    first_index = first_col.index('Sample Characteristics')
    
    # overview sheet attributes
    overview_expected_attributes = ['Tumor Type','Tumor Location','Sample Type',
                                    'Pathological Tumor Purity','Mutation based Tumor Purity',
                                    'Source of normal DNA', '' ,'Analysis Characteristics',
                                    'Analysis type','Enrichment approach','Genome regions analyzed ',
                                    'Bases sequenced','Sequence Read Length']
    present_attributes = map(encode, overview_sheet.col_values(0))[first_index+1:last_index]

    if overview_expected_attributes != present_attributes:
        print >>sys.stderr, 'Missing expected attributes are:'
        print >>sys.stderr, set(overview_expected_attributes) - set(present_attributes)
        print >>sys.stderr, 'Unexpected attributes are:'
        print >>sys.stderr, set(present_attributes) - set(overview_expected_attributes) 
        
    values = map(encode, overview_sheet.col_values(1))[first_index+1:last_index]
    data = data + zip(present_attributes, values)
    data = filter(lambda x: x[0] != 'Analysis Characteristics', data)

    # summary sheet attributes
    summary_expected_attributes = ['Somatic (Tumor-Specific) Alterations',
                                   'Number of somatic sequence alterations identified', 
                                   'Number of somatic copy number alterations identified', 
                                   '',
                                   'Overall Statistics',
                                   'Sequenced Bases Mapped to Genome', 
                                   'Sequenced Bases Mapped to Target Regions', 
                                   'Fraction of Sequenced Bases Mapped to Target Regions', 
                                   'Bases in target regions with at least 10 reads', 
                                   'Fraction of bases in target regions with at least 10 reads',
                                   '', 
                                   'Sequence Reads at Each Base',
                                   'Average Number of Total High Quality Sequences at Each Base',
                                   'Average Number of Distinct High Quality Sequences at Each Base',
                                   '',
                                   'Tumor/normal Matching',
                                   'Germline SNPs present ',
                                   'Percent T/N Matching']
    
    present_attributes = map(encode, summary_sheet.col_values(0))[10:27]
    tumor_values = map(encode, summary_sheet.col_values(1))[10:27]
    normal_values = map(encode, summary_sheet.col_values(2))[10:27]
    
    data = data + zip(['Tumor ' + item for item in  present_attributes], tumor_values)
    data = data + zip(['Normal ' + item for item in  present_attributes], normal_values)
    
    empty_attributes = ['', 'Somatic (Tumor-Specific) Alterations', 'Overall Statistics',
                        'Sequence Reads at Each Base', 'Tumor/normal Matching',
                        'Somatic (Tumor-Specific) Alterations']

    blacklist_attributes = ['Normal ' + x for x in empty_attributes] + \
                          ['Tumor ' + x for x in empty_attributes]

    data = filter(lambda x: x[0] not in blacklist_attributes, data)
    data = filter(lambda x: x[0] != '' , data)
    return data
#-----------------------------------------------------------------------# 
def write_values(data, path):
    f_w = file(path, 'w')
    print >>f_w, '\n'.join(map('\t'.join, zip(*data)))
    f_w.close()
#-----------------------------------------------------------------------# 
def main():
    args = InputParser()
    data = read_report_data(args.input)
    write_values(data, args.output)

if __name__ == '__main__':
    main()

