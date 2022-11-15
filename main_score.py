#
#  Djamal Boukerroui, 08-2022
#
import json
import textwrap

# from scoresLib import polygon_plot
from scoresLib import score_autocontours_lib as scl
import argparse
import csv
import numpy as np
import time
import os.path
from os import path

def load_cases_to_process_csv(name):
    """
    Load a set of cases to process from a csv file. The file has the bellow header
    Patient ID , GT RTSS , DLC_RTSS

    :param name: full path the input csv file
    :return: 3 lists corresponding to: test number, Reference RTSS filenames and test RTSS filenames
    """
    test_nbs = []
    ref_cases = []  # empty list
    test_cases = []
    with open(name, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            #  skip testing for the header
            if len(ref_cases) == 0:
                test_nbs.append(row[0].strip())
                ref_cases.append(row[1].strip())
                test_cases.append(row[2].strip())
                continue
            if os.path.isfile(row[1].strip()) and os.path.isfile(row[2].strip()):
                test_nbs.append(row[0].strip())
                ref_cases.append(row[1].strip())
                test_cases.append(row[2].strip())
            else:
                print('\tFound a non valid line in input case file')
    return test_nbs[1:], ref_cases[1:], test_cases[1:]


def usage():
    help_text = textwrap.dedent('''
        Contour metrics scoring:  
        Computes contouring performance measures between DICOM RTSS files''')
    parser = argparse.ArgumentParser(description=help_text, formatter_class=argparse.RawTextHelpFormatter)

    help_text = textwrap.dedent('''
    Full path to csv file having the following structures: 
    Test Nb, Reference RTSS, Test RTSS
    1, c:/full/path/GT/RTSS/filename1.dcm, c:/full/path/DLC/RTSS/filename1.dcm
    2, c:/full/path/GT/RTSS/filename2.dcm, c:/full/path/DLC/RTSS/filename2.dcm''')
    parser.add_argument('cases_csv', type=str, help=help_text)

    help_text = 'Full path to output file (without extension) to save results.'
    parser.add_argument('output_csv', type=str, help=help_text)

    return parser.parse_args()


def main():
    args = usage()
    print('\nReading input cases file:', args.cases_csv)
    test_ids, ref_cases, test_cases = load_cases_to_process_csv(args.cases_csv)
    print('\tFound {0} cases to process \n'.format(len(ref_cases)))


    start_time = time.time()
    case_score = {}
    nb_success = 0
    for i in range(0, len(ref_cases)):  # Exclude the first row of the csv file
        try:
            print('\tNow processing case {:s}'.format(test_ids[i]))
            scoresfound = scl.score_a_case(ref_cases[i], test_cases[i])
            case_score['case{0}'.format(i)] = (test_ids[i], scoresfound)
            nb_success += 1
        except Exception as ME:
            print(ME.message, flush=True)
            print('\t\t Something went wrong for this case')

    elapsed_time = time.time() - start_time

    print('\nSuccessfully processed {0} out of {1} in  {2} seconds\n'.format(nb_success, len(ref_cases), int(elapsed_time)))

    #FP = TEST and (not REF)
    #FN  = REF and (not TEST)
    dictfilt = lambda x, y: dict([(i, x[i]) for i in x if i in set(y)])
    measures = ['RefOrgan', 'TestOrgan',
                'HD', 'HD95_ref2test', 'HD95_test2ref',
                'MD_ref2test', 'MD_test2ref',
                'AD_ref2test', 'AD_test2ref',
                'RefVol', 'TestVol', 'TPVol', 'FNVol', 'FPVol',
                'three_D_DSC', 'ref_length', 'test_length', 'APL', 'cent_dist']
    with open(args.output_csv, mode='w', newline='\n', encoding='utf-8') as out_file:
        result_writer = csv.writer(out_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        result_writer.writerow(['Test ID'] + measures)
        for case_id, result in case_score.items():  # loop over patient
            test_id = result[0]
            current_results = result[1]
            for score in current_results:  # loop over all organs in results
                values = [test_id, score[measures[0]], score[measures[1]]]
                num_values = ['{:8.6f}'.format(score[x]) for x in measures[2:]]
                values = values + num_values
                #print(score['HD'], score['TestVol'])
                result_writer.writerow(values)

if __name__ == '__main__':
    main()
