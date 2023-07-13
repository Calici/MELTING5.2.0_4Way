import subprocess
import csv
import re
import math
import argparse


def strandIdentifierL(j):

    def divide_string(s, reverse=False):
        special_chars = ['L', '*', 'X_C', 'X_T']

        # Remove special characters for half length calculation
        s_no_special = s
        for char in special_chars:
            s_no_special = s_no_special.replace(char, '')
        half_length = len(s_no_special) // 2

        # Initialize variables
        half = ''
        count = 0
        cursor = len(s) - 1 if reverse else 0
        step = -1 if reverse else 1

        # Iterate over the string
        while 0 <= cursor < len(s):
            # Check if the substring from cursor is a special character
            is_special_char = any(s[cursor: cursor+len(char)].upper() == char for char in special_chars)
            if is_special_char:
                # If it's a special character, find which one and add it to the result
                for char in special_chars:
                    if s[cursor: cursor+len(char)].upper() == char:
                        if (count == half_length) and reverse:
                            return half
                        half = char + half if reverse else half + char
                        cursor += len(char)
            else:
                # If it's not a special character, add it to the result
                if (count == half_length):
                    return half
                count += 1
                half = s[cursor] + half if reverse else half + s[cursor]
                cursor += step
        return half

    input1 = row[j[0]].replace(" ", "")
    input2 = row[j[1]].replace(" ", "")
    
    half1 = divide_string(input1)
    half2 = divide_string(input2, reverse=True)
    # print(half1, half2)
    return half1, half2[::-1]






def strandIdentifier(j):
 
    input1 = row[j[0]].replace(" ", "")
    input1 = input1.replace("L", "").replace("*", "").replace("X_T", "").replace("X_C", "")
    length1 = len(input1)
    input1 = input1[0:len(input1)//2]
 
    input2 = row[j[1]].replace(" ", "")
    input2 = input2.replace("L", "").replace("*", "").replace("X_T", "").replace("X_C", "")
    input2 = input2[length1//2: ][::-1]
 
    return input1, input2
# Read rows from the CSV file
def melting_params(x):
    match_enthalpy = re.search(r'Enthalpy.*?([\-\d.]+)', x)
    match_entropy = re.search(r'Entropy.*?([\-\d.]+)', x)
    match_melting_temperature = re.search(r'Melting temperature.*?([\-\d.]+)', x)
    enthalpy = float(match_enthalpy.group(1) if match_enthalpy else 0)
    entropy = float(match_entropy.group(1) if match_entropy else 0)
    melting_temperature = float(match_melting_temperature.group(1) if match_melting_temperature else 0)

    gibbs_free_energy = enthalpy - (melting_temperature+273.15)* entropy
    return enthalpy, entropy, gibbs_free_energy
 
def mismatch_identifier(strand, complement_strand, mismatch_indices, strand_iteration):
    mismatches = 0
    for index, (s, c) in enumerate(zip(strand, complement_strand)):
        # if s == 'L' or c == 'L':
        #     continue
        if s == 'A' and c != 'T' or s == 'T' and c != 'A':
            mismatches += 1
            mismatch_indices[strand_iteration].append(index+1)
        elif s == 'C' and c != 'G' or s == 'G' and c != 'C':
            mismatches += 1
            mismatch_indices[strand_iteration].append(index+1)
    return mismatches, mismatch_indices
 
def write_to_csv(rows, out):
    with open(out, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Strand1', 'Strand2', 'Strand3', 'Strand4', 'Final Melting Temperature(degree Celsius)', 'Number of Mismatches', 'Mismatch Loc[Strand:Location on Strand]'])
        writer.writerows(rows)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Description of your program')
    # new options
    parser.add_argument('-out', '--output_file',  default='output.csv', help='Output file name')
    parser.add_argument('-in', '--input_file',  default='input.csv', help='Input file name')
    parser.add_argument('-sort', '--sort',  action='store_true', help='Sort the output file based on melting temperature')

    #existing options
    parser.add_argument('-E', '--electrolyte_concentration',  default='Na=1:Cl=1', help='Electrolyte concentration')          
    parser.add_argument('-P', '--ion_concentration', default='0.00000015', help='Ion concentration')
    parser.add_argument('-H', '--hybridization', default='dnadna', help='Hytbridization type')
    # parser.add_argument('-v', '--verbose', action='store_true', help='Verbose mode')
        # verbose option isn't useless for our case
    parser.add_argument('-T', '--threshold_value', default='60', help='Threshold for approximative computation')
    # parser.add_argument('-nnpath', '--folder_pathway', default='Data', help='Change the default pathway where to find the default calorimetric tables')
     # parser O option isn't implemented because our program doesn't use the option -O instead it uses -out and -in
    parser.add_argument('-self', '--self_complementary', action='store_true', help='To precise that the sequence entered is self complementary')
    parser.add_argument('-F', '--factor_value', default='4', help='Correction for the concentration of nucleic acid')
    parser.add_argument('-am', '--approximative_formula', default='wetdna91', help='Specific approximative formula')
    parser.add_argument('-nn', '--nearest_neighbor_model', default='all97', help='Specific nearest neighbor model')
    parser.add_argument('-sinMM', '--sinMM_model', default='allsanpey', help='Specific nearest neighbor model for single mismatch(es)')
    parser.add_argument('-GU', '--GU_model', default='ser12', help='Specific nearest neighbor model for GU base pairs')
    # parser.add_argument('-tanMM', '--tanMM_model', default='allsanpey', help='Specific nearest neighbor model for tandem mismatches')
    parser.add_argument('-intLP', '--intLP_model', default='san04', help='Specific nearest neighbor model for internal loop')
    parser.add_argument('-sinDE', '--sinDE_model', default='bom00', help='Specific nearest neighbor model for single dangling end(s)')
    parser.add_argument('-secDE', '--secDE_model', default='sugdna02', help='Specific nearest neighbor model for second dangling end(s)')
    parser.add_argument('-lonDE', '--lonDE_model', default='sugdna02', help='Specific nearest neighbor model for long dangling end(s)')
    parser.add_argument('-sinBU', '--sinBU_model', default='tan04', help='Specific nearest neighbor model for single bulge loop(s)')
    parser.add_argument('-lonBU', '--lonBU_model', default='san04', help='Specific nearest neighbor model for long bulge loop(s)')
    parser.add_argument('-CNG', '--CNG_model', default='bro05', help='Specific nearest neighbor model for RNA sequences composed of CNG repeats')
    parser.add_argument('-ino', '--ino_model', default='san05', help='Specific nearest neighbor model for inosine base')
    parser.add_argument('-ha', '--ha_model', default='sug01', help='Specific nearest neighbor model for hydroxyadenine base')
    parser.add_argument('-azo', '--azo_model', default='asa05', help='Specific nearest neighbor model for azobenzene')
    parser.add_argument('-lck', '--lck_model', default='owc11', help='Specific nearest neighbor model for locked nucleic acid')
    parser.add_argument('-tanLck', '--tanLck_model', default='owc11', help='Specific nearest neighbor model for consecutive locked nucleic acid')
    
    parser.add_argument('-sinMMLck', '--sinMMLck_model', default='owc11', help='Specific nearest neighbor model for consecutive locked nucleic acid with one single mismatch')
    parser.add_argument('-ion', '--ion_correction', default='owc2204', help='Specific ion correction')
    parser.add_argument('-naeq', '--naeq_correction', default='ahs01', help='Specific ion correction which gives a sodium equivalent concentration')
    parser.add_argument('-DMSO', '--DMSO_correction', default='ahs01', help='Specific DMSO correction')
    parser.add_argument('-for', '--formamide_correction', default='bla96', help='Specific formamide correction')
    # parser.add_argument('-help', '--help', action='help', help='Show this help message and exit')
    args = parser.parse_args()
    rows = []
    with open('input.csv', mode='r') as file:
        reader = csv.reader(file)
        next(reader)
        for row in reader:
            total_enthalpy=total_gibbs_free_energy=total_entropy=0
            Strands=[(0,3), (1,0), (2,1), (3,2)] #[[top-left], [bottom-right]]
            no_mismatches=0
            calc_temp=0
            strand_iteration=1
            mismatch_indices = {1:[], 2:[], 3:[], 4:[]}
            for one_side in Strands: #j is either top-left or bottom-right strand
                strandi, complement_strandi = strandIdentifier(one_side)
                mismatches, mismatch_indices = mismatch_identifier(strandi, complement_strandi, mismatch_indices, strand_iteration)
                no_mismatches+=mismatches

                strandL, complement_strandL= strandIdentifierL(one_side)


                command = ['./melting', 
                                        # '-v' if args.verbose else '', 
                                        "-S", strandL, 
                                        "-C", complement_strandL, 
                                        "-H", args.hybridization,
                                        '-P', args.ion_concentration, '-E', args.electrolyte_concentration, 
                                        '-T', args.threshold_value,
                                        # '-nnpath', args.folder_pathway,
                                        '-self' if args.self_complementary else '',
                                        '-F', args.factor_value,
                                        '-am', args.approximative_formula,
                                        '-nn', args.nearest_neighbor_model,
                                        '-sinMM', args.sinMM_model,
                                        '-GU', args.GU_model,
                                        # '-tanMM', args.tanMM_model,
                                        '-intLP', args.intLP_model,
                                        '-sinDE', args.sinDE_model,
                                        '-secDE', args.secDE_model,
                                        '-lonDE', args.lonDE_model,
                                        '-sinBU', args.sinBU_model,
                                        '-lonBU', args.lonBU_model,
                                        '-CNG', args.CNG_model,
                                        '-ino', args.ino_model,
                                        '-ha', args.ha_model,
                                        '-azo', args.azo_model,
                                        '-lck', args.lck_model,
                                        '-tanLck', args.tanLck_model,
                                        '-sinMMLck', args.sinMMLck_model,
                                        '-ion', args.ion_correction,
                                        '-naeq', args.naeq_correction,
                                        '-DMSO', args.DMSO_correction,
                                        '-for', args.formamide_correction
                                        ]
                # -nn sug96 gave out 35.1
                #'-ion','tanna06' gave 37.13
                #'-ion','ahs01' gave 41.86
                #'-ion','kam71' gave 35
                #'-ion','owc1904' gave 43
                #'-ion','owc2004' gave 45
                #'-ion','owc2104' gave 30
                #'-ion','owc2204' gave 30
                #'-ion','san96' gave 38.45
                # '-ion','sug04' gave 41.5
                # '-ion','schlif' gave 29.02
                


                # if(args.verbose):
                #     command.append('-v')
 
                # print(command)
                command_line = subprocess.run(command, 
                                            capture_output=True, 
                                            text=True)
                output = command_line.stdout
                # print(output)
                max_temp=0
                if mismatches>0:
                    calc_temp=float("nan")
                else:    
                    enthalpy, entropy, gibbs_free_energy = melting_params(output)
                    # if melting_temperature>max_temp:
                    #     max_temp=melting_temperature
                    total_enthalpy+=enthalpy
                    total_gibbs_free_energy+=gibbs_free_energy
                    total_entropy+=entropy
                strand_iteration+=1
            if(math.isnan(calc_temp)==False):
                calc_temp = ((total_enthalpy-total_gibbs_free_energy)/total_entropy)-273.15
            # row.extend([calc_temp, no_mismatches, mismatch_indices])
            row.extend([calc_temp, no_mismatches, mismatch_indices])
            # print((row))
            rows.append(row)
    if(args.sort):
        rows = sorted(rows, key=lambda x: (math.isnan(x[4]), x[4]), reverse=False) #if reverse is true then descending order
    print(rows)
    write_to_csv(rows, args.output_file)