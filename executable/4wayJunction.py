import subprocess
import csv
import re
import math
import argparse

def strandIdentifier(j):
 
    input1 = row[j[0]].replace(" ", "")
    length1 = len(input1)
    input1 = input1[0:len(input1)//2]
 
    input2 = row[j[1]].replace(" ", "")
    input2 = input2[length1//2: length1][::-1]
 
 
 
    return input1, input2
# Read rows from the CSV file
def melting_params(x):
    match_enthalpy = re.search(r'Enthalpy.*?([\-\d.]+)', x)
    match_entropy = re.search(r'Entropy.*?([\-\d.]+)', x)
    match_melting_temperature = re.search(r'Melting temperature.*?([\-\d.]+)', x)
    enthalpy = float(match_enthalpy.group(1))
    entropy = float(match_entropy.group(1))
    melting_temperature = float(match_melting_temperature.group(1))
    gibbs_free_energy = enthalpy - (melting_temperature+273.15)* entropy
    return enthalpy, entropy, gibbs_free_energy
 
def mismatch_identifier(strand, complement_strand, mismatch_indices, strand_iteration):
    mismatches = 0
    for index, (s, c) in enumerate(zip(strand, complement_strand)):
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
    parser.add_argument('-P', '--ion_concentration', default='0.00000015', help='Ion concentration')
    parser.add_argument('-E', '--electrolyte_concentration',  default='Na=0.005:Cl=0.005', help='Electrolyte concentration')
    parser.add_argument('-out', '--output_file',  default='output.csv', help='Output file name')
    parser.add_argument('-in', '--input_file',  default='input.csv', help='Input file name')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose mode')
    parser.add_argument('-sort', '--sort',  action='store_true', help='Sort the output file based on melting temperature')
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
                strand, complement_strand = strandIdentifier(one_side)
                mismatches, mismatch_indices = mismatch_identifier(strand, complement_strand, mismatch_indices, strand_iteration)
                no_mismatches+=mismatches


                command = ['./melting',"-S", strand, 
                                        "-C", complement_strand, "-H", "dnadna",
                                        '-P', args.ion_concentration, '-E', args.electrolyte_concentration]
 
 
                if(args.verbose):
                    command.append('-v')
 
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