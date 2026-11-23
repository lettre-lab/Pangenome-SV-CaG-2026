# Update ids for a SVs vcf (chr-pos-type-len)
#python update_ids.py input.vcf output.vcf
import sys

input_name = sys.argv[1]
output_name = sys.argv[2]
output = open(output_name, "w")

with open(input_name, "r") as input:
    prev_pos = 0
    allele_number = 2

    line = input.readline()

    while line:
        if line[0] == '#':
            output.write(line)
            line = input.readline()
      
        else:  
            
            
            next_line = input.readline()
            line_ = line.split('\t')
            info = line_[7].split(';')
            if line_[2].split('-')[2] == "DEL":
              id = line_[2].split('-')[0] + '-' + line_[2].split('-')[1] + '-' + line_[2].split('-')[2] + '-' + str(info[7].split('=')[1])
            else:
              id = line_[2].split('-')[0] + '-' + line_[2].split('-')[1] + '-' + line_[2].split('-')[2] + '-' + str(info[8].split('=')[1])
            #print(id)
            # new id has the form: chr-pos-type-len
            
            # cases of multiallelic
            if prev_pos == line_[1]:
                id = id + '-' + str(allele_number)
                allele_number += 1

            # check if it's the first multiallelic variant
            elif next_line:
                if line_[1] == next_line.split('\t')[1]:
                    id = id + '-' + '1'

                else:
                    allele_number = 2

            #break

            line_[2] = id
            new = ""
            prev_pos = line_[1]
            for i in line_:
                if i is line_[-1]:
                    new += str(i.strip())
                else:
                    new += str(i) + "\t"
      
            output.write(new + "\n")
            line = next_line

