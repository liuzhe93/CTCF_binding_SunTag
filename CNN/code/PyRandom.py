from __future__ import division
import os,sys
import random

for i in range(len(sys.argv)):
    "sys.argv[%d] = %s" % (i, sys.argv[i])
if len(sys.argv) != 4:
    print('''python PyRandom.py input output num_sample num_sum''')
    exit(0)

fin = sys.argv[1]
fout = sys.argv[2]
number_to_sample = int(sys.argv[3])                           


with open(fin) as input:
    num_lines = sum([1 for line in input])
total_records = int(num_lines / 2)    

print("sampling " + str(number_to_sample) + " out of " + str(total_records) + " records")

records_to_keep = set(random.sample(xrange(total_records + 1), number_to_sample))
record_number = 0

with open(fin) as input:
    with open(fout, "w") as output:
        for line1 in input:
            line2 = input.next()
            if record_number in records_to_keep:
                output.write(line1)
                output.write(line2)
            record_number += 1







