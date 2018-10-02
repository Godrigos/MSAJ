#!/usr/bin/env python3

"""
This file is part of MSAJ.

MSAJ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MSAJ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MSAJ.  If not, see <https://www.gnu.org/licenses/>.

Copyright 2018 Rodrigo Aluizio
"""

from Bio import AlignIO
from glob import glob
from sys import exit
from os.path import splitext, basename
from seq_check import seq_check

# Define fastA file type extensions
types = ('*.fas', '*.fa', '*.faa', '*.fasta')

# Create fastA files list and populate it
files = []

for file in types:
    files.extend(glob(file))

files = sorted(files, key=str.lower)

# Read fastA files and store data in a list
msa = []
i = 0
for i in range(len(files)):
    msa.append(AlignIO.read(files[i], "fasta"))

# Create combined sequences
try:
    diff_ids = seq_check(msa)

    if diff_ids == set():
        pass
    else:
        print("Strains:")
        for value in diff_ids:
            print(str(value))
        print("do not have matches in all files!\n")
        print("Not merging them!")
        print("Press Enter to exit!")
        input()
        exit(1)

    msa[0].sort()
    combined = msa[0]
    i = 1
    for i in range(1, len(msa)):
        msa[i].sort()
        combined = combined + msa[i]

    # Save multiocus multiple sequence alignment
    with open("MLMSA.fas", "w") as f:
        AlignIO.write(combined, f, "fasta")
    # Print some simple but useful information
    print("Multilocus Multiple Sequence Alignment created (MLMSA.fas)!")
    print("\n----- Useful Information -----")
    print("Number of Loci: " + str(len(msa)))
    print("Number of Sequences: " + str(len(msa[0])))
    count = 1
    i = 0
    for i in range(len(files)):
        print("Locus " + splitext(basename(files[i]))[0] + ": " + str(count) +
              " - " + str(msa[i].get_alignment_length()+count-1))
        count = count + msa[i].get_alignment_length()
    print("------------------------------")
        
    with open('MLMSA.txt', 'w') as f:
        print("Multilocus Multiple Sequence Alignment created (MLMSA.fas)!",
              file=f)
        print("\n----- Useful Information -----", file=f)
        print("Number of Loci: " + str(len(msa)), file=f)
        print("Number of Sequences: " + str(len(msa[0])), file=f)
        count = 1
        i = 0
        for i in range(len(files)):
            print("Locus " + splitext(basename(files[i]))[0] + ": " +
                  str(count) + " - " +
                  str(msa[i].get_alignment_length()+count-1), file=f)
            count = count + msa[i].get_alignment_length()
        print("------------------------------\n", file=f)
    
    print("\nPress Enter to exit ...")
    input()
    exit(0)

except IndexError:
    print("No alignments files in the directory!")
    print("Press Enter to exit ...")
    input()
    exit(1)

