import random
import os
import pandas as pd
import numpy as np
import math
import csv
from Bio.Blast.Applications import NcbiblastnCommandline
import subprocess
from subprocess import PIPE
import matplotlib.pyplot as plt

os.chdir("D:")
# generate random sequence
def generateRandomSequence(length):
    return ''.join(random.choice('ACGT') for _ in range(length))
# apply mutations to a sequence
def applyMutations(sequence, mutationRate):
    mutatedSequence = list(sequence)  
    for i in range(len(mutatedSequence)):
        if random.random() < mutationRate:  
            mutatedSequence[i] = random.choice('ACGT'.replace(mutatedSequence[i], ''))  # randomly select a base to replace
    return ''.join(mutatedSequence)
# impliment reverse complement
def reverseComplement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))
# insert repeated sequences at specified positions
def insertRepeatedSequenceAtPositions(genome, repeatSeqs, positions, mutation_rate):
    repeatInfo = {}  # record repeat information
    shiftPositions = []  # record shift positions
    positions.sort(key=lambda x: x[1]) 
    for repeat_type, pos, reverse, repeat_seq in positions: 
        original_repeat_seq = repeat_seq  
        # apply mutations
        repeat_seq = applyMutations(repeat_seq, mutation_rate)
        if reverse:
            repeat_seq = reverseComplement(repeat_seq)  
        genome = genome[:pos] + repeat_seq + genome[pos:]
        if repeat_type not in repeatInfo:
            repeatInfo[repeat_type] = {
                'sequence': original_repeat_seq,
                'mutatedSequence': repeat_seq,
                'positions': []
            }
        repeatInfo[repeat_type]['positions'].append((pos, pos + len(repeat_seq) - 1, reverse))
        shiftPositions.append((pos, len(repeat_seq)))
    genomeWithRepeats = genome
    for repeat_type, pos, reverse, repeat_seq in positions:
        for shiftPos, shiftLen in shiftPositions:
            # update positions
            if pos > shiftPos:
                pos += shiftLen  

    return genomeWithRepeats, repeatInfo
def generateFastaFileWithRepeats(total_length, repeatSeqs, repeat_count_per_seq, min_gap, mutation_rate,
                                     output_file, info_file):
    genome = generateRandomSequence(total_length)  # generate original genome
    # generate repeats in random positions
    positions = []
    for i, repeat_seq in enumerate(repeatSeqs):
        repeat_type = i + 1  
        for _ in range(repeat_count_per_seq):
            min_pos = random.randint(0, total_length - len(repeat_seq) - min_gap)
            max_pos = min_pos + min_gap  
            pos = random.randint(min_pos, max_pos)
            reverse = random.choice([True, False])  # reverse or not
            positions.append((repeat_type, pos, reverse, repeat_seq))
    # insert repeated sequences
    genome_with_repeats, repeatInfo = insertRepeatedSequenceAtPositions(genome, repeatSeqs, positions,
                                                                             mutation_rate)
    # generate fasta file
    with open(output_file, 'w') as f, open(info_file, 'w') as info_f:
        f.write(">contig_1\n")
        f.write(genome_with_repeats + "\n")  

        for _, start, _, seq in positions:
                end = start + len(seq) - 1
                # record start and end position
                info_f.write(f"{start}, {end}\n")
# parameters that can be changed
repeatSeqs = [generateRandomSequence(random.randint(40000, 60000)) for _ in range(8)]  
totalLength = 50000000  
repeatCountPerSeq = 20  
minGap = 1000  
mutationRate = 1e-4  

inputSequenceFile = "sequence.fasta"
def createBlastDb(fastaFile, dbType, outputPath):
    # makeblastdb command
    makeblastdbCommand = [
        "makeblastdb",
        "-in", fastaFile,
        "-dbtype", dbType,
        "-out", outputPath
    ]
    subprocess.run(makeblastdbCommand, check=True)
fastaFile = inputSequenceFile  
dbType = 'nucl'  
outputPath = inputSequenceFile  
def txtToCsv(file_path):
    with open(file_path, 'r') as txt_file, open('result.csv', 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        lines = txt_file.readlines()
        for line in lines:
            data = line.strip().split('\t')
            writer.writerow(data)
def calculateOverlapLength(seg1, seg2):
    start_max = max(seg1[6], seg2[0])
    end_min = min(seg1[7], seg2[1])
    return max(0, end_min - start_max)
def calculatePtr(groundTruth, sequencing_result):
    # precision = tp/(tp+fp)
    # recall = tp/(tp+fn)
    tp = 0
    total_length = sum(row[3] for row in df)
    total_length = total_length / 20
    gtl = sum([(gt[1] - gt[0]) for gt in groundTruth])
    for result in sequencing_result:
        for gt in groundTruth:
            overlap = calculateOverlapLength(result, gt)
            tp += overlap
    tp = tp / 20
    precision = tp / total_length if total_length > 0 else 0
    recall = tp / gtl
    return precision, recall
for loop in range(5):
    generateFastaFileWithRepeats(totalLength, repeatSeqs, repeatCountPerSeq, minGap, mutationRate,
                                 "sequence.fasta",
                                 "groundTruth.csv")
    createBlastDb(fastaFile, dbType, outputPath)
    blastnCline = NcbiblastnCommandline(query=inputSequenceFile, db=inputSequenceFile, evalue=0.00001, outfmt=6, out="result.txt", word_size=11)
    stdout, stderr = blastnCline()
    print('BLASTn completed.')
    txtToCsv("result.txt")
    groundTruth = np.array(
        pd.read_csv('groundTruth.csv'))
    df=pd.read_csv('result.csv')
    df = df[df.iloc[:, 3] >= 20000]
    df.to_csv('result1.csv', index=False)
    sequencing_result = np.array(
    pd.read_csv('result1.csv'))
    sequencing_result=sequencing_result[1:]
    df = np.array(
        pd.read_csv('result.csv'))
    df=df[1:]
    precision, recall = calculatePtr(groundTruth, sequencing_result)
    f1 = 2 * (precision * recall) / (precision + recall) 
    gscore = math.sqrt(precision * recall)  
    #save the scores in a csv file
    with open("record.csv",'a',newline='') as f:
        writer = csv.writer(f)
        writer.writerows([[f"{precision:.4f}",f"{recall:.4f}",f"{f1:.4f}",f"{gscore:.4f}"]])
    df = pd.read_csv('record.csv')
column_name = df.columns[2] 
last_rows_data = df[-5:][column_name]
last_rows_data.plot(kind='box', title=f'F1-score with {repeatCountPerSeq} repeat counts and {mutationRate} mutation rate')
plt.xlabel('mutation rate') 
plt.ylabel('F1-score')
xlabels = [f'{mutationRate}']
plt.show()