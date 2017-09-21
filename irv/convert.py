from __future__ import print_function
import csv
import collections

FileDetails = collections.namedtuple('FileDetails', ['ballotCount', 'candidateNames','candidateArrayLen'])

def convertPrefs(prefs, candidates, candArrayLen):
    out = []
    for i in range(0, len(prefs)):
        row = [0] * candArrayLen
        for j in range(0, len(prefs)):
            if prefs[j] == i + 1:
                row[j] = 1
        out.append(row)

    #print(' '.join(str(x) for x in out))
    #print("ot")
    #print(out)
    outstring = ""
    for inner in out:
        outstring = outstring + ' '.join(str(x) for x in inner) + '\n'
    return outstring

def next_greater_power_of_2(x):
    return 2**(x-1).bit_length()


def createSampleFile(fileDetails, outFile):
    print('# generated ballot file',file=outFile)
    print('# CandidateNames:' + str(fileDetails.candidateNames),file=outFile)
    print(str(fileDetails.ballotCount),file=outFile)
    print(str(len(fileDetails.candidateNames)),file=outFile)
    print(str(len(fileDetails.candidateNames)),file=outFile)
    print(str(fileDetails.candidateArrayLen),file=outFile)

def analyse_file(file_path):
    current = -1
    f = open(file_path, 'rb')
    name = set()
    reader = csv.reader(f, dialect='excel-tab')
    next(reader, None)
    ballot_counter = 0
    for row in reader:
        if row[3] == 'Formal':
            if row[2] != current:
                if current != -1:
                    ballot_counter = ballot_counter + 1
                current = row[2]
            name.add(row[4])
    f.close()
    ballot_counter = ballot_counter + 1
    #candArrayLen =next_greater_power_of_2(len(list(name))+1)
    return FileDetails(ballot_counter, list(name), len(list(name)))

def convert_file(file_details, infile, outfile):
    reader = csv.reader(infile, dialect='excel-tab')
    current = -1
    candidates = len(file_details.candidateNames)
    cand_name = file_details.candidateNames
    prefs = [0] * candidates
    next(reader, None)
    for row in reader:
        if row[3] == 'Formal':
            if row[2] != current:
                if current != -1:
                    print(convertPrefs(prefs, candidates, file_details.candidateArrayLen), file=outfile)
                current = row[2]
                prefs = [0] * candidates
            if row[6] != '':
                prefs[cand_name.index(row[4])] = int(row[6])
    print(convertPrefs(prefs, candidates, file_details.candidateArrayLen), file=outfile)
    infile.close()
    outfile.close()

file_details = analyse_file('../data/SGE2015 LA Pref Data_NA_Auburn.txt')
out_file = open('../data/auburn.txt', 'w')
#createSampleFile(fileDetails,outFile)
in_file = open('../data/SGE2015 LA Pref Data_NA_Auburn.txt','rb')
convert_file(file_details, in_file, out_file)
