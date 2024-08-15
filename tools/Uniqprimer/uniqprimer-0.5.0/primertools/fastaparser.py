'''
Created on Jan 1, 2011

@author: John L. Herndon
@contact: herndon@cs.colostate.edu
@organization: Colorado State University
@group: Computer Science Department, Asa Ben-Hur's laboratory 
'''


from . import utils
from . import primersequence

from Bio import SeqIO
from Bio import Seq
from functools import reduce

def parseFastaFileAsPrimerSequence(fileName):
    utils.logMessage("fastaparser::parseFastaFileAsPrimerSequence( )", f"parsing fasta file {fileName}")
    returnValue = {}

    sequences = SeqIO.parse(open(fileName), "fasta")

    for sequence in sequences:
        seqdata = primersequence.PrimerSequence(sequence.id, len(sequence), sequence.seq)
        returnValue[sequence.id] = seqdata

    utils.logMessage("fastaparser::parseFastaFileAsPrimerSequence( )", f"read {len(list(returnValue.keys()))} sequences")

    return returnValue

def parseFastaFile(fileName):
    '''
    parse a fasta file and return a list of Bio.Seq
    '''
    utils.logMessage("fastaparser::parseFastaFile( )", f"parsing fasta file {fileName}")

    sequences = SeqIO.parse(open(fileName), "fasta")

    return sequences

def writeFastaFile(sequences, fileName):
    '''
    write a set of sequences to a fasta file.
    returns the name of the new file
    ''' 
    primerSequenceIdent = "primer_sequences"
    utils.logMessage("PrimerManager::writeFastaFile( )", f"Writing {len(sequences)} sequences to fasta file")
    seqRecords = []
    i = 0
    for sequence in sequences:
        seqStr = str(reduce(lambda x, y: str(x) + str(y), sequence))
        seqRecord = SeqIO.SeqRecord(Seq.Seq(seqStr), id=f"seq_{i}")
        seqRecord.annotations["molecule_type"] = "DNA"
        seqRecords.append(seqRecord)
        i += 1

    SeqIO.write(seqRecords, open(fileName, "w"), "fasta")

    utils.logMessage("PrimerManager::writeFastaFile( )", "writing fasta file complete")    
    return fileName
