'''
Created on Jan 1, 2011

@author: John L. Herndon
@contact: herndon@cs.colostate.edu
@organization: Colorado State University
@group: Computer Science Department, Asa Ben-Hur's laboratory 
'''

from . import utils
import tempfile
from . import programs
from . import eprimerparser
from . import primersearchutils
from . import fastaparser
from Bio import Seq

# Custom exception class
class PrimerManagerError(Exception):
    pass

class PrimerManager(object):
    '''
    A class used to find primers given a set of sequences.
    '''

    def __init__(self, eprimerOptions):
        self.eprimer = programs.Eprimer(eprimerOptions)
        self.primersearch = programs.PrimerSearch()
    
    def findPrimers(self, sequences, outputFile, primerpairs=20, returnPrimers=False):
        '''
        A method to find a set of primers based on the given sequences
        '''
        
        utils.logMessage("PrimerManager::findPrimers", "writing sequences to a fasta file")
        
        # Eliminate all sequences that are less than the desired amplification size...
        if len(sequences) == 4:
            print(sequences)
        sequences = [x for x in sequences if len(x) >= 200]
        
        primerFastaFile = utils.getTemporaryDirectory() + "/sequenceForEprimer.fasta"
        fastaparser.writeFastaFile(sequences, primerFastaFile)

        utils.logMessage("PrimerManager::findPrimers", "executing eprimer3 program")
        self.eprimer.execute([primerFastaFile, outputFile])
        utils.logMessage("PrimerManager::findPrimers", "eprimer3 file {0} created. Parsing for primers.".format(outputFile))
        
        primers = eprimerparser.parsePrimerSequences(outputFile)
        
        utils.logMessage("PrimerManager::findPrimers", "parsing for sequences complete")
        
        if returnPrimers:
            return primers
    
    def getPrimers(self, sequences):
        utils.logMessage("PrimerManager::getPrimers", "finding primers that are common to all include files")
        
        if not sequences:
            raise PrimerManagerError("No primers exist.")
        
        referenceEPrimerFile = utils.getTemporaryDirectory() + "/referenceprimers.ep3"
        
        # Run eprimer to find primers in the reference file
        primers = self.findPrimers(sequences, referenceEPrimerFile, 20, True)
        
        if not primers:
            raise PrimerManagerError("No primers found.")
        
        return primers
    
    def crossValidatePrimers2(self, primers, includeFile, j):    
        includeSequences = fastaparser.parseFastaFile(includeFile)
        # Write a primer search input file using the primers argument
        primerInputFileName = utils.getTemporaryDirectory() + "/tmpinputprimers2.ps" + str(j)
        primerOutputFileName = utils.getTemporaryDirectory() + "/tmpoutputprimers2.ps" + str(j)
        primersearchutils.writePrimerSearchInputFile(primers, primerInputFileName)

        utils.logMessage("PrimerManager::crossValidatePrimers2", "finding primers that are in the supplied include file")
        # Run primer search to identify the primers
        self.primersearch.execute([includeFile, primerInputFileName, primerOutputFileName, "0"])

        # Read the found primers from the file
        commonPrimers = primersearchutils.parsePrimerSearchFile(primerOutputFileName)

        # Compose a list of primers that are not found in the exclude file...
        returnValue = [primer for primer in primers if primer.id in commonPrimers]

        utils.logMessage("PrimerManager::crossValidatePrimers2", "{0} unique primers identified out of {1}".format(len(returnValue), len(primers)))

        if not returnValue:
            raise PrimerManagerError("No unique primers found.")
        
        return returnValue
    
    def crossValidatePrimers(self, primers, excludeFile):
        excludeSequences = fastaparser.parseFastaFile(excludeFile)
        
        # Write a primer search input file using the primers argument
        primerInputFileName = utils.getTemporaryDirectory() + "/tmpinputprimers.ps"
        primerOutputFileName = utils.getTemporaryDirectory() + "/tmpoutputprimers.ps"
        primersearchutils.writePrimerSearchInputFile(primers, primerInputFileName)
        
        utils.logMessage("PrimerManager::crossValidatePrimers", "finding primers that are not in the supplied exclude file")
        # Run primer search to identify the primers
        self.primersearch.execute([excludeFile, primerInputFileName, primerOutputFileName, "10"])
        
        # Read the found primers from the file
        commonPrimers = primersearchutils.parsePrimerSearchFile(primerOutputFileName)
        
        # Compose a list of primers that are not found in the exclude file...
        returnValue = [primer for primer in primers if primer.id not in commonPrimers]
        
        utils.logMessage("PrimerManager::crossValidatePrimers", "{0} unique primers identified out of {1}".format(len(returnValue), len(primers)))
        
        if not returnValue:
            raise PrimerManagerError("No unique primers found.")
        
        return returnValue
