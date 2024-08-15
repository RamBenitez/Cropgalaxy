#!/usr/bin/env python3

'''
Created on Jan 1, 2011

@author: John L. Herndon
@contact: herndon@cs.colostate.edu
@organization: Colorado State University
@group: Computer Science Department, Asa Ben-Hur's laboratory
'''

import sys
import time
import os
import getopt
from primertools import *

version = "0.5.0"

class UniqPrimerFinder(object):

    def __init__(self, includeFiles, excludeFiles, crossValidate, eprimerOptions):
        utils.logMessage("UniqPrimerFinder::__init__()", "Initializing UniqPrimerFinder")
        self.includeFiles = includeFiles
        self.includeFileManager = includefilemanager.IncludeFileManager()
        self.excludeFiles = excludeFiles
        self.excludeFileManager = excludefilemanager.ExcludeFileManager()
        self.primerManager = primermanager.PrimerManager(eprimerOptions)
        self.crossValidate = crossValidate
        utils.logMessage("UniqPrimerFinder::__init__()", "Initializing UniqPrimerFinder - complete")

    def writeOutputFile(self, primers, outputFileName, maxresults=100):
        with open(outputFileName, 'w') as outputFile:
            for i, primer in enumerate(primers, 1):
                outputFile.write("{0}\t{1}\t{2}\t{3}\n".format(i, primer.forwardPrimer, primer.reversePrimer, primer.productSize))
                if i > maxresults:
                    break
        utils.logMessage("UniqPrimerFinder::writeOutputFile()", "output file written.")

    def findPrimers(self, outputFile="uPrimer.txt"):
        outputFile = uPrimer
        utils.logMessage("UniqPrimerFinder::findPrimers()", "Finding primers for include files")
        startTime = time.time()
        
        utils.printProgressMessage("*** Creating Combined Fasta File for Exclude Files ***")
        for excludeFile in self.excludeFiles:
            self.excludeFileManager.addExcludeFile(excludeFile)
        self.excludeFileManager.exportSequences()
        
        self.includeFileManager.setExcludeFile(self.excludeFileManager.getOutputFileName())
        utils.printProgressMessage("*** Finding Sequences Unique to Target Genome ***")
        for includeFile in self.includeFiles:
            self.includeFileManager.processIncludeFile(includeFile)
                
        uniqueSequences = self.includeFileManager.getUniqueSequences()
        utils.printProgressMessage("*** Finding Primers ***")
        
        primers = self.primerManager.getPrimers(uniqueSequences)
         
        if self.crossValidate:
            utils.printProgressMessage("*** Cross Validating Primers ***")
            primers = self.primerManager.crossValidatePrimers(primers, self.excludeFileManager.getOutputFileName())
            for j, includeFile in enumerate(self.includeFiles, 1):
                primers = self.primerManager.crossValidatePrimers2(primers, includeFile, j)
        
        utils.logMessage("UniqPrimerFinder::findPrimers()", "found {0} unique sequences".format(len(primers)))
        self.writeOutputFile(primers, outputFile)
        utils.logMessage("UniqPrimerFinder::findPrimers()", "Finished finding primers")
        
        endTime = time.time()
        elapsedMinutes = int((endTime - startTime) / 60)
        elapsedSeconds = int((endTime - startTime) % 60)
        print("*** Time Elapsed: {0} minutes, {1} seconds ***".format(elapsedMinutes, elapsedSeconds))
        print("*** Output Written to {0} ***".format(outputFile))

def printUsageAndQuit():
    global version
    print("uniqprimer - finds primers unique to a genome")
    print("Version: " + str(version))
    print("Summary of Options.")
    print("Required Arguments:")
    print(" -i <filename>: use <filename> as an include file. Primers will be identified for this genome")
    print(" -x <filename>: use <filename> as an exclude file. Primers for this genome will be excluded")
    print(" -o <filename>: specify the name of the unique primer output file (default is uPrimer.txt)")
    print(" -l <filename>: specify the name of the log output file")
    print(" -f <filename>: specify the name of the Fasta of differential sequences")
    print("\nOptional Arguments:")
    print(" --productsizerage: set a range for the desired size of PCR product (default=200-250). Example: ./uniqprimer -productsizerage 100-150")
    print(" --primersize: set the desired primer size (default=20)")
    print(" --minprimersize: set the minimum primer size (default=27)")
    print(" --maxprimersize: set the maximum primer size (default=18)")
    print(" --crossvalidate: force the program to cross validate primers against exclude files for extra certainty")
    print(" --keeptempfiles: force the program to keep temporary files")
    print("\n\nExample:")
    print("uniqprimer -i <includefile1> -i <includefile2> ... -i <includefileN> -x <excludefile1> -x <excludefile2> ... -x <excludefileN> -o primers.txt -l logfile.txt -f seqForPrimer3.fa")
    utils.shutdownLogging()
    sys.exit()

opts = 'i:x:h:o:l:f:'
longopts = ["productsizerange=", "primersize=", "minprimersize=", "maxprimersize=", "crossvalidate", "keeptempfiles"]

def parseArgs(args):
    global uPrimer
    global lf
    global fastaDiff

    crossValidate = False
    cleanup = True
    optlist, args = getopt.getopt(args, opts, longopts)
    
    includeFiles = []
    excludeFiles = []
    eprimerOptions = utils.EPrimerOptions()
    
    verbose = False
    for opt in optlist:
        if opt[0] == '-i':
            includeFiles.append(opt[1])
        elif opt[0] == '-x':
            excludeFiles.append(opt[1])
        elif opt[0] == '-v':
            verbose = True 
        elif opt[0] == '-o':
            uPrimer = str(opt[1])
        elif opt[0] == '-l':
            lf = str(opt[1])
        elif opt[0] == '-f':
            fastaDiff = str(opt[1])
        elif opt[0] == '--productsizerange':
            eprimerOptions.setProductRange(opt[1])
        elif opt[0] == '--primersize':
            eprimerOptions.setPrimerSize(opt[1])
        elif opt[0] == '--minprimersize':
            eprimerOptions.setMinPrimerSize(opt[1])
        elif opt[0] == '--maxprimersize':
            eprimerOptions.setMaxPrimerSize(opt[1])
        elif opt[0] == '--crossvalidate':
            crossValidate = True
        elif opt[0] == '--keeptempfiles':
            cleanup = False
        elif opt[0] == '-h':
            printUsageAndQuit()
        else:
            print("Unknown option: " + str(opt[0]))
            printUsageAndQuit()
    
    if len(includeFiles) == 0 or len(excludeFiles) == 0:
        print("You must specify at least one include file and at least one exclude file")
        printUsageAndQuit()

    return includeFiles, excludeFiles, crossValidate, cleanup, verbose, eprimerOptions, lf, uPrimer, fastaDiff

def main(args, debug=False):
    includeFiles, excludeFiles, crossValidate, cleanup, verbose, eprimerOptions, lf, uPrimer, fastaDiff = parseArgs(args)
    utils.initialize(True, cleanup, lf)
    
    tmpdir = utils.getTemporaryDirectory()
    command = "cp -rf " + tmpdir + "/sequenceForEprimer.fasta" + " " + fastaDiff
 
    try:
        utils.logMessage("uniqprimer::Main()", "Logging include files: ")
        utils.logList("uniqprimer::Main()", includeFiles)
        utils.logMessage("uniqprimer::Main()", "Logging exclude files: ")
        utils.logList("uniqprimer::Main()", excludeFiles)
        print("*** Finding Primers ***")
        uniqPrimer = UniqPrimerFinder(includeFiles, excludeFiles, crossValidate, eprimerOptions)
        uniqPrimer.findPrimers()
    except utils.NoFileFoundException as nfe:
        print("File not found: " + str(nfe.filename))
        printUsageAndQuit()
    except utils.ProgramNotFoundException as pnfe:
        print(str(pnfe.programName) + ": program is not installed or is not in your path.")
        print(str(pnfe.details))
    except utils.NoPrimersExistException:
        print("Failure: No unique primers exist for this combination")
    except BaseException as e:
        print("An unexpected error occurred. Please check 'log_uniqprimer.txt' and report the issue.")
        print("Details:")
        print(e)
    
    os.system(command)
    utils.shutdown()

    print("*** Finished ***")
    
if __name__ == '__main__':
    if len(sys.argv) == 1:
        printUsageAndQuit()
    main(sys.argv[1:], debug=True)
