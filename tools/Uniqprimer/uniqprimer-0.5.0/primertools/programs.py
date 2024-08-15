'''
Created on Jan 1, 2011

@author: John L. Herndon
@contact: herndon@cs.colostate.edu
@organization: Colorado State University
@group: Computer Science Department, Asa Ben-Hur's laboratory 
'''

from . import utils
import subprocess


class ProgramBase:
    
    def __init__(self):
        self.programName = None
        self.proc = None
    
    def getProcessArgs(self, args):
        raise NotImplementedError("This method should be overridden by subclasses")
        
    def execute(self, args, async_=False):
        '''
        Run the program with given arguments
        '''
        
        utils.logMessage("ProgramBase::Execute()", f"Running the {self.programName} program.")
        
        args, outputFile = self.getProcessArgs(args)
        
        print(f"*** Running {self.programName} ***")
        
        utils.logList("ProgramBase::Execute()", args)
        
        proc = subprocess.Popen(args)
        
        if not async_:
            proc.wait()
            print(f"*** Running {self.programName} Complete ***")
        
        return outputFile
    
    def isFinished(self):
        return self.proc.poll() is not None
        
class Nucmer(ProgramBase):
    
    def __init__(self):
        nucmerPath = utils.search_file('nucmer')
        super().__init__()
        if nucmerPath is None:
            raise utils.ProgramNotFoundException('nucmer', "Please ensure that the MUMmer package is installed and configured on your system.")
        
        self.nucmer = nucmerPath
        self.programName = "nucmer"
        self.outputExtension = ".coords"
        
    def getProcessArgs(self, inputArgs):
        timestamp = utils.getTimeStamp()
        identifier = "nucmer_alignments"
        args = [self.nucmer, '-p', identifier, '-o', '--minmatch', '300', '--maxgap', '1']
        args.extend(inputArgs)
        outputFile = f"{identifier}.coords"
        
        return args, outputFile

class Eprimer(ProgramBase):
    
    def __init__(self, eprimerOptions):
        super().__init__()
        self.programName = "EPrimer3"
        self.options = eprimerOptions
        
        primer3corePath = utils.search_file("primer3_core")
        if primer3corePath is None:
            raise utils.ProgramNotFoundException("primer3_core", "Please ensure that the primer3 package is installed on your system. It can be obtained from http://primer3.sourceforge.net/")
        
        eprimerPath = utils.search_file("eprimer3")
        if eprimerPath is None:
            raise utils.ProgramNotFoundException('eprimer3', "Please ensure that the EMBOSS package is installed and configured on your system.")
        
        self.primer3core = primer3corePath
        self.eprimer3 = eprimerPath
        
    def getProcessArgs(self, inputArgs):
        inputFasta = inputArgs[0]
        outputFile = inputArgs[1]
        args = [
            self.eprimer3, inputFasta, outputFile, '-numreturn', '2', 
            '-prange', self.options.getProductRange(), '-osize', str(self.options.getPrimerSize()),
            '-minsize', str(self.options.getMinPrimerSize()), '-maxsize', str(self.options.getMaxPrimerSize())
        ]
        
        return args, outputFile
    
class PrimerSearch(ProgramBase):
    def __init__(self):
        super().__init__()
        self.programName = "PrimerSearch"
        primerSearchPath = utils.search_file("primersearch")
        if primerSearchPath is None:
            raise utils.ProgramNotFoundException("primersearch", "Please ensure that the EMBOSS package is installed on your system.")
    
        self.primerSearch = primerSearchPath
        
    def getProcessArgs(self, inputArgs):
        args = [self.primerSearch]
        args.extend(['-seqall', inputArgs[0]])
        args.extend(['-infile', inputArgs[1]])
        args.extend(['-mismatchpercent', inputArgs[3]])
        args.extend(['-outfile', inputArgs[2]])
    
        return args, inputArgs[2]
