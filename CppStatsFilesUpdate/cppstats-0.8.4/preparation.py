#!/usr/bin/env python
# -*- coding: utf-8 -*-

# #################################################
# imports from the std-library

import os, sys, platform
import shutil  # for copying files and folders
import errno  # for error/exception handling
import threading  # for parallelism
import subprocess  # for calling other commands
import re  # for regular expressions
from abc import ABCMeta, abstractmethod  # abstract classes
from argparse import ArgumentParser, RawTextHelpFormatter  # for parameters to this script
from collections import OrderedDict

# #################################################
# paths

__preparation_scripts_subfolder = "preparations"
__preparation_lib_subfolder = "lib"
__preparation_lib_srcml_subfolder = "srcml"

def getPreparationScript(filename):
    return os.path.join(__preparation_scripts_subfolder, filename)

def getLib(path):
    return os.path.abspath(os.path.join(__preparation_lib_subfolder, path))


# #################################################
# platform specific preliminaries

# cf. https://docs.python.org/2/library/sys.html#sys.platform
__platform = sys.platform.lower()

__iscygwin = False
if (__platform.startswith("linux")):
    __s2sml_executable = "linux/src2srcml"
    __sml2s_executable = "linux/srcml2src"
elif (__platform.startswith("darwin")):
    __s2sml_executable = "darwin/src2srcml"
    __sml2s_executable = "darwin/srcml2src"
elif (__platform.startswith("cygwin")) :
    __s2sml_executable = "win/src2srcml.exe"
    __sml2s_executable = "win/srcml2src.exe"
    __iscygwin = True
else:
    print("Your system '" + __platform + "' is not supported by SrcML right now.")

_s2sml = getLib(os.path.join(__preparation_lib_srcml_subfolder, __s2sml_executable))
_sml2s = getLib(os.path.join(__preparation_lib_srcml_subfolder, __sml2s_executable))


# #################################################
# imports from subfolders

import cppstats, cli

# for rewriting of #ifdefs to "if defined(..)"
# for turning multiline macros to oneliners
# for deletion of include guards in H files
from preparations import rewriteIfdefs, rewriteMultilineMacros, deleteIncludeGuards

from lib.cpplib import cpplib

# #################################################
# global constants

_filepattern_c = ('.c', '.C')
_filepattern_h = ('.h', '.H')
_filepattern = _filepattern_c + _filepattern_h


# #################################################
# helper functions

def notify(message):
    if (__iscygwin):
        return

    # import pynotify  # for system notifications
    #
    # pynotify.init("cppstats")
    # notice = pynotify.Notification(message)
    # notice.show()


# function for ignore pattern
def filterForFiles(dirpath, contents, pattern=_filepattern):
    mylist = [filename for filename in contents if
              not filename.endswith(pattern) and
              not os.path.isdir(os.path.join(dirpath, filename))
    ]
    return mylist


def runBashCommand(command, shell=False, stdin=None, stdout=None):
    # split command if not a list/tuple is given already
    if type(command) is str:
        command = command.split()

    process = subprocess.Popen(command, shell=shell, stdin=stdin, stdout=stdout, stderr=stdout)
    out, err = process.communicate() # TODO do something with the output
    process.wait()

    # FIXME do something with return value of process.wait()!
    # if ret is not 0:
    #     print "#### " + " ".join(command) + " returned " + str(ret)


def replaceMultiplePatterns(replacements, infile, outfile):
    with open(infile, "rb") as source:
        with open(outfile, "w") as target:
            data = source.read()
            for pattern, replacement in replacements.items():
                data = re.sub(pattern, replacement, data, flags=re.MULTILINE)
            target.write(data)


def stripEmptyLinesFromFile(infile, outfile):
    with open(infile, "rb") as source:
        with open(outfile, "w") as target:
            for line in source:
                if line.strip():
                    target.write(line)


def silentlyRemoveFile(filename):
    try:
        os.remove(filename)
    except OSError as e:  # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occured


def getCygwinPath(filename):
    return subprocess.check_output(['cygpath', '-m', filename]).strip()

def src2srcml(src, srcml):
    global _s2sml

    if (__iscygwin):
        src = getCygwinPath(src)
        #srcml = getCygwinPath(srcml)
        _s2sml = getCygwinPath(_s2sml)

    runBashCommand([_s2sml, src, "--language=C"], stdout = open(srcml, 'w+'))# + " -o " + srcml)
    # FIXME incorporate "|| rm ${f}.xml" from bash


def srcml2src(srcml, src):

    if (__iscygwin) :
        global _sml2s
        src = getCygwinPath(src)
        srcml = getCygwinPath(srcml)
        _sml2s = getCygwinPath(_sml2s)

    runBashCommand([_sml2s, srcml], stdout = open(src, 'w+'))# + " -o " + src)


# #################################################
# abstract preparation thread

class AbstractPreparationThread(object, metaclass=ABCMeta):
    '''This class prepares a single folder according to the given kind of preparations in an independent thread.'''
    sourcefolder = "source"

    def __init__(self, options, inputfolder=None, inputfile=None):
        self.options = options
        self.notrunnable = False

        if (inputfolder):
            self.file = None
            self.folder = inputfolder
            self.source = os.path.join(self.folder, self.sourcefolder)

            self.project = os.path.basename(self.folder)

            # get full path of subfolder "_cppstats"
            self.subfolder = os.path.join(self.folder, self.getSubfolder())

        elif (inputfile):
            self.file = inputfile
            self.outfile = self.options.outfile
            self.folder = os.path.dirname(self.file)

            self.project = os.path.basename(self.file)

            # get full path of temp folder for
            import tempfile
            self.subfolder = tempfile.mkdtemp(suffix=self.getSubfolder())


        else:
            self.notrunnable = True

    def startup(self):
        # LOGGING
        notify("starting '" + self.getPreparationName() + "' preparations:\n " + self.project)
        print("# starting '" + self.getPreparationName() + "' preparations: " + self.project)

    def teardown(self):

        # delete temp folder for file-based preparation
        if (self.file):
            shutil.rmtree(self.subfolder)

        # LOGGING
        notify("finished '" + self.getPreparationName() + "' preparations:\n " + self.project)
        print("# finished '" + self.getPreparationName() + "' preparations: " + self.project)

    def run(self):

        if (self.notrunnable):
            print("ERROR: No single file or input list of projects given!")
            return

        self.startup()

        if (self.file):

            self.currentFile = os.path.join(self.subfolder, self.project)
            shutil.copyfile(self.file, self.currentFile)

            self.backupCounter = 0
            self.prepareFile()

            shutil.copyfile(self.currentFile + ".xml", self.outfile)
        else:
            # copy C and H files to self.subfolder
            self.copyToSubfolder()
            # preparation for all files in the self.subfolder (only C and H files)
            for root, subFolders, files in os.walk(self.subfolder):
                for file in files:
                    f = os.path.join(root, file)
                    self.currentFile = f

                    self.backupCounter = 0
                    self.prepareFile()

        self.teardown()

    def copyToSubfolder(self):

        # TODO debug
        # echo '### preparing sources ...'
        # echo '### copying all-files to one folder ...'

        # delete folder if already existing
        if os.path.isdir(self.subfolder):
            shutil.rmtree(self.subfolder)

        # copy all C and H files recursively to the subfolder
        shutil.copytree(self.source, self.subfolder, ignore=filterForFiles)

    def backupCurrentFile(self):
        '''# backup file'''
        if (not self.options.nobak):
            bak = self.currentFile + ".bak" + str(self.backupCounter)
            shutil.copyfile(self.currentFile, bak)
            self.backupCounter += 1

    @classmethod
    @abstractmethod
    def getPreparationName(cls):
        pass

    @abstractmethod
    def getSubfolder(self):
        pass

    @abstractmethod
    def prepareFile(self):
        pass

    # TODO refactor such that file has not be opened several times! (__currentfile)
    def rewriteMultilineMacros(self):
        tmp = self.currentFile + "tmp.txt"

        self.backupCurrentFile()  # backup file

        # turn multiline macros to oneliners
        shutil.move(self.currentFile, tmp)  # move for script
        rewriteMultilineMacros.translate(tmp, self.currentFile)  # call function

        os.remove(tmp)  # remove temp file

    def formatCode(self):
        tmp = self.currentFile + "tmp.txt"

        self.backupCurrentFile()  # backup file

        # call astyle to format file in Java-style
        shutil.move(self.currentFile, tmp)  # move for script
        runBashCommand(["astyle", "--style=java"], stdin=open(tmp, 'r'), stdout=open(self.currentFile, 'w+'))

        os.remove(tmp)  # remove temp file

    def deleteComments(self):
        tmp = self.currentFile + "tmp.xml"
        tmp_out = self.currentFile + "tmp_out.xml"

        self.backupCurrentFile()  # backup file

        # call src2srcml to transform code to xml
        src2srcml(self.currentFile, tmp)

        # delete all comments in the xml and write to another file
        runBashCommand(["xsltproc", getPreparationScript("deleteComments.xsl"), tmp], stdout=open(tmp_out, 'w+'))

        # re-transform the xml to a normal source file
        srcml2src(tmp_out, self.currentFile)

        # delete temp files
        silentlyRemoveFile(tmp)
        silentlyRemoveFile(tmp_out)

    def deleteWhitespace(self):
        """deletes leading, trailing and inter (# ... if) whitespaces,
        replaces multiple whitespace with a single space"""
        tmp = self.currentFile + "tmp.txt"

        self.backupCurrentFile()  # backup file

        # replace patterns with replacements
        replacements = {
            '^[ \t]+': '',  # leading whitespaces
            '[ \t]+$': '',  # trailing whitespaces
            '#[ \t]+': '#',  # inter (# ... if) whitespaces # TODO '^#[ \t]+' or '#[ \t]+'
            '\t': ' ',  # tab to space
            '[ \t]{2,}': ' '  # multiple whitespace to one space

        }
        replaceMultiplePatterns(replacements, self.currentFile, tmp)

        # move temp file to output file
        shutil.move(tmp, self.currentFile)

    def rewriteIfdefsAndIfndefs(self):
        tmp = self.currentFile + "tmp.txt"

        self.backupCurrentFile()  # backup file

        # rewrite #if(n)def ... to #if (!)defined(...)
        d = rewriteIfdefs.rewriteFile(self.currentFile, open(tmp, 'w'))

        # move temp file to output file
        shutil.move(tmp, self.currentFile)

    def removeIncludeGuards(self):
        # include guards only exist in H files, otherwise return
        _, extension = os.path.splitext(self.currentFile)
        if (extension not in _filepattern_h):
            return

        tmp = self.currentFile + "tmp.txt"

        self.backupCurrentFile()  # backup file

        # delete include guards
        deleteIncludeGuards.apply(self.currentFile, open(tmp, 'w'))

        # move temp file to output file
        shutil.move(tmp, self.currentFile)

    def removeOtherPreprocessor(self):
        tmp = self.currentFile + "tmp.txt"

        self.backupCurrentFile()  # backup file

        # delete other preprocessor statements than #ifdefs
        cpplib._filterAnnotatedIfdefs(self.currentFile, tmp)

        # move temp file to output file
        shutil.copyfile(tmp, self.currentFile)

    def deleteEmptyLines(self):
        tmp = self.currentFile + "tmp.txt"

        self.backupCurrentFile()  # backup file

        # remove empty lines
        stripEmptyLinesFromFile(self.currentFile, tmp)

        # move temp file to output file
        shutil.move(tmp, self.currentFile)

    def transformFileToSrcml(self):
        source = self.currentFile
        dest = self.currentFile + ".xml"

        # transform to srcml
        src2srcml(source, dest)


# #################################################
# preparation-thread implementations

class GeneralPreparationThread(AbstractPreparationThread):
    @classmethod
    def getPreparationName(cls):
        return "general"

    def getSubfolder(self):
        return "_cppstats"

    def prepareFile(self):
        # multiline macros
        self.rewriteMultilineMacros()

        # delete comments
        self.deleteComments()

        # delete leading, trailing and inter (# ... if) whitespaces
        self.deleteWhitespace()

        # rewrite #if(n)def ... to #if (!)defined(...)
        self.rewriteIfdefsAndIfndefs()

        # removes include guards from H files
        self.removeIncludeGuards()

        # delete empty lines
        self.deleteEmptyLines()

        # transform file to srcml
        self.transformFileToSrcml()


class DisciplinePreparationThread(AbstractPreparationThread):
    @classmethod
    def getPreparationName(cls):
        return "discipline"

    def getSubfolder(self):
        return "_cppstats_discipline"

    def prepareFile(self):
        # multiline macros
        self.rewriteMultilineMacros()

        # delete comments
        self.deleteComments()

        # delete leading, trailing and inter (# ... if) whitespaces
        self.deleteWhitespace()

        # rewrite #if(n)def ... to #if (!)defined(...)
        self.rewriteIfdefsAndIfndefs()

        # removes include guards from H files
        self.removeIncludeGuards()

        # removes other preprocessor than #ifdefs
        self.removeOtherPreprocessor()

        # delete empty lines
        self.deleteEmptyLines()

        # transform file to srcml
        self.transformFileToSrcml()


class FeatureLocationsPreparationThread(AbstractPreparationThread):
    @classmethod
    def getPreparationName(cls):
        return "featurelocations"

    def getSubfolder(self):
        return "_cppstats_featurelocations"

    def prepareFile(self):
        # multiline macros
        self.rewriteMultilineMacros()

        # delete comments
        self.deleteComments()

        # delete leading, trailing and inter (# ... if) whitespaces
        self.deleteWhitespace()

        # rewrite #if(n)def ... to #if (!)defined(...)
        self.rewriteIfdefsAndIfndefs()

        # transform file to srcml
        self.transformFileToSrcml()


class PrettyPreparationThread(AbstractPreparationThread):
    @classmethod
    def getPreparationName(cls):
        return "pretty"

    def getSubfolder(self):
        return "_cppstats_pretty"

    def prepareFile(self):
        # multiline macros
        self.rewriteMultilineMacros()

        # format the code
        self.formatCode()

        # # delete comments
        # self.deleteComments()
        #
        # # delete empty lines
        # self.deleteEmptyLines()


# #################################################
# collection of preparation threads

# add all subclass of AbstractPreparationThread as available preparation kinds
__preparationkinds = []
for cls in AbstractPreparationThread.__subclasses__():
    entry = (cls.getPreparationName(), cls)
    __preparationkinds.append(entry)

# exit, if there are no preparation threads available
if (len(__preparationkinds) == 0):
    print("ERROR: No preparation tasks found! Revert your changes or call the maintainer.")
    print("Exiting now...")
    sys.exit(1)
__preparationkinds = OrderedDict(__preparationkinds)

def getKinds():
    return __preparationkinds


# #################################################
# main method


def applyFile(kind, inputfile, options):
    kinds = getKinds()

    # get proper preparation thread and call it
    threadClass = kinds[kind]
    thread = threadClass(options, inputfile=inputfile)
    thread.run()


def getFoldersFromInputListFile(inputlist):
    ''' This method reads the given inputfile line-wise and returns the read lines without line breaks.'''

    file = open(inputlist, 'r')  # open input file
    folders = file.read().splitlines()  # read lines from file without line breaks
    file.close()  # close file

    folders = [f for f in folders if not f.startswith("#")]  # remove commented lines
    folders = list(filter(os.path.isdir, folders))  # remove all non-directories
    folders = list(map(os.path.normpath, folders)) # normalize paths for easier transformations

    return folders


def applyFolders(kind, inputlist, options):
    kinds = getKinds()

    # get the list of projects/folders to process
    folders = getFoldersFromInputListFile(inputlist)

    # for each folder:
    for folder in folders:
        # start preparations for this single folder

        # get proper preparation thread and call it
        threadClass = kinds[kind]
        thread = threadClass(options, inputfolder=folder)
        thread.run()


def applyFoldersAll(inputlist, options):
    kinds = getKinds()
    for kind in list(kinds.keys()):
        applyFolders(kind, inputlist, options)


if __name__ == '__main__':
    kinds = getKinds()

    # #################################################
    # options parsing

    options = cli.getOptions(kinds, step = cli.steps.PREPARATION)

    # #################################################
    # main

    if (options.inputfile):

        # split --file argument
        options.infile = os.path.normpath(os.path.abspath(options.inputfile[0])) # IN
        options.outfile = os.path.normpath(os.path.abspath(options.inputfile[1])) # OUT

        # check if inputfile exists
        if (not os.path.isfile(options.infile)):
            print("ERROR: input file '{}' cannot be found!".format(options.infile))
            sys.exit(1)

        applyFile(options.kind, options.infile, options)

    elif (options.inputlist):
        # handle --list argument
        options.inputlist = os.path.normpath(os.path.abspath(options.inputlist)) # LIST

        # check if list file exists
        if (not os.path.isfile(options.inputlist)):
            print("ERROR: input file '{}' cannot be found!".format(options.inputlist))
            sys.exit(1)

        if (options.allkinds):
            applyFoldersAll(options.inputlist, options)
        else:
            applyFolders(options.kind, options.inputlist, options)

    else:
        print("This should not happen! No input file or list of projects given!")
        sys.exit(1)
