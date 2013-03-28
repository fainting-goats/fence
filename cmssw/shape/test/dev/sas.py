#!/bin/env python
import optparse
import sys
import ROOT
import numpy
import re
import warnings
import os.path

# for trigger efficiency fits
from HWWAnalysis.ShapeAnalysis.hwwtools import confirm
from HWWAnalysis.ShapeAnalysis.ROOTtree import Event,Leaves

#   _______                 
#  / ___/ /__  ___  ___ ____
# / /__/ / _ \/ _ \/ -_) __/
# \___/_/\___/_//_/\__/_/   
#                           

class TreeCloner(object):
    def __init__(self):
        self.ifile = None
        self.itree = None
        self.ofile = None
        self.otree = None
        self.label = None
        self.event = None
    
    def _openRootFile(self,path, option=''):
        f =  ROOT.TFile.Open(path,option)
        if not f.__nonzero__() or not f.IsOpen():
            raise NameError('File '+path+' not open')
        return f

    def _getRootObj(self,d,name):
        o = d.Get(name)
        if not o.__nonzero__():
            raise NameError('Object '+name+' doesn\'t exist in '+d.GetName())
        return o

    def connect(self, tree, input):
        self.ifile = self._openRootFile(input)
        self.itree = self._getRootObj(self.ifile,tree)

    def clone(self,output,branches=[]):
        
        #reset event
        self.event = None

        self.ofile = self._openRootFile(output, 'recreate')

        for b in self.itree.GetListOfBranches():
            if b.GetName() not in branches: continue
            b.SetStatus(0)

        self.otree = self.itree.CloneTree(0)

        ## BUT keep all branches "active" in the old tree
        self.itree.SetBranchStatus('*'  ,1)

        # create and event on the output tree (only cloned branches)
        self.event(self.otree)
        self.event.link(self.itree)



    def disconnect(self):
        self.otree.Write()
        self.ofile.Close()
        self.ifile.Close()

        self.ifile = None
        self.itree = None
        self.ofile = None
        self.otree = None
        self.event = None

#    ___                    __   _____         _____         
#   / _ )_______ ____  ____/ /  / ___/______ _/ _/ /____ ____
#  / _  / __/ _ `/ _ \/ __/ _ \/ (_ / __/ _ `/ _/ __/ -_) __/
# /____/_/  \_,_/_//_/\__/_//_/\___/_/  \_,_/_/ \__/\__/_/   
#                                                            

class Grafter(TreeCloner):
    '''Adds or replace variables to the tree '''

    def __init__(self):
        self.variables = {}
        self.regex = re.compile('([a-zA-Z0-9]*)/(['+''.join(Leaves.flags)+'])=(.*)')

    def help(self):
        return '''Makes a copy of the original tree adding or replacing variables. (TODO: Add warning for order dependent effects)'''

    def addOptions(self,parser):
        description = self.help()
        group = optparse.OptionGroup(parser,self.label, description)
        group.add_option('-v','--var',dest='variables',action='append',default=[])
        parser.add_option_group(group)
        return group

    def configure(self,opts):
        if not opts.variables:
            raise ValueError('No variables defined?!?')

        for s in opts.variables:
            #TODO  make the type optional, F by default
            r = self.regex.match(s)
            if not r:
                raise RuntimeError('Malformed option '+s)
            name=r.group(1)
            type=r.group(2)
            formula=r.group(3)
#             if type=='F':
#                 numtype = numpy.float32
#             elif type=='D':
#                 numtype = numpy.float64
#             elif type=='I':
#                 numtype = numpy.int32
#             else:
#                 RuntimeError('Type '+type+' not supported')

            self.variables[name] = (formula, type)

    def process(self,**kwargs):

        tree  = kwargs['tree']
        input = kwargs['input']
        output = kwargs['output']

        self.connect(tree,input)

        formulas = [ ( ROOT.TTreeFormula(name,formula, self.itree), type) for name, (formula,type) in self.variables.iteritems() ]


        print 'Adding/replacing the following branches'
        template=' {0:10} | {1:^3} | "{2}"'
        for name  in sorted(self.variables):
            (formula, type) = self.variables[name]
            print template.format(name,type,formula)
        print

#         #TODO Add the capability to completely force a replacement of the variables
#         oldbranches = [ b.GetName() for b in self.itree.GetListOfBranches() ]
#         hasMutation = False
#         for bname in self.variables:
#             # not there, continue
#             if bname not in oldbranches: continue
#             # found, check for consistency
#             branch = self.itree.GetBranch(bname)
#             btype = self.variables[bname][1]
#             newtitle = bname+'/'+btype
#             if ( branch.GetTitle() != newtitle ):
#                 print 'WARNING: Branch mutation detected: from',branch.GetTitle(),'to',newtitle
#                 hasMutation = True

#         if hasMutation:
#             confirm('Mutation detected. Do you _really_ want to continue?') or sys.exit(0)

        # there is a risk for entangled variables. if y = f(x) and x is changed first, y will take the new value. Is it desired?
        # if not the solution is to separate the variables to be overwritten between the old and new tree
        self.clone(output,self.variables.keys())

        self.clone(output)
        e = self.event
        for (formula,type) in formulas:
            leaves = e.leaves()
            name = formula.GetName()
            if name not in leaves:
                self.event.add( (name,type) )
            elif type and type != e.leaftype(name) :
                # if type is defined, check it maches
                raise AttributeError('Overwrite requested to type \''+type+'\' for variable '+name+' originally \''+e.leaftype+'\'')
                

        nentries = self.itree.GetEntries()
        print 'Entries:',nentries

        # avoid dots in the loop
        itree = self.itree
        otree = self.otree

        step = 5000
        for i in xrange(0,nentries):
            itree.GetEntry(i)

            if i > 0 and i%step == 0:
                print str(i)+' events processed.'

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                for (formula,type) in formulas:
                    # event.name = formula
                    getattr(e,formula.GetName()) = formula.EvalInstance()

            otree.Fill()
        
        self.disconnect()
        print '- Eventloop completed'

def execute(module,tree,iofiles):
    nfiles=len(iofiles)
    for i,(ifile,ofile) in enumerate(iofiles):
        print '-'*80
        print 'Entry {0}/{1} | {2}'.format(i+1,nfiles,ifile)
        print '-'*80
        odir  = os.path.dirname(ofile)
        if odir and not os.path.exists(odir):
           os.system('mkdir -p '+odir)
           print file,ofile

        print 'Input: ',ifile
        print 'Output:',ofile
        print '-'*80
    
        module.process( input=ifile, output=ofile, tree=tree )



def cli_test():
    usage = '''
    Usage:
        %prog <command> <options> filein.root fileout.root
        %prog <command> <options> file1.root file2.root ... dirout
        %prog -r <command> <options> dirin dirout

    In the latter case the directory tree in dirin is rebuilt in dirout

    Valid commands:
        '''+', '.join(modules.keys()+['help'])+'''

    Type %prog <command> -h for the command specific help
    '''

    parser = optparse.OptionParser(usage)
    parser.add_option('-t','--tree',        dest='tree',                            default='latino')
    parser.add_option('-r','--recursive',   dest='recursive',action='store_true',   default=False)
    parser.add_option('-F','--force',       dest='force',action='store_true',       default=False)

    # some boring argument handling
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    modname = sys.argv[1]
    if modname.startswith('-'):
        parser.print_help()
        sys.exit(0)

    if modname == 'help':
        if len(sys.argv) == 2:
            parser.print_help()
            sys.exit(0)
        
        module = sys.argv[2]
        if module not in modules:
            print 'Help: command',module,'not known'
            print 'The available commands are',', '.join(modules.keys())
            sys.exit(0)
        
        print 'Help for module',module+':'
        modules[module].help()
        modules[module].addOptions(parser)
        parser.print_help()
        sys.exit(0)


    if modname not in modules:
        print 'Command',modname,'unknown'
        print 'The available commands are',modules.keys()
        sys.exit(0)

    module = modules[modname]
    group = module.addOptions(parser)

    sys.argv.remove(modname)

    (opt, args) = parser.parse_args()

    print opt,args

    sys.argv.append('-b')

    try:
        module.configure(opt)
    except Exception as e:
        print 'Error in modutle',module.label
        print e
        sys.exit(1)

    tree = opt.tree

    nargs = len(args)
    if nargs < 2:
        parser.error('Input and/or output files missing')


    # file1 file2 file3 > outdir
    elif nargs > 2:
        inputs = args[:-1]
        output = args[-1]

        print inputs
        
        iofiles = [ (f,os.path.join(output,os.path.basename(f))) for f in inputs ]

        execute( module, tree, iofiles )
    
    elif nargs == 2:
        input = args[0]
        output = args[1]

        # recursiveness here
        if os.path.isdir(input):
            if not opt.recursive:
                print input,'is a directory. Use -r to go recursive'
                sys.exit(0)

            if os.path.exists(output) and not os.path.isdir(output):
                print output,'exists and is not a directory!'
                sys.exit(0)

            fileList = []
            for root, subFolders, files in os.walk(input):
                for file in files:
                    if not file.endswith('.root'): continue
                    fileList.append(os.path.join(root,file))

            print 'The directory tree',input,'will be gardened and copied to',output
            print 'The following files will be copied:'
            print '\n'.join(fileList)
            print 'for a grand total of',len(fileList),'files'
            opt.force or ( confirm('Do you want to continue?') or sys.exit(0) )

            output = output if output[-1]=='/' else output+'/'
            iofiles = [ (f,os.path.join(output,os.path.basename(f))) for f in fileList ]

            execute( module, tree, iofiles )

        else:
            if os.path.exists(output) and os.path.isdir(output):
                output = os.path.join( output , os.path.basename(input) )
            execute(module,tree,[(input,output)])

if __name__ == '__main__':
    
    banner='''
     (                                      (                        
 )\ )          (       (           (    )\ )                     
(()/(        ) )\  (   )\          )\ )(()/(   )     (    ) (    
 /(_)) (  ( /(((_)))((((_)(  (    (()/( /(_)) (     ))\( /( )(   
(_))   )\ )(_))_ /((_)\ _ )\ )\ )  ((_)|_))   )\  '/((_)(_)|()\  
/ __| ((_|(_)_| (_)) (_)_\(_)(_/(  _| |/ __|_((_))(_))((_)_ ((_) 
\__ \/ _|/ _` | / -_) / _ \| ' \)) _` |\__ \ '  \() -_) _` | '_| 
|___/\__|\__,_|_\___|/_/ \_\_||_|\__,_||___/_|_|_|\___\__,_|_|   
                                                                 
'''
    print banner
    cli_test()
