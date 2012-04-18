#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#USEFUL PACKAGES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import os,sys,commands,tempfile,time
import numpy as np
import commands
from pylab import *

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ALIASES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#plt=matplotlib.pyplot
exit=sys.exit
argv=sys.argv
system=os.system
sleep=time.sleep
tmpf=tempfile.NamedTemporaryFile
isfile=os.path.isfile
isdir=os.path.isdir

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#GLOBAL VARIABLES AND CONFIGURATION
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Do you want logs of the scripts
VERBOSE_LOG=True

#Numerical constants
MAXREAL=1E38

#!!!FROM HERE DO NOT TOUCH!!!
LOGSTART=False

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#CLASSES AND ROUTINES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class dictobj(object):
    """
    Class that allows the conversion from dictionary to class-like
    objects
    """
    def __init__(self,dic={}):self.__dict__.update(dic)
    def __add__(self,other):
        self.__dict__.update(other.__dict__)
        return self

def shiftList(a):
    b=a[::-1]
    try:
        val=list.pop(b)
    except IndexError:
        return False,[]
    a=[v for v in b[::-1]]
    return val,a

def shellExec(cmd,out=True,sim=False):
    """
    Execute a command
    """
    global FLOG
    if sim:return ""
    if not out:
        output=system(cmd)
        status=0
    else:
        status,output=commands.getstatusoutput(cmd)
        if VERBOSE_LOG and LOGSTART:
            FLOG.write("CMD:\n\t%s"%cmd)
            FLOG.write("OUTPUT:\n")
            FLOG.write(output)
            FLOG.write("\n")
            FLOG.write("STATUS:%d\n\n"%status)

    if status:
        print "Error:\n"
        print output
        exitScript(status)

    return output,status

def errorMsg(msg,usage="Usage:\n\t./script\n"):
    print "Error: ",msg
    print usage
    exitScript(0)

def exitScript(status):
    if VERBOSE_LOG:
        FLOG.close()
        print "Check the logfile '%s'"%LOGFILE
    exit(status)
    return 0

def getHeader(datafile,field):
    f=open(datafile,"r")
    for line in f:
        if "#" in line:
            if field in line:
                datos=line.strip().split(":")
                values=" ".join(datos[1:])
        else:
            break
    f.close()
    return values.strip().split(" ")

def gnuPlot(cmd):
    plotfile="plot.gpl"
    fp=open(plotfile,"w")
    fp.write(cmd);
    fp.close()
    shellExec("gnuplot %s"%plotfile)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ROUTINES RELATED TO DYNAMO MODULES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def filterFile(simfile,filfile,col,minval,maxval,equal='false'):
    pass

def transFile(simfile,trnfile,
              poscm="0,0,0",velcm="0,0,0",angmom="1E-6,1E-6,1",
              rflabel=0):
    pass

def scalarMap(simfile,mapfile,scalarmap,ngrid,
              coords='car',ini='-,-,-',end='-,-,-'):
    pass

def updateSim(simfile):
    pass

def statsFile(simfile,statfile,col,numbins,
              typebin='l',npb=100):
    pass

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#INFORMATION ABOUT THE PROGRAM
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM,argv=shiftList(argv)
RUNDIR,status=shellExec("dirname %s"%PROGRAM)
PROGRAM,status=shellExec("basename %s"%PROGRAM)
RUN="%s/../run "%RUNDIR
if VERBOSE_LOG:
    LOGSTART=True
    LOGFILE="%s.log"%PROGRAM
    FLOG=open(LOGFILE,"w")
