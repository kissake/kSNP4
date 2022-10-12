#!/usr/bin/python3

import sys

def clearLine():
    # Clear whatever might be already written on the current line 
    # to the left using backspaces.

    # Not enough backspaces means we don't erase all the text.  Too many
    # doesn't really hurt?

    clearstring = " \b\b" * 20
    print(clearstring * 10,end="")


def writeElapsed(seconds):
    # Convert number of seconds elapsed to hours, minutes and seconds,
    # and display this info formatted to the screen.

    hoursElapsed = int ( seconds / 60 / 60 )
    minutesElapsed = int( ( seconds - ( hoursElapsed * 60 * 60 ) ) / 60 ) 
    secondsElapsed = int( seconds - (hoursElapsed * 60 * 60) - ( minutesElapsed * 60 ) )

    print("[%02i:%02i:%02i]" % (hoursElapsed, minutesElapsed, secondsElapsed), end="")

def writeSeparator():
    print(" ", end="")

def writePhase(phaseNum, phaseTotal, phaseLabel):

    print("(%i/%i) %s"%(phaseNum, phaseTotal, phaseLabel), end="")


def updateLine(seconds, phaseNum, phaseTotal, phaseLabel, subphaseNum=None, subphaseTotal=None, subphaseLabel=None, overwrite=True):
    if overwrite:
        clearLine()
    writeElapsed(seconds)
    writeSeparator()
    writePhase(phaseNum, phaseTotal, phaseLabel)
    if subphaseTotal != None:
        writeSeparator()
        writePhase(subphaseNum, subphaseTotal, subphaseLabel)

    if not overwrite:
        print()



if __name__ == "__main__":

    overwrite = False
    sys.stderr.write(" ".join(sys.argv) + "\n")
    secondsDuration = int( sys.argv[1] )
    phaseTotal = int(sys.argv[2])
    phaseNumber = int(sys.argv[3])
    phaseName = sys.argv[4]
    if len(sys.argv) > 6:
        subphaseTotal = int(sys.argv[5])
        subphase = int(sys.argv[6])
        subphaseLabel = sys.argv[7]
        
        updateLine(secondsDuration, phaseNumber, phaseTotal, phaseName, subphase, subphaseTotal, subphaseLabel, overwrite=overwrite)

    else:
        
        updateLine(secondsDuration, phaseNumber, phaseTotal, phaseName, overwrite=overwrite)
