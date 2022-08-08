#!/usr/bin/python3

import sys

def clearLine():
    # Clear whatever might be already written on the current line 
    # to the left using backspaces.

    # Not enough backspaces means we don't erase all the text.  Too many
    # doesn't really hurt?
    print("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",end="")
    print("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",end="")
    print("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",end="")
    print("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",end="")
    print("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",end="")


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


def updateLine(seconds, phaseNum, phaseTotal, phaseLabel, subphaseNum, subphaseTotal, subphaseLabel):
    clearLine()
    writeElapsed(seconds)
    writeSeparator()
    writePhase(phaseNum, phaseTotal, phaseLabel)
    writeSeparator()
    writePhase(subphaseNum, subphaseTotal, subphaseLabel)



if __name__ == "__main__":

    secondsDuration = int( sys.argv[1] )
    phaseNumber = int(sys.argv[2])
    phaseName = sys.argv[3]
    totalPhases = 10 # Move this to an argument the first time it needs to change.
    subphase = int(sys.argv[4])
    subphaseTotal = int(sys.argv[5])
    subphaseLabel = sys.argv[6]

    updateLine(secondsDuration, phaseNumber, totalPhases, phaseName, subphase, subphaseTotal, subphaseLabel)
