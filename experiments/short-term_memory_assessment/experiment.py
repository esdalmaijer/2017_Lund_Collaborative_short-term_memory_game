import copy
import math
import random

from constants import *
from pygaze.display import Display
from pygaze.screen import Screen
from pygaze.keyboard import Keyboard
from pygaze.mouse import Mouse
from pygaze.logfile import Logfile
import pygaze.libtime as timer

from custom import calc_angular_dist, generate_oris, generate_locs, pol2car, StimScreen

import numpy


# # # # #
# PREPARATION

# Intialise the Display.
disp = Display()

# Initialise the basic input devices.
kb = Keyboard(keylist=None, timeout=None)
mouse = Mouse(mousebuttonlist=None, timeout=None)

# Initialise a log.
log = Logfile()
header = ['trialnr', 'nstim', 'fixonset', 'stimonset', 'maintenanceonset', \
    'probeonset', 'RT', 'response']
for i in range(max(NSTIM)):
    header.extend(['stimx%d' % (i), 'stimy%d' % (i), 'stimori%d' % (i), \
        'stimerror%d' % (i)])
header.extend(['E', 'X', 'T'])
for i in range(max(NSTIM)-1):
    header.append('NT%d' % i)
log.write(header)

# Initialise a blank Screen for ad-hoc drawing.
scr = Screen()

# Initialise a blank Screen.
blankscr = Screen()

# Initialise a fixation Screen.
fixscr = Screen()
fixscr.draw_fixation(fixtype=FIXTYPE, diameter=FIXSIZE, pw=FIXPW)

# Initialise stimulus and probe Screens.
stimscr = {}
probescr = {}
for nstim in NSTIM:
    locs = nstim * [DISPCENTRE]
    oris = nstim * [0]
    stimscr[nstim] = StimScreen(nstim, locs, oris, \
        linewidth=STIMLINEWIDTH, stimtypes='gabor', showcolourwheel=False)
    probescr[nstim] = StimScreen(nstim, locs, oris, \
        linewidth=STIMLINEWIDTH, stimtypes='noise', showcolourwheel=False)

# Randomise trials.
trials = []
for nstim in NSTIM:
    t = { 
        'nstim':nstim, \
        }
    trials.extend(TRIALSPERCELL * [t])
random.shuffle(trials)


# # # # #
# RUN TRIALS

# Run through all trials.
for trialnr, trial in enumerate(trials):
    
    # PREPARE
    
    # Randomly choose a probe.
    probed = random.randint(0, trial['nstim']-1)
    # Generate stimulus locations.
    stimlocs = generate_locs(trial['nstim'], MINLOCDIST)
    # Generate stimulus orientations.
    stimoris = generate_oris(trial['nstim'], MINORIDIST, zero=MAXORI)
    
    # Update the stimulus Screen.
    stimtypes = trial['nstim'] * ['gabor']
    stimscr[trial['nstim']].update(stimlocs, stimoris, \
        linewidth=None, stimtypes=stimtypes)
    
    # Update the probe Screen.
    probelw = trial['nstim'] * [STIMLINEWIDTH]
    probelw[probed] = PROBELINEWIDTH
    probeoris = trial['nstim'] * [0]
    probestimtypes = trial['nstim'] * ['noise']
    probescr[trial['nstim']].update(stimlocs, probeoris, \
        linewidth=probelw, stimtypes=probestimtypes)
    

    # RUN
    
    # Show the fixation Screen.
    disp.fill(fixscr)
    fixonset = disp.show()
    timer.pause(random.randint(FIXTIME[0], FIXTIME[1]))
    
    # Show the stimulus Screen.
    disp.fill(stimscr[trial['nstim']])
    stimonset = disp.show()
    timer.pause(STIMTIME)
    
    # Show a blank Screen.
    disp.fill(blankscr)
    maintenanceonset = disp.show()
    timer.pause(MAINTENANCETIME)
    
    # Show the probe Screen.
    disp.fill(probescr[trial['nstim']])
    probeonset = disp.show()
    
    # Get a response.
    # flush keyboard
    kb.get_key(keylist=None, timeout=1, flush=True)
    # show mouse position
    mouse.set_visible(True)
    # set mouse position to centre
    #mouse.set_pos(DISPCENTRE)
    # interactive dial
    resp = False
    respori = None
    while not resp and probed >= 0:
        
        # Update dial
        if respori != None:
            probeoris[probed] = respori
            probestimtypes[probed] = 'gabor'
            probescr[trial['nstim']].update(stimlocs, probeoris, \
                linewidth=probelw, stimtypes=probestimtypes)

        # show display
        disp.fill(probescr[trial['nstim']])
        t1 = disp.show()
        
        # check for keyboard input
        key, presstime = kb.get_key(timeout=1, flush=False)
    
        # check mouse state
        ms = mouse.get_pressed()
        
        # check mouse position
        mpos = mouse.get_pos()
        
        # use mouse position to update orientation, but only if the left
        # mouse button is pressed
        if ms[0] == 1:

            # calculate distance between click and probed stimulus
            dist = numpy.sqrt((stimlocs[probed][0] - mpos[0])**2 + \
                (stimlocs[probed][1] - mpos[1])**2)
            
            # only continue if the click is on a stimulus or the surrounding colour wheel
            # Also make sure that the distance is not 0.
            if 0 < dist <= STIMSIZE:
                # correct the mouse position to a central location
                mpos = (\
                    mpos[0] + (DISPCENTRE[0]-stimlocs[probed][0]), \
                    mpos[1] + (DISPCENTRE[1]-stimlocs[probed][1]) \
                    )
                # calculate the position of a neutral point (at 0 degrees) at
                # the same distance
                npos = pol2car((dist, 0))
                # calculate the angle between the two points
                pdistsquare = float(mpos[0] - npos[0])**2 + \
                    float(mpos[0] - npos[0])**2
                # law of cosines
                mang = math.degrees(math.acos((2 * dist**2 - pdistsquare) \
                    / (2 * dist**2)))
                # adjust value to match mouse position
                # (cosine rule does not work >90 degrees in this setup)
                if mpos[1] > DISPCENTRE[1]:
                    mang = 90.0 + (90.0 - mang)
                if mpos[0] < DISPCENTRE[0]:
                    mang = 180.0 + (180.0 - mang)

                # Correct if angle is over MAXORI degrees
                mang = mang % MAXORI
                
                # Store the current response in variables
                respori = int(mang)

        if key == 'space' and respori != None:
            resp = True
        elif key == 'escape' and ESCKILL:
            log.write(["DEBUGKILL"])
            log.close()
            tracker.close()
            disp.close()
            raise Exception("DEBUG KILL")
        
        # auto response if needed
        if AUTORESP:
            resp = True
            respori = random.randint(0, MAXORI-1)
            
    # hide mouse again
    mouse.set_visible(False)
    
    # Log the response.
    line = [trialnr, trial['nstim'], fixonset, stimonset, maintenanceonset, \
        probeonset, presstime-probeonset, respori]
    for i in range(max(NSTIM)):
        if i < trial['nstim']:
            line.extend([stimlocs[i][0], stimlocs[i][1], stimoris[i], \
                calc_angular_dist(stimoris[i], respori, zero=MAXORI)[0]])
        else:
            line.extend([numpy.NaN, numpy.NaN, numpy.NaN, \
                numpy.NaN])
    nts = copy.deepcopy(stimoris)
    t = nts.pop(probed)
    line.extend([calc_angular_dist(t, respori, zero=MAXORI)[0], \
        respori, t])
    for i in range(max(NSTIM) - 1):
        if i < trial['nstim'] - 1:
            line.append(nts[i])
        else:
            line.append(numpy.NaN)
    log.write(line)


# # # # #
# CLOSE

# Close the log file.
log.close()

# Present an exit message.
scr.draw_text(text="This is the end of the experiment.", fontsize=FONTSIZE)
disp.fill(scr)
disp.show()
kb.get_key(keylist=None, timeout=None, flush=True)

# Close the display.
disp.close()
