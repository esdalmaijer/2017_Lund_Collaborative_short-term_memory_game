import sys
import copy
import math
import random

from constants import *
from custom import calc_angular_dist, pol2car, StimScreen
from libclient import Client

from pygaze.display import Display
from pygaze.screen import Screen
from pygaze.sound import Sound
from pygaze.keyboard import Keyboard
from pygaze.mouse import Mouse
from pygaze.logfile import Logfile
from pygaze.eyetracker import EyeTracker
from pygaze import libtime as timer

from libcolour import create_colourwheel

import numpy


# # # # #
# INITIALISE

if DEBUG:
    print("Opening new DEBUG file.")
    debugfile = open(os.path.join(DEBUGDIR, 'experiment_%s.txt' % (STARTTIME)), 'w')

def _print(message):
    
    if DEBUG:
        debugfile.write("(%s %d): %s\n" % (time.strftime("%H:%M:%S"), time.time()*1000, message))

# Check if an IP was passed
if len(sys.argv) > 1:
    multicast_ip = sys.argv[1]
    multicast_ip = int(multicast_ip[1:])
    _print("Multicast IP = %s" % multicast_ip)
else:
    multicast_ip = 0
    _print("No multicast IP found.")

# GENERAL
# Total number of points.
total = 0.0

# Create a colour wheel to select colours from for the stimuli.
cw = create_colourwheel(STIML, STIMR, savefile='colourwheel.png')

# PYGAZE
# Initialise a new Display instance
disp = Display()

# Initialise a Screen instance for arbitrary drawing.
scr = Screen()
scr.draw_text(text="Initialising the experiment...", fontsize=FONTSIZE)
disp.fill(scr)
disp.show()

# Initialise the ka-ching sound.
kaching = Sound(soundfile=KACHING)

# Initialise a Keyboard and a Mouse instance for participant interaction.
kb = Keyboard()
mouse = Mouse()

# COMMUNICATIONS
timer.pause(5000)
_print("Initialising Client.")
# Initialise a new Client instance.
client = Client(multicast_ip)

# Establish a connection with the server.
scr.clear()
scr.draw_text(text="Connecting to the server...", fontsize=FONTSIZE)
disp.fill(scr)
disp.show()
_print("Connecting to the Server.")
success = client.contact_server(timeout=CONTIMEOUT)

# Get the experiment parameters.
scr.clear()
scr.draw_text(text="Getting experiment parameters...", fontsize=FONTSIZE)
disp.fill(scr)
disp.show()
_print("Getting experiment parameters.")
clientnr, clientcolour = client.get_parameters(timeout=CONTIMEOUT)

# Create a fixation screen.
fixscr = Screen()
fixscr.draw_fixation(fixtype=FIXTYPE, diameter=FIXSIZE)

# Create stimulus screens to accommodate all stimulus numbers.
stimscr = {}
respscr = {}
for nstim in NSTIM:
    locs = nstim * [DISPCENTRE]
    oris = nstim * [0]
    stimscr[nstim] = StimScreen(nstim, locs, oris, \
        linewidth=STIMLINEWIDTH, stimtypes='gabor', showcolourwheel=False)
    respscr[nstim] = StimScreen(nstim, locs, oris, \
        linewidth=STIMLINEWIDTH, stimtypes='noise', showcolourwheel=False)

# Create a feedback Screen for showing how much people earned.
fbscr = Screen()
fbscr.draw_rect(x=DISPCENTRE[0]-FBBARW//2, y=DISPCENTRE[1]-FBBARH//2, \
    w=FBBARW, h=FBBARH, pw=FBBARPW, fill=False)

# Create Screens for the post-run questions and the frame-by-frame updated
# fill of a bar.
qscr = Screen()
barqscr = Screen()

# LOGGING
# Create a logfile now that we know the client number
logpath = os.path.join(DATADIR, '%d_%s' % (clientnr, LOGFILENAME))
log = Logfile(filename=logpath+'_beh')
# Write a header to the log
header = ['trialnr', \
    'nstim', 'rewtype', 'rewamount', 'reward', 'total', \
    'fixonset', 'stimonset', 'maintenanceonset', 'probeonset', 't1', 'rewardonset', \
    'RT', 'nfixself', 'nfixother', 'probedstim', 'response']
for stimnr in range(max(NSTIM)):
    header.extend(['self_fixated%d' % stimnr, 'other_fixated%s' % stimnr, \
        'rewarded%d' % stimnr, 'ori%d' % stimnr, 'XYloc%d' % stimnr, \
        'error%d' % stimnr])
header.extend(['E', 'X', 'T'])
for i in range(max(NSTIM)-1):
    header.append('NT%d' % i)
log.write(header)

# Initialise the eye tracker.
_print("Initialising EyeTracker.")
tracker = EyeTracker(disp, logfile=logpath+'_eye')
#tracker.calibrate()


# # # # #
# RUN TRIALS

_print("Starting experiment. Clientnr=%d.\n" % (clientnr))

# Run until the server stops the experiment.
ntrials = 0
running = True
while running:
    
    # PREPARE TRIAL
    # Show the fixation screen.
    disp.fill(fixscr)
    fixonset = disp.show()

    # Update the trial counter.
    ntrials += 1
    # Get this trials' parameters.
    trial = client.get_trial(timeout=CONTIMEOUT)
    # Get the horizontal and vertical stimulus positions.
    stimx = []
    stimy = []
    for (x, y) in trial['locs']:
        stimx.append(x)
        stimy.append(y)
    stimx = numpy.array(stimx)
    stimy = numpy.array(stimy)
    # Construct the background colours.
    bgcols = trial['nstim'] * [BGC]
    
    # Update the stimulus screen.
    stimscr[trial['nstim']].update(trial['locs'], trial['oris'], \
        bgcols, linewidth=None)
    
    # Wait for the inter-trial interval minus the preparation time.
    timer.pause(ITI - (timer.get_time()-fixonset))

    
    # STIMULUS ARRAY
    # Start recording eye movements.
    tracker.start_recording()
    tracker.log("TRIALSTART")
    _print("Trial start.")

    # Show the stimulus screen.
    disp.fill(stimscr[trial['nstim']])
    stimonset = disp.show()
    tracker.log("stimulus_onset")
    _print("Stimulus onset.")
    
    # WAIT FOR FIXATIONS (both locally and from server)
    fixstim = None
    fixonset = None
    fixlist = {'self':[], 'other':[], 'all':[]}
    trialrunning = True
    fliptime = copy.deepcopy(stimonset)
    while trialrunning:

        # Get the current gaze position.
        gazepos = tracker.sample()
        if AUTORESP:
            if (fixonset != None) and (fixstim not in fixlist['all']):
                gazepos = (stimx[fixstim], stimy[fixstim])
            else:
                i = random.choice(range(len(trial['locs'])))
                gazepos = tuple(trial['locs'][i])
        _print("gazepos=(%s,%s)" % (gazepos))

        # Calculate the distance between the current gaze position and
        # all stimuli. But only if the sample is valid.
        if None not in gazepos:
            d = numpy.sqrt((gazepos[0]-stimx)**2 + (gazepos[1]-stimy)**2)
            _print("Distance to stimuli = %s" % (list(d)))
            # Do nothing if the current fixation is still on the same stimulus.
            if (numpy.min(d) < MAXFIXDIST) and (numpy.argmin(d) == fixstim):
                _print("Continuing fixation on stimulus %d" % (fixstim))
                pass
            # If the current fixation is on a different stimulus, update variables.
            elif (numpy.min(d) < MAXFIXDIST) and (numpy.argmin(d) != fixstim):
                fixstim = numpy.argmin(d)
                fixonset = timer.get_time()
                _print("New fixonset at %d on stimulus %d" % (fixstim, fixonset))
        # Check whether the fixation has lasted long enough to be counted, and
        # whether the currently fixated stimulus hasn't been fixated before.
        if (fixonset is not None) \
            and (timer.get_time() - fixonset >= MINFIXDUR) \
            and (fixstim not in fixlist['all']):
            _print("Fixated %d" % (fixstim))
            # Notify the server that a new fixation was made.
            ok = client.set_fixation(fixstim, timeout=CONTIMEOUT)
            # Add the stimulus to the fixation list, and set the background
            # colour to the right one.
            if ok:
                tracker.log("fixated_%d" % (fixstim))
                _print("Fixated stimulus %d accepted by server" % (fixstim))
                fixlist['self'].append(copy.deepcopy(fixstim))
                fixlist['all'].append(copy.deepcopy(fixstim))
                bgcols[fixstim] = STIMBGCOL['self']
        # If no fixation occurred, tell the server that nothing happened.
        else:
            _print("Fixated nothing")
            client.set_fixation(-1, timeout=CONTIMEOUT)
        # Ask the server for an update on the current fixations.
        timer.pause(5)
        _print("Asking for new fixations.")
        newfixations, trialrunning = client.get_fixations(timeout=CONTIMEOUT)
        _print("New fixations = %s" % newfixations)
        # Update the fixation lists and the stimulus background colours.
        if len(newfixations) > 0:
            for fixnr in newfixations:
                if fixnr not in fixlist['self']:
                    _print("Other fixated %d" % (fixnr))
                    fixlist['other'].append(copy.deepcopy(fixnr))
                    fixlist['all'].append(copy.deepcopy(fixnr))
                    bgcols[fixnr] = STIMBGCOL['other']

        # Update the stimulus screen every FRAMETIME ms.
#        if timer.get_time() - fliptime >= FRAMETIME:
        stimscr[trial['nstim']].update(trial['locs'], trial['oris'], \
            bgcols, linewidth=None)
        disp.fill(stimscr[trial['nstim']])
        fliptime = disp.show()

    # Post stimulus time.
    timer.pause(POSTSTIMTIME)
    
    # Maintenance phase
    disp.fill(fixscr)
    maintenanceonset = disp.show()
    tracker.log("stimulus_offset")
    _print("Stimulus offset")
    timer.pause(MAINTENANCETIME)
    
    # Stop recording eye movements.
    tracker.log("TRIALEND")
    timer.pause(5)
    tracker.stop_recording()
    
    
    # RESPONSE

    # Get the probed stimulus.
    probed = client.get_probe(timeout=CONTIMEOUT)
    _print("Probenr = %d" % (probed))

    # If no stimulus was fixated, the server will set probed to -1.
    if probed == -1:
        scr.clear()
        scr.draw_text(text="You didn't look at any stimulus.", \
            fontsize=FONTSIZE)
        disp.fill(scr)
        probeonset = disp.show()
        t1 = timer.get_time()
    
    # Update and show the response screen.
    else:
        probelw = trial['nstim'] * [STIMLINEWIDTH]
        probelw[probed] = PROBELINEWIDTH
        probeoris = trial['nstim'] * [0]
        probestimtypes = trial['nstim'] * ['noise']
        respscr[trial['nstim']].update(trial['locs'], probeoris, \
            bgcols, linewidth=probelw, stimtypes=probestimtypes)
        disp.fill(respscr[trial['nstim']])
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
        
        # Don't collect any responses when the participant didn't fixate a
        # stimulus.
        if probed == -1:
            break
        
        # update dial
        if respori != None:
            probeoris[probed] = respori
            probestimtypes[probed] = 'gabor'
            respscr[trial['nstim']].update(trial['locs'], probeoris, \
                bgcols, linewidth=probelw, stimtypes=probestimtypes)

        # show display
        disp.fill(respscr[trial['nstim']])
        t1 = disp.show()
        
        # check for keyboard input
        key, presstime = kb.get_key(timeout=1, flush=False)
        
        # Check if there is a timeout.
        if t1 - probeonset > RESPTIMEOUT:
            # If no response was given yet, choose the most wrong answer.
            if respori is None:
                respori = (trial['oris'][probed] + MAXORI/2.0) % MAXORI
            break
    
        # check mouse state
        ms = mouse.get_pressed()
        
        # check mouse position
        mpos = mouse.get_pos()
        
        # use mouse position to update orientation, but only if the left
        # mouse button is pressed
        if ms[0] == 1:

            # calculate distance between click and probed stimulus
            dist = numpy.sqrt((trial['locs'][probed][0] - mpos[0])**2 + \
                (trial['locs'][probed][1] - mpos[1])**2)
            
            # only continue if the click is on a stimulus or the surrounding colour wheel
            # Also make sure that the distance is not 0.
            if 0 < dist <= STIMSIZE:
                # correct the mouse position to a central location
                mpos = (\
                    mpos[0] + (DISPCENTRE[0]-trial['locs'][probed][0]), \
                    mpos[1] + (DISPCENTRE[1]-trial['locs'][probed][1]) \
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
            respori = random.randint(0, MAXORI-1)
            resp = True
            
    # hide mouse again
    mouse.set_visible(False)
    _print("Response = %s" % respori)
    
    # Send the response to the server, but only do so if a stimulus was probed.
    if probed >= 0:
        client.set_response(respori, timeout=CONTIMEOUT)
    
    # Show a fixation screen so that the participant knows their response has
    # registered.
    disp.fill(fixscr)
    disp.show()
    
    
    # TRIAL END
    
    # Get the trial outcome (error=int, reward=float, others=dict) from the 
    # server.
    error, reward, others = client.get_outcome(timeout=CONTIMEOUT)
    if probed == -1:
        error = numpy.NaN
        reward = 0.0
    _print("Error = %s, reward = %.2f" % (error, reward))
    _print("Others = %s" % (others))

    # Generate the reward text.
    groupsum = 0.0
    players = others.keys()
    players.sort()
    rewtxt = {}
    for k in players:
        if others[k]['probe'] == -1:
            groupsum += 0.0
            rewtxt[k] = "nothing"
        else:
            groupsum += others[k]['reward']
            rewtxt[k] = round(others[k]['reward'])
    players.insert(0, 'You')
    if probed == -1:
        rewtxt['You'] = "nothing"
    else:
        rewtxt['You'] = round(reward)
    groupsum += reward

    # Add the reward to the total.
    if trial['rewtype'] == 'collaborate':
        total += groupsum / (1.0 + len(others.keys()))
    else:
        total += reward

    # Play the reward sound.
    if reward >= REWSOUNDTHRESHOLD:
        kaching.play()
    
    # Show the rewards for all participants on the screen.
    if trial['rewtype'] == 'collaborate':
        # In the collaboration experiment, we only draw one bar.
        fbscr.clear()
        fbscr.draw_rect(x=DISPCENTRE[0]-FBBARW//2, y=DISPCENTRE[1]-FBBARH//2, \
            w=FBBARW, h=FBBARH, pw=FBBARPW, fill=False)
        # For every participant, add an extra bit of the bar.
        y = DISPCENTRE[1] + FBBARH//2 - FBBARPW
        for k in players:
            # Draw the bar fill.
            if k == 'You':
                col = STIMBGCOL['self']
            else:
                col = STIMBGCOL['other']
            if rewtxt[k] == "nothing":
                h = 1
            else:
                h = FBBARH * (rewtxt[k] / (len(players) * MAXREW))
            y -= h
            fbscr.draw_rect(x=DISPCENTRE[0]-FBBARW//2+FBBARPW, y=y, \
                w=FBBARW-FBBARPW*2, h=h, pw=0, fill=True, colour=col)
            # Draw the associated text.
            fbscr.draw_text(text="Player %s: %s" % (k, rewtxt[k]), \
                pos=(DISPCENTRE[0]+FBBARW, y+h//2), fontsize=FONTSIZE)
            rewardonset = disp.fill(fbscr)
            disp.show()
            timer.pause(500)
        # Show group total and share.
        s = "Group: %d; your share %d" % \
            (round(groupsum), round(groupsum / float(len(players))))
        fbscr.draw_text(text=s, fontsize=FONTSIZE, \
            pos=(DISPCENTRE[0], DISPCENTRE[1] - FBBARH//2 - 3*FONTSIZE))
        rewardonset = disp.fill(fbscr)
        disp.show()
        timer.pause(500)
    # TODO: In the competition experiment, show one bar per participant.
    elif trial['rewtype'] == 'compete':
        pass
    # Keep the reward screen on for a bit.
    timer.pause(REWTIME)
    
    # Log stuff
    outline = [ntrials, \
        trial['nstim'], trial['rewtype'], trial['rewamount'], reward, total, \
        fixonset, stimonset, maintenanceonset, probeonset, t1, rewardonset, \
        t1-probeonset, len(fixlist['self']), len(fixlist['other']), probed, respori]
    for stimnr in range(max(NSTIM)):
        if stimnr < trial['nstim']:
            if probed == -1:
                err = numpy.NaN
            else:
                err = calc_angular_dist(trial['oris'][stimnr], respori, zero=MAXORI)[0]
            outline.extend([int(stimnr in fixlist['self']), int(stimnr in fixlist['other']), \
                1, trial['oris'][stimnr], trial['locs'][stimnr], \
                err])
        else:
            outline.extend([numpy.NaN, numpy.NaN, \
                numpy.NaN, numpy.NaN, numpy.NaN, \
                numpy.NaN])
    nt = copy.deepcopy(trial['oris'])
    t = nt.pop(probed)
    outline.append(error)
    outline.append(respori)
    outline.append(t)
    for i in range(max(NSTIM)-1):
        if i < len(nt):
            outline.append(nt[i])
        else:
            outline.append(numpy.NaN)
    log.write(outline)
    
    # Wait for the server to continue the experiment.
    _print("Waiting for experiment to continue.")
    running = client.get_continue(timeout=BREAKTIMEOUT)
    _print("Continuing experiment.\n")
    
    # Allow for re-calibration, and send server 'ready' message when done.
    if running:
        scr.clear()
        scr.draw_text(text="Press Space to continue.", fontsize=FONTSIZE)
        disp.fill(scr)
        disp.show()
        if not AUTORESP:
            key, presstime = kb.get_key(keylist=['space','q'], timeout=5000, \
                flush=True)
        if key == 'q':
            tracker.calibrate()
        else:
            disp.fill(fixscr)
            disp.show()
        _print("Telling Server I'm ready.")
        client.set_ready()

# Close the client connection.
client.close()


# # # # #
# QUESTION TIME

# Ask about the collaboration, the other's likability, and their fairness.
questions = [ \
    "How good was the collaboration between you and player %s?", \
    "How likable did you find player %s?", \
    "How fair did you think player %s was?", \
    ]

# Get all the player names.
players = others.keys()
players.sort()

# Create a log file for the post-run questions.
qlog = Logfile(filename=logpath+'_que')
header = ['player', 'question', 'proportion']
qlog.write(header)

# Loop through all the players.
for name in players:
    
    # Loop through all the questions.
    for question in questions:

        # Create a Screen for drawing the post-game questions on.
        qscr.clear()
        qscr.draw_text(question % (name), pos=POSTQPOS, fontsize=POSTQFONTSIZE)
        qscr.draw_rect( \
            x=DISPCENTRE[0]-POSTQBARW//2, \
            y=DISPCENTRE[1]-POSTQBARH//2, \
            w=POSTQBARW, \
            h=POSTQBARH, \
            pw=POSTQBARPW, \
            fill=False)
        
        # Run until the mouse is pressed.
        pressed = False
        while not pressed:
            # Get the current mouse position.
            mpos = mouse.get_pos()
            # Calculate the bar fill.
            if mpos[0] < DISPCENTRE[0] - POSTQBARW//2 + POSTQBARPW//2:
                fillp = 0.0
            elif mpos[0] > DISPCENTRE[0] + POSTQBARW//2 - POSTQBARPW//2:
                fillp = 1.0
            else:
                fillp = (mpos[0] - (DISPCENTRE[0] - POSTQBARW//2 + POSTQBARPW//2)) \
                    / float(POSTQBARW-POSTQBARPW)
            # Draw the new bar fill.
            barqscr.clear()
            barqscr.copy(qscr)
            barqscr.draw_rect( \
                colour=POSTQBARCOL, \
                x=DISPCENTRE[0] - POSTQBARW//2 + POSTQBARPW//2, \
                y=DISPCENTRE[1]-POSTQBARH//2, \
                w=max(1, int((POSTQBARW-POSTQBARPW)*fillp)), \
                h=POSTQBARH-POSTQBARPW, \
                fill=True)
            disp.fill(barqscr)
            disp.show()
            # Check if the left mouse button is pressed.
            ms = mouse.get_pressed()
            # Stop if the left mouse button was pressed.
            if ms[0]:
                pressed = True
            if AUTORESP:
                fillp = random.random()
                pressed = True
        
        # Log the response.
        qlog.write([name, question, fillp])
        
        # Short break between the questions.
        disp.fill()
        disp.show()
        timer.pause(1000)


# # # # #
# CLOSE DOWN

_print("Closing experiment.")

# Calculate payment.
maxpay = (ntrials/3) * MAXREW
minpay = maxpay * MINPERFORMANCE
payrate = MAXADDPAY / float(maxpay - minpay) # pounds / credits
if total > minpay:
    pay = BASEPAY + (total-minpay) * payrate
else:
    pay = BASEPAY
_print("\n\nPAY = %.2f pounds per hour\n\n" % pay)

# On-screen message to indicate stuff is happening.
scr.clear()
scr.draw_text(text="Closing the experiment...", fontsize=FONTSIZE)
disp.fill(scr)
disp.show()

# Close the log file.
log.close()

# Close the connection to the eye tracker.
tracker.close()

# On-screen ending message.
scr.clear()
scr.draw_text(text="You earned %d credits!\n\nThanks for participating!" \
    % (total), fontsize=FONTSIZE)
disp.fill(scr)
disp.show()
if not AUTORESP:
    kb.get_key(keylist=None, timeout=None, flush=True)

# Close the Display.
disp.close()
