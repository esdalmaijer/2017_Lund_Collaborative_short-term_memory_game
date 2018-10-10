import copy
import time
import random

from constants import *
from custom import calc_angular_dist, generate_locs, generate_oris
from libserver import Server

from libcolour import create_colourwheel

import numpy
from scipy import stats


# # # # #
# INITIALISE

# COMMUNICATIONS
# Create a new Server instance.
print("Starting a new Server.")
server = Server()

# Establish a connection with all the clients we're expecting.
print("\nWaiting for clients...")
success = server.wait_for_clients(NPARTICIPANTS, timeout=None)

# All clients expect their own number, and their own colour.
for i in range(NPARTICIPANTS):
    col = (numpy.random.randint(0, 256), numpy.random.randint(0, 256), \
        numpy.random.randint(0, 256))
    print("\tClient %d, colour %s." % (i, col))
    server.set_parameters(i, col)

# Create a normal distribution or scaling rewards to
x = range(-180, 181)
y = stats.norm.pdf(x, loc=0, scale=REWSD)
y *= 1.0 / y.max()
rewardscale = y[180:]


# # # # #
# RANDOMISATION

print("\nRandomising trials")

# Each trial is a dict that has the following keys:
#   'rewtype'   Indicates the reward type. Should be one of the following
#               strings: 'collaborate', 'compete'.
#   'rewamount' Indicates how the amount of rewarded credits should be
#               calculated. Should be one of these strings: 'incidental',
#               'contingent'.
#   'nstim'     The number of stimuli in this trial's stimulus array. Should
#               be an integer.
#   'locs'      The locations of this trial's stimuli. Should be a list with
#               len=nstim, and an (x,y) tuple at each index.
#   'oris'      The orientations of this trial's stimuli. Should be a list
#               with len=nstim, with an int between 0 and MAXORI at every index.

# Create a list of all unique trials.
utrials = []
# Loop through all possible combinations of conditions.
for nstim in NSTIM:
    for rewtype in REWTYPES:
        t = {  \
            'nstim':nstim, \
            'rewtype':rewtype, \
            'rewamount':REWAMOUNT, \
            'oris':generate_oris(nstim, MINORIDIST, zero=MAXORI), \
            'locs':generate_locs(nstim, MINLOCDIST, edgeprop=0.1)
            }
        utrials.append(t)
# Get the number of unique trials
nutrials = len(utrials)
# Calculate the number of times each unique trial should be repeated
nreps = NTRIALS / nutrials
# Generate a list of trials
trials = nreps * utrials
# Shuffle the list of all trials
random.shuffle(trials)


# # # # #
# RUN EXPERIMENT

print("\nStarting trials!")

# Loop through all trials.
for trialnr, trial in enumerate(trials):
    
    
    # PREPARE TRIAL
    
    # Send the new trial parameters.
    print("\tBroadcasting parameters for trial %d" % (trialnr))
    server.set_trial(trial)
    time.sleep(0.05)

    
    # STIMULUS ARRAY

    print("\nPresenting stimuli")

    # Process incoming fixations.
    t0 = time.time() * 1000.0
    t1 = time.time() * 1000.0
    fixlist = {'all':[]}
    for clientnr in range(NPARTICIPANTS):
        fixlist[clientnr] = []
    trialrunning = True
    while trialrunning:
        # Check for new incoming fixations from all clients.
        newfixlist = []
        for clientnr in range(NPARTICIPANTS):
            # Get this client's new fixations, and confirm any new ones.
            # (New fixations are checked against the list of all previous
            # fixations, and clients receive a confirmation if the fixated
            # stimulus they reported is new.)
            fixstim = server.get_fixation(clientnr, fixlist['all'], \
                timeout=CONTIMEOUT)
            # Store any new fixations.
            if fixstim >= 0:
                print("\tClient %d fixated stimulus %d" % (clientnr, fixstim))
                fixlist['all'].append(copy.deepcopy(fixstim))
                fixlist[clientnr].append(copy.deepcopy(fixstim))
                newfixlist.append(copy.deepcopy(fixstim))
                    
        # Check whether the trial should still continue.
        if STIMTIME == None:
            # Check if all stimuli are fixated yet.
            if len(fixlist['all']) >= trial['nstim']:
                trialrunning = False
        else:
            # Check the time.
            if t1 - t0 > STIMTIME:
                trialrunning = False
        # Broadcast any new fixations, and whether the trial is still running.
        server.set_fixations(newfixlist, trialrunning, timeout=CONTIMEOUT)

    
    # RESPONSE

    print("\nChoosing probes")

    # Choose a stimulus to be probed for every client.
    time.sleep(0.05)
    probes = []
    for clientnr in range(NPARTICIPANTS):
        # If the participant didn't fixate anything, we can't probe anything!
        if len(fixlist[clientnr]) == 0:
            probed = -1

        else:
            # Choose a random stimulus to be probed.
            probed = random.choice(fixlist[clientnr])
            
        # Send the probed stimulus number to the client.
        print("\tClient %d will be probed on stimulus %d" % (clientnr, probed))
        server.set_probe(clientnr, probed, timeout=CONTIMEOUT)
        probes.append(copy.deepcopy(probed))
    time.sleep(0.05)
    
    # Wait for all clients to report the response orientations.
    print("\nAwaiting responses")
    responses = []
    for clientnr in range(NPARTICIPANTS):
        # Only clients that had a probed stimulus will report a response.
        if probes[clientnr] >= 0:
            respori = server.get_response(clientnr, timeout=CONTIMEOUT)
            print("\tClient %d responded %d" % (clientnr, respori))
        # Clients that did not have a probed stimulus will convert the
        # response and the error into NaNs, and the reward to 0.0
        else:
            respori = 0
        responses.append(copy.deepcopy(respori))
    time.sleep(0.05)
    
    
    # TRIAL END
    
    # Compute all the errors and rewards.
    errors = []
    rewards = []
    others = {}
    for clientnr in range(NPARTICIPANTS):
        # Calculate the error.
        error = calc_angular_dist(trial['oris'][probed], responses[clientnr], zero=MAXORI)
        errors.append(error)
        # Calculate the reward.
        reward  = 0
        if trial['rewamount'] == 'incidental':
            reward = 10
        elif trial['rewamount'] == 'contingent':
            reward = (rewardscale[int(abs(error))] * (MAXREW - MINREW)) + MINREW
        rewards.append(copy.deepcopy(reward))
        # Store the rewards of all participants in the 'others' dict. We will
        # use this to select specific participants from a bit later.
        others[clientnr+1] = copy.deepcopy(reward)
    
    # Send the info to the clients.
    print("\nCalculating outcomes")
    for clientnr in range(NPARTICIPANTS):
        # Copy and modify the others dict.
        o = copy.deepcopy(others)
        o.pop(clientnr+1)
        # Send the info to the client.
        print("\tClient %d: error=%d, reward=%.2f" % (clientnr, errors[clientnr], rewards[clientnr]))
        server.set_outcome(clientnr, errors[clientnr], rewards[clientnr], o)
    time.sleep(0.05)
    
    # Tell all clients whether the experiment should continue or not.
    print("\nTelling clients continue=%s" % (trialnr < len(trials)-1))
    for clientnr in range(NPARTICIPANTS):
        server.set_continue(clientnr, trialnr < len(trials)-1)
    time.sleep(0.05)
    
    # Wait for all clients to signal that they are ready.
    if trialnr < len(trials)-1:
        print("\nWaiting for clients to get ready for the next trial.")
        for clientnr in range(NPARTICIPANTS):
            server.get_ready(clientnr)
            print("\tClient %d is ready!" % (clientnr))
        time.sleep(0.05)


# # # # #
# CLOSE DOWN

# Close the connection with all clients.
print("\nEnd of experiment. Closing the connections with all clients.")
server.close()
