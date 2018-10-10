import copy
import time
import random
import socket

from constants import *
from custom import calc_angular_dist, generate_locs, generate_oris
from libserver import Server

import numpy
from scipy import stats


class PairServer:
    
    """Server for a pair of clients.
    """
    
    def __init__(self, ip_addresses):

        # # # # #
        # INITIALISE
        
        # COMMUNICATIONS
        # Create a new Server instance.
        print("Starting a new Server.")
        self.server = Server(ip_addresses)
        
        # Start the clients' experiment scripts.
        if not MANUALSTART:
            print("Starting client experiment scripts")
            sock = socket.socket(socket.AF_INET, # Internet
                                 socket.SOCK_DGRAM) # UDP    
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            for clientnr in ip_addresses:
                print("Starting client %d" % clientnr)
                sock.sendto(STARTCMD % (ip_addresses[0]), \
                    (CLIENIPBASE % (clientnr), UDP_TIM_CLIENT_PORT))
                time.sleep(5.0)
            sock.close()
        
        # Establish a connection with all the clients we're expecting.
        print("\nWaiting for clients...")
        success = self.server.wait_for_clients(timeout=None)
        
        # Throw an exception if we failed to find all clients.
        if not success:
            raise Exception("Could not connect to all expected clients (%s)" \
                % ip_addresses)
        # Store the IP addresses if we did manage.
        self._clients = copy.deepcopy(self.server._clientlist)
        
        # All clients expect their own number, and their own colour.
        for clientnr in self._clients:
            col = (numpy.random.randint(0, 256), numpy.random.randint(0, 256), \
                numpy.random.randint(0, 256))
            print("\tClient %d, colour %s." % (clientnr, col))
            self.server.set_parameters(clientnr, col)
        
        # Create a normal distribution or scaling rewards to
        x = range(-180, 181)
        y = stats.norm.pdf(x, loc=0, scale=REWSD)
        y *= 1.0 / y.max()
        self._rewardscale = y[180:]


    def prepare_block(self, nstim, rewtypes, ntrials):

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
        for ns in nstim:
            for rt in rewtypes:
                t = {  \
                    'nstim':ns, \
                    'rewtype':rt, \
                    'rewamount':REWAMOUNT, \
                    'oris':generate_oris(ns, MINORIDIST, zero=MAXORI), \
                    'locs':generate_locs(ns, MINLOCDIST, edgeprop=0.1)
                    }
                utrials.append(t)
        # Get the number of unique trials
        nutrials = len(utrials)
        # Calculate the number of times each unique trial should be repeated
        nreps = NTRIALS / nutrials
        # Generate a list of trials
        self._trials = nreps * utrials
        # Shuffle the list of all trials
        random.shuffle(self._trials)


    def run_block(self):

        # # # # #
        # RUN EXPERIMENT
        
        print("\nStarting trials!")
        
        # Loop through all trials.
        for trialnr, trial in enumerate(self._trials):
            
            
            # PREPARE TRIAL
            
            # Randomise this trials' stimulus orientations and locations.
            trial['oris'] = generate_oris(trial['nstim'], MINORIDIST, zero=MAXORI)
            trial['locs'] = generate_locs(trial['nstim'], MINLOCDIST, edgeprop=0.1)
            
            # Send the new trial parameters.
            print("\tBroadcasting parameters for trial %d" % (trialnr))
            self.server.set_trial(trial)
            time.sleep(0.05)
        
            
            # STIMULUS ARRAY
        
            print("\nPresenting stimuli")

            # Process incoming fixations.
            t0 = time.time() * 1000.0
            t1 = time.time() * 1000.0
            fixlist = {'all':[]}
            for clientnr in self._clients:
                fixlist[clientnr] = []
            trialrunning = True
            while trialrunning:
                # Check for new incoming fixations from all clients.
                newfixlist = []
                for clientnr in self._clients:
                    # Get this client's new fixations, and confirm any new ones.
                    # (New fixations are checked against the list of all previous
                    # fixations, and clients receive a confirmation if the fixated
                    # stimulus they reported is new.)
                    fixstim = self.server.get_fixation(clientnr, fixlist['all'], \
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
                self.server.set_fixations(newfixlist, trialrunning, timeout=CONTIMEOUT)
        
            
            # RESPONSE
        
            print("\nChoosing probes")
        
            # Choose a stimulus to be probed for every client.
            time.sleep(0.05)
            probes = {}
            for clientnr in self._clients:
                # If the participant didn't fixate anything, we can't probe anything!
                if len(fixlist[clientnr]) == 0:
                    probed = -1
        
                else:
                    # Choose a random stimulus to be probed.
                    probed = random.choice(fixlist[clientnr])
                    
                # Send the probed stimulus number to the client.
                print("\tClient %d will be probed on stimulus %d" % (clientnr, probed))
                self.server.set_probe(clientnr, probed, timeout=CONTIMEOUT)
                probes[clientnr] = copy.deepcopy(probed)
            time.sleep(0.05)
            
            # Wait for all clients to report the response orientations.
            print("\nAwaiting responses")
            responses = {}
            for clientnr in self._clients:
                # Only clients that had a probed stimulus will report a response.
                if probes[clientnr] >= 0:
                    respori = self.server.get_response(clientnr, timeout=CONTIMEOUT)
                    print("\tClient %d responded %d" % (clientnr, respori))
                # Clients that did not have a probed stimulus will convert the
                # response and the error into NaNs, and the reward to 0.0
                else:
                    respori = 0
                responses[clientnr] = copy.deepcopy(respori)
            time.sleep(0.05)
            
            
            # TRIAL END
            
            # Compute all the errors and rewards.
            errors = {}
            rewards = {}
            others = {}
            for clientnr in self._clients:
                # Calculate the error.
                if probes[clientnr] == -1:
                    error = (MAXORI / 2)
                else:
                    error = calc_angular_dist(trial['oris'][probes[clientnr]], \
                        responses[clientnr], zero=MAXORI)
                errors[clientnr] = copy.deepcopy(error)
                # Calculate the reward.
                if probes[clientnr] == -1:
                    reward  = 0.0
                else:
                    if trial['rewamount'] == 'incidental':
                        reward = 10.0
                    elif trial['rewamount'] == 'contingent':
                        reward = (self._rewardscale[int(abs(error))] * (MAXREW - MINREW)) + MINREW
                rewards[clientnr] = copy.deepcopy(reward)
                # Store the rewards of all participants in the 'others' dict. We will
                # use this to select specific participants from a bit later.
                others[clientnr] = {'probe':probes[clientnr], \
                    'error':int(error), 'reward':reward}
            
            # Send the info to the clients.
            print("\nCalculating outcomes")
            for clientnr in self._clients:
                # Copy and modify the others dict by removing this client's
                # own score form the temporary copy.
                o = copy.deepcopy(others)
                o.pop(clientnr)
                # Send the info to the client.
                print("\tClient %d: error=%d, reward=%.2f" % \
                    (clientnr, errors[clientnr], rewards[clientnr]))
                self.server.set_outcome(clientnr, errors[clientnr], \
                    rewards[clientnr], o)
            time.sleep(0.05)
            
            # Tell all clients whether the experiment should continue or not.
            print("\nTelling clients to continue=%s" % (trialnr < len(self._trials)-1))
            for clientnr in self._clients:
                self.server.set_continue(clientnr, trialnr < len(self._trials)-1)
            time.sleep(0.05)
            
            # Wait for all clients to signal that they are ready.
            if trialnr < len(self._trials)-1:
                print("\nWaiting for clients to get ready for the next trial.")
                for clientnr in self._clients:
                    self.server.get_ready(clientnr)
                    print("\tClient %d is ready!" % (clientnr))
                time.sleep(0.05)


    def close(self):

        # # # # #
        # CLOSE DOWN
        
        # Close the connection with all clients.
        print("\nClosing the connections with all clients.")
        self.server.close()
