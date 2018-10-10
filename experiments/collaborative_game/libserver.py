from constants import *
import codebook as cb

import UDPClient

import os
import copy
import json
import time
from threading import Lock, Thread

import numpy


class Server:
    
    def __init__(self, ip_addresses):
        
        if DEBUG:
            if not os.path.isdir('DEBUG'):
                os.mkdir('DEBUG')
            self._debugfile = open(os.path.join(DEBUGDIR, \
                'server_debug_%s_%s.txt' % \
                (ip_addresses[0], STARTTIME)), 'w')
        
        # Set how many messages should be sent each time (Copies can be
        # included to make it more likely that all messages arrive. This is
        # a very lazy implementation. Don't judge me.)
        self._message_reps = 0
        
        # Who am I?
        self._servernr = int(SERVERIP.split('.')[-1])
        
        # Set up a dict for keeping track of the last sent messages of each
        # type, for each client.
        self._lastmessage = {'unknown':{}}
        
        # Start with an empty client list.
        self._unconfirmed_ip_addresses = ip_addresses
        self._clientlist = []
        # Start with empty outgoing and incoming dicts.
        self._outgoing = {}
        self._incoming = {'unknown':[]}
        # Set up a dict to store time differences in, and set the maximum
        # length (this will determine the recency and accuracy of the running
        # average of time differences). Each client will be assigned a list in
        # this dict, so we can keep track of all timing discrepancies.
        self._timelist = {}
        self._timelistlen = 5
        # And we need locks to prevent simultaneous access to said dicts.
        self._udplock = Lock()
        self._outlock = Lock()
        self._inlock = Lock()
        self._timelock = Lock()
        
        # Select the right group address.
        if MANUALSTART:
            group_address = 0
        else:
            group_address = ip_addresses[0]

        # Open a new UDP client
        self._print("Opening a new UDPClient")
        self.udp = UDPClient.UDPClient()
        self.udp.port = 4444
        self.udp.groupAddress = '225.0.0.%s' % (group_address)
        self.udp.loopBack = False
        self.udp.reuseSocket = 1
        self.udp.numReceiverThreads = 2
        self.udp.init()
        
        # Set the waiting time between checking for new incoming and outgoing
        # messages, to allow the other processes some time to run.
        # SECONDS!
        self._refresh_delay = 0.0
        # Set the grace period for the time difference between outgoing
        # messages and incoming responses.
        # CURRENTLY UNUSED
        self._graceperiod = 0
        # Set the time that samples should be kept for.
        # MILISECONDS
        self._cleantime = 1000
        # Time spent waiting for an expected message until a request is sent
        # to repeat it. IN SECONDS!
        self._reptimeout = 10.0
        
        # Start the message processing Threads.
        self._print("Starting new message processing Threads")
        self._outgoing_thread = Thread(target=self._outgoing_messages)
        self._outgoing_thread.name = "outgoing_messages"
        self._outgoing_thread.daemon = True
        self._outgoing_thread.start()
        self._incoming_thread = Thread(target=self._incoming_messages)
        self._incoming_thread.name = "incoming_messages"
        self._incoming_thread.daemon = True
        self._incoming_thread.start()
    
    
    def _print(self, message):
        
        if DEBUG:
            self._debugfile.write("(%d) %s\n" % (time.time()*1000, message))
            self._debugfile.flush() # internal buffer to RAM
            os.fsync(self._debugfile.fileno()) # RAM file cache to disk

    
    def _clean_incoming(self, clientnr, t0):
        
        """For internal use. Cleans the incoming queue by throwing out all
        messages that were received before the passed time.
        """
        
        # Obtain the incoming Lock, so that no concurrent access can occur.
        self._inlock.acquire()

        # Loop until all messages are removed.
        cleaning = True
        while cleaning:
            # Only go on if there are messages in the queue.
            if len(self._incoming[clientnr]) > 0:
                # Check if the next message was before the requested removal time.
                if self._incoming[clientnr][0][0] < t0:
                    t, m = self._incoming[clientnr].pop(0)
                    self._print("Removed message '%s' (%.3f) from the incoming queue" \
                        % (m, t))
                else:
                    cleaning = False
            else:
                cleaning = False
        
        # Release the incoming Lock, so that other Threads can do their thang.
        self._inlock.release()
    
    
    def _incoming_messages(self):
        
        """For internal use. Listens for new messages, and adds them to the
        cue with a timestamp of when they came in.
        """
        
        # Set up special cases in which a direct reply is sent back from this
        # thread.
        expected = {}
        expected['waitwhat'] = cb.WAITWHATSERVER[:cb.WAITWHATSERVER.find('_')]
        
        # Run indefinitively.
        while True:
            
            # Pause a bit, we don't want to overdo it.
            time.sleep(self._refresh_delay)
            
            # Get new incoming commands.
            self._udplock.acquire()
            cmds = self.udp.getCommands()
            self._udplock.release()
            
            # Add new commands to the queue.
            for c in cmds:
                # Parse the message.
                target, message, clienttime = c.text.split('|')
                self._print("Found message (%s to %s, t=%s) '%s'" % \
                    (c.ip, target, clienttime, message))
                # Only process the messages that were directed at the server.
                if int(target) == self._servernr:
                    # Check if the sender has an incoming queue already.
                    if c.ip in self._clientlist:
                        # Check if the message is a 'waitwhat' message.
                        if expected['waitwhat'] in message:
                            # Parse the waitwhat message, which looks like this:
                            # 'waitwhatserver_expected=%s'
                            msg, xpctd = message.split('_')
                            xpctd = xpctd[xpctd.find('=')+1:]
                            # Re-send the last version of the expected message.
                            if xpctd in self._lastmessage[c.ip].keys():
                                self._print("Resending the last version of expected message '%s' to client %s: '%s'" % \
                                    (expected, c.ip, self._lastmessage[c.ip][xpctd]))
                                self._msg_client(self._lastmessage[c.ip][xpctd], c.ip)
                            else:
                                self._print("Do not have a last version of expected message '%s' for client %s" % \
                                    (xpctd, c.ip))
                        # Process the regular message.
                        else:
                            # Use the timestamps to sync timing.
                            clienttime = int(clienttime) / 1000.0
                            servertime = c.timestamp / 1000.0
                            td = self._time_sync(c.ip, clienttime, servertime)
                            # Store the message in the correct queue.
                            self._print("Adding message '%s' (t=%d) to incoming queue for client %d" \
                                % (message, clienttime+td, c.ip))
                            self._inlock.acquire()
                            self._incoming[c.ip].append((clienttime+td, message))
                            self._inlock.release()
                    # Otherwise, add the message to the 'unknown' queue.
                    else:
                        self._print("Adding message '%s' to incoming queue for unregistered clients" \
                            % (message))
                        self._inlock.acquire()
                        self._incoming['unknown'].append(c)
                        self._inlock.release()
                else:
                    self._print("Ignoring message '%s', as it wasn't for server %d" \
                        % (message, self._servernr))
    
    
    def _outgoing_messages(self):
        
        """For internal use. Empties out the list of outgoing messages. Call
        from a Thread.
        """
        
        # Run indefinitively.
        while True:
            
            # Pause a bit, we don't want to overdo it.
            time.sleep(self._refresh_delay)
            
            # Loop through all clients.
            for clientnr in self._clientlist:
                
                # Loop through all outgoing messages.
                while len(self._outgoing[clientnr]) > 0:
                    
                    # Get the next message.
                    self._outlock.acquire()
                    message = self._outgoing[clientnr].pop(0)
                    self._outlock.release()
                    
                    # Send dat phat message!
                    self._print("Sending '%s' to client %d." % (message, clientnr))
                    self._udplock.acquire()
                    msg = 'cmd,%d|%s' % (clientnr, message)
                    self.udp.sendWithTimeStamp(msg, '|')
                    for i in range(self._message_reps):
                        self.udp.sendWithTimeStamp(msg, '|')
                    self._udplock.release()
                    
                    # Update the last-message-sent dict.
                    if clientnr not in self._lastmessage.keys():
                        self._lastmessage[clientnr] = {}
                    if '_' in message:
                        m = message[:message.find('_')]
                    else:
                        m = message
                    self._lastmessage[clientnr][m] = message


    def _msg_all_clients(self, message):
        
        """Convenience funtion for sending messages to all clients at once.
        """
        
        for clientnr in self._clientlist:
            self._msg_client(message, clientnr)


    def _msg_client(self, message, clientnr):
        
        """Convenience funtion for sending messages to a client.
        """
        
        # Send a message to a specific client.
        self._print("Adding '%s' to outgoing queue for client %d." % \
            (message, clientnr))
        self._outlock.acquire()
        self._outgoing[clientnr].append(message)
        self._outlock.release()
    
    
    def _time_sync(self, clientnr, clienttime, servertime):
        
        """Compates the times of the server and the client, computes the
        difference, adds this to a list, and computers the average time
        difference between server and this particular client
        """
        
        # Add the new timestamp to the list.
        self._timelock.acquire()
        self._timelist[clientnr].append(servertime - clienttime)
        # Chuck out the oldest timestamp.
        while len(self._timelist[clientnr]) > self._timelistlen:
            self._timelist[clientnr].pop(0)
        # Compute the average time difference.
        td = numpy.mean(self._timelist[clientnr])
        self._timelock.release()
        
        return td


    def _wait_for_message(self, clientnr, expectedmessage, timeout, \
        mintime=None):
        
        """Convenience funtion for waiting for specific messages from a client.
        """
        
        # Parse the first part (the message) from the expected reply. We need
        # to do this, because parts of some replies will contain parameters
        # that differ between replies.
        if '_' in expectedmessage:
            expected = copy.copy(expectedmessage[:expectedmessage.find('_')])
        else:
            expected = copy.copy(expectedmessage)
        
        self._print("Waiting for message '%s' from %s" % (expected, clientnr))

        # Wait for a message or a timeout.
        t0 = time.time()
        last_attempt = time.time()
        no_message = True
        no_timeout = True
        while no_message and no_timeout:
            # Get the current message queue.
            self._inlock.acquire()
            cmds = copy.deepcopy(self._incoming[clientnr])
            self._inlock.release()
            # Loop through the queue.
            for t, c in cmds:
                self._print("Examining message '%s' from %s (t=%d)" % (c, clientnr, t))
                # Check if the message fits the expected time (with a grace
                # period of several milliseconds).
                if mintime != None:
                    if t < mintime:
                        # If the timestamp was before the minimum time,
                        # skip this particular message.
                        self._print("Message ('%s') was too early (t=%d < mintime=%d)" \
                            % (c, t, mintime))
                        continue

                # Check if the message is the expected message.
                if expected in c:
                    # Remove the message from the queue.
                    self._inlock.acquire()
                    i = self._incoming[clientnr].index((t,c))
                    self._incoming[clientnr].pop(i)
                    self._inlock.release()
                    # Stop the while loop.
                    no_message = False
                    break
                else:
                    self._print("Message ('%s') was not expected ('%s')" \
                        % (c, expectedmessage))

                # Check if there is a timeout.
                if timeout != None:
                    if time.time() - t0 > timeout:
                        no_timeout = False
                        break
            
            # Check if we should re-send the message.
            if time.time() - last_attempt > self._reptimeout:
                # Let the client know what we're expecting from them.
                self._wait_what(clientnr, expected)
                # Update the last attempt time.
                last_attempt = time.time()
        
        # Clean up the incoming queue.
        self._clean_incoming(clientnr, t0*1000 - self._cleantime)

        # Return a success Boolean and the message/fault.
        if no_message == False:
            return (True, c)
        if no_timeout == False:
            return (False, 'timeout')
        return (False, 'unknown')


    def _wait_for_reply(self, message, clientnr, expectedreply, timeout):
        
        """Convenience funtion for sending messages to the server, and
        to wait for a specific response.
        """
        
        # Parse the first part (the message) from the expected reply. We need
        # to do this, because parts of some replies will contain parameters
        # that differ between replies.
        if '_' in expectedreply:
            expected = copy.copy(expectedreply[:expectedreply.find('_')])
        else:
            expected = copy.copy(expectedreply)
        
        # Send the message to the client, and get a rough timestamp for the
        # moment of sending.
        sendtime = time.time() * 1000
        self._msg_client(message, clientnr)
        
        # Wait for the expected reply.
        success, reply = self._wait_for_message(clientnr, expected, timeout, \
            mintime=sendtime)

        return (success, reply)
    
    
    def _wait_what(self, clientnr, expected):
        
        """Sends a message to the client, effectively asking to repeat the
        last message of a specific type.
        """
        
        self._msg_client(cb.WAITWHATCLIENT % (expected), clientnr)


    def close(self):
        
        self._clientlist = []
        self.udp.deInit()

    
    def get_fixation(self, clientnr, fixlist, timeout=None):
        
        """Checks whether a particular client has reported any new fixations,
        checks them against the existing list of fixations, and sends a 1 back
        if the fixated stimulus is new, or a 0 if a different client has
        already fixated the newly reported stimulus.
        """
        
        # Wait until the client reports a fixation.
        success, message = self._wait_for_message(clientnr, cb.SETFIXATION, \
            timeout)
        
        # Parse the message, which looks like this:
        # 'setfix_stimnr=%d'
        msg, stimnr = message.split('_')
        # Get the number of the fixated stimulus.
        stimnr = int(stimnr[stimnr.find('=')+1:])
        
        # Send back an OK if the fixated stimulus wasn't fixated before.
        # (A fixnr of -1 means there was no new fixation.)
        ok = int( (stimnr not in fixlist) and (stimnr >= 0) )
        self._msg_client(cb.FIXATIONRECEIVED % (ok), clientnr)
        
        if ok:
            return stimnr
        else:
            return None
    
    
    def get_ready(self, clientnr, timeout=None):
        
        """Waits until a particular client says it's ready for the next trial.
        """
        
        # Wait for a message from the client.
        self._wait_for_message(clientnr, cb.CLIENTREADY, timeout)
        
        # Reply to let the client know that the server is ready for the next
        # trial too.
        self._msg_client(cb.CLIENTGOGOGO, clientnr)
    
    
    def get_response(self, clientnr, timeout=None):
        
        """Waits for a client to report the participant's response.
        """
        
        # Wait for the message from the client.
        success, message = self._wait_for_message(clientnr, cb.RESPONSE, \
            timeout)
        
        # Tell the client we got the response.
        self._msg_client(cb.GOTRESPONSE, clientnr)
        
        # Parse the message, which looks like this:
        # 'response_ori=%d'
        msg, respori = message.split('_')
        # Get the number of the fixated stimulus.
        respori = int(respori[respori.find('=')+1:])
        
        return respori
    
    
    def set_continue(self, clientnr, running, timeout=None):
        
        """Tells all clients that the experiment is still running (or not).
        """
        
        # Wait for the message from the client.
        success, message = self._wait_for_message(clientnr, cb.DOESEXPCONTINUE, \
            timeout)
        
        # Tell the client whether the experiment continues.
        self._msg_client(cb.EXPCONTINUES % (int(running)), clientnr)
    
    
    def set_fixations(self, newfixlist, trialrunning, timeout=None):
        
        """Broadcasts all newly fixated stimulus numbers.
        """
        
        # Send the message to all clients. The message looks like this:
        # 'newfixations_fixlist=%s_trialrunning=%d'
        self._msg_all_clients(cb.NEWFIXATIONS % (newfixlist, int(trialrunning)))
    
    
    def set_outcome(self, clientnr, error, reward, others):
        
        """Tells a specific client what their trial outcome is.
        """
        
        # Send a message to the client.
        self._msg_client(cb.TRIALOUTCOME % (error, reward, others), clientnr)
    
    
    def set_parameters(self, clientnr, clientcolour):
        
        """Sends a client's parameters.
        """
        
        # Message the client.
        self._msg_client(cb.PARAMETERS % (clientnr, clientcolour), clientnr)
    
    
    def set_probe(self, clientnr, probenr, timeout=None):
        
        """Tells a client what stimulus should be probed.
        """
        
        # Wait until the client asks for the probe.
        self._wait_for_message(clientnr, cb.PROBEPLZ, timeout)

        # Message the client.
        self._msg_client(cb.PROBENR % (probenr), clientnr)
        
    
    def set_trial(self, trialdict):
        
        """Sends the trial parameters to all clients.
        """
        
        # Convert the trial parameters into a JSON-fomatted string.
        paramstring = json.dumps(trialdict)
        
        # Send the message.
        self._msg_all_clients(cb.TRIALPARAMETERS % (paramstring))
    
    
    def wait_for_clients(self, timeout=None):
        
        """Waits for clients to contact, then establishes a connection with
        each client, and then waits for more clients until nclients are
        connected.
        """
        
        self._print("Waiting for clients %s to join" % (self._unconfirmed_ip_addresses))
        
        # Look out for running clients.
        success = False
        timed_out = False
        t0 = time.time()
        while not success and not timed_out:
            
            # Get all broadcasted commands.
            self._inlock.acquire()
            messages = copy.deepcopy(self._incoming['unknown'])
            self._inlock.release()
            
            # Loop through received commands.
            for m in messages:
                # Check whether this is a starting message.
                if cb.START in m.text:
                    # Check whether the IP address is expected, but not in
                    # the list of already accepted addresses.
                    if (m.ip in self._unconfirmed_ip_addresses) and (m.ip not in self._clientlist):
                        self._print("Adding %d to the list of clients" % (m.ip))
                        # Add the IP address to the list of clients.
                        self._outlock.acquire()
                        self._outgoing[m.ip] = []
                        self._outlock.release()
                        self._inlock.acquire()
                        self._incoming[m.ip] = []
                        self._inlock.release()
                        self._timelock.acquire()
                        self._timelist[m.ip] = []
                        self._timelock.release()
                        self._clientlist.append(m.ip)
                        # Remove the IP address from the list of expected
                        # clients, and add it to the accepted clients.
                        self._unconfirmed_ip_addresses.pop(self._unconfirmed_ip_addresses.index(m.ip))
                        # Send a message back to tell the client it's
                        # accepted.
                        self._msg_client(cb.CLIENTCONFIRM % (m.ip), m.ip)
           
                # Check if we have all clients yet.
                if len(self._unconfirmed_ip_addresses) == 0:
                    success = True
                
                # Check if a timeout occurred.
                if timeout != None:
                    if time.time() - t0 > timeout:
                        return False
        
        # Set the new clientlist as a filter to the UDP Socket.
        self.udp.setComputerFilter(self._clientlist)
        
        return True

