from constants import *
import codebook as cb

import UDPClient

import os
import ast
import copy
import json
import time
import socket
from threading import Lock, Thread


class Client:
    
    def __init__(self, multicast_ip):
        
        if DEBUG:
            self._debugfile = open(os.path.join(DEBUGDIR, 'client_debug_%s.txt' \
                % (STARTTIME)), 'w')
        
        # Start with an unknown client number. The server tells each client
        # who they are.
        self._clientnr = None
        
        # Set how many messages should be sent each time (Copies can be
        # included to make it more likely that all messages arrive. This is
        # a very lazy implementation. Don't judge me.)
        self._message_reps = 0
        
        # Set up a dict for keeping track of the last sent messages of each
        # type.
        self._lastmessage = {}
        
        # Time spent waiting for an expected message until a request is sent
        # to repeat it. IN SECONDS!
        self._reptimeout = 10.0
        
        # Set up a list with incoming and outgoing messages.
        self._incoming = []
        self._outgoing = []
        self._maxincominglen = 10
        self._incominglock = Lock()
        self._outgoinglock = Lock()
        
        # Open a new UDP client
        self._print("%s: Opening a new UDPClient (multicast IP %s)" % \
            (self._clientnr, multicast_ip))
        self.udp = UDPClient.UDPClient()
        self.udp.port = 4444
        self.udp.groupAddress = '225.0.0.%s' % (multicast_ip)
        self.udp.loopBack = False
        self.udp.reuseSocket = 1
        self.udp.numReceiverThreads = 2
        self.udp.init()
        
        # Only listen to the Server.
        self._servernr = int(SERVERIP.split('.')[-1])
        self.udp.setComputerFilter([self._servernr])
        
        # Start message processing Thread.
        self._print("%s: Starting the message processing Thread" % \
            (self._clientnr))
        self._message_processer = Thread(target=self._process_messages)
        self._message_processer.name = "client_message_processor"
        self._message_processer.daemon = True
        self._message_processer.start()
    
    
    def _process_messages(self):
        
        """For internal use. Listens in on the communications that come in,
        and sends off messages that should go out. Also makes sure that the
        queues for incoming and outgoing stimuli don't get too long.
        """
        
        self._print("%s: Starting _process messages, looking out for special messages:" \
            % (self._clientnr))
        
        # Set some expected messages.
        expected = {}
        expected['clientconfirm'] = cb.CLIENTCONFIRM[:cb.CLIENTCONFIRM.find('_')]
        expected['waitwhat'] = cb.WAITWHATCLIENT[:cb.WAITWHATCLIENT.find('_')]
        
        for key in expected.keys():
            self._print("%s: Special message '%s': '%s'" % \
                (self._clientnr, key, expected[key]))
        
        # Run idefinitively
        while True:
            
            # Get new incoming commands.
            cmds = self.udp.getCommands()
            self._print("%s: Found %d new UDP commands." % \
                (self._clientnr, len(cmds)))
            # Add new commands to the queue.
            for c in cmds:
                # Parse the message.
                target, message, clienttime = c.text.split('|')
                self._print("%s: Found message (%s to %s, t=%s) '%s'" % \
                    (self._clientnr, c.ip, target, clienttime, message))
                # Only process messages from the server.
                if c.ip == self._servernr:
                    # Check if this is a client confirmation message.
                    if expected['clientconfirm'] in message:
                        self._print("%s: Adding message '%s' (t=%s) to the incoming queue" \
                            % (self._clientnr, message, clienttime))
                        self._incominglock.acquire()
                        self._incoming.append(message)
                        self._incominglock.release()
                    # Only process the messages that were directed at this client.
                    elif target in ['None', str(self._clientnr)]:
                        # Check if this is a confused message to find out what
                        # the client is waiting for.
                        if expected['waitwhat'] in message:
                            self._print("%s: Received '%s' from server" % \
                                (self._clientnr, message))
                            # Parse the waitwhat message, which looks like this:
                            # 'waitwhatclient_expected=%s'
                            msg, xpctd = message.split('_')
                            xpctd = xpctd[xpctd.find('=')+1:]
                            # Re-send the last version of the expected message.
                            if xpctd in self._lastmessage.keys():
                                self._outgoing.append(self._lastmessage[xpctd])
                                self._print("%s: Resending the last version of expected message '%s': '%s'" % \
                                    (self._clientnr, xpctd, self._lastmessage[xpctd]))
                            else:
                                self._print("%s: Do not have a last version of expected message '%s'" % \
                                    (self._clientnr, xpctd))
                        else:
                            # Add the message to the queue.
                            self._print("%s: Adding message '%s' (t=%s) to the incoming queue" \
                                % (self._clientnr, message, clienttime))
                            self._incominglock.acquire()
                            self._incoming.append(message)
                            self._incominglock.release()
                            # Chuck a message out if the queue is getting too long.
                            if len(self._incoming) > self._maxincominglen:
                                self._incominglock.acquire()
                                delmsg = self._incoming.pop(0)
                                self._incominglock.release()
                                self._print("%s: Removed message '%s' from the incoming queue" \
                                    % (self._clientnr, delmsg))
                    else:
                        self._print("%s: Ignoring message '%s', as it wasn't for me (%s)" \
                            % (self._clientnr, message, self._clientnr))
                else:
                    self._print("%s: Ignoring message '%s', as it wasn't from the server (%s)" \
                        % (self._clientnr, message, self._servernr))
            
            # Process outgoing commands.
            while len(self._outgoing) > 0:
                # Send a message to the server.
                self._outgoinglock.acquire()
                message = self._outgoing.pop(0)
                self._outgoinglock.release()
                self._print("%s: Sending '%s' to %s" % \
                    (self._clientnr, message, self._servernr))
                msg = 'cmd,%s|%s' % (self._servernr, message)
                self.udp.sendWithTimeStamp(msg, '|')
                for i in range(self._message_reps):
                    self.udp.sendWithTimeStamp(msg, '|')
                # Store the message in the 'last sent' dict.
                if '_' in message:
                    m = message[:message.find('_')]
                else:
                    m = message
                self._lastmessage[m] = message
    
    
    def _print(self, message):
        
        if DEBUG:
            self._debugfile.write("(%d) %s\n" % (time.time()*1000, message))
            self._debugfile.flush() # internal buffer to RAM
            os.fsync(self._debugfile.fileno()) # RAM file cache to disk


    def _msg_server(self, message):
        
        """Convenience funtion for sending messages to the server.
        """
        
        # Add message to the outgoing queue.
        self._outgoinglock.acquire()
        self._outgoing.append(message)
        self._outgoinglock.release()
        
#        # Send a message to the server.
#        self._print("%s: Sending '%s' to %s" % (self._clientnr, message, self._servernr))
#        msg = 'cmd,%d|%s' % (self._servernr, message)
#        self.udp.sendWithTimeStamp(msg, '|')
#        for i in range(self._message_reps):
#            self.udp.sendWithTimeStamp(msg, '|')


    def _wait_for_message(self, expectedmessage, timeout):
        
        """Convenience funtion for waiting for specific messages from the
        server.
        """
        
        # Parse the first part (the message) from the expected reply. We need
        # to do this, because parts of some replies will contain parameters
        # that differ between replies.
        if '_' in expectedmessage:
            expected = copy.copy(expectedmessage[:expectedmessage.find('_')])
        else:
            expected = copy.copy(expectedmessage)
        
        self._print("%s: Waiting for message '%s'" % (self._clientnr, expected))
        
        # Wait for a message or a timeout.
        t0 = time.time()
        last_attempt = time.time()
        no_message = True
        no_timeout = True
        while no_message and no_timeout:
            # Get the current message queue.
            cmds = copy.deepcopy(self._incoming)
            # Loop through the queue.
            for msg in cmds:
                self._print("%s: Examining message '%s'" % \
                    (self._clientnr, msg))

                # Check if the message is the expected message.
                if expected in msg:
                    # Remove the message from the queue.
                    self._incominglock.acquire()
                    self._incoming.pop(self._incoming.index(msg))
                    self._incominglock.release()
                    no_message = False
                    break
                else:
                    self._print("Message ('%s') was not expected ('%s')" \
                        % (msg, expectedmessage))

                # Check if there is a timeout.
                if timeout != None:
                    if time.time() - t0 > timeout:
                        no_timeout = False
                        break
            
            # Check if we should re-send the message.
            if time.time() - last_attempt > self._reptimeout:
                # Let the server know what we're expecting from them.
                self._wait_what(expected)
                # Update the last attempt time.
                last_attempt = time.time()

        # Return a success Boolean and the message/fault.
        if no_message == False:
            return (True, msg)
        if no_timeout == False:
            return (False, 'timeout')
        return (False, 'unknown')


    def _wait_for_reply(self, message, expectedreply, timeout):
        
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
        
        # Send the message to the server.
        self._msg_server(message)

        # Wait for the expected reply.
        success, reply = self._wait_for_message(expected, timeout)

        # Return a success Boolean and the reply/fault.
        return (success, reply)
    
    
    def _wait_what(self, expected):
        
        """Sends a message to the server, effectively asking to repeat the
        last message of a specific type.
        """
        
        self._msg_server(cb.WAITWHATSERVER % (expected))


    def close(self):
        
        self._print("%s: Closing connection to %s" % (self._clientnr, self._servernr))
        self._servernr = None
        self._clientnr = None
        self.udp.deInit()


    def contact_server(self, timeout=None):
        
        """Let's the server know that this client is ready to go.
        """
        
        # TODO: Generate UUID, send UUID to Server; wait for server to reply
        # with same UUID.
        # uid = uuid.uuid1()
        
        # Send a message to the Server that indicates we're OK to start,
        # and wait for the Server to acknowledge this.
        success, reply = self._wait_for_reply(cb.START, cb.CLIENTCONFIRM, \
            timeout)
        
        # Parse the reply, which looks like this:
        # 'clientoktojoin_nr=%d'
        msg, clientnr = reply.split('_')
        # Set the number for this client.
        self._clientnr = int(clientnr[clientnr.find('=')+1:])

        # Return the successfullness.
        return success


    def get_continue(self, timeout=None):
        
        """Check with the server whether the experiment should be continued.
        """
        
        # Send a request to the server.
        success, reply = self._wait_for_reply(cb.DOESEXPCONTINUE, \
            cb.EXPCONTINUES, timeout)
        
        # Parse the reply, which looks like this:
        # 'expcontinue_running=%d'
        msg, running = reply.split('_')
        # Convert the 'running=0' or 'running=1' string into a bool.
        running = bool(int(running[running.find('=')+1:]))
        
        return running
    
    
    def get_fixations(self, timeout=None):
        
        """Get a list of the numbers of stimuli that were fixated since this
        function was last called.
        """
        
        # Send a request to the server.
        success, reply = self._wait_for_message(cb.NEWFIXATIONS, \
            timeout=timeout)
        
        # Parse the reply, which looks like this:
        # newfixations_fixlist=%s_trialrunning=%d
        msg, fixlist, trialrunning = reply.split('_')
        # Convert the fixlist from str to list
        fixlist = ast.literal_eval(fixlist[fixlist.find('=')+1:])
        # Convert the 'trialrunning=0' or 'trialrunning=1' string into a bool.
        trialrunning = bool(int(trialrunning[trialrunning.find('=')+1:]))
        
        return fixlist, trialrunning


    def get_outcome(self, timeout=None):
        
        """Get the outcome of the trial.
        """
        
        # Wait for the server to send a message.
        success, reply = self._wait_for_message(cb.TRIALOUTCOME, \
            timeout=timeout)
        
        # Parse the reply, which looks like this:
        # 'trialoutcome_error=%d_reward=%.2f_others=%s'
        msg, error, reward, others = reply.split('_')
        # Get the response error.
        error = int(error[error.find('=')+1:])
        # Get the reward for this trial.
        reward = float(reward[reward.find('=')+1:])
        # Get the reward for all others on this trial.
        others = ast.literal_eval(others[others.find('=')+1:])
        
        return error, reward, others


    def get_parameters(self, timeout=None):
        
        """Get the experiment parameters from the server.
        """
        
        # Send a request to the server.
        success, reply = self._wait_for_message(cb.PARAMETERS, timeout)
        
        # Parse the reply, which looks like this:
        # 'expparameters_nr=%d_col=%s'
        msg, clientnr, colour = reply.split('_')
        # Get the number of stimuli.
        clientnr = int(clientnr[clientnr.find('=')+1:])
        # Get the colour for this participant.
        colour = ast.literal_eval(colour[colour.find('=')+1:])
        
        return clientnr, colour


    def get_probe(self, timeout=None):
        
        """Gets the number of the probed stimulus.
        """

        # Send a request to the server.
        success, reply = self._wait_for_reply(cb.PROBEPLZ,  \
            cb.PROBENR, timeout=timeout)
        
        # Parse the reply, which looks like this:
        # 'probenr_nr=%d'
        msg, probenr = reply.split('_')
        # Get the number of the probed stimulus.
        probenr = int(probenr[probenr.find('=')+1:])
        
        return probenr


    def get_trial(self, timeout=None):
        
        """Get the parameters for the next trial from the server.
        """
        
        # Send a request to the server.
        success, reply = self._wait_for_message(cb.TRIALPARAMETERS, \
            timeout=timeout)
        
        # Parse the reply (formatted as json).
        msg, jsondict = reply.split('_')
        # Parse the json into a dict.
        trial = json.loads(jsondict)
        
        return trial


    def set_fixation(self, stimnr, timeout=None):
        
        """Tells the server that a new stimulus has been fixated.
        """

        # Send a message to the server, and wait for confirmation of receipt.
        success, reply = self._wait_for_reply(cb.SETFIXATION % (stimnr), \
            cb.FIXATIONRECEIVED, timeout=timeout)
        
        # Parse the reply, which looks like this:
        # 'newfixreceived_ok=%d'
        msg, ok = reply.split('_')
        # Convert the 'ok=0' or 'ok=1' string into a bool.
        ok = bool(int(ok[ok.find('=')+1:]))
        
        return ok


    def set_ready(self, timeout=None):
        
        """Tells the server that the client is ready for the next trial.
        """

        # Send a message to the server, and wait for confirmation of receipt.
        success, reply = self._wait_for_reply(cb.CLIENTREADY, \
            cb.CLIENTGOGOGO, timeout=timeout)
        
        return success


    def set_response(self, respori, timeout=None):
        
        """Sends the response (the reported stimulus orientation) to the
        server.
        """
        
        # Send a message to the server, and wait for confirmation of receipt.
        success, reply = self._wait_for_reply(cb.RESPONSE % (int(respori)), \
            cb.GOTRESPONSE, timeout=timeout)
        
        return success

