# -*- coding: utf-8 -*-

# The maximum length of a message
MAXLEN = 1024

# START indicates a client would like to join the experiment.
START = 'startexp'
# CLIENTCONFIRM is sent by the server to confirm that the client has been
# registered.
CLIENTCONFIRM = 'clientoktojoin_nr=%d'

## PARAMETERSPLZ is sent by a client to ask for the experimental parameters.
#PARAMETERSPLZ = 'parametersplz'
# PARAMETERS is sent by the server to tell the client what this experiments'
# parameters are going to be.
PARAMETERS = 'expparams_nr=%d_col=%s'

## TRIALPLZ is sent by the client to ask for the parameters of the next trial.
#TRIALPLZ = 'trialplz'
# TRIALPARAMETERS is sent by the server to provide the client with all trial
# information. The wildcard is replaced by a JSON-formatted dict.
TRIALPARAMETERS = 'trialparams_%s'

## GETFIXATIONS is sent by the client to ask for updates on which stimuli have
## been fixated since the last time GETFIXATIONS was sent.
#GETFIXATIONS = 'newfixplz'
# NEWFIXATIONS is sent by the server to update the client on which stimuli
# have been fixated since the client last asked, and on whether the trial
# should continue to run. fixlist is a list (as a string), and trialrunning
# is either 1 or 0.
NEWFIXATIONS = 'newfixations_fixlist=%s_trialrunning=%d'

# SETFIXATION is sent by the client to let the server know that a new
# stimulus has been fixated. The %d will be converted to the number of the
# fixated stimulus.
SETFIXATION = 'setfix_stimnr=%d'
# FIXATIONRECEIVED is sent by the server to acknowledge that a new fixation
# was detected by the client, and also tells the client whether this new
# fixation has been excepted. %d will be converted to 1 if it was accepted,
# or 0 if the new fixation was not accepted.
FIXATIONRECEIVED = 'newfixreceived_ok=%d'

## PROBEPLZ is sent by the client to ask which stimulus should be probed.
PROBEPLZ = 'probeplz'
# PROBENR is sent by the server to provide the client with the number of the
# stimulus that should be probed.
PROBENR = 'probenr_nr=%d'

# RESPONSE is sent by the client to report what a participant's response was.
RESPONSE = 'response_ori=%d'
# GOTRESPONSE is sent by the server to confirm to the client that their
# participant's response was received.
GOTRESPONSE = 'gotrespok'

## OUTCOMEPLZ is sent by the client to ask for the outcome of a trial, namely
## the response error and the earned reward.
#OUTCOMEPLZ = 'outcomeplz'
# TRIALOUTCOME is sent by the server to report the outcome of a trial, namely
# the response error and the earned reward.
TRIALOUTCOME = 'trialoutcome_error=%d_reward=%.2f_others=%s'

## DOESEXPCONTINUE is sent by the client to ask whether the experiment is still
## running.
DOESEXPCONTINUE = 'doesexpcontinue'
# EXPCONTINUES is sent by the server to let the client know whether the
# experiment is still on.
EXPCONTINUES = 'expcontinues_running=%d'

# CLIENTREADY is sent by the client to confirm that the client is ready for
# the next trial to start.
CLIENTREADY = 'clientisready'
# CLIENTGOGOGO is sent by the server to confirm that the server knows that the
# client is ready for the next trial.
CLIENTGOGOGO = 'clientgogogo'

# WAITWHATCLIENT is sent by the server to ask the client what the last version
# of a particular expected message was.
WAITWHATCLIENT = 'waitwhatclient_expected=%s'
# WAITWHATSERVER is sent by the client to ask the server what the last version
# of a particular expected message was.
WAITWHATSERVER = 'waitwhatserver_expected=%s'

# CLOSE indicates the experiment is terminated.
CLOSE = 'stopexp'

