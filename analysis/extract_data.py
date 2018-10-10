# -*- coding: utf-8 -*-

import os
import copy
import pickle

import numpy

import vonmises_toolbox as vt
from entropy import KL_divergence_continuous


# # # # #
# CONSTANTS

# GENERAL
# Sort participants within a pair by short-term memory capacity or by apathy.
# (Or by any other key in the idata dict.)
SORTBY = 'capacity'
# Included runs for the analysis.
RUNS = [ \
    #'2017-04-13_pm', \ # pilot experiment for tech setup, lacks questions
    '2017-04-19_13-30', \
    '2017-04-19_15-00', \
    #'2017-04-20_13-30', \ # excluded for only having one pair, thereby ruining anonymisation.
    '2017-04-21_09-30', \
    '2017-04-21_13-30', \
    '2017-04-21_15-00']

# SHORT-TERM MEMORY TEST
# Number of stimuli used.
WM_NSTIM = [2, 4]
# Number of stimuli used in the collaborative game.
NSTIM = 8
# Number of trials in the collaborative game.
NTRIALS = 25
# Fit the non-targets?
FITNONTAR = True

# Files and folders.
# Folder that contains this script file.
CODEDIR = os.path.dirname(os.path.abspath(__file__))
# Highest-level directory that we need.
DIR = os.path.dirname(CODEDIR)
# Directory that contains all data we collected.
DATADIR = os.path.join(DIR, 'runs')
if not os.path.isdir(DATADIR):
    raise Exception("Could not find data folder at '%s'" % (DATADIR))
# Directory for pre-processed data (intermediate steps with numbers that are
# generated through different analysis steps.)
PPDIR = os.path.join(DIR, 'processed_data_sorted-by-%s' % (SORTBY))
if not os.path.isdir(PPDIR):
    os.mkdir(PPDIR)
# Directory for generated output, such as graphs.
OUTDIR = os.path.join(DIR, 'output_sorted-by-%s' % (SORTBY))
if not os.path.isdir(OUTDIR):
    os.mkdir(OUTDIR)


# # # # #
# EXTRACT DATA

print("\nExtracting data from text files.")

participants = []
wmdata = {}
qdata = {}
pairs = []
sesdata = {}


# APATHY DATA (AMI QUESTIONNAIRE)
print("\tProcessing apathy questionnaire.")
# Count the number of files (one for each participant) we have.
ami_file = os.path.join(DATADIR, 'questionnaires', 'AMI.csv')
qdata['ami'] = {}
# Read the data from the current file. This results in a series of
# NumPy arrays: One for each column in the data file. This vector is
# all string values, and the first index contains the name of the
# column.
raw = numpy.loadtxt(ami_file, unpack=True, dtype=str, delimiter=',')
# To be able to process the data better, we will stick it in a dict
# where each of the columns in the data file will be associated with
# a key (the variable name from the header in the data file), and the
# values in a numerical data type where appropriate. A conversion to
# floating point numbers will be attempted for all variables, and will
# be given up on when it fails.
qd = {}
for j in range(len(raw)):
    var = raw[j][0]
    try:
        val = raw[j][1:].astype(float)
    except:
        val = raw[j][1:]
    qd[var] = copy.deepcopy(val)
# Now we can go through all trials and sort the data by runs and IPs.
n = len(qd['u'])
for i in range(n):
    # Get the session number and the IP number for the computer.
    run = qd['session'][i]
    ip = qd['ip'][i]
    # Only questionnaire data from people that have partaken were assigned
    # a session string. Other will be an empty string.
    if run != '':
        # Ignore if the current participant's run is not in the RUNS we're
        # currently analysing.
        if run not in RUNS:
            continue
        # Add a new entry in the dict for the current participant.
        ppname = '%d_%d' % (RUNS.index(run), int(ip))
        qdata['ami'][ppname] = {}
        # Start with 0 as total apathy score, and the sub-scores.
        t = 0.0
        es = 0.0
        sm = 0.0
        ba = 0.0
        # Loop through all the questions, Q1-Q18.
        for j in range(1,19):
            # Construct the question name.
            qname = 'Q%d' % (j)
            # Update the total and the relevant sub-score with the question's
            # response value.
            t += float(qd[qname][i])
            # Behavioural activation sub-scale.
            if j in [5, 9, 10, 11, 12, 15]:
                ba += float(qd[qname][i])
            # Social motivation sub-scale.
            elif j in [2, 3, 4, 8, 14, 17]:
                sm += float(qd[qname][i])
            # Emotional sensitivity sub-scale.
            elif j in [1, 6, 7, 13, 16, 18]:
                es += float(qd[qname][i])
            # Unaccounted for questions (there shouldn't be any).
            else:
                raise Exception("ERROR: The stupid programmer writes ugly code, and also missed a question on the AMI.")
        # Store the data.
        qdata['ami'][ppname]['total'] = copy.deepcopy(t)
        qdata['ami'][ppname]['ba'] = copy.deepcopy(ba)
        qdata['ami'][ppname]['sm'] = copy.deepcopy(sm)
        qdata['ami'][ppname]['es'] = copy.deepcopy(es)

# BIG FIVE DATA
print("\tProcessing big five questionnaire.")
# Count the number of files (one for each participant) we have.
bigfive_file = os.path.join(DATADIR, 'questionnaires', 'BigFive.csv')
qdata['bigfive'] = {}
# Read the data from the current file. This results in a series of
# NumPy arrays: One for each column in the data file. This vector is
# all string values, and the first index contains the name of the
# column.
raw = numpy.loadtxt(bigfive_file, unpack=True, dtype=str, delimiter=',')
# To be able to process the data better, we will stick it in a dict
# where each of the columns in the data file will be associated with
# a key (the variable name from the header in the data file), and the
# values in a numerical data type where appropriate. A conversion to
# floating point numbers will be attempted for all variables, and will
# be given up on when it fails.
qd = {}
for j in range(len(raw)):
    var = raw[j][0]
    try:
        val = raw[j][1:].astype(float)
    except:
        val = raw[j][1:]
    qd[var] = copy.deepcopy(val)
# Now we can go through all trials and sort the data by runs and IPs.
n = len(qd['u'])
for i in range(n):
    # Get the session number and the IP number for the computer.
    run = qd['session'][i]
    ip = qd['ip'][i]
    # Only questionnaire data from people that have partaken were assigned
    # a session string. Other will be an empty string.
    if run != '':
        # Ignore if the current participant's run is not in the RUNS we're
        # currently analysing.
        if run not in RUNS:
            continue
        # Add a new entry in the dict for the current participant.
        ppname = '%d_%d' % (RUNS.index(run), int(ip))
        qdata['bigfive'][ppname] = {}
        # Start with 0s as the sub-scores.
        extravert = 0.0
        agreeable = 0.0
        conscient = 0.0
        emostable = 0.0
        opentoexp = 0.0
        # Loop through all the questions, Q119_1-Q119_10.
        for j in range(1,11):
            # Construct the question name.
            qname = 'Q119_%d' % (j)
            # Convert and flip the score from 8-14 to 1-7
            #    Original: 8 (strongly agree) - 14 (strongly disagree)
            #    Intended: 1 (disagree strongly) - 7 (strongly agree)
            score = 15 - float(qd[qname][i])
            # Check if the question is reverese-scored.
            if j in [1, 3, 4, 7, 9]:
                score = 8 - score
            # Extraversion subscale.
            if j in [3,10]:
                extravert += score
            # Agreeableness subscale.
            elif j in [2,4]:
                agreeable += score
            # Conscientiousness subscale.
            elif j in [5,7]:
                conscient += score
            # Emotional stability subscale.
            elif j in [1,8]:
                emostable += score
            # Openness to experience subscale.
            elif j in [6,9]:
                opentoexp += score
            # Unaccounted for questions (there shouldn't be any).
            else:
                raise Exception("ERROR: The stupid programmer writes ugly code, and also missed a question on the BigFive.")
        # Store the data.
        qdata['bigfive'][ppname]['extraversion'] = copy.deepcopy(extravert)
        qdata['bigfive'][ppname]['agreeableness'] = copy.deepcopy(agreeable)
        qdata['bigfive'][ppname]['conscientiousness'] = copy.deepcopy(conscient)
        qdata['bigfive'][ppname]['emotional_stability'] = copy.deepcopy(emostable)
        qdata['bigfive'][ppname]['openness'] = copy.deepcopy(opentoexp)


# Loop through all runs, extract the relevant data, and store it in a pickle
# file for further processing.
for i, run in enumerate(RUNS):
    
    # WM CAPACITY
    print("\tProcessing short-term memory capacity test for run '%s' (%d/%d)." \
        % (run, i+1, len(RUNS)))
    # Count the number of files (one for each participant) we have.
    wm_files = os.listdir(os.path.join(DATADIR, run, 'wm_test'))
    wm_files.sort()
    n_participants = len(wm_files)
    # Go through all the files to extract data.
    for fname in wm_files:
        # Parse the file name to find the participant code. This will be a
        # combination of the session number and the computer's IP number.
        name, ext = os.path.splitext(fname)
        ppname = '%d_%d' % (RUNS.index(run), int(name.replace('ip','')))
        participants.append(ppname)
        # Read the data from the current file. This results in a series of
        # NumPy arrays: One for each column in the data file. This vector is
        # all string values, and the first index contains the name of the
        # column.
        raw = numpy.loadtxt(os.path.join(DATADIR, run, 'wm_test', fname), \
            unpack=True, dtype=str, delimiter='\t')
        # To be able to process the data better, we will stick it in a dict
        # where each of the columns in the data file will be associated with
        # a key (the variable name from the header in the data file), and the
        # values in a numerical data type where appropriate. A conversion to
        # floating point numbers will be attempted for all variables, and will
        # be given up on when it fails.
        wmdata[ppname] = {}
        for j in range(len(raw)):
            var = raw[j][0]
            try:
                val = raw[j][1:].astype(float)
            except:
                val = raw[j][1:]
            wmdata[ppname][var] = copy.deepcopy(val)
    
    
    # SESSION DATA
    print("\tProcessing pair data for run '%s' (%d/%d)." \
        % (run, i+1, len(RUNS)))
    # Count the number of files (one for each participant) we have.
    pair_files = os.listdir(os.path.join(DATADIR, run, 'collaborate'))
    pair_files.sort()
    # Go through the files and parse them to find all unique files in this run.
    fnames = []
    for fname in pair_files:
        name, ext = os.path.splitext(fname)
        if ext == '.txt' and name[-3:] == 'beh':
            # Only store files for which a matching collaboration-question
            # file is available, as these are the ones that actually finished
            # the experiment (some runs in 2017-04-19_13-30 crashed).
            if os.path.isfile(os.path.join(DATADIR, run, 'collaborate', \
                "%s.txt" % (name.replace('_beh', '_que')))):
                fnames.append(name.replace('_beh',''))
            else:
                print("EXCLUSION: '%s' did not complete experiment." % \
                    (name.replace('_beh','')))
    # Count the number of files we need to process.
    n_pair_files = len(fnames)
    # Go through all pair files, and store their data.
    for fname in fnames:
        # Parse the file name to find the current participant number.
        ip = int(fname[:fname.find('_')])
        # Load the question data.
        qraw = numpy.loadtxt( \
            os.path.join(DATADIR, run, 'collaborate', '%s_que.txt' % fname), \
            unpack=True, dtype=str, delimiter='\t')
        # Put the question data in a reasonable format.
        q = {}
        for j in range(len(qraw)):
            # The stored questions are really long, so we're going to shorten
            # them a bit to 'good' for "How good was the collaboration between
            # you and player %d?", 'likeable' for "How likeable did you find
            # player %d?", and 'fair' for "How fair did you think player %s
            # was?". This will greatly help in later data processing.
            if qraw[j][0] == 'question':
                for k in range(1, len(qraw[j])):
                    for var in ['good', 'likable', 'fair']:
                        if var in qraw[j][k]:
                            if var == 'likable':
                                var = 'likeable'
                            qraw[j][k] = var
                            break
            q[qraw[j][0]] = qraw[j][1:]
        # Find the IP of the other participant in this pair.
        ip_other = int(q['player'][0])
        # Generate a new key for the current pair in the session data dict.
        pname = '%d_%d-%d' % (i, ip, ip_other)
        sesdata[pname] = {}
        # Save this participant's responses to the questions.
        for j, var in enumerate(q['question']):
            sesdata[pname][var] = q['proportion'][j]
        # Load the data from this pair's game.
        raw = numpy.loadtxt( \
            os.path.join(DATADIR, run, 'collaborate', '%s_beh.txt' % fname), \
            unpack=True, dtype=str, delimiter='\t')
        # Parse the loaded raw data.
        for j in range(len(raw)):
            # Get the variable name (from the header: first row in the data
            # file, which is now the first index of every vector in raw).
            var = raw[j][0]
            # Replace all occurences of None by NaN, which can be converted to
            # floats automatically.
            raw[j][1:][raw[j][1:]=='None'] = 'nan'
            # XYloc values have been stored as strings, e.g. '[518, 724]'.
            # These need to be converted to NumPy arrays so that we can process
            # them more easily.
            if 'XYloc' in var:
                xy = numpy.zeros((len(raw[j][1:]), 2), dtype='|S5')
                for k in range(len(raw[j][1:])):
                    xy[k,:] = numpy.array(raw[j][1:][k].replace('[', '').replace(']', '').split(', '))
                sesdata[pname][var] = xy.astype(float)
                continue
            # Get the values. Try to convert the str-type to float; and just
            # use them as-is if the conversion fails (they're probably not
            # numbers if the float conversion fails).
            try:
                val = raw[j][1:].astype(float)
            except:
                val = raw[j][1:]
            sesdata[pname][var] = copy.deepcopy(val)
        # Store this pair.
        p = '%d_%d-%d' % (i, min([ip, ip_other]), max([ip, ip_other]))
        if p not in pairs:
            pairs.append(p)

# Store the extracted data.
print("Storing extracted data in a pickle file.")
with open(os.path.join(PPDIR, 'qdata.pickle'), 'wb') as f:
    pickle.dump(qdata, f)
with open(os.path.join(PPDIR, 'wmdata.pickle'), 'wb') as f:
    pickle.dump(wmdata, f)
with open(os.path.join(PPDIR, 'sesdata.pickle'), 'wb') as f:
    pickle.dump(sesdata, f)


# # # # #
# PROCESS WM DATA

print("\nComputing individuals' short-term memory capacity.")

# Create a new dict to hold all the individual data in.
idata = {}
for nstim in WM_NSTIM:
    for var in ['abserr_%d', 'rawsd_%d', 'kT_%d', 'pT_%d', 'pNT_%d', 'pU_%d', \
        'sd_%d', 'dkl_%d', 'capacity_%d']:
        idata[var % (nstim)] = numpy.ones(len(participants)) * numpy.NaN
idata['ppnames'] = copy.deepcopy(participants)

# Loop through all participants.
for i, ppname in enumerate(participants):
    
    print("\tProcessing participant '%s' (%d/%d)" % \
        (ppname, i+1, len(participants)))
    
    # Keep track of all estimates of the working memory capacity for this
    # participant.
    c = numpy.zeros(len(WM_NSTIM))
    
    # Split the data into conditions.
    for nstim in WM_NSTIM:
        # Create a Boolean vector to select the data from this condition.
        sel = wmdata[ppname]['nstim'] == nstim
        # Compute the average error.
        abserr = numpy.nanmean(numpy.abs(wmdata[ppname]['E'][sel]))
        # Compute the standard deviation in the raw data.
        rawsd = numpy.nanstd(numpy.abs(wmdata[ppname]['E'][sel]))
        # NOTE: Our models require data between 0 (0 deg) and 2*pi (360 deg),
        # while our data is between 0 and 180 degrees. To convert, we multiply
        # we multiply the data by 2 (dirty, but works!), and then convert
        # the result to radians.
        responses = numpy.deg2rad(wmdata[ppname]['X'][sel] * 2.0)
        target_orientations = numpy.deg2rad(wmdata[ppname]['T'][sel] * 2.0)
        # Collect all non-target orientations in one list.
        if FITNONTAR:
            non_target_orientations = []
            for j in range(nstim - 1):
                non_target_orientations.append( \
                    numpy.deg2rad(wmdata[ppname]['NT%d' % j][sel] * 2.0))
        else:
            non_target_orientations = None
        # Fit mixture model to extract guessing and misbinding.
        params, ll = vt.fullspacestimate( \
            responses, \
            target_orientations, \
            non_target_orientations, \
            fituniform=True)
        kT, kNT, pT, pNT, pU = params
        # Compute Kullback-Leibler divergence to quantify short-term memory
        # capacity in bits. In order to do so, we need the probability density
        # function that fits best with the observed data, and one that reflects
        # complete guessing (i.e. a uniform one).
        dx = 0.01
        x = numpy.arange(-numpy.pi, numpy.pi, dx)
        uniform = numpy.ones(len(x)) / (2.0 * numpy.pi)
        vonmises = (pT+pNT) * vt.vonmisespdf(x, 0, kT) + pU * uniform
        dkl = KL_divergence_continuous(vonmises, uniform, dx, mode='bits')
        # Store the computed values.
        idata['abserr_%d' % (nstim)][i] = copy.deepcopy(abserr)
        idata['rawsd_%d' % (nstim)][i] = copy.deepcopy(abserr)
        idata['kT_%d' % (nstim)][i] = copy.deepcopy(kT)
        idata['sd_%d' % (nstim)][i] = copy.deepcopy(vt.kappa2sd(kT))
        idata['pT_%d' % (nstim)][i] = copy.deepcopy(pT)
        idata['pNT_%d' % (nstim)][i] = copy.deepcopy(pNT)
        idata['pU_%d' % (nstim)][i] = copy.deepcopy(pU)
        idata['dkl_%d' % (nstim)][i] = copy.deepcopy(dkl)
        idata['capacity_%d' % (nstim)][i] = copy.deepcopy(dkl * nstim)
    
    # Compute the average working memory capacity for this participant.
    c_avg = (idata['capacity_2'][i] + idata['capacity_4'][i]) / 2.0
    wmdata[ppname]['capacity'] = copy.deepcopy(c_avg)
    if 'capacity' not in idata.keys():
        idata['capacity'] = numpy.ones(len(participants)) * numpy.NaN
    idata['capacity'][i] = copy.deepcopy(c_avg)
    
    # Add the questionnaire data to the individual data.
    if 'apathy' not in idata.keys():
        idata['apathy'] = numpy.ones(len(participants)) * numpy.NaN
        idata['apathy_ba'] = numpy.ones(len(participants)) * numpy.NaN
        idata['apathy_sm'] = numpy.ones(len(participants)) * numpy.NaN
        idata['apathy_es'] = numpy.ones(len(participants)) * numpy.NaN
    idata['apathy'][i] = copy.deepcopy(qdata['ami'][ppname]['total'])
    idata['apathy_ba'][i] = copy.deepcopy(qdata['ami'][ppname]['ba'])
    idata['apathy_sm'][i] = copy.deepcopy(qdata['ami'][ppname]['sm'])
    idata['apathy_es'][i] = copy.deepcopy(qdata['ami'][ppname]['es'])

    if 'B5_extraversion' not in idata.keys():
        idata['B5_extraversion'] = numpy.ones(len(participants)) * numpy.NaN
        idata['B5_agreeableness'] = numpy.ones(len(participants)) * numpy.NaN
        idata['B5_conscientiousness'] = numpy.ones(len(participants)) * numpy.NaN
        idata['B5_emotional_stability'] = numpy.ones(len(participants)) * numpy.NaN
        idata['B5_openness'] = numpy.ones(len(participants)) * numpy.NaN
    idata['B5_extraversion'][i] = copy.deepcopy(qdata['bigfive'][ppname]['extraversion'])
    idata['B5_agreeableness'][i] = copy.deepcopy(qdata['bigfive'][ppname]['agreeableness'])
    idata['B5_conscientiousness'][i] = copy.deepcopy(qdata['bigfive'][ppname]['conscientiousness'])
    idata['B5_emotional_stability'][i] = copy.deepcopy(qdata['bigfive'][ppname]['emotional_stability'])
    idata['B5_openness'][i] = copy.deepcopy(qdata['bigfive'][ppname]['openness'])

    # Peak at the pair data to compute average horizontal and vertical
    # location of stimuli claimed by this individual in all runs, as well as
    # their rankings of the collaboration.
    this_run, this_ip = ppname.split('_')
    n = 0
    x = []
    y = []
    collab = { \
        'fair':[], \
        'good':[], \
        'likeable':[], \
        }
    # Loop through all sessions.
    for sesname in sesdata.keys():
        # Get the current session details.
        run_nr, ips = sesname.split('_')
        players = ips.split('-')
        # Only process sessions where this participant was the first player.
        if (run_nr == this_run) and (this_ip == players[0]):
            # Store the collaboration values.
            for var in collab.keys():
                collab[var].append(float(sesdata[sesname][var]))
            # Find the locations of all items that this player claimed.
            n = 0
            x = []
            y = []
            # Loop through all stimuli (stimulus number = k).
            for k in range(NSTIM):
                # Find all trials in which this participant fixated stimulus k.
                sel = sesdata[sesname]['self_fixated%d' % k].astype(bool)
                # Sum all stimulus locations, and count the number of stimuli for
                # this player.
                n += numpy.sum(sel.astype(int))
                x.append(numpy.mean(sesdata[sesname]['XYloc%d' % k][sel][:,0]))
                y.append(numpy.mean(sesdata[sesname]['XYloc%d' % k][sel][:,1]))
    # Add the collaboration perception data to the individual data.
    if 'fair' not in idata.keys():
        idata['fair'] = numpy.ones(len(participants)) * numpy.NaN
        idata['good'] = numpy.ones(len(participants)) * numpy.NaN
        idata['likeable'] = numpy.ones(len(participants)) * numpy.NaN
        idata['collab'] = numpy.ones(len(participants)) * numpy.NaN
    for var in collab.keys():
        idata[var][i] = numpy.mean(collab[var])
    idata['collab'][i] = (idata['fair'][i] + idata['good'][i] + idata['likeable'][i]) / 3.0
    # Add the location data to the individual data.
    if 'avgx' not in idata.keys():
        idata['nclaimed'] = numpy.ones(len(participants)) * numpy.NaN
        idata['avgx'] = numpy.ones(len(participants)) * numpy.NaN
        idata['avgy'] = numpy.ones(len(participants)) * numpy.NaN
    idata['nclaimed'][i] = float(n)
    idata['avgx'][i] = numpy.nanmean(x)
    idata['avgy'][i] = numpy.nanmean(y)

# Store the pre-processed data.
print("Storing individuals' data in a pickle file.")
with open(os.path.join(PPDIR, 'idata.pickle'), 'wb') as f:
    pickle.dump(idata, f)


# # # # #
# PROCESS PAIR DATA

print("\nComputing game performance for pairs.")

# Create an empty dict to hold data for all the pairs.
pdata = {'player_high':{}, 'player_low':{}}
varnames = idata.keys()
varnames.extend(['name', \
    'collab', 'good', 'likeable', 'fair', \
    'nclaimed', 'earned', 'rewprop', \
    'abserr', 'estcap', 'kT', 'pT', 'pNT', 'pU', 'sd', 'dkl', \
    'avgx', 'avgy', 'varx', 'vary', \
    'apathy', 'apathy_ba', 'apathy_sm', 'apathy_es', \
    'B5_extraversion', 'B5_agreeableness', 'B5_conscientiousness', 'B5_emotional_stability', 'B5_openness', \
    ])
for var in varnames:
    for player in pdata.keys():
        if var in ['name']:
            pdata[player][var] = numpy.ones(len(pairs), dtype='|S8')
        else:
            pdata[player][var] = numpy.ones(len(pairs)) * numpy.NaN
diffvars = ['d_capacity', 'd_nclaimed', 'd_earned', 'd_rewprop', \
    'd_abserr', 'd_estcap', 'd_collab', 'd_good', 'd_likeable', 'd_fair', \
    'd_avg', 'd_var', 'd_avgx', 'd_avgy', 'd_varx', 'd_vary', \
    'm_capacity', 'm_earned', 'm_abserr', \
    'm_collab', 'm_good', 'm_likeable', 'm_fair', \
    'd_apathy', 'd_apathy_ba', 'd_apathy_sm', 'd_apathy_es', \
    'd_B5_extraversion', 'd_B5_agreeableness', 'd_B5_conscientiousness', 'd_B5_emotional_stability', 'd_B5_openness', \
    ]
for var in diffvars:
    if var not in pdata.keys():
        pdata[var] = numpy.ones(len(pairs)) * numpy.NaN
pdata['cont'] = {}
pdata['player_high']['cont'] = {}
pdata['player_low']['cont'] = {}
icontvars = ['nclaimed', 'abserr', 'earned', 'avgx', 'avgy', 'varx', 'vary']
for var in icontvars:
    pdata['player_high']['cont'][var] = numpy.ones((len(pairs), NTRIALS)) * numpy.NaN
    pdata['player_low']['cont'][var] = numpy.ones((len(pairs), NTRIALS)) * numpy.NaN
contvars = ['d_nclaimed', 'd_abserr', 'd_earned', 'd_avg', 'd_avgx', 'd_avgy']
for var in contvars:
    pdata['cont'][var] = numpy.ones((len(pairs), NTRIALS)) * numpy.NaN

# Loop through all pairs.
for i, pname in enumerate(pairs):
    
    print("\tProcessing pair %s (%d/%d)" % (pname, i+1, len(pairs)))
    
    # Parse the pair name to get the individual player names.
    run_nr, players = pname.split('_')
    players = players.split('-')
    
    # Check which of the players had the highest SORTBY value.
    ppname = ['%s_%s' % (run_nr, players[0]), '%s_%s' % (run_nr, players[1])]
    ppi = [participants.index(ppname[0]), participants.index(ppname[1])]
    if idata[SORTBY][ppi[0]] > idata[SORTBY][ppi[1]]:
        ph = copy.deepcopy(players[0])
        pl = copy.deepcopy(players[1])
    else:
        ph = copy.deepcopy(players[1])
        pl = copy.deepcopy(players[0])
    # Re-order the players by their capacity.
    players = [ph, pl]

    # Get the data from both players.
    for j, ip in enumerate(players):
        if j == 0:
            ptype = 'high'
        elif j == 1:
            ptype = 'low'
        # Compose the participant name, and the session name.
        ppname = '%s_%s' % (run_nr, ip)
        ppi = participants.index(ppname)
        sesname = '%s_%s-%s' % (run_nr, ip, players[1-j])
        # Store the individual data (e.g. short-term memory capacity).
        for var in idata.keys():
            if var != 'ppnames':
                pdata['player_%s' % (ptype)][var][i] = copy.deepcopy(idata[var][ppi])
        # Store the participant's name.
        pdata['player_%s' % (ptype)]['name'][i] = copy.deepcopy(ppname)
        # Find the locations of all items that this player claimed.
        n = 0
        x = []
        y = []
        xy = numpy.zeros((NSTIM, NTRIALS, 2)) * numpy.NaN
        # Loop through all stimuli (stimulus number = k).
        for k in range(NSTIM):
            # Find all trials in which this participant fixated stimulus k.
            sel = sesdata[sesname]['self_fixated%d' % k].astype(bool)
            # Sum all stimulus locations, and count the number of stimuli for
            # this player.
            n += numpy.sum(sel.astype(int))
            x.append(numpy.mean(sesdata[sesname]['XYloc%d' % k][sel][:,0]))
            y.append(numpy.mean(sesdata[sesname]['XYloc%d' % k][sel][:,1]))
            # Store these locations.
            xy[k,sel,:] = sesdata[sesname]['XYloc%d' % k][sel]
        # Note that we store these values per participant here, in a bit of a
        # crude way (averaged across all items and trials). Further in this
        # script, the per-trial difference in average location is compared
        # between the two individuals in a pair. This is a much more precise
        # way to look at inter-individual differences and variability in the
        # locations of claimed items.
        pdata['player_%s' % (ptype)]['avgx'][i] = numpy.nanmean(x)
        pdata['player_%s' % (ptype)]['avgy'][i] = numpy.nanmean(y)
        pdata['player_%s' % (ptype)]['varx'][i] = numpy.nanstd(x)
        pdata['player_%s' % (ptype)]['vary'][i] = numpy.nanstd(y)

        # Compute the continuous variables (= one value for every trial).
        pdata['player_%s' % (ptype)]['cont']['avgx'][i,:] = numpy.nanmean(xy[:,:,0], axis=0)
        pdata['player_%s' % (ptype)]['cont']['avgy'][i,:] = numpy.nanmean(xy[:,:,1], axis=0)
        pdata['player_%s' % (ptype)]['cont']['varx'][i,:] = numpy.nanstd(xy[:,:,0], axis=0)
        pdata['player_%s' % (ptype)]['cont']['vary'][i,:] = numpy.nanstd(xy[:,:,1], axis=0)
        pdata['player_%s' % (ptype)]['cont']['nclaimed'][i,:] = \
            numpy.copy(sesdata[sesname]['nfixself'])
        pdata['player_%s' % (ptype)]['cont']['earned'][i,:] = \
            numpy.copy(sesdata[sesname]['reward'])
        pdata['player_%s' % (ptype)]['cont']['abserr'][i,:] = \
            numpy.abs(sesdata[sesname]['X'])
        # Store the data for this player's game performance: Several variables
        # will be stored for each individual player, and then compaired. The
        # difference between these variables is taken as the data for one pair.
        # Pairs are treated as independent data points.
        # Grab the ratings of this player on the other.
        collab = 0.0
        for var in ['good', 'likeable', 'fair']:
            pdata['player_%s' % (ptype)][var][i] = float(sesdata[sesname][var])
            collab += float(sesdata[sesname][var])
        # Compute the average of all three questions as an aggregated
        # collaboration score.
        pdata['player_%s' % (ptype)]['collab'][i] = collab / 3.0
        # The average number of claimed stimuli in this game.
        pdata['player_%s' % (ptype)]['nclaimed'][i] = \
            numpy.nanmean(sesdata[sesname]['nfixself'])
        # The total earned reward by this player.
        pdata['player_%s' % (ptype)]['earned'][i] = \
            numpy.nansum(sesdata[sesname]['reward'])
        # The reward earned by this player as a proportion of the total amount
        # of reward earned by this pair (which is the total for this individual
        # multiplied by two, as pair earnings were summed and divided by two).
        pdata['player_%s' % (ptype)]['rewprop'][i] = \
            numpy.nansum(sesdata[sesname]['reward']) \
            / (2*sesdata[sesname]['total'][-1])
        # The average absolute error by this player.
        pdata['player_%s' % (ptype)]['abserr'][i] = \
            numpy.nanmean(numpy.abs(sesdata[sesname]['E']))
        # Fit this player's performance to a mixture model, to extract
        # guessing and misbinding.
        not_nan = numpy.isnan(sesdata[sesname]['X']) == False
        responses = numpy.deg2rad(sesdata[sesname]['X'][not_nan] * 2.0)
        target_orientations = numpy.deg2rad(sesdata[sesname]['T'][not_nan] * 2.0)
        non_target_orientations = []
        for k in range(NSTIM-1):
            non_target_orientations.append( \
                numpy.deg2rad(sesdata[sesname]['NT%d' % k][not_nan] * 2.0))
        params, ll = vt.fullspacestimate( \
            responses, \
            target_orientations, \
            non_target_orientations, \
            fituniform=True)
        kT, kNT, pT, pNT, pU = params
        pdata['player_%s' % (ptype)]['kT'][i] = kT
        pdata['player_%s' % (ptype)]['sd'][i] = vt.kappa2sd(kT)
        pdata['player_%s' % (ptype)]['pT'][i] = pT
        pdata['player_%s' % (ptype)]['pNT'][i] = pNT
        pdata['player_%s' % (ptype)]['pU'][i] = pU
        # Estimate this player's short-term memory capacity based on the
        # fitted propability density function.
        dx = 0.01
        x = numpy.arange(-numpy.pi, numpy.pi, dx)
        uniform = numpy.ones(len(x)) / (2.0 * numpy.pi)
        vonmises = (pT+pNT) * vt.vonmisespdf(x, 0, kT) + pU * uniform
        pdata['player_%s' % (ptype)]['dkl'][i] = \
            KL_divergence_continuous(vonmises, uniform, dx, mode='bits')
        pdata['player_%s' % (ptype)]['estcap'][i] = pdata['player_%s' % (ptype)]['dkl'][i] \
            * pdata['player_%s' % (ptype)]['nclaimed'][i]
        
    # Compute the differences between the players.
    for var in diffvars:
        # If the variable is d_avg, calculate it.
        if var == 'd_avg':
            # Compute per-trial difference, and then average (good way, in
            # case players don't stick to consistent sides).
            dx = pdata['player_high']['cont']['avgx'][i,:] - pdata['player_low']['cont']['avgx'][i,:]
            dy = pdata['player_high']['cont']['avgy'][i,:] - pdata['player_low']['cont']['avgy'][i,:]
            pdata[var][i] = numpy.nanmean(numpy.sqrt(dx**2 + dy**2))
            continue
        elif var == 'd_var':
            # Compute per-trial difference, and then compute the variance.
            dx = pdata['player_high']['cont']['avgx'][i,:] - pdata['player_low']['cont']['avgx'][i,:]
            dy = pdata['player_high']['cont']['avgy'][i,:] - pdata['player_low']['cont']['avgy'][i,:]
            pdata[var][i] = numpy.nanstd(numpy.sqrt(dx**2 + dy**2))
            continue
        elif var == 'd_varx':
            # Compute per-trial difference, and then compute the variance.
            dx = pdata['player_high']['cont']['avgx'][i,:] - pdata['player_low']['cont']['avgx'][i,:]
            pdata[var][i] = numpy.nanstd(dx)
            continue
        elif var == 'd_vary':
            # Compute per-trial difference, and then compute the variance.
            dy = pdata['player_high']['cont']['avgy'][i,:] - pdata['player_low']['cont']['avgy'][i,:]
            pdata[var][i] = numpy.nanstd(dy)
            continue
        # If the variable relates to anything else, calculate it.
        if var[0] == 'd':
            ivar = var.replace('d_', '')
            pdata[var][i] = pdata['player_high'][ivar][i] - pdata['player_low'][ivar][i]
        elif var[0] == 'm':
            ivar = var.replace('m_', '')
            pdata[var][i] = (pdata['player_high'][ivar][i] + pdata['player_low'][ivar][i]) / 2.0
        else:
            raise Exception("ERROR: Did not recognise variable '%s'." % var)
    for var in contvars:
        if var == 'd_avg':
            dx = pdata['player_high']['cont']['avgx'] - pdata['player_low']['cont']['avgx']
            dy = pdata['player_high']['cont']['avgy'] - pdata['player_low']['cont']['avgy']
            pdata['cont'][var] = numpy.sqrt(dx**2 + dy**2)
            continue
        ivar = var.replace('d_', '')
        pdata['cont'][var] = pdata['player_high']['cont'][ivar] - pdata['player_low']['cont'][ivar]

# Store the pre-processed data.
print("Storing pairs' data in a pickle file.")
with open(os.path.join(PPDIR, 'pdata.pickle'), 'wb') as f:
    pickle.dump(pdata, f)

print("\nAll done!")
print("\a")
