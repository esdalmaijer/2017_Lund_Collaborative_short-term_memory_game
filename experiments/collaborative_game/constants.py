# -*- coding: utf-8 -*-

import os
import sys
import time

from helpers import ang2cm, cm2ang

# Communications
SERVERIP = '192.168.1.28'
CLIENIPBASE = '192.168.1.%d'
#CLIENTIPS = [2, 5, 8, 9, 12, 15, 17, 22, 23, 24]
#CLIENTIPS = [2, 5, 8, 9, 12, 15, 17, 24]
#CLIENTIPS = [2, 5, 8, 12, 15, 17]
#CLIENTIPS = [2, 5, 8, 12]
#CLIENTIPS = [2, 5]
#CLIENTIPS = [22, 23]
#CLIENTIPS = [5, 9]
#CLIENTIPS = [2, 5, 17, 24]
#CLIENTIPS = [5, 12, 17, 24]
CLIENTIPS = [5, 9, 17, 22]
SOCKTIMEOUT = None
CONTIMEOUT = None
BREAKTIMEOUT = None
# IP, PORTNR, AND NPARTICIPANTS AREN'T USED IN MULTICAST SETUP
MYIP = '192.168.1.17'
PORTNR = 7777
NPARTICIPANTS = 2

# Debugging options.
MANUALSTART = False
DEBUG = False
ESCKILL = True
AUTORESP = False

# Number of trials.
NTRIALS = 25
# Number of trials between each break
BREAKFREQ = 60
# Number of stimuli. ONLY USE EVEN NUMBERS!
NSTIM = [8]
# Reward types
REWTYPES = ['collaborate'] #, 'compete']
# Reward amount
# 'incidental' for same reward regardless of performance
# 'contingent' for higher rewards for better performance
REWAMOUNT = 'contingent'
# Maximum and minimum reward in contingent trials
MAXREW = 100.0 # credits
MINREW = 1.0 # credits
# Reward from which upwards a sound is played
REWSOUNDTHRESHOLD = 110
# Standard deviation of the reward-scaling function
REWSD = 17 # degrees

# Payment details
BASEPAY = 6.0 # pounds per hour
MAXADDPAY = 2.0 # pounds per hour
MINPERFORMANCE = 0.5 # proportion of best possible performance

# Display properties
DISPTYPE = 'psychopy'
DISPSIZE = (1680, 1050)
DISPCENTRE = (int(DISPSIZE[0]/2.0), int(DISPSIZE[1]/2.0))
SCREENSIZE = (40.5,30.5)
SCREENDIST = 62.0
PIXPERCM = (DISPSIZE[0]/SCREENSIZE[0] + DISPSIZE[1]/SCREENSIZE[1]) / 2.0
FGC = (0, 0, 0)
BGC = (128, 128, 128)

# Stimulus properties
# Maximum stimulus orientation in degrees
MAXORI = 180
# Stimulus size in visual angle
STIMSIZE = 1.5
# Distance between display centre and stimulus centre.
STIMDIST = 5
# Stimulus outline thickness in pixels.
STIMLINEWIDTH = 3
PROBELINEWIDTH = 7
# Stimulus background colours once they are fixated.
STIMBGCOL = {'self':(200,200,200), 'other':(50,50,50)}
# Gabor properties
STIMSF = 5 # cycles per stimulus
STIMALPHA = 1.0
STIMCONTRAST = 1.0
STIMNOISERES = 32 # resolution of noise patch
# Minimal angular distance between stimulus orientations.
MINORIDIST = 10
# Minimal Cartesian distance between the stimulus positions (in pixels).
MINLOCDIST = 150
# Stimulus colour properties (lightness and radius around whitepoint)
STIML = 50
STIMR = 22
# Type and size (=diameter) of the fixation mark
FIXTYPE = 'cross'
FIXSIZE = 12
# Size of the font on screens with text.
FONTSIZE = 24
# Colour wheel properties
CWL = 50 # colour wheel L
CWR = 22 # colour wheel eccentricity
CWIMGSIZE = 300.0 # pixels; half of the images width/height
# Feedback bar.
FBBARW = DISPSIZE[0] * 0.15
FBBARH = DISPSIZE[1] * 0.65
FBBARPW = 5
# Post-block question properties.
POSTQPOS = (DISPCENTRE[0], DISPSIZE[1]*0.25)
POSTQFONTSIZE = 24
POSTQBARW = DISPSIZE[0] * 0.5
POSTQBARH = DISPSIZE[1] * 0.1
POSTQBARPW = 5
POSTQBARCOL = (255,0,0)

# Timing
# FRAMETIME is the duration (in milliseconds) of each frame in the stimulus array.
FRAMETIME = 50
# STIMTIME is the duration (in milliseconds) of the stimulus presentation. If
# set to None, the trial will last until one second after all stimuli have
# been fixated.
STIMTIME = None
# POSTSTIMTIME is how long the stimulus presentation lingers after all the
# stimuli have been fixated.
POSTSTIMTIME = 1000
# REWTIME is the duration (in milliseconds) of the reward screen.
REWTIME = 1000
# MAINTENANCETIME is the dureation of the maintenance phase (between stimulus
# presentation and the response probe.)
MAINTENANCETIME = 3000
# RESPTIMEOUT is the maximum time participants are allowed to take to give a
# response.
RESPTIMEOUT = 7000
# Inter trial interval (time during which a central fixation mark is visible).
ITI = 1250

# Eyetracker properties
MINFIXDUR = 150 # ms
MAXFIXDIST = 2.0 # degrees of visual angle
TRACKERTYPE = 'smi'
DRIFTCHECKFREQ = 20
DUMMYMODE = False

# Files and folders.
# path to the client files
CLIENTDIR = r"\\LAPTOP%d\Share\WinPython-PyGaze-0.5.1\Lund\reward_joint_wm_v2"
# converts the file adress to the proper adress (C:\ ect)
DIR = os.path.dirname(os.path.abspath(__file__))
# folder for all resources
RESDIR = os.path.join(DIR, 'resources')
if not os.path.isdir(RESDIR):
	raise Exception("Could not find resources directory ('%s')" % (RESDIR))
# constructs a new path a folder that containts data
DATADIR = os.path.join(DIR, 'data') 
# checks if the folder already exists
if not os.path.isdir(DATADIR):
    #creates the folder if not
    os.mkdir(DATADIR)
# Construct a path to the DEBUG directory, and create one if it doesn't exist
DEBUGDIR = os.path.join(DIR, 'DEBUG')
if not os.path.isdir(DEBUGDIR):
    os.mkdir(DEBUGDIR)
# Path to the image file for the colour wheel.
CWIMG = os.path.join(RESDIR, 'colour_wheel_L%d_r%d.png' % (CWL, CWR))
# Path to the sound file for reward feedback.
KACHING = os.path.join(RESDIR, 'kaching.wav')
 
# Ask for a name for the new file which is not empty and under 8 characters   
STARTTIME = time.strftime("%Y-%m-%d_%H-%M-%S")
LOGFILENAME = STARTTIME[:]
LOGFILE = os.path.join(DATADIR, LOGFILENAME)
# Set a debug log file.
DEBUGLOG = os.path.join(DEBUGDIR, "stdout_%s.txt" % \
    time.strftime("%Y-%m-%d_%H-%M-%S"))
ERRORLOG = os.path.join(DEBUGDIR, "stderr_%s.txt" % \
    time.strftime("%Y-%m-%d_%H-%M-%S"))

# Redirect STOUT to a file.
sys.stdout = open(DEBUGLOG, 'w')
sys.stderr = open(ERRORLOG, 'w')

# Convert distances and sizes from visual angle to pixels.
STIMSIZE = int(ang2cm(STIMSIZE, SCREENDIST) * PIXPERCM)
STIMDIST = int(ang2cm(STIMDIST, SCREENDIST) * PIXPERCM)
MAXFIXDIST = ang2cm(MAXFIXDIST, SCREENDIST) * PIXPERCM
CWIMGSCALE = 1.0 / (CWIMGSIZE/STIMSIZE)

# tegen Tim's client praten om clients op te starten
STARTCMD = r'C:\Share\WinPython-PyGaze-0.5.1\Lund\reward_joint_wm_v2\RUNME.bat -%s'
UDP_TIM_CLIENT_PORT = 9090