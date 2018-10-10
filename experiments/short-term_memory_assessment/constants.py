import os
from helpers import ang2cm

# Experimental propeties
NSTIM = [2, 4]
TRIALSPERCELL = 40
AUTORESP = False

# Stimulus properties.
# Maximum stimulus orientation in degrees
MAXORI = 180
# Stimulus size in visual angle
STIMSIZE = 1.5
# Distance between display centre and stimulus centre.
STIMDIST = 5
# Stimulus outline thickness in pixels.
STIMLINEWIDTH = 3
PROBELINEWIDTH = 7
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
# Fixation properties.
FIXPW = 3
FIXTYPE = 'cross'
FIXSIZE = 12

# Timing.
FIXTIME = [1000, 1500]
STIMTIME = 5000
MAINTENANCETIME = 3000

# Display properties
DISPTYPE = 'psychopy'
DISPSIZE = (1920, 1080) #(1366, 768)
DISPCENTRE = (int(DISPSIZE[0]/2.0), int(DISPSIZE[1]/2.0))
SCREENSIZE = (40.5,30.5)
SCREENDIST = 62.0
PIXPERCM = (DISPSIZE[0]/SCREENSIZE[0] + DISPSIZE[1]/SCREENSIZE[1]) / 2.0
FGC = (0, 0, 0)
BGC = (128, 128, 128)
FONTSIZE = 24

# path to the actual file
# converts the file adress to the proper adress (C:\ ect)
DIR = os.path.dirname(os.path.abspath(__file__))
# constructs a new path a folder that containts data
DATADIR = os.path.join(DIR, 'data') 
# checks if the folder already exists
if not os.path.isdir(DATADIR):
    #creates the folder if not
    os.mkdir(DATADIR)
 
# Ask for a name for the new file which is not empty and under 8 characters   
LOGFILENAME = ''
while LOGFILENAME == '' or len(LOGFILENAME) > 8:
    LOGFILENAME = raw_input("What's the file name? ")
LOGFILE = os.path.join(DATADIR, LOGFILENAME)

# Convert visual angle into pixels.
STIMSIZE = ang2cm(STIMSIZE, SCREENDIST) * PIXPERCM
STIMDIST = ang2cm(STIMDIST, SCREENDIST) * PIXPERCM

# Set times to lower numbers during AUTORESP.
if AUTORESP:
    FIXTIME = [10, 50]
    STIMTIME = 50
    MAINTENANCETIME = 50
