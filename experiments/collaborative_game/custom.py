# -*- coding: utf-8 -*-

import math
import random

from constants import *
import pygaze
from pygaze.screen import Screen
from pygaze._misc.misc import pos2psychopos, rgb2psychorgb

import numpy
from psychopy.visual import Circle, GratingStim, ImageStim, Rect


# # # # #
# FUNCTIONS

def pol2car(cor, mode='degrees'):
    
    """Returns the cartesian coordinates corresponding to given polar coordinates
    (see circle below)
    
            0
        315        45
        
    270                90
    
        225        135
            180
    
    arguments
    cor        -- a single or a list of (distance, angle) polar TUPLES
    
    keyword arguments
    mode        -- inicates whether angles are in degrees or radians (default
               = 'degrees')
    
    returns
    cor        -- a single of a list of (x, y) Cartesian tuples
    """
    
    # polar coordinate
    if type(cor) == tuple:
        cor = [cor]
    
    # Cartesian coordinates
    newcor = []
    for dist, ang in cor:
        if mode == 'degrees':
            ang = math.radians(ang-90)
        else:
            ang = ang - math.radians(90)
        x = (math.cos(ang) * float(dist)) + DISPSIZE[0]/2.0
        y = (math.sin(ang) * float(dist))  + DISPSIZE[1]/2.0
        newcor.append((int(x),int(y)))
    
    # single tuple or list
    if len(newcor) == 1:
        newcor = newcor[0]
    
    return newcor

def calc_angular_dist(ori1, ori2, zero=360):
    
    """Calculates the angular distance between two angles (in DEGREES!),
    each angle representing a point on an imaginary circle (default). Note
    that the origin of the round distribution can be determined by setting the
    'zero' keyword.
    
    arguments

    ori1    -    first orientation(s) (the target orientation), in degrees
    ori2    -    second orientation (the response), in degrees
    
    keyword arguments
    
    zero    -    origin of the distribution (360 for a full circle, 180 for lines
            along the diameter of a circle)
    """
    
    if type(ori1) in [int,float, numpy.int, numpy.int16, numpy.int32, numpy.int64, numpy.float, numpy.float16, numpy.float32, numpy.float64]:
        ori1 = [ori1]
    elif type(ori1) == numpy.ndarray and ori1.shape == ():
        ori1 = [ori1]
    ori1 = numpy.array(ori1)
    if type(ori2) in [int,float, numpy.int, numpy.int16, numpy.int32, numpy.int64, numpy.float, numpy.float16, numpy.float32, numpy.float64]:
        ori2 = [ori2]
    elif type(ori2) == numpy.ndarray and ori2.shape == ():
        ori2 = [ori2]
    ori2 = numpy.array(ori2)
    
    # moduli after division
    ori1 = ori1 % zero
    ori2 = ori2 % zero
    
    # calculate the difference
    error = ori2 - ori1
    # where the difference is larger then a clockwise 90 degree rotation, it
    # should be counted as a counter-clockwise (negative) rotation
    error[error>zero/2] = -1 * (zero - error[error>zero/2])
    error[error<=-zero/2] = (zero + error[error<=-zero/2])
    
    return error

def generate_locs(nlocs, mindist, edgeprop=0.1):
    
    """Generates a list of (x,y) Cartesian coordinates in pixels that are
    randomly distributed over the display, with a minimum distance between
    them.
    
    Arguments
    
    nlocs           -   Integer that determines the number of orientations
                        to be generated.

    mindist         -   Integer that determines the minimal Cartesian
                        distance between locations.

    Keyword Arguments
    
    edgeprop        -   Float that determines the proportion of the screen
                        around the edges that cannot be used to place stimulus
                        locations in.
    """
    
    # Calculate the minimum and maximum coordinates.
    xmin = int(0 + edgeprop * DISPSIZE[0])
    xmax = int(DISPSIZE[0] - edgeprop * DISPSIZE[0])
    ymin = int(0 + edgeprop * DISPSIZE[1])
    ymax = int(DISPSIZE[1] - edgeprop * DISPSIZE[1])
    
    # Generate an empty list of locations.
    stimx = numpy.zeros(nlocs) * numpy.NaN
    stimy = numpy.zeros(nlocs) * numpy.NaN
    locs = []
    
    # Keep trying to add a random location until the list is of the requested
    # size, or until the number of attempts is OVER 9000!
    attempt = 0
    max_attempts = 9000
    while len(locs) < nlocs:
        # Update the attempt counter.
        attempt += 1
        if attempt > max_attempts:
            raise Exception("ERROR in custom.generate_locs: Used to many attempts to generate %d locations with a minimal distance of %d" \
                % (nlocs, mindist))
        # Generate a random coordinate that is not near the centre.
        tooclose = True
        while tooclose:
            x = random.randint(xmin, xmax)
            y = random.randint(ymin, ymax)
            if numpy.sqrt((DISPCENTRE[0]-x)**2 + (DISPCENTRE[1]-y)**2) > 2*MAXFIXDIST:
                tooclose = False
        # If there is no location yet, add the current one.
        if len(locs) == 0:
            d = numpy.inf
        # Calculate the distance between the new coordinate and all other
        # coordinates.
        else:
            d = (stimx-x)**2 + (stimy-y)**2
            d = numpy.sqrt(d[numpy.isnan(d)==False])
        # Check whether the distance between the new coordinate and all other
        # coordinates is large enough.
        if numpy.min(d) > mindist:
            i = len(locs)
            stimx[i] = x
            stimy[i] = y
            locs.append((x,y))

    return locs

def generate_oris(noris, mindist, zero=360):
    
    """Generates orientations with a certain minimal distance between them. The
    generated orientations will be a list of integers.
    
    Arguments
    
    noris            -      Integer that determines the number of orientations
                            to be generated.

    mindist          -      Integer that determines the minimal circular
                            distance between orientations.

    Keyword Arguments
    
    zero             -      Integer that determines the maximum value in the
                            circular space. If you're using a round space,
                            0 = 360, for example. Default value is 360.

    Returns
    
    oris             -      A list of integers between 0 and the maximum value
                            of your space (set by the zero keyword).
    """
    
    # Generate an empty list to hold the orientations.
    oris = []
    
    # Keep trying to add a new random number, until the list is of the
    # requested size.
    attempt = 0
    max_attempts = 9000
    while len(oris) < noris:
        # Generate a random number.
        o = random.randint(0, zero-1)
        # If there are no orientations yet, add the new one.
        if oris == []:
            oris.append(o)
        # Check if the number is distant enough from all other numbers.
        else:
            # Calculate the distance between the new orientation and all
            # existing orientations.
            d = numpy.abs(calc_angular_dist(numpy.array(oris), o, zero=zero))
            # Check whether the lowest distance is higher than the minimal
            # distance.
            if d.min() >= mindist:
                oris.append(o)
            else:
                attempt += 1
        # Break if we used to many attempts.
        if attempt > max_attempts:
            raise Exception("ERROR in custom.generate_oris: Used to many attempts to generate %d orientations with a minimal distance of %d" \
                % (noris, mindist))
    
    return oris


# # # # #
# STIMULUS CLASS
class StimScreen:
    
    """Custom class for stimulus screens in this experiment."""
    
    def __init__(self, nstim, locs, oris, linewidth=3, \
        stimtypes='gabor', showcolourwheel=False):
        
        """Initialises a new StimScreen instance.
        
        Arguments
        
        nstim        -      Integer that determines the number of stimuli on
                            this screen.
        
        locs         -      List of (x,y) tuples that determine the positions
                            of all stimuli. The list's length should be equal
                            to the number of stimuli as defined by nstim.
        
        oris         -      List of integers that determine the orientations
                            of all stimuli. The list's length should be equal
                            to the number of stimuli as defined by nstim.
        
        Keyword Arguments
        
        linewidth    -      Integer or a list of integers that determines the
                            width of the lines around stimuli (in pixels).
                            Default value is 3.
        
        stimtypes    -      String or a list of strings that determines the
                            type of stimulus. Options are 'gabor' and 'noise'.
                            Default is 'gabor'.
        
        showcolourwheel-    Boolean that determines whether a central colour
                            wheel should be drawn of not. Default is False.
        """
        
        # Settings from the constants.
        self._sf = float(STIMSF) / float(STIMSIZE)
        self._alpha = STIMALPHA
        self._contrast = STIMCONTRAST
        
        # Noise texture.
        self._noisetex = numpy.random.rand(STIMNOISERES,STIMNOISERES)
        
        # Convert the linewidth to a list (if necessary).
        if type(linewidth) in [int, float]:
            linewidth = nstim * [int(linewidth)]

        # Convert the stimulus types to a list (if necessary).
        if type(stimtypes) in [str, unicode]:
            stimtypes = nstim * [stimtypes]
        
        # Create a Screen to use its wonderful drawing functions.
        self.scr = Screen()

        # Draw the fixation cross.
        self.scr.draw_fixation(fixtype=FIXTYPE, diameter=FIXSIZE)
        
        # Draw the colouw wheel
        if showcolourwheel:
            # Load the image.
            self.scr.draw_image(CWIMG, scale=CWIMGSCALE)
        
        # Create an internal list of stimuli (=PsychoPy stimulus instances) by
        # copying it from the internal screen.
        self.screen = self.scr.screen
        
        # Draw the backgrounds to the stimuli, which will appear in the
        # background colour initially, but will turn to a different colour to
        # mark that they are fixated by a participant.
        # Keep a list of the index numbers of all stimuli. The indices refer
        # to positions within the self.scr.screen list of PsychoPy stimuli.
        self._bgindexnrs = []
        # Draw the stimuli.
        for i in range(nstim):
            # Add the stimulus' index number to the list of indices.
            self._bgindexnrs.append(len(self.screen))
            # Create a new Circle stimulus instance.
            stim = Circle(pygaze.expdisplay, \
                radius=STIMSIZE, \
                edges=64, \
                pos=pos2psychopos(locs[i]), \
                fillColor=rgb2psychorgb(BGC), \
                lineColor=rgb2psychorgb(BGC), \
                lineWidth=0)
            # Add the new stimulus to our list of stimuli.
            self.screen.append(stim)

        # Keep a list of the index numbers of all stimuli. The indices refer
        # to positions within the self.scr.screen list of PsychoPy stimuli.
        self._stimindexnrs = []
        self._outlineindexnrs = []
        # Draw the stimuli.
        for i in range(nstim):
            # Add the stimulus' index number to the list of indices.
            self._stimindexnrs.append(len(self.screen))
#            # Create a new Rect stimulus instance.
#            stim = Rect(pygaze.expdisplay, \
#                pos=pos2psychopos(locs[i]), \
#                fillColor=rgb2psychorgb(list(oris[i])), \
#                lineColor=rgb2psychorgb(list(linecols[i])), \
#                lineWidth=linewidth[i], \
#                width=STIMSIZE, \
#                height=STIMSIZE)
            # Create a Gabor-ish GratingStim.
            if stimtypes[i] == 'gabor':
                tex = 'sin'
            else:
                tex = self._noisetex
            stim = GratingStim(pygaze.expdisplay, \
                pos=pos2psychopos(locs[i]), \
                ori=oris[i], \
                size=STIMSIZE, \
                sf=self._sf, \
                opacity=self._alpha, \
                contrast=self._contrast, \
                tex=tex, \
                mask='circle', \
                color=(1,1,1)
                )
            # Add the new stimulus to our list of stimuli.
            self.screen.append(stim)

            # Add an outline for the stimulus.
            self._outlineindexnrs.append(len(self.screen))
            stim = Circle(pygaze.expdisplay, \
                pos=pos2psychopos(locs[i]), \
                lineWidth=linewidth[i], \
                radius=STIMSIZE//2, \
                edges=32, \
                closeShape=False, \
                fillColor=None, \
                lineColor=(0,0,0)
                )
            # Add the new stimulus to our list of stimuli.
            self.screen.append(stim)

    
    def update(self, locs, oris, bgcols, linewidth=None, stimtypes=None):

        """Updates the locations, colours, line colours, and line widths of
        this stimulus array.
        
        Arguments
        
        locs         -      List of (x,y) tuples that determine the positions
                            of all stimuli. The list's length should be equal
                            to the number of stimuli as defined by nstim.
        
        oris         -      List of integers that determine the orientations
                            of all stimuli. The list's length should be equal
                            to the number of stimuli as defined by nstim.
        
        bgcols     -        List of (r,g,b) tuples that determine the colours
                            of the background to all stimuli. The list's length
                            should be equal to the number of stimuli as
                            defined by nstim.

        Keyword Arguments
        
        linewidth    -      Integer or a list of integers that determines the
                            width of the lines around stimuli (in pixels), or
                            None to leave the width as it is. Default value is
                            None.
        
        stimtypes    -      String or a list of strings that determines the
                            type of stimulus. Options are 'gabor' and 'noise',
                            or None to not update.
                            Default is None.
        """
        
        # Convert the linewidth to a list (if necessary).
        if type(linewidth) in [int, float]:
            linewidth = len(self._stimindexnrs) * [int(linewidth)]

        # Convert the stimulus types to a list (if necessary).
        if type(stimtypes) in [str, unicode]:
            stimtypes = len(self._stimindexnrs) * [stimtypes]

        # Loop through all stimuli.
        #     stimnr is a number between 0 and nstim
        #     stimindexnr refers to the index of a stimulus in self.screen
        for stimnr, stimindexnr in enumerate(self._stimindexnrs):
            
            # Update the stimulus location.
            self.screen[stimindexnr].pos = pos2psychopos(locs[stimnr])
            self.screen[self._bgindexnrs[stimnr]].pos = pos2psychopos(locs[stimnr])
            self.screen[self._outlineindexnrs[stimnr]].pos = pos2psychopos(locs[stimnr])

#            # Update the stimulus colour.
#            self.screen[stimindexnr].fillColor = rgb2psychorgb(list(oris[stimnr]))

            # Update the stimulus orientation.
            self.screen[stimindexnr].setOri(oris[stimnr])

            # Update the stimulus background colour
            self.screen[self._bgindexnrs[stimnr]].fillColor = rgb2psychorgb(list(bgcols[stimnr]))
            self.screen[self._bgindexnrs[stimnr]].lineColor = rgb2psychorgb(list(bgcols[stimnr]))
        
            # Update the stimulus line width and colour.
            if linewidth != None:
                self.screen[self._outlineindexnrs[stimnr]].lineWidth = linewidth[stimnr]
                if linewidth[stimnr] == PROBELINEWIDTH:
                    self.screen[self._outlineindexnrs[stimnr]].lineColor = \
                        (1,-1,-1)
                else:
                    self.screen[self._outlineindexnrs[stimnr]].lineColor = \
                        (0,0,0)
            # Update the stimulus texture.
            if stimtypes != None:
                if stimtypes[stimnr] == 'gabor':
                    self.screen[stimindexnr].setTex('sin')
                else:
                    self.screen[stimindexnr].setTex(self._noisetex)
