from constants import *
from pygaze.display import Display
from pygaze.screen import Screen
from pygaze.eyetracker import EyeTracker

disp = Display()
tracker = EyeTracker(disp)

scr.draw_text(text="Resetting connection to %s tracker" % (TRACKERTYPE), \
    fontsize=24)
disp.fill(scr)
disp.show()

try:
    tracker.stop_recording()
except:
    print("Could not stop recording")

tracker.close()
disp.close()
