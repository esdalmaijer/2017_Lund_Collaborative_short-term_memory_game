import math

def cm2ang(cm, screendist):
    
    """Converts a distance in centimeters to a distance in degrees of visual
    angle.
    
    Arguments
    
    cm              -       Float or interger of a size in centimeters.

    screendist      -       Float or interger of the distance between the
                            monitor and the participant in centimeters.
    
    Returns
    
    ang             -       Float; size in visual angle.
    """
    
    ang = 2.0 * math.degrees(math.atan( (cm/100.0) / (2.0 * (screendist/100.0)) ))
    
    return ang

def ang2cm(ang, screendist):
    
    """Converts a distance in centimeters to a distance in degrees of visual
    angle.
    
    Arguments
    
    ang             -       Float; size in visual angle.

    screendist      -       Float or interger of the distance between the
                            monitor and the participant in centimeters.
    
    Returns
    
    cm              -       Float or interger of a size in centimeters.
    """
    
    cm = math.tan(math.radians(ang) / 2.0) * (2*screendist)
    
    return cm