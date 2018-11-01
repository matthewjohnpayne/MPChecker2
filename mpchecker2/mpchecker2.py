# -*- coding: utf-8 -*-
# mpchecker2/mpchecker2/mpchecker2.py

'''
--------------------------------------------------------------
Main MPChecker2 module.

Oct 2018
Matt Payne

To contain all (high-level) functions required to predict
minor planets within the field-of-view of a pointing

As currently sketched-out 
(i)     Pre-Calculations are handled in mpchecker2/mpchecker2/precalc.py
(ii)    MPCHECKER class takes POINTINGs as input and then
        calculates the SUSPECTS within the FoV
(iii)   IDENTIFY class takes SUSPECTS and compares them to data
        (either EXPOSURE or TRACKLETS) and returns IDENTIFICATIONS

--------------------------------------------------------------
'''


# Import third-party packages
# --------------------------------------------------------------
import sys
import numpy as np


# Import neighboring packages
# --------------------------------------------------------------
import mpchecker2.precalc


# Data classes for MPChecker2
# --------------------------------------------------------------

class POINTING:
    '''
    To contain the pointing (exposure) which is the input to MPChecker

    Will have to contain information like ...
    ... date, obscode, (RA,Dec)/Unit-Vector, radius, ...

    '''

    # Add some content ...
    def __init__(self, something ):
        '''

        Parameters
        ----------

        '''



class SUSPECTS:
    '''
    To contain the MPCHECKER results regarding the minor planets which are
    potentially in the field of view of a POINTING.
    
    For each name/desig, will need to hold some specification of the orbit
    and its associated uncertainty 
     - *** SHOULD DEFINE SOME SIMPLE CLASS TO REFER TO HERE ***
    
    '''

    # Add some content ...
    def __init__(self ):
        '''

        Parameters
        ----------

        '''
        pass


class IDENTIFICATIONS:
    '''
    To contain the IDENTIFY results which map known SUSPECTS to
    data from TRACKLETS / EXPOSURES 
    
    - *** SHOULD DEFINE SOME SIMPLE CLASS TO REFER TO HERE ***
    
    '''
    
    # Add some content ...
    def __init__(self ):
        '''

        Parameters
        ----------

        '''
        pass




# Main functional classes for MPChecker2
# --------------------------------------------------------------

class MPCHECKER:
    '''
    Provide all methods associated with MPChecker
     - Get "shortlists" for night
     - Do detailed orbital advance
     - Get "refined" list of actual overlap with FoV
     - Associated refined list with detections/tracklets
        
    '''
    
    
    def __init__(self, pointing ):
        '''
        Must initialize with a valid POINTING object

        Parameters
        ----------
        pointing    :   POINTING-object or iterable of POINTING-objects
        
        Examples
        --------
        >>> import mpchecker2
        >>> p = POINTING(something)
        >>> m = MPCHECKER(p)
        >>> s = m.suspects

        '''
        
        # My working assumpting is that by default we want to
        # calculate the objects in the FoV of the pointing(s)
        if isinstance(pointing, POINTING):
            self.suspects = get_objects_which_overlap_FoV(pointing)
        elif isinstance(pointing, list) and isinstance(pointing[0], POINTING):
            self.suspects = [ get_objects_which_overlap_FoV(pointing) for pointing in pointings ]
        else:
            sys.exit("Incompatible data types input to MPCHECKER class")


    def get_objects_which_overlap_FoV(pointing):
        '''
        We want to calculate which objects are in (or have some chance of being in) 
        the Fov of a single pointing
        We load precalculated quantities to help us understand which objects have 
        any chance at all of being close to the FoV
        We
        
        Parameters
        ----------
        
        Returns
        -------
        suspects    :   SUSPECTS-object
        
        Examples
        --------
        >>> ...
        
        '''


        # Get precalculated nightly datasets
        # - hp2object contains GEOCENTRIC healpix (i.e., only an approximation to what we need)
        hp2object , objects2func = PRECALC().get_nightly_precalcs( pointing.JDutc )
        
        # Advance all objects of interest to specific time of pointing
        for HP in pointing.HParray :
            for objectName in hp2object[HP]:
                helioPosns = objects2func[objectName][pointing.JDutc]

                # Calculate TOPOCENTRIC sky-plane coordinates of objects at time of pointing
                pass # skyUnitVecs     = ...  # <<-- Something to do with MPCSky
    
                # Calculate whether objects are in the FoV
                # (more specifically, have some part of the their sky-plane-uncert in the FoV)
                # ****** CAREFUL ABOUT *MULTIPLE* skyUnitVecs & SINGLE pointing.UnitVec *******
                angleRad = np.arccos(np.clip( np.dot(skyUnitVecs,pointing.UnitVec) , -1, 1))

        return SUSPECTS(something)








class IDENTIFY:   ## <<-- What is the difference between identification and attribution ?
    '''
    To match SUSPECTS from POINTINGs to
    a set of detections and/or tracklets
    '''
    
    # Add some content ...
    def __init__(self , suspects, data):
        '''
            
        Parameters
        ----------
        suspects    :   single SUSPECTS-object or list of SUSPECTS-objects
            These are the results of running MPChecker
        data        :   single EXPOSURE-object or list of TRACKLET-objects
            These are the data we wish to identify (attribute?) 
            with previously known object-orbits
            
        Returns
        -------
        ids         :   IDENTIFICATIONS-object
        
        
        Examples
        --------
        >>> ...
        
        '''
        # Single detections ...
        if isinstance(suspects, SUSPECTS) and isinstance(data, exposure):
            ids = identify_SUSPECTS_from_multiple_pointings_with_tracklets_from_multiple_exposures(suspects, data)
        # Multiple detections in tracklets ...
        elif isinstance(suspects, list) and isinstance(data, trackletList):
            ids = identify_SUSPECTS_from_multiple_pointings_with_tracklets_from_multiple_exposures(suspects, data)
        # Erroroneous input ...
        else:
            sys.exit("Incompatible data types input to IDENTIFY class")

        # Does any additional reformatting or proocessing need to be done ?
        
        return ids


    def identify_SUSPECTS_from_a_single_pointing_with_detections_from_a_single_exposure(suspects, exposure):
        '''
        
        Parameters
        ----------
        suspects    :   SUSPECTS-object
        exposure    :   EXPOSURE-object
        
        Returns
        -------
        
        Examples
        --------
        >>> ...
        
        '''
            
        return IDENTIFICATIONS()

    def identify_SUSPECTS_from_multiple_pointings_with_tracklets_from_multiple_exposures(suspectsList, trackletList):
        '''
        Identify SUSPECTS from multiple POINTINGS with a list of TRACKLETS
        Dates of the pointings need to correspond to the dates of the tracklet detections

        Parameters
        ----------
        suspectsList    :   list of SUSPECTS-objects
        trackletList    :   list of TRACKLET-objects

        Returns
        -------

        Examples
        --------
        >>> ...

        '''

        return IDENTIFICATIONS()





