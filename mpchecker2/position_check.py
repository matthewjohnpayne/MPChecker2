# -*- coding: utf-8 -*-
# mpchecker2/mpchecker2/position_check.py

'''
    --------------------------------------------------------------
    Main module for position_check function.
    
    Feb 2020
    Matt Payne
    
    To contain all core functionality required to evaluate the equivalent of the old "pcheck"
    functionality. I.e. comparing reported observations of an object to the
    predicted ephemeris position for the same object and evaluating the residuals. 
    
    The basic goal is (I think) to have some estimate of whether the observastions
    really are of the claimed object.
    
    Is part of the overall cheby/check/mpchecker infrastructure:
     - See https://drive.google.com/open?id=1F86lnaHL01xPAACX2SYVdOgvvLLnxfYD

    --------------------------------------------------------------
'''


# Import third-party packages
# --------------------------------------------------------------
import sys
import numpy as np


# Import neighboring packages
# --------------------------------------------------------------
import sql


# Data classes for Ephemeris
# --------------------------------------------------------------

class PCHECK():
    '''
        Provide functionalities that wrap around orbit_cheby / ephemeris_cheby
        Allows the user to compare the observed (RA,Dec) with the expected 
        (RA,Dec) from orbital ephemeris
    '''

    def __init__(self, *args ):
        '''
            ...
        '''
        # (i) initialize multi-sector-cheby-dict
        if   FROM_PROVID in args and provid in args and FROM_PROVID:
            self.mscd = self._get_cheby(args)
        elif FROM_CHEBY in args  and mscd in args   and FROM_CHEBY :
            self.mscd = args['mscd']
        else:
            self.mscd = None
            self._generate_empty():

        # (ii) if possible, do main position-residual-check
        if obs in args and self.mscd != None:
            self.calculate_residuals( mscd , args )


    def _get_cheby():
        '''
            We might need to have a wrapper around some sql call to the 
            database of pre-calculated Multi-Sector Cheby Dicts
            
            inputs:
            -------
            provisional_id : string 
             - name of the object for which you want to retrieve an orbit-cheby 
             
            returns:
            --------
            MSCD : Multi-Sector Cheby Dict [see orbit_cheby]
            
        '''
        pass

    def calculate_residuals( mscd , args ):
        '''
            ...
        '''
        pass 
