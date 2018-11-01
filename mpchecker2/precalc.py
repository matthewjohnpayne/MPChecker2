# -*- coding: utf-8 -*-
# mpcadvancer/mpcadvancer/mpcadvancer.py

'''
--------------------------------------------------------------
Precalculation routines for MPChecker2

Oct 2018
Matt Payne

We will want to facilitate preprocessing of various orbit-related
quantities to speed up calculations for MPChecker2

Brain-dump:
(i) There needs to exist some kind of master-list of orbit-objects
(ii) They need to have some kind of fitted state at epoch (or many epochs)
(iii) We want to be able to propagate states to new epochs
(iv) We want to produce a local (~nightly?) poly-func (cheby?) 
     that will allow for rapid positional evaluation during the night
     for each object
(v) We want to calculate the list of HP touched by each object during
    the night (needs to include uncert. region)
(vi) We want to invert to produce a list of objects-per-HP

--------------------------------------------------------------
'''


# Import third-party packages
# --------------------------------------------------------------
import sys
import numpy as np


# Import neighboring packages
# --------------------------------------------------------------



# Class(es) / Routines for MPChecker2
# --------------------------------------------------------------


class PRECALC:
    '''
    MPChecker2 will work faster if precalculations can be done
    We want to calculate (and save) nightly versions of:
    (a) poly-func (cheby?) for calculating heliocentric positions
    (b) lists of objects-per-HP (geocentric?)
    
    Once I have decide where funcstions wil live ...
    ... *** "self" needs to be added in a lot of places below ***
    
    
    '''

    # Add some content ...
    def __init__(self):
        '''

        Parameters
        ----------

        '''


    def get_nightly_precalcs(JDutc):
        '''
        
        Parameters
        ----------
        JDutc           : float or int
        JD Specification of the night. Gets rounded to ...
        
        Returns
        -------
        hp2object       : dictionary (or perhaps list) of lists
        For an input healpix-integer, it contains a list of the objects
        which touched / may-have-touched the healpix in the stated night
        objects2func    : dictionary (or perhaps list) of functions
        For an input object, a function is returned to provide
        rapid evaluation of heliocentric position for this night.
        E.g. cheby/poly-func
        
        Examples
        --------
        >>> ...
        
        '''
        # Ensure JDutc is in the desired form
        JDutc = self.JDutcForNight(JDutc)

        # Don't know whether I really need this check-and-create function
        try:
            objects2func = read_nightlyPolynomialFunctions(JDutc)
        except:
            objects2func = create_polynomialfunctions_to_allow_position_calculations_per_object(JDutc)
            write_nightlyPolynomialFunctions(JDutc, objects2func)


        # Don't know whether I really need this check-and-create function
        try:
            hp2object   =   read_nightlyListsOfObjectsPerHP(JDutc)
        except:
            hp2object   = produce_a_list_of_HP_per_object_and_then_invert_to_produce_a_list_of_objects_per_HP(JDutc , objects2func)
            write_nightlyListsOfObjectsPerHP(JDutc, hp2object)

        return objects2func , hp2object




    def produce_a_list_of_HP_per_object_and_then_invert_to_produce_a_list_of_objects_per_HP(JDutc, objects2func ):
        '''
        For each object ...
        ... advance the object as a function of time
        ... calculate the HP it touches 
        
        The "touch" needs to incorporate the uncertainty, 
        not just the nominal orbit
        
        Need to refine my thinking regarding 
        ... what time-span (12? 24? hours) 
        ... what time grid is necessary ?
        ... do I want to allow this to NOT be the geocenter, and instead be for a specific observatory? 
        
        Is there any need to treat NEOs differently?
        
        Initial sketch-out below seems likely to have a lot 
        of list-loops: need to make these more efficient if possible
            
        Parameters
        ----------
        JDutc           :   float or int
            JD to specify the night of interest
        objects2func    :
            Polynomialfunction for each object
            Allows calculations of heliocentric position
        
        Returns
        -------
        hp2object    :   dictionary (or perhaps list) of lists
            For each healpix-integer, it contains a list of the objects
            which touched / may-have-touched the healpix in the stated night
        
        Examples
        --------
        >>> ...
        
        '''
        # Ensure JDutc is in the desired form
        JDutc = self.JDutcForNight(JDutc)

        # Initialize the list-of-lists we aim to populate
        N = 768 # <<-- to be linked to standard healpix spec param for MPC stuff, e.g. hp.nside2npix(8)
        listObjectsPerHP = [ [] for i in range(N) ]
        
        #(ii) Do we want to evaluate a fixed grid of times ?
        timeArray       = ...
        sys.exit(0)
        
        for objectName in objectNames:
            #(iii) Evaluate heliocentric positions across the night
            helioPosns  = PosnFunctions[objectName][timeArray]
            #(iv) Evaluate Geo/Observatory-Centric healpix across the night
            skyVecs     = ...  # <<-- Something to do with MPCSky
            HParray     = hp.vec2piv(skyVecs)
            #(v) Invert relation and store as a function of healpix
            for HP in HParray:
                listObjectsPerHP[HP].append(objectName)
        
        return listObjectsPerHP



    def write_nightlyPolynomialFunctions(JDutc, objects2func):
        '''
        This may want to live elsewhere ...
        '''
        pass


    def read_nightlyPolynomialFunctions(JDutc):
        '''
        This may want to live elsewhere ...
        '''
        pass


    def write_nightlyListsOfObjectsPerHP(JDutc, hp2object):
        '''
        This may want to live elsewhere ...
        '''
        pass


    def read_nightlyListsOfObjectsPerHP(JDutc):
        '''
        This may want to live elsewhere ...
        '''
        pass


    def function_to_propagate_extant_states_to_new_epoch():
        '''
        This may want to live in MPCAdvancer
        '''
        pass
        
        
        def create_polynomialfunctions_to_allow_position_calculations_per_object(JDutc):
        '''
        This may want to live in MPCAdvancer
        Should this be passed a list of objects ...
        ... or do we want it to default to working on all known/available objects ?
        '''
        pass

