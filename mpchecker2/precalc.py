# -*- coding: utf-8 -*-
# mpchecker2/mpchecker2/precalc.py

'''
    --------------------------------------------------------------
    MPChecker2's precalculation module.
    
    Jan 2020
    Matt Payne
    
    This module provides functionalities to
    (a) create & save new precalculations
    (b) load extant precalculations
    
    
    Note that at present *no* attempt is made to ***ALWAYS RETURN NEOs***
     - This needs to be added-to / improved
     - We could (perhaps) return everything closer than 0.X au (on a given JD)
     - Or return everything faster than 0.Y au/day (on a given JD)
     
    --------------------------------------------------------------
    '''


# Import third-party packages
# --------------------------------------------------------------
import sys, os
import numpy as np
import operator
from collections import OrderedDict, defaultdict
from astropy_healpix import HEALPix
from functools import lru_cache
import json


# Import neighboring packages
# --------------------------------------------------------------



# Default for caching stuff using lru_cache
# -------------------------------------------------------------
cache_size_default = 16

# Data classes/methods
# -------------------------------------------------------------

class Base():
    '''
       Parent class to hold some file/directory definitions & methods
    '''
    def __init__(self):
        
        # Filing ...
        self.Names  = 'names.txt'
        self.Coeffs = 'coefficients.json'
        self.HP     = '_healpix.txt'

        # Healpix ...
        self.HP_nside=16
        self.HP_order='nested'
        self.hp = HEALPix(nside=self.HP_nside, order=self.HP_order)
        self.npix = self.hp.npix
    
        # Default list of Julian Dates to use (2440000 ==> 1968, 2460000.0 ==> 2023)
        # - Could/Should change this so that the end date is checked to be X-days in the future, or something similar
        self.JDlist = np.arange(2440000., 2460000.0, 1.0)

    def _fetch_data_directory(self):
        '''
            Returns the default path to the directory where data will be downloaded.
            
            By default, this method will return ~/.mpchecker2/data
            and create this directory if it does not exist.
            
            If the directory cannot be accessed or created, then it returns the local directory (".")
            
            Returns
            -------
            data_dir : str
            Path to location of `data_dir` where data (FITs files) will be downloaded
        '''
        
        data_dir = os.path.join(os.path.expanduser('~'), '.mpchecker_data')
        
        # If it doesn't exist, make a new data directory
        if not os.path.isdir(data_dir):
            
            try:
                os.mkdir(data_dir)
        
            # downloads locally if OS error occurs
            except OSError:
                warnings.warn('Warning: unable to create {}. '
                              'Download directory set to be the current working directory instead.'.format(data_dir))
                data_dir = '.'

        return data_dir

class Names(Base):
    '''
        List of names/designation (strings) 
        Allows implicit mapping of names <--> integers (posn in list)
        Allows easier lookups later-on
        
        Methods
        -------
        
        insert()
        get_number()
        get_name() 
        
        _read_file()
        _write_file()
        _get_filepath()
        
    '''
    def __init__(self):
        # Give access to "Base" methods
        super().__init__()
        
        # Read file by default
        self.names_dict = self._read_names_file()


    # --------------------------------------------------------
    # Read extant Name/Number data
    # --------------------------------------------------------
    def get_number(self, name):
        ''' 
            Map name to a number (the position in the file)
            If name NOT known, returns False
            
            input
            -------
            name: string
             - name/designation of an object

            return
            -------
            number: int
             - position of name in file / ordered-dict

        '''
        if name not in self.names_dict:
            return None
        else:
            return self.names_dict[name]

    def get_name(self, number):
        ''' 
            Map number to a name 
            If number not in file, returns False
        
            input
            -------
            number: int
            - position of name in file / ordered-dict
            
            return
            -------
            name: string
            - name/designation of an object

        '''
        if number < 0 or number > len(self.names_dict):
            return False
        else:
            return list(self.names_dict.keys())[number]

    def _get_names_filepath(self,):
        ''' Get the filepath to the Name file '''
        return os.path.join(self._fetch_data_directory() , self.Names )

    def _read_names_file(self,):
        ''' 
            Read the contents of the name file
             - File is just a list of names
            
            return
            -------
            Ordered Dict
             - key : str; name/designation
             - value : integer; enumerating position in file/dict
        '''
        if os.path.isfile(self._get_names_filepath()):
            with open(self._get_names_filepath(), 'r') as fh :
                return OrderedDict( ((_.strip() , i) for i,_ in enumerate(fh.readlines()) ) )
        else:
                return OrderedDict()

    # --------------------------------------------------------
    # Write new Name/Number data
    # --------------------------------------------------------
    def insert(self, name):
        ''' Insert name into dict of names
            and/or
            Get position of name in name-file
        '''
        
        # insert if required
        if self.get_number(name) == None:
            self.names_dict[name] = len(self.names_dict)
        
        # save the data (do we really want to do this every time, or do we want to push this to the end)
        self._write_names_file()
        
        return self.names_dict[name]

    def _write_names_file(self,):
        ''' 
            Write to the name file
        '''
        with open(self._get_names_filepath(), 'w') as fh :
            for key,value in self.names_dict.items() :
                # Ensure we have a line-break
                item = key + '\n' if key[-1] != '\n' else key
                fh.write(item)
    
class LRU(OrderedDict):
    '''
        Implementing variant of functools.lru_cache():
        Limit size, evicting the least recently looked-up key when full
        https://docs.python.org/3/library/collections.html#collections.OrderedDict
    '''
    def __init__(self, maxsize=cache_size_default, /, *args, **kwds):
        self.maxsize = maxsize
        super().__init__(*args, **kwds)
    
    def __getitem__(self, key):
        ''' Any time we get an item, we move it to the end/right/last/most-recent'''
        value = super().__getitem__(key)
        self.move_to_end(key)
        return value
    
    def __setitem__(self, key, value):
        ''' Any time we set an item, we check to see if we have too many items.
            If so then we take the oldest/least-recently-used and
            (i) save it to disk *** NOT YET IMPLEMENTED ***
            (ii) remove it from the dict
            
        '''
        super().__setitem__(key, value)
        if len(self) > self.maxsize:
            oldest = next(iter(self))
            del self[oldest]

class HPData(Base):
    '''
        Precalculated HP location of objects
        - Probably from the geocenter
        - Probably on a nightly basis
        
        
        Methods
        -------
        
        '''
    def __init__(self):
        # Give access to "Base" methods
        super().__init__()

        # Variables to facilitate reading from "block files"
        self.logBlockSize = 2
        self.block_dict   = LRU()
    
    
    # --------------------------------------------------------
    # Read extant HP data
    # --------------------------------------------------------
    def _query_HP(JD, HP):
        ''' 
            For a given JD, find the object-integers in a given HP
            
            inputs:
            -------
            JD: float or int
             - julian date
            HP: integer
             - healpix
            
            returns:
            --------
            list of integers
            - object integers contained within the HP
            - can be empty
        '''
        
        # Get filing "block" for supplied JD
        JD_int, JD_block, block_posn = _round_JD(JD)
        
        # Return data for specific HP on specific night
        return self._get_block_data(JD_block)[block_posn * self.npix + HP]
    
    def _get_block_data(JD_block):
        '''
            Get block data, either from file or from block_dict
            
            input:
            ------
            JD_block : integer
            - see *_round_JD()* for derivation of "JD_block" from "JD"
            
            return
            -------
            contents of file
        '''
        # If we do not already have the data in memory, read from file
        if JD_block not in block_dict:
            
            # N.B. block_dict == LRU / OrderedDict, so this puts it to the "end/right/last/most-recent"
            self.block_dict[JD_block] = self._read_HP_block_file(JD_block)
        
        # Return block entry
        return self.block_dict[JD_block]

    @lru_cache(maxsize=cache_size_default)
    def _read_HP_block_file(self, JD_block):
        '''
            Read the contents of the HEALPIX "block file"
            
            input:
            ------
            JD_block : integer
            - see *_round_JD()* for derivation of "JD_block" from "JD"
            
            return
            -------
            contents of file
        '''
        # Check existence of blockfile
        blockfile = self._get_HP_filepath(JD_block)
        assert os.path.isfile(blockfile, ' blockfile [%r] does not exist' % blockfile
        
        # Read from file
        with open(blockfile, 'r') as fh :
            return fh.readlines()
    

    # --------------------------------------------------------
    # Write new HP data
    # --------------------------------------------------------
    def update_object_HP(self, number, JDlist, old_HPlist, new_HPlist):
        '''
            Update the HP information for the location of object-number for a list of JDs
            This includes deleting the old information if necessary (see update_object_JD)
        '''
        
        # check that the supplied lists are of equal lengths
        assert len(JDlist)==len(old_HPlist)==len(new_HPlist), \
            'lengths not equal: update_object() : len(JDlist), len(old_HPlist), len(new_HPlist)  = [%r], [%r], [%r]' % \
            (len(JDlist),len(old_HPlist),len(new_HPlist))
        
        # update the data for each JD
        for JD, oldHP, newHP in zip( JDlist, old_HPlist, new_HPlist ):
            self._update_object_JD(number, JD, oldHP, newHP)


        # save to file
        # - unclear whether I want this here or to defer to a higher authority ...
        #
        # Hmmmmmmm : JDlist can have multiple JDs, and these could span multiple JD_blocks, hence involving multiple files
        # => might want to save for every night, or save whenever a different block is encountered ...
        #
        #self._write_HP_file()


    def _update_object_JD(self, number, JD, oldHP, newHP):
        '''
            Update the HP information for the location of object-number for a specific JD
            This includes deleting the old information if necessary
            
            At present this only deals with single-integer HP, *not* lists
            
        '''

        # if oldHP == newHP: don't need to do anything
        if oldHP == newHP:
            pass
        else:
        
        
            # if oldHP exists, need to remove
            if isinstance(oldHP , int):
                
                # read extant data (I.e. get object-numbers within given HP on given JD)
                numbers_list = _query_HP(JD, oldHP)
                
                # remove the number from the old position
                # - the assumption is that number is the the oldHP you said it is in: if not, emit error ...
                try :
                    numbers_list.remove(number)
                except:
                    sys.exit('Attempted to remove number [%r] from numbers_list, but failed. Checking number in numbers_list = %r' % \
                             (number, number in numbers_list) )
        
            # if oldHP is None, don't need to remove
            elif oldHP is None:
                pass
            # if oldHP is anything else, don't know how to process
            else:
                sys.exit(' *_update_object_JD()* does not know how to process oldHP of type: %r' % type(oldHP))


            # Now add the number to the object-numbers in the list for the newHP ...
            numbers_list = _query_HP(JD, newHP)
            numbers_list.append(number)

            # Use the updated numbers_list to update the JD_dict
            self.JD_dict[JD][newHP] = numbers_list


    # --------------------------------------------------------
    # Assorted HP convenience functions
    # --------------------------------------------------------
    @lru_cache(maxsize=cache_size_default)
    def _round_JD(self, JD):
        '''
        Convenience function to round supplied JD to
        (a) nearest integer
        (b) nearest "block" of days (for filing purposes)
        - default to blocks of 10^logBlockSize == 10^2
        (c) posn of integer_JD within block
        -  0 <= integer < 10^logBlockSize
        '''
        integer_JD = round(JD)
        block_JD   = int(round(round(JD, -self.logBlockSize )))
        block_posn = integer_JD - block_JD + round( (10**self.logBlockSize)/2 )
        return integer_JD , block_JD, block_posn

    @lru_cache(maxsize=cache_size_default)
    def _get_HP_filepath(self, JD_block):
        ''' Get the filepath to the healpix block-file '''
        return os.path.join(self._fetch_data_directory() , 'HP_' + str(JD_block) + '.dat' )


    def _update_block_list(self, block_list):
        '''
        ...
        '''
        pass

    def _update_HP_file(self, JD):
        '''
        ...
        '''
        # Get filing "block" for supplied JD
        integer_JD , block_JD, block_posn = _round_JD(JD)

        # Open the correct file to write to ...
        with open(self._get_HP_filepath(JD_block), 'w') as fh :
            
            # Write the block list to file
            return fh.readlines()

    def _create_empty_block_list(self, ):
        '''
        Returns a list of empty lists
        List is of length N_healpix * N_days_in_block
        Convenience function: Probably only used at initialization
        '''
        return [ [] for i in range( self.npix * 10**self.logBlockSize ) ]

    def _create_empty_block_file(self, block_JD):
        '''
        Creates an empty file (populated with an empty block list)
        Convenience function: Probably only used at initialization
        '''
        with open(self._get_HP_filepath(block_JD), 'w') as fh :
            for empty_list in self._create_empty_block_list():
                fh.write(','.join(empty_list)+'\n')

    def _create_all_required_empty_block_files(self, ):
        '''
        Creates all_required_empty_block_files (each populated with an empty block list)
        Convenience function: Probably only used at initialization
        '''
        # Create a list containing all of the required block_JD's
        block_JDs = list( set( [self._round_JD(JD)[1] for JD in self.JDlist] ) )

        # Create an empty file (populated with an empty block list) for each block_JD
        for block_JD in block_JDs:
            self._create_empty_block_file( block_JD )

class ObjectCoefficients(Names, HPData):
    '''
        Fitted coefficients for each object 
        
        Basic structure is a list of dictionaries
        
        WIP to define what these will be
         - At least positional information 
         - Probably velocity information
         - Some kind of Uncert/Covar information
         
        
        Methods
        -------
        
        get_coeff()
        
        _read_coeff_file()
        _write_coeff_file()
        
    '''
    def __init__(self):
        
        # Give access to "Names" methods & attributes
        Names.__init__(self)
        HPData.__init__(self)

        # List of coefficient-dictionaries: populated from file
        self.coeff_list = self._read_coeff_file()



    # --------------------------------------------------------
    # Read extant Coeff data
    # --------------------------------------------------------
    def get_coeff(self,number):
        '''
            Read co-efficients from list
        '''
        
        # Get coefficients from position in list
        if number >= len(self.coeff_list):
            return False
        else:
            return self.coeff_list[number]

    def _get_coeff_filepath(self,):
        ''' Returns the filepath to the coefficient file '''
        return os.path.join(self._fetch_data_directory() , Base().Coeffs )

    def _read_coeff_file(self,):
        '''
            Read the contents of the coefficients file
             - File is a list of dictionaries
            
            return
            -------
            Ordered Dict
            - key : str; name/designation
            - value : integer; enumerating position in file/dict
        '''
        if os.path.isfile(self._get_coeff_filepath()):
            with open(self._get_coeff_filepath(), 'r') as fh :
                return json.load(fh)
        else:
            return []
    
    # --------------------------------------------------------
    # Write new Coeff data
    # --------------------------------------------------------
    def upsert_coeff(self, number, new_coeff_dict ):
        '''
            Insert/Update coefficients for an object
            
            Inputs:
            -------
            number: integer 
             - position in the list into which we will insert/update
            coeff_dict: dictionary
            - all coefficient data to be saved for an object
            - must contain "name" as a key
            - must contain/enable-calculation-of geocentric-healpix data
            
            Returns:
            --------
            ????:
            -
            
        '''

        # update list of coefficients
        # - N.B. it is useful for later HP-updates to keep the old-coefficients
        
        # If there is an extant entry, update it ...
        if number < len(self.coeff_list):
            old_coeff_dict = self.coeff_list[number]
            self.coeff_list[number] = new_coeff_dict
        # Add an entry on the end
        #(but only if the position-number suggested is the same as the current length) ...
        elif number == len(self.coeff_list):
            old_coeff_dict = None
            self.coeff_list.append(new_coeff_dict)
        # Otherwise, crap out ...
        else:
            sys.exit('ObjectCoefficients.insert: number [%r] > len(self.coeff_list) [%r]' % \
                     (number, len(self.coeff_list))
                     )


        # Not sure whether I ultimately want this here: is doing file-updates every time ...
        self._write_coeff_file()

        return old_coeff_dict


    def _write_coeff_file(self,):
        '''
            Write list-of-dicts to file 
             - N.B. have put in hack to turn any np-arrays to lists before dumping to json
        '''
        with open(self._get_coeff_filepath(), 'w') as fh :
            # N.B. json hates arrays, so this converts them to lists before dumping to file
            json.dump( [self._arrays_to_lists(dict) for dict in self.coeff_list], fh)


    def _arrays_to_lists(self, dict):
        ''' 
            json hates arrays
            so this converts them to lists
        '''
        return { k:v.tolist() if isinstance(v, np.ndarray ) else v for k,v in dict.items() }


    # --------------------------------------------------------
    # Crucially important routines !!!
    # --------------------------------------------------------
    """
    def getObservationVector(JD, coeff_dict):
        '''
            #Given a set of objects and a time, this returns a set of vectors from Pan-STARRs to the objects as they appear in
            #the sky at that time (speed-of-light delay factored in).

            *** A METHOD NEEDS TO EXIST TO GET THE ObservationVector FROM THE COEFFICIENT-DICTIONARY ***

            This function calculates the ObservationVector of a SINGLE-OBJECT on a SINGLE NIGHT

        '''
        c = MPCL.Constants.speed_of_light
        baseTime = 2454466.5
        offset = 100 #100 days before for speed-of-light buffer
        t0Base = baseTime + offset
        t0=MPCL.EOP.jdTDB(t0Base)

        t1Base = t0Base + atTime
        t = MPCL.EOP.jdTDB(t1Base)
        obsPosition = np.array(MPCL.Observatories.getObservatoryPosition('F51', t))

        totalTime = 0
        #for objNum in range(len(chebyshevs)):
        lightDelay = 0.0
        for i in range(2):
            delayedTime = t - lightDelay
            centaurPosition = evaluate3DPiecewiseCheb(delayedTime-t0+offset, chebyshevs[objNum]) # Get from Margaret
            obsToCentaurVector = centaurPosition - obsPosition
            d = np.linalg.norm(obsToCentaurVector)
            lightDelay = d/c
        
        x, y, z = obsToCentaurVector
        return x,y,z
    """


    def _generate_HP_from_CoeffDict(self, JDlist, coeff_dict):
        '''
            This function calculates the HP location of a SINGLE-OBJECT over a series of MANY dates
            
            inputs:
            -------
            JDlist: list of integers
             - list of integer Julian Days for which calculations should be performed
             
            coeff_dict: dictionary 
             - generalized coefficient dictionary 
             - must contain sufficient information to allow HP to be extracted/calculated
             
            returns:
            --------
            HPlist: list
             - will be list of integers if coeff-dict properly specified
             - in the event of improper specificiation and/or coeff_dict == None, then a list of "None" will be returned
             - len(HPlist) == len(JDlist)
        '''
        # set default to be list of "None"
        HPlist = [None]*len(JDlist)
        
        # calculate the HP for each JD
        try:
            if coeff_dict != None:
                print("\n\t *** WARNING: RANDOM HEALPIX !!! \n ")
                for i, JD in enumerate(JDlist):
                    #x,y,z = getObservationVector(JD, coeff_dict)
                    #HPlist[i]  = hp.vec2pix(self.HP_nside, x, y, z, order=self.HP_order)
                    HPlist[i]=np.random.randint(self.npix)
        except:
            print("\n\t *** WARNING: FAILURE IN *_generate_HP_from_CoeffDict()* !!! \n ")
        return HPlist

# --------------------------------------------------------
# Primary External Class for Access to MPChecker2's
#  Pre-calculated Data
# --------------------------------------------------------

class PreCalc(ObjectCoefficients):
    '''
        Primary External Class for accessing MPChecker2's 
        pre-calculated data

        Methods
        -------
        
        upsert()
        query_HPlist()
        
    '''
    def __init__(self):
        
        # Give access to "ObjectCoefficients" methods & attributes
        # - Under-the-hood also allows access to "Names" and "HPData" and "Base" methods & attributes
        ObjectCoefficients.__init__(self)

    def upsert(self, new_coeff_dict_variable ):
        '''
            Main method used to insert/update coefficients for an object
            Also handles the healpix calculations used for efficiency-of-read
            
            Inputs:
            -------
            new_coeff_dict_variable: dictionary (or iterable of dictionaries)
            - all coefficient data to be saved for an object
            - must contain "name" as a key
            - must contain/enable-calculation-of geocentric-healpix data
            
            Returns:
            --------
            ????:
            -
            
        '''
        # ensure that the supplied variable is formatted correctly
        new_coeff_dict_list = self._rectify_inputs(new_coeff_dict_variable)
        
        # check contents of dictionary(ies)
        self._check_dictionary_contents(new_coeff_dict_list)
        
        # iterate over pairs in lists ...
        for new_coeff_dict in new_coeff_dict_list:
        
            # map name to number
            number = self.insert(new_coeff_dict["name"])

            # update list of coefficients
            # - N.B. it is useful for later HP-updates to keep the old-coefficients
            old_coeff_dict = self.upsert_coeff( number , new_coeff_dict )
            
            # Use the coefficient-dictionary(ies) to get the HP for each JD in JDlist
            old_HPlist = self._generate_HP_from_CoeffDict(self.JDlist, old_coeff_dict)
            new_HPlist = self._generate_HP_from_CoeffDict(self.JDlist, new_coeff_dict)

            # update HP data
            self.update_object_HP(number, self.JDlist, old_HPlist, new_HPlist)

    def get_nightly_precalcs(JD, HPlist):
        '''
            Main method used to query the precalculated healpix
            For a given JD, finds the object-integers for each HP in a list of HPs
            
            Returns the object-names and coefficient-dictionaries for
            the objects present in each healpix on the specified julian date
            
            
            Note that at present *no* attempt is made to ***ALWAYS RETURN NEOs***
            - This needs to be added-to / improved
            - We could (perhaps) return everything closer than 0.X au (on a given JD)

            inputs
            ------
            JD: float or int
             - julian date of the night. If <float> will be silently converted to <int>
             
            HPlist: list-of-integers 
             - healpix to be queried 
             - if integer (single healpix) supplied, is silently converted to list
             
            returns
            -------
            dictionary
            - key   = HP
            - value = list of coeff-dictionaries
        '''
        HPlist = [HPlist] if isinstance(HPlist, int) else HPlist
        return {HP: [ self.get_coeff(number) for number in self._query_HP(JD, HP)] for HP in HPlist}


    def _rectify_inputs(self,  new_coeff_dict_variable ):
        '''
            Private method called to rectify the input to upsert()
            
            inputs:
            -------
            new_coeff_dict: dictionary or list-of-dictionaries
            
            returns: 
            --------
            name_list: list-of-strings 
             - names of each object being "upcerted" 
             
            new_coeff_dict_list: list-of-dictionaries
             - 
            
        '''
        # If singular quantity, make into lists
        if  isinstance(new_coeff_dict_variable, dict):
            new_coeff_dict_list = [new_coeff_dict_variable]

        # If non-singular, check plausibly formatted
        # (N.B. Not doing detailed content checks at this point)
        elif    isinstance(new_coeff_dict_variable, (list, np.ndarray)) and \
                np.all([ isinstance(_, dict) for _ in new_coeff_dict_variable]):
            new_coeff_dict_list = new_coeff_dict_variable
        else:
            sys.exit('Cannot process inputs of type ... %r , %r' % (type(name), type(new_coeff_dict_variable)) )


        # return everything in list form
        return new_coeff_dict_list


    def _check_dictionary_contents( self , new_coeff_dict_list):
        '''
            Implement various sanity checks on the contents of the coefficient dictionaries
            We'll want to check for 
             - "name" variable
             -
        '''
        for n, d in enumerate(new_coeff_dict_list) :
            assert "name" in d, ' the variable *name* is missing from a/the [%d]^th supplied dictionary' % (n)







