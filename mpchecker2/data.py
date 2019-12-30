# -*- coding: utf-8 -*-
# mpchecker2/mpchecker2/data.py

'''
    --------------------------------------------------------------
    MPChecker2's data module.
    
    Dec 2019
    Matt Payne
    
    I think that the following data containers are necessary ...
    (i) Mapping of unique names/designations (string) to integers
    (ii) Cheby Polys to represent orbital positions, uncertainties, etc. 
         Will need to be able to be able to be "segmented"
    (iii) "Nightly" HealPix Dictionaries to store approximate sky-plane
        location from the geo-center
        
    Note that at present *no* attempt is made to ***ALWAYS RETURN NEOs***
     - This needs to be added-to / improved
     - We could (perhaps) return everything closer than 0.X au (on a given JD)
        
    --------------------------------------------------------------
    '''


# Import third-party packages
# --------------------------------------------------------------
import sys, os
import numpy as np
import operator
from collections import OrderedDict, defaultdict
from astropy_healpix import HEALPix as hp
from functools import lru_cache
#import json

# Import neighboring packages
# --------------------------------------------------------------
#   import mpchecker2.precalc


# --------------------------------------------------------
# Primary External Class for MPChecker2-Data Access
# --------------------------------------------------------

class Data(ObjectCoefficients):
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

    def upsert(self,name,new_coeff_dict):
        '''
            Main method used to insert/update coefficients for an object
            Also handles the healpix calculations used for efficiency-of-read
            
            Inputs:
            -------
            name: string (or iterable of strings)
            - this is the provisional-designation of the object
            
            coeff_dict: dictionary (or iterable of dictionaries)
            - all coefficient data to be saved for an object
            - must contain/enable-calculation-of geocentric-healpix data
            
            Returns:
            --------
            ????:
            -
            
        '''
        # ensure that the supplied variables are formatted as required
        name_list , new_coeff_dict_list = _rectify_inputs(name,new_coeff_dict)
        
        # iterate over pairs in lists ...
        for name, new_coeff_dict in zip(name_list , new_coeff_dict_list):
        
            # map name to number
            number = self.insert(name)

            # update list of coefficients
            # - N.B. it is useful for later HP-updates to keep the old-coefficients
            if number < len(self.coeff_list):    # update extant entry
                old_coeff_dict = self.coeff_list[number]
                self.coeff_list[number] = new_coeff_dict
            elif number == len(self.coeff_list): # append new entry to end
                old_coeff_dict = None
                self.coeff_list.append(new_coeff_dict)
            else:
                sys.exit('ObjectCoefficients.insert: number [%r] > len(self.coeff_list) [%r]' % \
                         (number, len(self.coeff_list))
                         )

            # Use the coefficient-dictionary(ies) to get the HP for each JD in JDlist
            old_HPlist = self._generate_HP_from_CoeffDict(JDlist, old_coeff_dict)
            new_HPlist = self._generate_HP_from_CoeffDict(JDlist, new_coeff_dict)

            # update HP data
            self.update_object_HP(number, JDlist, old_HPlist, new_HPlist)

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
            - value = list of object-integers
        '''
        HPlist = [HPlist] if isinstance(HPlist, int) else HPlist
        return {HP:self._query_HP(JD, HP) for HP in HPlist}


    def _rectify_inputs(name,new_coeff_dict):
        '''
            Private method called to rectify the inputs to upsert()
            
            inputs:
            -------
            name: string or list-of-strings
            
            new_coeff_dict: dictionary or list-of-dictionaries
            
            returns: 
            --------
            name_list: list-of-strings 
            new_coeff_dict_list: list-of-dictionaries
            
        '''
        # If singular quantities, make into lists
        if isinstance(name, str) and isinstance(new_coeff_dict, dict):
            name_list = [name]
            new_coeff_dict_list = [new_coeff_dict]

        # If non-singular, check that they are plausibly formatted
        # (N.B. Not doing detailed content checks at this point)
        elif:
            isinstance(name, (list, np.ndarray))           and np.all([ isinstance(_, str)  for _ in name]) and \
            isinstance(new_coeff_dict, (list, np.ndarray)) and np.all([ isinstance(_, dict) for _ in new_coeff_dict]):
                name_list = name
                new_coeff_dict_list = new_coeff_dict
        else:
            sys.exit('Cannot process inputs of type ... %r , %r' % (type(name), type(new_coeff_dict)) )

        # return everything in list form
        return name_list, new_coeff_dict_list






# Data classes/methods used under-the-hood
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
        self.npix = hp(nside=self.HP_nside, order=self.HP_order)
    
        # Caching
        self.cache_size_default = 16
    
        # Default list of Julian Dates to use (2440000 ==> 1968, 2460000.0 ==> 2023)
        # - Should probably change this so that the end date is checked to be X-days in the future, or something similar
        self.JDlist = np.arange(2440000., 2460000.0, 1.0)

    def _fetch_data_directory(self):
        '''
            Returns the default path to the directory where data will be downloaded.
            
            By default, this method will return ~/.shifty/data
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
        if number < 0 or number > len(self.names):
            return False
        else:
            return self.names[name]

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
        if self.get_number(name) == None:
            self.names_dict[name] = len(self.names_dict)
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
        
        upsert()
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
        # might be nice to try to convert any str to a number using Names ....
        
        # Get coefficients from position in list
        if number >= len(self.coeff_list):
            return False
        else:
            return coeff_list[number]

    def _get_coeff_filepath(self,):
        ''' Get the filepath to the Name file '''
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
    # Crucially important routine !!!
    # --------------------------------------------------------
    def _generate_HP_from_CoeffDict(self, JDlist, coeff_dict):
        '''
            *** A METHOD NEEDS TO EXIST TO GET THE HEALPIX FROM THE COEFFICIENT-DICTIONARY ***
            *** THIS COULD BE EITHER A SIMPLE GET/READ, OR IT COULD BE A CALCULATION ***
            
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
                for i, JD in enumerate(JDlist):
                    ### *** INSERT METHOD HERE *** ###
                    print("The current HP calculation in _generate_HP_from_CoeffDict() is WRONG/TEMPORARY/RANDOM")
                    HPlist[i] = np.random.randint(self.npix)
        except:
            pass

        return HPlist


class HPData(Base):
    '''
        Precalculated HP location of objects
        - Probably from the geocenter
        - Probably on a nightly basis
        
        
        Methods
        -------
        
        '''
    def __init__(self):
        self.logBlockSize = 2
        self.JD_dict      = defaultdict()
    
    
    # --------------------------------------------------------
    # Read extant HP data
    # --------------------------------------------------------
    @lru_cache(maxsize=Base.cache_size_default)
    def _query_HP(JD, HP):
        ''' For a given JD, find the object-integers in a given HP '''
        return _query_JD(JD)[HP]
   
    @lru_cache(maxsize=Base.cache_size_default)
    def _query_JD(JD):
        ''' For a given JD, return all HP data from file '''
    
        # Get filing "block" for supplied JD
        JD_int, JD_block, block_posn = _round_JD(JD)
        
        # Read data from block-file
        data_list = _read_HP_file(JD)
        
        # Return data for specific night
        self.defaultdict[JD] = data_list[block_posn * Base.npix : (block_posn + 1) * Base.npix ]
        return self.defaultdict[JD]
    
    @lru_cache(maxsize=Base.cache_size_default)
    def _read_HP_file(JD):
        '''
            (i) Decide which HEALPIX file needs to be opened
            - based on value of JD
            (ii) Read the contents of the HEALPIX file
            - file is a ...
            
            return
            -------
        '''
        # Get filing "block" for supplied JD
        JD_int, JD_block, block_posn = _round_JD(JD)
        
        # Read from file
        with open(self._get_HP_filepath(JD_block), 'r') as fh :
            return fh.readlines()
    

    # --------------------------------------------------------
    # Write new HP data
    # --------------------------------------------------------
    def update_object_HP(number, JDlist, old_HPlist, new_HPlist):
        '''
            ...
        '''
        
        # check that the supplied lists are of equal lengths
        assert len(JPlist)==len(HPlist), \
            'lengths not equal: update_object() : len(JDlist), len(HPlist) = [%r], [%r]' % \
            (len(JPlist),len(HPlist))
        
        # update the data for each JD
        for JD, HP in zip( JDlist, HPlist ):
            update_object_JD(JD, HP)

    def update_object_JD(number, JD, HP):
        '''
            ...
        '''
            
        if isinstance(HP , int):

            # read extant data
            numbers_list = _query_HP(JD, HP)

            # update list of object-numbers
            if HP not in numbers_list:
                numbers_list.append(HP)


        else:
            sys.exit('No method in place to process HP of type [%r]' % type(HP))


    # --------------------------------------------------------
    # Assorted HP convenience functions
    # --------------------------------------------------------
    @lru_cache(maxsize=Base.cache_size_default)
    def _round_JD(JD):
        '''
        Convenience function to round supplied JD to
        (a) nearest integer
        (b) nearest "block" of days (for filing purposes)
        - default to blocks of 10^logBlockSize == 10^2
        (c) posn of integer_JD within block
        -  0 <= integer < 10^logBlockSize
        '''
        integer_JD = round(JD)
        block_JD   = round(round(JD, -self.logBlockSize ))
        block_posn = integer_JD - block_JD + round( (10**self.logBlockSize)/2 )
        return integer_JD , block_JD, block_posn

    @lru_cache(maxsize=Base.cache_size_default)
    def _get_HP_filepath(JD_block):
        ''' Get the filepath to the healpix block-file '''
        return os.path.join(self._fetch_data_directory() , str(JD_block) , Base.HP )


    def _create_empty_block_list():
        '''
        Returns an list of empty lists
        List is of length N_healpix * N_days_in_block
        '''
        return [ [] for i in range( Base.npix * 10**self.logBlockSize ) ]

    def _update_block_list(block_list):
        '''
        ...
        '''
        pass

    def _update_HP_file(JD):
        '''
        ...
        '''
        # Get filing "block" for supplied JD
        JD_int, JD_block = _round_JD(JD)



