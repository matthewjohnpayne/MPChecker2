# -*- coding: utf-8 -*-
# mpchecker2/mpchecker2/sql

'''
    --------------------------------------------------------------
    mpchecker2's sqlite module.
    
    Feb 2020
    Matt Payne
    
    This module provides functionalities to
    ...
    
    *THIS IS STILL PRIMARILY A STRAIGHT COPY OF THE STUFF FROM 'SIFTER' IT HAS NOT YET BEEN FULLY CHANGED TO MPCHECKER2 REQUIREMENTS (ONLY PARTLY DONE)*
    
    --------------------------------------------------------------
    '''


# Import third-party packages
# --------------------------------------------------------------
import sys, os
import numpy as np
import sqlite3
from sqlite3 import Error
import pickle

# Import neighboring packages
# --------------------------------------------------------------
print(os.path.dirname(os.path.realpath(__file__)))
sys.path.append( os.path.dirname(os.path.realpath(__file__)) )
import precalc




# Data classes/methods
#
# N.B. Many of these sqlite functions are copied straight from
# https://www.sqlitetutorial.net/sqlite-python/
# E.g.
# https://www.sqlitetutorial.net/sqlite-python/creating-database/
# https://www.sqlitetutorial.net/sqlite-python/create-tables/
# ...
#
# -------------------------------------------------------------


def fetch_db_filepath():
    '''
    '''
    B = precalc.Base()
    db_dir = B._fetch_data_directory()
    return os.path.join(db_dir , B.db_filename)

def create_connection(db_file):
    """ Create a database connection to the SQLite database
        specified by db_file
        
        inputs:
        -------
        db_file: database file
        
        return: 
        -------
        Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)
    
    return conn


def create_table(conn, create_table_sql):
    """ Create a table from the create_table_sql statement
        
        inputs:
        -------
        conn: Connection object
        
        create_table_sql: a CREATE TABLE statement
        
        return:
        -------
    """
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Error as e:
        print(e)



def create_specific_tables(conn):
    """ Create the specific table(s) that we need for *mpchecker2*
        Currently creates:
        (i) objects_by_jdhp
        (ii) object_coefficients

        inputs:
        -------
        conn: Connection object

        return:
        -------

    """
    
    
    # Create tables ...
    # (i) simple table w/ 3 user-supplied fields
    sql_create_objects_by_jdhp_table = """ CREATE TABLE IF NOT EXISTS objects_by_jdhp (
        id integer PRIMARY KEY,
        jd integer NOT NULL,
        hp integer NOT NULL,
        object_integer_list blob
        ); """
    
    # (ii) larger more complex table that needs many fields, one per coefficient-set
    #  - Dynamically generate the field-specs that will be required for the coeffs-by-sector ...
    sector_names = "".join( [ 'sector_%d_%d blob, ' % (i, jd) for i, jd in precalc.Base().required_sector_dict.items() ] )
    sector_names = sector_names[:-2] + " " # just getting rid of the final ","
    
    sql_create_object_coefficients_table = """ CREATE TABLE IF NOT EXISTS object_coefficients (
        id integer PRIMARY KEY,
        object_designation TEXT UNIQUE, """ + \
        sector_names + "); "

    # create table(s)
    if conn is not None:
        
        # create tables
        create_table(conn, sql_create_objects_by_jdhp_table)
        create_table(conn, sql_create_object_coefficients_table)

        # Create compound/combined indicees
        createSecondaryIndex =  "CREATE INDEX index_jdhp ON objects_by_jdhp (jd, hp);"
        conn.cursor().execute(createSecondaryIndex)



def upsert_multi_sector_cheby_dict(conn, multi_sector_cheby_dict):
    """
        insert/update multi_sector_cheby_dict for a single object
        
        N.B ...
        https://stackoverflow.com/questions/198692/can-i-pickle-a-python-dictionary-into-a-sqlite3-text-field
        pdata = cPickle.dumps(data, cPickle.HIGHEST_PROTOCOL)
        curr.execute("insert into table (data) values (:data)", sqlite3.Binary(pdata))


            To insert multiple rows into a table, you use the following form of the INSERT statement:

            INSERT INTO table1 (column1,column2 ,..)
            VALUES
            (value1,value2 ,...),
            (value1,value2 ,...),
            ...
            (value1,value2 ,...);


        inputs:
        -------
        conn: Connection object
        
        multi_sector_cheby_dict : dictionary
         - see orbit_cheby module for detailed specification
        
        return:
        -------
        


    """
    # Extract necessary info from dict
    unpacked_designation = multi_sector_cheby_dict['unpacked_designation']
    print( multi_sector_cheby_dict.keys() )

    # Explit loop over sectors to generate insert statement
    for sector_dict in multi_sector_cheby_dict['sectors']:
        

    #
    pdata = pickle.dumps(tracklet_dict, pickle.HIGHEST_PROTOCOL)

    sql =  """ INSERT OR REPLACE INTO tracklets(jd,hp,tracklet_name,tracklet)
        VALUES(?,?,?,?)
        """
    cur = conn.cursor()
    cur.execute(sql, (jd, hp, tracklet_name, sqlite3.Binary(pdata),))
    conn.commit()


'''


def delete_tracklet(conn, tracklet_name):
    """
        delete tracklet data
        
        inputs:
        -------
        tracklet_name: string
        
        return:
        -------
        
        
    """
    sql = 'DELETE FROM tracklets WHERE tracklet_name=?'
    cur = conn.cursor()
    cur.execute(sql, (tracklet_name,))
    conn.commit()

"""

def delete_tracklets(conn, tracklet_name_list):
    """
        delete list of tracklet data
        
        inputs:
        -------
        tracklet_name: string
        
        return:
        -------
        
        
        """
    sql = 'DELETE FROM tracklets WHERE tracklet_name IN (?)'
    cur = conn.cursor()
    cur.execute(sql, (tracklet_name_list,))
    conn.commit()
"""





def query_tracklets_jdhp(conn, JD, HP):
    """
       Define standard query used to find all tracklets for which jdhp == input
       
       inputs:
       -------
       JD: integer
       HP: integer
       
       return:
       -------
       list of tracklet_names

    """
    cur = conn.cursor()
    cur.execute("SELECT tracklet_name, tracklet FROM tracklets WHERE jd=? AND hp=?", ( int(JD), int(HP) , ))
    
    # return a dictionary-of-dictionaries
    return { row[0]: pickle.loads( row[1] ) for row in cur.fetchall() }

'''
