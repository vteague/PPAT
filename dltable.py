from __future__ import print_function

import subprocess
import time
import mmap
import hashlib
import mathTools.otosEC as oEC
import gmpy2 as gmpy
from gmpy2 import mpz




class DLTable:
    """
    Discrete Log lookup table that can handle very large files. Provides methods to add
    to the table during writing. Includes sort method to call Unix/Linux sort on the
    generated file. Uses a memory mapped file for lookup to allow very large files
    to be used. Performs a binary search for indexes which are the SHA256 hash of the
    string representation of the input element.

    Each line is a fixed size (64 chars for Hex SHA256) + value string (fixed) + newline

    A table should be created (written to via add_row) before closed and sorted. Once sorted
    the table can be used to lookup values via a binary search.
    """

    @classmethod
    def __init__(cls, cryptoGroup, table_file_path, line_length):
        cls.group = cryptoGroup
        cls.table_path = table_file_path
        cls.line_length = line_length
        cls.lookupcount = 0
        cls.totaltime = 0

    @staticmethod
    def calculate_line_length(maxval):
        """
        Calculates the necessary line length for storing values that are at most 2^maxval
        """
        return len(str((2**maxval)+1)) + 64 + 1

    @classmethod
    def open(cls, for_writing=False):
        """
        Opens the table file for reading or writing. Set for_writing to True to overwrite or
        build a new table
        """
        if for_writing:
            cls.table_file = open(cls.table_path, 'w')
            cls.open_for_write = True
        else:
            # open for reading, therefore create memory map
            cls.table_file = open(cls.table_path, 'r')
            cls.open_for_write = False
            cls.map = mmap.mmap(cls.table_file.fileno(), 0, prot=mmap.PROT_READ)
            cls.lines = cls.map.size()/cls.line_length

    @classmethod
    def get_lookup_time(cls):
        """
        Debug method that keeps track of how long lookups take on average
        """
        return cls.totaltime / cls.lookupcount

    @classmethod
    def lookup(cls, elem):
        """
        Looks up elem in the table. First it converts elem to its string representation
        before SHA256'ing that value and performing a binary search of the memory mapped
        file. If it finds an entry it will return the string value (including padding)
        associated with that key. Otherwise it returns None
        """
        if cls.open_for_write:
            raise Exception('Cannot lookup from a table that is being written to')
        start = time.time()

        # Create hash representation of element
        hashval = cls._create_hash_string(elem)

        # Find middle point of file to start at
        index = cls.lines//2

        lowerindex = 0
        upperindex = cls.lines + 1

        # Loop until we either find the key or run out of search space
        while upperindex - lowerindex > 1:
            index = lowerindex + ((upperindex - lowerindex)//2)
            if index == 0:
                index = 1
            cls.map.seek((index-1) * cls.line_length)
            currentline = cls.map.readline()
            hashlookup = currentline[:64]
            if hashlookup == hashval:
                end = time.time()
                cls.totaltime = cls.totaltime + (end-start)
                cls.lookupcount = cls.lookupcount + 1
                return currentline[64:]
            elif hashlookup < hashval:
                lowerindex = index
            else:
                upperindex = index
        return None


    @classmethod
    def close(cls):
        """
        Close the file and any open memory map
        """
        if not cls.open_for_write:
            cls.map.close()
        cls.table_file.close()

    @classmethod
    def add_row(cls, element, value):
        """
        Adds a row to the table, converting the element to its hash representation
        and padding appropriately
        """
        if cls.open_for_write != True:
            raise Exception('Table not openned for writing')
        hashval = cls._create_hash_string(element)
        outstring = hashval + str(value)
        while len(outstring) < cls.line_length-1:
            outstring = outstring + ' '
        outstring = outstring + '\n'
        cls.table_file.write(outstring)

    @classmethod
    def _create_hash_string(cls, element):
        """
        Private utility method to contruct a padded hash value
        """
        hashval = hashlib.sha256(element.__str__()).hexdigest()
        while len(hashval) < 64:
            hashval = hashval + ' '
        return hashval

    @classmethod
    def sort(cls, outfile):
        """
        Sort the underlying DLTable file and save to the specified outfile. This must
        be called once a table has been constructed and before trying to lookup a value.
        """
        proc = subprocess.Popen(['sort', cls.table_path], stdout=subprocess.PIPE)
        stdout, _ = proc.communicate()
        with open(outfile, 'w') as out:
            out.write(stdout)

