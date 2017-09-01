"""
# Copyright 2017 Chris Culnane
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
from __future__ import print_function

import subprocess
import time
import mmap
import hashlib
from table import Table



class MapTable(Table):
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

    def __init__(self, table_file_path, line_length):
        self.table_path = table_file_path
        self.line_length = line_length
        self.lookupcount = 0
        self.totaltime = 0
        self.table_file = None
        self.open_for_write = False
        self.map = None
        self.lines = 0

    @staticmethod
    def calculate_line_length(maxval):
        """
        Calculates the necessary line length for storing values that are at most 2^maxval
        """
        return len(str((2**maxval)+1)) + 64 + 1

    def open(self, for_writing=False):
        """
        Opens the table file for reading or writing. Set for_writing to True to overwrite or
        build a new table
        """
        if for_writing:
            self.table_file = open(self.table_path, 'w')
            self.open_for_write = True
        else:
            # open for reading, therefore create memory map
            self.table_file = open(self.table_path, 'r')
            self.open_for_write = False
            self.map = mmap.mmap(self.table_file.fileno(), 0, prot=mmap.PROT_READ)
            self.lines = self.map.size()/self.line_length

    def get_lookup_time(self):
        """
        Debug method that keeps track of how long lookups take on average
        """
        return self.totaltime / self.lookupcount

    def lookup(self, elem):
        """
        Looks up elem in the table. First it converts elem to its string representation
        before SHA256'ing that value and performing a binary search of the memory mapped
        file. If it finds an entry it will return the string value (including padding)
        associated with that key. Otherwise it returns None
        """
        if self.open_for_write:
            raise Exception('Cannot lookup from a table that is being written to')
        start = time.time()

        # Create hash representation of element
        hashval = self._create_hash_string(elem)

        # Find middle point of file to start at
        index = self.lines//2

        lowerindex = 0
        upperindex = self.lines + 1

        # Loop until we either find the key or run out of search space
        while upperindex - lowerindex > 1:
            index = lowerindex + ((upperindex - lowerindex)//2)
            if index == 0:
                index = 1
            self.map.seek((index-1) * self.line_length)
            currentline = self.map.readline()
            hashlookup = currentline[:64]
            if hashlookup == hashval:
                end = time.time()
                self.totaltime = self.totaltime + (end-start)
                self.lookupcount = self.lookupcount + 1
                return currentline[64:]
            elif hashlookup < hashval:
                lowerindex = index
            else:
                upperindex = index
        return None


    def close(self):
        """
        Close the file and any open memory map
        """
        if not self.open_for_write:
            self.map.close()
        self.table_file.close()

    def add_row(self, element, value):
        """
        Adds a row to the table, converting the element to its hash representation
        and padding appropriately
        """
        if self.open_for_write != True:
            raise Exception('Table not openned for writing')
        hashval = self._create_hash_string(element)
        outstring = hashval + str(value)
        while len(outstring) < self.line_length-1:
            outstring = outstring + ' '
        outstring = outstring + '\n'
        self.table_file.write(outstring)

    def _create_hash_string(self, element):
        """
        Private utility method to contruct a padded hash value
        """
        hashval = hashlib.sha256(element.__str__()).hexdigest()
        while len(hashval) < 64:
            hashval = hashval + ' '
        return hashval


    def sort(self, outfile):
        """
        Sort the underlying DLTable file and save to the specified outfile. This must
        be called once a table has been constructed and before trying to lookup a value.
        """
        proc = subprocess.Popen(['sort', self.table_path], stdout=subprocess.PIPE)
        stdout, _ = proc.communicate()
        with open(outfile, 'w') as out:
            out.write(stdout)

