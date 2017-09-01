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

import time
from table import Table

class MemTable(Table):
    """
    In memory dictionary that acts as a lookup table
    """

    def __init__(self):
        self.table = {}
        self.open_for_write = False
        self.lookupcount = 0
        self.totaltime = 0

    def open(self, for_writing=False):
        """
        Opens the table file for reading or writing. Set for_writing to True to overwrite or
        build a new table
        """
        self.open_for_write = for_writing

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
        start = time.time()
        ret = None
        if elem in self.table:
            ret = self.table[elem]
        end = time.time()
        self.totaltime = self.totaltime + (end-start)
        self.lookupcount = self.lookupcount + 1
        return ret

    def close(self):
        """
        Close the file and any open memory map
        """

    def is_empty(self):
        if not self.table:
            return True
        return False

    def add_row(self, element, value):
        """
        Adds a row to the table, converting the element to its hash representation
        and padding appropriately
        """
        self.table[element] = value

    
    