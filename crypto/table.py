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
import abc

class Table(object):
    """
    Abstract table class representing both memory and file based tables
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def open(self, for_writing=False):
        """
        Open the table for either reading or writing
        """
        raise NotImplementedError('users must define __open__ to use this base class')

    @abc.abstractmethod
    def close(self, group, base, max_dl=2 ** 32, max_search=2 ** 12):
        """
        Closes the table
        """
        raise NotImplementedError('users must define __build__ to use this base class')

    @abc.abstractmethod
    def add_row(self, element, value):
        """
        Add a row to the table with the specified index of element and value
        """
        raise NotImplementedError('users must define __add_row__ to use this base class')

    @abc.abstractmethod
    def lookup(self, element):
        """
        Lookup an element in the table returning value or None
        """
        raise NotImplementedError('users must define __lookup__ to use this base class')

    @abc.abstractmethod
    def is_empty(self):
        """
        Checks if the table is empty
        """
        raise NotImplementedError('users must define __is_empty__ to use this base class')

        