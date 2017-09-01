import abc

class Table(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def open(self, for_writing=False):
        raise NotImplementedError('users must define __open__ to use this base class')

    @abc.abstractmethod
    def close(self, group, base, max_dl=2 ** 32, max_search=2 ** 12):
        raise NotImplementedError('users must define __build__ to use this base class')

    @abc.abstractmethod
    def add_row(self, element, value):
        raise NotImplementedError('users must define __add_row__ to use this base class')

    @abc.abstractmethod
    def lookup(self, elem):
        raise NotImplementedError('users must define __lookup__ to use this base class')

    @abc.abstractmethod
    def is_empty(self):
        raise NotImplementedError('users must define __is_empty__ to use this base class')

        