import abc

class DLTable(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, table, group):
        self.table = table
        self.group = group

    @abc.abstractmethod
    def build(self, base, max_dl=2 ** 32, max_search=2 ** 12):
        raise NotImplementedError('users must define build to use this base class')

    @abc.abstractmethod
    def extract(self, a, b):
        raise NotImplementedError('users must define extract to use this base class')

