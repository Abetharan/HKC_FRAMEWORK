import abc
class CouplingMethod(metaclass = abc.ABCMeta):
    @abc.abstractmethod
    def method(self):
        """
        Coupling method 
        """
        raise NotImplementedError