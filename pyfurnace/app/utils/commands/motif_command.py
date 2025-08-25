from abc import ABC, abstractmethod


class MotifCommand(ABC):

    @abstractmethod
    def execute(self):
        pass
