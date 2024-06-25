from abc import ABC



class Module(ABC):
    def __init__(self, method : str, **kwargs):
        
        pass
    
    @abstractmethod
    def call(self, adata : sc.AnnData) -> sc.AnnData:
        
        pass