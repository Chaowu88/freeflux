'''Define the Context class.'''


__author__ = 'Chao Wu'


class Context():
    
    def __init__(self):
        
        self.operations = []
        
    
    def add_undo(self, op):
        '''
        Parameters
        ----------
        op: callable
            Operation to reset the model.
        '''
        
        self.operations.append(op)


    def undo(self):
        
        while self.operations:
            op = self.operations.pop()
            op()
