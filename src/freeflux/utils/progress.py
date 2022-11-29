#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '04/15/2022'




from time import time, sleep
from datetime import timedelta
from threading import Thread




class Progress():
    '''
    Progress instances can be either used in a WITH statement or assigned to a variable. 
    In the second case, the start and stop method should be called explicitely.
    '''
    
    def __init__(self, descp, silent = False):
        '''
        Parameters
        ----------
        descp: str
            Description.
        silent: bool
            Whether to show progress bar.
        '''
        
        self.descp = descp
        self.count = 0
        self.is_running = False
        self.silent = silent
        
    
    def __enter__(self):
        
        self.start()
        
        
    def __exit__(self, type, value, traceback):
        
        self.stop()
    
        
    def _sleeptime(self):
        
        if self.descp == 'optimization':
            sleeptime = 0.1    
        elif self.descp == 'fitting with CIs':
            sleeptime = 10
        else:
            sleeptime = 2
            
        return sleeptime    
    
        
    def _show_progress(self, descp):
        
        t1 = time()
        while self.is_running:
            sleep(1)
            t2 = time()
            print('%s [elapsed: %s]' % (descp, timedelta(seconds = round(t2 - t1))), end = '\r')
        else:
            print('\n')
        
        
    def start(self):
        
        if not self.silent:
            print()
            self.is_running = True
            self.thread = Thread(target = self._show_progress, args = (self.descp,))
            self.thread.start()
        else:
            pass

    
    def stop(self):
        
        if not self.silent:
            self.is_running = False
            self.thread.join()
        else:
            pass
        
        