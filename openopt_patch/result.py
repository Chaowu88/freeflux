from oologfcn import OpenOptException
from numpy import ndarray, matrix, asscalar, asarray, asfarray
from time import time, perf_counter, strftime

class OpenOptResult: 
    # TODO: implement it
    #extras = EmptyClass() # used for some optional output
    def __call__(self, *args):
        if not self.isFDmodel:
            raise OpenOptException('Is callable for FuncDesigner models only')
        r = []
        for arg in args:
            tmp = [(self._xf[elem] if isinstance(elem,  str) \
                    else self.xf[elem]) for elem in (arg.tolist() if isinstance(arg, ndarray) \
                    else arg if type(arg) in (tuple, list) else [arg])]
            tmp = [asscalar(item) if type(item) in (ndarray, matrix) and item.size == 1 \
                   #else item[0] if type(item) in (list, tuple) and len(item) == 0 \
                   else item for item in tmp]
            r.append(tmp if type(tmp) not in (list, tuple) or len(tmp)!=1 else tmp[0])
        r = r[0] if len(args) == 1 else r
        if len(args) == 1 and type(r) in (list, tuple) and len(r) >1: 
            r = asfarray(r, dtype = type(r[0]))
        return r
        
    def __init__(self, p):
        self.rf = asscalar(asarray(p.rf))
        self.ff = asscalar(asarray(p.ff))
        self.isFDmodel = p.isFDmodel
        self.probType = p.probType
        self.probName = p.name
        from openopt import oosolver
        self.solverName = oosolver(p.solver).__name__
        
        self._misc = {}
        Arr = ['_localtime_started', '_localtime_finished', 
                     'isObjFunValueASingleNumber', 'isFDmodel', 'err', '_bg_color']
        if p.isFDmodel:
            Arr += ['_init_ooarrays',
                     '_init_fixed_ooarrays', 'freeVarsSet', 'fixedVarsSet',
                     '_init_fixed_ooarrays']
        for attr in Arr:
                        self._misc[attr] = getattr(p, attr)
        
        
#        self.p = p
#        self.name = p.name
        if p.probType == 'EIG':
            self.eigenvalues, self.eigenvectors = p.eigenvalues, p.eigenvectors
        
        if p.isFDmodel:
            from FuncDesigner import oopoint
            self.xf = dict((v, asscalar(val) if isinstance(val, ndarray) and val.size ==1 \
                            else dict((field, v.domain[int(val)][j]) for j, field in enumerate(v.fields)) if v.fields != ()\
                            else v.reverse_aux_domain[int(val)] if 'reverse_aux_domain' in v.__dict__\
                            else v.aux_domain[val] if 'aux_domain' in v.__dict__ else val) \
                            for v, val in p.xf.items())
            if not hasattr(self, '_xf'):
                #self._xf = dict([(v.name, asscalar(val) if isinstance(val, ndarray) and val.size ==1 else val) for v, val in p.xf.items()])
                self._xf = dict((v.name, val) for v, val in self.xf.items())
            self.xf = oopoint(self.xf, maxDistributionSize = p.maxDistributionSize, skipArrayCast = True)
        else:
            self.xf = p.xf
        
        # TODO: mb use the fields in MOP
        if p.probType == 'MOP':
            for attr in ('xf', 'ff', 'rf', '_xf'):
                delattr(self, attr)

        
        #if len(p.solutions) == 0 and p.isFeas(p.xk): p.solutions = [p.xk]
        
        # TODO: mb perform check on each solution for more safety?
        # although it can slow down calculations for huge solutions number
        #self.solutions = p.solutions 

        self.elapsed = dict()
        self.elapsed['solver_time'] = round(100.0*(time() - p.timeStart))/100.0
        self.elapsed['solver_cputime'] = round(100.0*(perf_counter() - p.cpuTimeStart))/100.0
        self.elapsed['initialization_time'] = round(100.0*p.initTime)/100.0
        self.elapsed['initialization_cputime'] = round(100.0*p.initCPUTime)/100.0

        for fn in ('ff', 'istop', 'duals', 'isFeasible', 'msg', 'stopcase', 'iterValues',  'special', 'extras', 'solutions'):
            if hasattr(p, fn):  setattr(self, fn, getattr(p, fn))

        #TODO: mention in doc; mb somehow remove for large sequence of subprobs to clear memory
        if hasattr(p.solver, 'innerState'):
            self.extras['innerState'] = p.solver.innerState

        self.solverInfo = dict()
        for fn in ('homepage',  'alg',  'authors',  'license',  'info', 'name'):
            self.solverInfo[fn] =  getattr(p.solver,  '__' + fn + '__')

            # note - it doesn't work for len(args)>1 for current Python ver  2.6
            #self.__getitem__ = c # = self.__call__
            
        if p.plot:
            #for df in p.graphics.drawFuncs: df(p)    #TODO: include time spent here to (/cpu)timeElapsedForPlotting
            self.elapsed['plot_time'] = round(100*p.timeElapsedForPlotting[-1])/100 # seconds
            self.elapsed['plot_cputime'] = p.cpuTimeElapsedForPlotting[-1]
        else:
            self.elapsed['plot_time'] = 0
            self.elapsed['plot_cputime'] = 0

        self.elapsed['solver_time'] -= self.elapsed['plot_time']
        self.elapsed['solver_cputime'] -= self.elapsed['plot_cputime']

        self.evals = dict((key, val if type(val) == int else round(val *10) /10.0) for key, val in p.nEvals.items())
        self.evals['iter'] = p.iter
        
    _encode = lambda self: None # overloaded in some subclasses

    def _tkview(self, w=None):
        r = self
#        p = r.p
        try:
            from GUI import Text, Scrollbar, Tk, Button
        except:
            self._misc['err']('''
            Tkinter is not installed. 
            If you have Linux you could try using "apt-get install python-tk"
            ''')
        useRoot = w is None
        if useRoot:
            from openopt import __version__ as ooversion
            w = Tk()
            w.wm_title('OpenOpt  ' + ooversion)
            C = Button(w, text = 'Close', command = w.destroy)
            C.pack(side='bottom', fill='both', expand=True)
        if not getattr(w, 'made', False): 
            scrollbar = Scrollbar(w)
            R = Text(w,  yscrollcommand=scrollbar.set, bg=self._misc['_bg_color'])
            R.pack(side = 'left', expand=True, fill='both')
            scrollbar.pack(side='right', fill='y')
            scrollbar.config(command=R.yview)
            
            rr = _txtrepr(r)
            for elem in rr:
                R.insert('end', elem + '\n')
            
            w.made = True
        else:
            w.pack()#is it reachable?
        if useRoot:
            w.mainloop()
    
    def __getattr__(self, attr):
        if attr == 'view':
            print(''.join([elem + '\n' for elem in _txtrepr(self)]))
        elif attr == 'tk':
            self._tkview()
        else:
            raise AttributeError('''
            you are trying to obtain incorrect attribute "%s" 
            from OpenOpt result struct''' % attr
            )
        
def _txtrepr(r):
    rr = []
    if r is None:
        rr.append('result is None, no other info is available')
        return
        
    rr.append('Name: ' +  r.probName)
    rr.append('Type: ' +  r.probType)
    
    M = r._misc
    _localtime_started, _localtime_finished, \
     isObjFunValueASingleNumber, isFDmodel  =     M['_localtime_started'], M['_localtime_finished'], \
     M['isObjFunValueASingleNumber'], M['isFDmodel']
    
    rr.append('Solver: ' +  r.solverName)
#        R.highlight_pattern("word", foreground="red")
    rr.append(strftime("Started: \t%a, %d %b %Y %H:%M:%S", _localtime_started))
    rr.append(strftime("Finished:\t%a, %d %b %Y %H:%M:%S", _localtime_finished))
    
    if isObjFunValueASingleNumber:
        rr.append('Objective value: %0.9g' %  r.ff)
        if 'extremumBounds' in r.extras:
            t = r.extras['extremumBounds']
            rr.append('Estimated |f - f*|: %0.2g' % (t[1]-t[0]))
    
    for fn in ('isFeasible', 'stopcase','istop', 'msg'):
        rr.append(fn +': ' + str(getattr(r, fn)))
    
    rr.append('-'*45)
    # TODO: MCP etc
    if isFDmodel:
        _init_ooarrays,\
        _init_fixed_ooarrays, freeVarsSet, fixedVarsSet,\
        _init_fixed_ooarrays = M['_init_ooarrays'],\
     M['_init_fixed_ooarrays'], M['freeVarsSet'], M['fixedVarsSet'],\
     M['_init_fixed_ooarrays']
    if isFDmodel and r.probType not in ('KSP', 'BPP', 'TSP'):
        from FuncDesigner import _Stochastic

        '''                                           Free vars                                                         '''
        usedVars = set()
        rr2 = []
        free_ooarrays = _init_ooarrays.difference(_init_fixed_ooarrays)
        for v in free_ooarrays:
            if v.name.startswith('unnamed') and not v[0].name.startswith('unnamed'):
                for _v in v.view(ndarray):
                    if isinstance(_v(r), _Stochastic): continue
                    rr2.append(_v.name +': ' + str(_v(r)))
            else:
                if isinstance(v[0](r), _Stochastic): continue #TODO: check other array elements
                rr2.append(v.name +': ' + str(v(r).tolist()))
            usedVars.update(v.view(ndarray).tolist())
        for v in freeVarsSet.difference(usedVars):
            V = v(r)
            if isinstance(V, _Stochastic): continue
            rr2.append(v.name +': ' + str(V if type(V)!=ndarray else V.tolist()))
        rr2.sort(key = lambda elem:len(elem))
        rr += rr2
        
        '''                                           Fixed vars                                                         '''
        if len(fixedVarsSet) != 0:
            usedVars = set()
            rr2 = []
            fixed_ooarrays = _init_fixed_ooarrays
            for v in fixed_ooarrays:
                if v.name.startswith('unnamed') and not v[0].name.startswith('unnamed'):
                    for _v in v.view(ndarray):
                        if isinstance(_v(r), _Stochastic): continue
                        rr2.append(_v.name +': ' + str(_v(r)))
                else:
                    if isinstance(v[0](r), _Stochastic): continue#TODO: check other array elements
                    rr2.append(v.name +': ' + str(v(r).tolist()))
                usedVars.update(v.view(ndarray).tolist())
            for v in fixedVarsSet.difference(usedVars):
                V = v(r)
                if isinstance(V, _Stochastic): continue
                rr2.append(v.name +': ' + str(V if type(V)!=ndarray else V.tolist()))
            
            if len(rr2):
                rr2.sort(key = lambda elem:len(elem))
                rr.append('\n'+'-'*13+'  fixed variables  ' + '-'*13)
                rr += rr2

    elif r.probType == 'TSP' and hasattr(r, 'Edges'): # MOP
        for elem in [r.nodes, r.edges, r.Edges]:
            rr.append('\n' + str(elem))
            if elem is not r.Edges:
                rr.append('\n'+'-'*45)
    else:
        rr.append(str(r.xf if r.probType not in ('MCP', 'STAB', 'DSP') else r.solution))
    return rr
