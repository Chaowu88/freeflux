__docformat__ = "restructuredtext en"
from time import time, perf_counter, localtime
from numpy import asfarray, nan, ones, all, atleast_1d, any, isnan, \
array_equal, asscalar, asarray, ndarray, isscalar, seterr, isinf,inf#, where
import numpy as np
from setDefaultIterFuncs import stopcase,  SMALL_DELTA_X,  SMALL_DELTA_F, IS_MAX_ITER_REACHED
from check import check
from oologfcn import OpenOptException, errSet
from openopt import __version__ as version
from result import OpenOptResult


######################
# don't change to mere ooMisc! 
from openopt.kernel.ooMisc import isSolved 
######################

#from baseProblem import ProbDefaults
from nonOptMisc import getSolverFromStringName, EmptyClass
# for PyPy
from openopt.kernel.nonOptMisc import where

try:
    import setproctitle
    hasSetproctitleModule = True
except ImportError:
    hasSetproctitleModule = False
    
try:
    from FuncDesigner import FuncDesignerException
except:
    FuncDesignerException = OpenOptException
    
#from openopt.kernel.ooMisc import __solverPaths__
ConTolMultiplier = 0.8

#if __solverPaths__ is None:
#    __solverPaths__ = {}
#    file = string.join(__file__.split(os.sep)[:-1], os.sep)
#    for root, dirs, files in os.walk(os.path.dirname(file)+os.sep+'solvers'):
#        rd = root.split(os.sep)
#        if '.svn' in rd: continue
#        rd = rd[rd.index('solvers')+1:]
#        for file in files:
#            print file
#            if len(file)>6 and file[-6:] == '_oo.py':
#                __solverPaths__[file[:-6]] = 'openopt.solvers.' + string.join(rd,'.') + '.'+file[:-3]

#import pickle
#f = open('solverPaths.py', 'w')
#solverPaths = pickle.load(f)


def runProbSolver(p_, solver_str_or_instance=None, *args, **kwargs):
    #p = copy.deepcopy(p_, memo=None, _nil=[])
    errSet.clear()
    p = p_
    p._localtime_started = localtime()
    if len(args) != 0: p.err('unexpected args for p.solve()')
    if hasattr(p, 'was_involved'): p.err("""You can't run same prob instance for twice. 
    Please reassign prob struct. 
    You can avoid it via using FuncDesigner oosystem.""")
    else: p.was_involved = True

    if solver_str_or_instance is None:
        if hasattr(p, 'solver'): solver_str_or_instance = p.solver
        elif 'solver' in kwargs: solver_str_or_instance = kwargs['solver']

    if type(solver_str_or_instance) is str and ':' in solver_str_or_instance:
        isConverter = True
        probTypeToConvert,  solverName = solver_str_or_instance.split(':', 1)
        p.solver = getSolverFromStringName(p, solverName)
        solver_params = {}
        #return converter(solverName, *args, **kwargs)
    else:
        isConverter = False
        if solver_str_or_instance is None:
            p.err('you should provide name of solver')
        elif type(solver_str_or_instance) is str:
            p.solver = getSolverFromStringName(p, solver_str_or_instance)
        else: # solver_str_or_instance is oosolver
            if not solver_str_or_instance.isInstalled:
                p.err('''
                solver %s seems to be uninstalled yet, 
                check http://openopt.org/%s for install instructions''' % (solver_str_or_instance.__name__, p.probType))
            p.solver = solver_str_or_instance
            for key, value  in solver_str_or_instance.fieldsForProbInstance.items():
                setattr(p, key, value)
    p.isConverterInvolved = isConverter
    if isConverter and p.probType not in ['DFP']: #, 'NLLSP']:
        p.err('converters are temporarily disabled, use native solvers')

    old_err = seterr(all= 'ignore')
    
    if 'debug' in kwargs.keys():
       p.debug =  kwargs['debug']

    probAttributes = set(p.__dict__)
    solverAttributes = set(p.solver.__dict__)
    intersection = list(probAttributes.intersection(solverAttributes))
    if len(intersection) != 0:
        if p.debug:
            p.warn('''
            attribute %s is present in both solver and prob 
            (probably you assigned solver parameter in prob constructor), 
            the attribute will be assigned to solver''' % intersection[0])
        for elem in intersection:
            setattr(p.solver, elem, getattr(p, elem))

    solver = p.solver.__solver__
    
    for key, value in kwargs.items():
        if hasattr(p.solver, key):
            if isConverter:
                solver_params[key] = value
            else:
                setattr(p.solver, key, value)
        elif hasattr(p, key) or key in ('lb', 'ub'):
            setattr(p, key, value)
        else: 
            p.warn('incorrect parameter for prob.solve(): "' + str(key) + '" - will be ignored (this one has been not found in neither prob nor ' + p.solver.__name__ + ' solver parameters)')
    if p.probType == 'EIG' and 'goal' in kwargs:
        p.err('for EIG parameter "goal" should be used only in class instance definition, not in "solve" method')
        
    setSomeDefaults(p)
    
    T = time()
    C = perf_counter()
    try:
        p._Prepare()
    except OpenOptException as ooe:
        p.err(ooe.msg)
    except FuncDesignerException as fde:
        p.err(fde.msg)
    
    p.initTime = time() - T
    p.initCPUTime = perf_counter() - C
    if p.initTime > 1 or p.initCPUTime > 1:
        p.disp('Initialization: Time = %0.1f CPUTime = %0.1f' % (p.initTime, p.initCPUTime))
        
    for fn in ['FunEvals', 'Iter', 'Time', 'CPUTime']:
        if hasattr(p,'min'+fn) and hasattr(p,'max'+fn) and getattr(p,'max'+fn) < getattr(p,'min'+fn):
            p.warn('min' + fn + ' (' + str(getattr(p,'min'+fn)) +') exceeds ' + 'max' + fn + '(' + str(getattr(p,'max'+fn)) +'), setting latter to former')
            setattr(p,'max'+fn, getattr(p,'min'+fn))

    for fn in ['maxFunEvals', 'maxIter']: setattr(p, fn, int(getattr(p, fn)))# to prevent warnings from numbers like 1e7

    if hasattr(p, 'x0'): 
        try:
            p.x0 = atleast_1d(asfarray(p.x0).copy())
        except NotImplementedError:
            p.x0 = asfarray(p.x0.tolist())
    for fn in ['lb', 'ub', 'b', 'beq']:
        if hasattr(p, fn):
            fv = getattr(p, fn)
            if fv is not None:# and fv != []:
                if str(type(fv)) == "<class 'map'>":
                    p.err("Python3 incompatibility with previous versions: you can't use 'map' here, use rendered value instead")
                setattr(p, fn, asfarray(fv).flatten())
            else:
                setattr(p, fn, asfarray([]))

    if p.solver._requiresFiniteBoxBounds:
        ind1, ind2 = isinf(p.lb), isinf(p.ub)
        if isscalar(p.implicitBounds): p.implicitBounds = (-p.implicitBounds, p.implicitBounds) # may be from lp2nlp converter, thus omit nlp init code
        p.lb[ind1] = p.implicitBounds[0] if asarray(p.implicitBounds[0]).size == 1 else p.implicitBounds[0][ind1]
        p.ub[ind2] = p.implicitBounds[1] if asarray(p.implicitBounds[1]).size == 1 else p.implicitBounds[0][ind2]


#    if p.lb.size == 0:
#        p.lb = -inf * ones(p.n)
#    if p.ub.size == 0:
#        p.ub = inf * ones(p.n)

    p.stopdict = {}

    for s in ['b','beq']:
        if hasattr(p, s): setattr(p, 'n'+s, len(getattr(p, s)))

    #if p.probType not in ['LP', 'QP', 'MILP', 'LLSP']: p.objFunc(p.x0)

    p.isUC = p._isUnconstrained()
    
    isIterPointAlwaysFeasible = p.solver.__isIterPointAlwaysFeasible__ if type(p.solver.__isIterPointAlwaysFeasible__) == bool \
        else p.solver.__isIterPointAlwaysFeasible__(p)
    if isIterPointAlwaysFeasible:
        #assert p.data4TextOutput[-1] == 'log10(maxResidual)'
        if p.data4TextOutput[-1] == 'log10(maxResidual)': 
            p.data4TextOutput = p.data4TextOutput[:-1]
#        else:
#            p.err('bug in runProbSolver.py')
    elif p.useScaledResidualOutput:
        p.data4TextOutput[-1] = 'log10(MaxResidual/ConTol)'

    if p.showFeas and p.data4TextOutput[-1] != 'isFeasible': 
        p.data4TextOutput.append('isFeasible')
    # TODO: rework
    if p.solver.__name__ == 'interalg' and p.isObjFunValueASingleNumber:
        p.data4TextOutput = p.data4TextOutput[:1] + \
        ['f*_distance_estim', 'f*_bound_estim'] + p.data4TextOutput[1:]
    if p.maxSolutions != 1:
        p._nObtainedSolutions = 0
        p.data4TextOutput.append('nSolutions')

    if not p.solver.iterfcnConnected:
        if SMALL_DELTA_X in p.kernelIterFuncs: p.kernelIterFuncs.pop(SMALL_DELTA_X)
        if SMALL_DELTA_F in p.kernelIterFuncs: p.kernelIterFuncs.pop(SMALL_DELTA_F)

    if not p.solver._canHandleScipySparse:
        if hasattr(p.A, 'toarray'): p.A = p.A.toarray()
        if hasattr(p.Aeq, 'toarray'): p.Aeq = p.Aeq.toarray()
    
    if isinstance(p.A, ndarray) and type(p.A) != ndarray: # numpy matrix
        p.A = p.A.A 
    if isinstance(p.Aeq, ndarray) and type(p.Aeq) != ndarray: # numpy matrix
        p.Aeq = p.Aeq.A 

    if hasattr(p, 'optVars'):
        p.err('"optVars" is deprecated, use "freeVars" instead ("optVars" is not appropriate for some prob types, e.g. systems of (non)linear equations)')

#    p.xf = nan * ones([p.n, 1])
#    p.ff = nan
    #todo : add scaling, etc
    p.primalConTol = p.contol
    if not p.solver.__name__.startswith('interalg'): 
        p.contol *= ConTolMultiplier

    p.timeStart = time()
    p.cpuTimeStart = perf_counter()
    
    # TODO: move it into solver parameters
    if p.probType not in ('MINLP', 'IP'):
        p.plotOnlyCurrentMinimum = p.__isNoMoreThanBoxBounded__()


    ############################
    # Start solving problem:

    if p.iprint >= 0:
        p.disp('\n' + '-'*25 + ' OpenOpt %s ' % version + '-'*25)
        pt = p.probType if p.probType != 'NLSP' else 'SNLE'
        s = 'problem: ' + p.name + '   type: %s' % pt
        if p.showGoal: 
            s += '    goal: %s' % p.goal
        p.disp(s)

        s = 'solver: ' +  p.solver.__name__ 
        
        if p.solver.__name__ == 'interalg':
            
            if p.isObjFunValueASingleNumber:
                if p.fTol is None and p.probType in ('IP', 'ODE'):
                    p.pWarn('''
                    please use p.fTol parameter for interalg with this problem;
                    p.ftol currently will be used instead
                    ''')
                    p.fTol = p.ftol
                if p.fTol is None:
                    p.fTol = 1e-7
                    p.warn('''
                    interalg requires p.fTol value (required objective function tolerance); 
                    10^-7 will be used'''
                    )
                s += '   fTol: %0.2g ' % p.fTol 
                if p.probType != 'IP':# TODO: implement for IP
                    s += '   rTol: %0.1g ' % p.rTol 
            elif p.probType in ('SNLE', 'NLSP'):
                s += '   maxSolutions: %s ' % p.maxSolutions 
        p.disp(s)

        s = 'variables: %d' % p.n
        #TODO: add num of smooth, bool, int, discrete, fixed, mb string & categorical variables
        p.disp(s)

        if not p.__isNoMoreThanBoxBounded__():
            s = 'constaints:  '
            hasLinearC = p.nb != 0 or p.nbeq != 0
            hasNonLinearC = p.nc != 0 or p.nh != 0
            if p.hasLogicalConstraints:
                s += '{%d logical}' % p.nLogicalConstraints
                if hasLinearC or hasNonLinearC:
                    s += '   '
            if hasLinearC:
                s += 'linear {%d ineq  %d eq}' % (p.nb, p.nbeq)
                if hasNonLinearC:
                    s += '   '
            if hasNonLinearC:
                s += 'non-linear {%d ineq  %d eq}' % (p.nc, p.nh-p.nLogicalConstraints)
            p.disp(s)

    p.extras = {}

#    solver(p)
    try:
#        if 0 and isConverter:
#            pass
#            # TODO: will R be somewhere used?
#            #R = converter(solverName, **solver_params)
#        else:
        nErr = check(p)
        if nErr: p.err("prob check results: " +str(nErr) + "ERRORS!")#however, I guess this line will be never reached.
        if p.probType not in ('IP', 'EIG'): p.iterfcn(p.x0)
        originalName = setProcTitle(p)
        
#        import os, cProfile
#        print '----'
#        print 'p.n, p.nb, p.nbeq, p.nc, p.nh:',  p.n, p.nb, p.nbeq, p.nc, p.nh
#        if p.solver.__name__ == 'interalg':
#            cProfile.runctx( 'solver(p)', globals(), locals(), filename="t.profile" )
#        #os.system('runsnake t.profile')
#        else:
#            solver(p)
            
        solver(p)
        
#        if p.solver.__name__ == 'interalg':
#            
#            print p.extras
#            os.system('pyprof2calltree -i t.profile -o myfile.prof.grind')
#            os.system('kcachegrind myfile.prof.grind')
#        print p.xk
#        print '===='

#    except killThread:
#        if p.plot:
#            print 'exiting pylab'
#            import pylab
#            if hasattr(p, 'figure'):
#                print 'closing figure'
#                #p.figure.canvas.draw_drawable = lambda: None
#                pylab.ioff()
#                pylab.close()
#                #pylab.draw()
#            #pylab.close()
#            print 'pylab exited'
#        return None
    except isSolved:
#        p.fk = p.f(p.xk)
#        p.xf = p.xk
#        p.ff = p.objFuncMultiple2Single(p.fk)

        if p.istop == 0: p.istop = 1000
    except OpenOptException as ooe:
        p.istop = -1
        p.err(ooe.msg)
    except FuncDesignerException as fde:
        p.err(fde.msg)
    finally:
        if p.isFDmodel:
            for v in p.freeVarsSet | p.fixedVarsSet:
                if v.fields != ():
                    v.domain, v.aux_domain = v.aux_domain, v.domain
        seterr(**old_err)

    if hasSetproctitleModule and originalName is not None:
        setproctitle.setproctitle(originalName)
    if len(errSet) != 0:
        finalGUIroutines(p)
        p._immutable = True
        return None
    ############################
    p.contol = p.primalConTol

    # Solving finished
    if hasattr(p, '_bestPoint') and not any(isnan(p._bestPoint.x)) and p.probType != 'ODE':
        try:
            p.iterfcn(p._bestPoint)
        except isSolved:
            pass
            
    if p.probType != 'EIG':
        if not hasattr(p, 'xf') and not hasattr(p, 'xk'): 
            p.xf = p.xk = ones(p.n)*nan
        if hasattr(p, 'xf') and (not hasattr(p, 'xk') or array_equal(p.xk, p.x0)): 
            p.xk = p.xf
        if not hasattr(p,  'xf') or all(isnan(p.xf)): 
            p.xf = p.xk
        if p.xf is nan: 
            p.xf = p.xk = ones(p.n)*nan
        
        if p.isFeas(p.xf) and (not p.probType=='MINLP' or p.discreteConstraintsAreSatisfied(p.xf)):
            p.isFeasible = True
        else: p.isFeasible = False
    else:
        p.isFeasible = True # check it!
    
    p.isFinished = True # After the feasibility check above!
    p._localtime_finished = localtime()
    
    if p.probType == 'MOP':
        p.isFeasible = True
    elif p.probType == 'IP':
        p.isFeasible = p.rk < p.ftol
    else:
        p.ff = p.fk = p.objFunc(p.xk)
        
        # walkaround for PyPy:
        if type(p.ff) == ndarray and p.ff.size == 1:
            p.ff = p.fk = asscalar(p.ff)
        
    if not hasattr(p,  'ff') or any(isnan(p.ff)): 
        p.iterfcn, tmp_iterfcn = lambda *args: None, p.iterfcn
        p.ff = p.fk
        p.iterfcn = tmp_iterfcn

    if p.invertObjFunc:  
        p.fk, p.ff = -p.fk, -p.ff

    if asfarray(p.ff).size > 1: 
        p.ff = p.objFuncMultiple2Single(p.fk)

    #p.ff = p.objFuncMultiple2Single(p.ff)
    #if not hasattr(p, 'xf'): p.xf = p.xk
    if type(p.xf) in (list, tuple) or isscalar(p.xf): 
        p.xf = asarray(p.xf)
    p.xf = p.xf.flatten()
    p.rf = p.getMaxResidual(p.xf) if p.probType not in ('IP', 'ODE') else p.rk

    if not p.isFeasible and p.istop > 0: p.istop = -100-p.istop/1000.0
    
    if p.istop == 0 and p.iter >= p.maxIter:
        p.istop, p.msg = IS_MAX_ITER_REACHED, 'Max Iter has been reached'
    
    p.stopcase = stopcase(p)

    p.xk, p.rk = p.xf, p.rf
    if p.invertObjFunc: 
        p.fk = -p.ff
        p.iterfcn(p.xf, -p.ff, p.rf)
    else: 
        p.fk = p.ff
        p.iterfcn(p.xf, p.ff, p.rf)
        
    p.lb, p.ub = p._lb, p._ub
    checkImplicitBounds(p)

    p.__finalize__()
    if not p.storeIterPoints: delattr(p.iterValues, 'x')

    r = OpenOptResult(p)

    #TODO: add scaling handling!!!!!!!
#    for fn in ('df', 'dc', 'dh', 'd2f', 'd2c', 'd2h'):
#        if hasattr(p, '_' + fn): setattr(r, fn, getattr(p, '_'+fn))

    p.invertObjFunc = False
    
    if p.isFDmodel:
        p.x0 = p._x0

    finalTextOutput(p, r)
    finalGUIroutines(p)
    
    finish = p.finish if type(p.finish) in (list, tuple, set) else [p.finish]
    for f in finish:
        f(r, p)
    p._immutable = True
    
    return r
    
##########################################################
def finalGUIroutines(p):
    if p.isManagerUsed: 
        p.GUI_items['Quit'].config(state='normal')
        p.GUI_items['time'].set('%d' % (time() - p.timeStart))
        p.GUI_items['cputime'].set('%d' % (perf_counter() - p.cpuTimeStart))
    elif p.plot:
        finalShow(p)
        
##########################################################
def finalTextOutput(p, r):
    if type(p.msg) == ndarray: # from fmincon etc
        msg = asarray(p.msg.flatten(),int).tolist()
        p.msg = ''.join(chr(c) for c in msg)
    if p.msg == '':
        msg = 'finished'
    if p.iprint >= 0:
        if len(p.msg):  
            p.disp("istop: " + str(r.istop) + ' (' + p.msg +')')
        else: 
            p.disp("istop: " + str(r.istop))

        p.disp('Solver:   Time Elapsed = %0.3g  \tCPU Time Elapsed = %0.3g' % (r.elapsed['solver_time'], r.elapsed['solver_cputime']))
        if p.plot:
            p.disp('Plotting: Time Elapsed = %0.3g  \tCPU Time Elapsed = %0.3g' % (r.elapsed['plot_time'], r.elapsed['plot_cputime']))
        
        if p.probType == 'MOP':
            msg = '%d solutions have been obtained' % len(p.solutions.coords)
            p.disp(msg)
            return
        
        # TODO: add output of NaNs number in constraints (if presernt)
        rMsg = 'max(residuals/requiredTolerances) = %g' % \
        (r.rf / p.contol) if p.useScaledResidualOutput else 'MaxResidual = %g' % r.rf

        if not p.isFeasible:
            nNaNs = (len(where(isnan(p.c(p.xf)))[0]) if p.c is not None and type(p.c)!=list else 0) \
            + (len(where(isnan(p.h(p.xf)))[0]) if p.h is not None and type(p.h)!=list else 0)
            if nNaNs == 0:
                nNaNsMsg = ''
            elif nNaNs == 1:
                nNaNsMsg = '1 constraint is equal to NaN, '
            else:
                nNaNsMsg = ('%d constraints are equal to NaN, ' % nNaNs)
            p.disp('NO FEASIBLE SOLUTION has been obtained (%s%s, objFunc = %0.8g)' % (nNaNsMsg,  rMsg, r.ff))
        else:
            if p.maxSolutions == 1:
                msg = "objFuncValue: " + (p.finalObjFunTextFormat % r.ff)
                if not p.isUC: msg += ' (feasible, %s)' % rMsg
            else:
                msg = '%d solutions have been obtained' % len(p.solutions)
            p.disp(msg)

##########################################################
def finalShow(p):
    if not p.plot: return
    pylab = __import__('pylab')
    pylab.ioff()
    if p.show:
#        import os
#        if os.fork():
            pylab.show()

##########################################################
def checkImplicitBounds(p):
    if p.solver._requiresFiniteBoxBounds and p.implicitBounds is not None:
        delta = 0.01 * (p.implicitBounds[1] - p.implicitBounds[0])
        tmp = where(isinf(p.lb), p.implicitBounds[0], -inf)
        delta_r = p.xk - tmp
        tmp = where(isinf(p.ub), p.implicitBounds[1], inf)
        delta_l = tmp - p.xk
        undefined_x = np.all(np.isnan(p.xk))# & np.isnan(p.xf))
        if not undefined_x and (any(delta_r<delta) or any(delta_l<delta)):
            Msg = '''
            some coordinates of solution are too close to implicitBounds
            (difference is less than 1%)
            '''
            p.warn(Msg)
            if p.istop > 0:
                p.istop = p.stopcase = 0 
                p.msg += ';' + Msg

##########################################################
def setProcTitle(p):
    originalName = None
    if hasSetproctitleModule:
        try:
            originalName = setproctitle.getproctitle()
            if originalName.startswith('OpenOpt-'):
                originalName = None
            else:
                s = 'OpenOpt-' + p.solver.__name__
                # if p.name != 'unnamed':
                s += '-' + p.name
                setproctitle.setproctitle(s)
        except:
            pass
    else:
        p.pWarn('''
        please install setproctitle module 
        (it's available via easy_install 
        and Linux software channels like apt-get)''')
    return originalName

##########################################################
def setSomeDefaults(p):
    p.iterValues = EmptyClass()
    p.iterCPUTime = []
    p.iterTime = []
    p.iterValues.x = [] # iter points
    p.iterValues.f = [] # iter ObjFunc Values
    p.iterValues.r = [] # iter MaxResidual
    p.iterValues.rt = [] # iter MaxResidual Type: 'c', 'h', 'lb' etc
    p.iterValues.ri = [] # iter MaxResidual Index
    p.solutions = [] # list of solutions, may contain several elements for interalg and mb other solvers
    if p._baseClassName == 'NonLin':
        p.iterValues.nNaNs = [] # number of constraints equal to numpy.nan

    if p.goal in ['max','maximum']: p.invertObjFunc = True

    #TODO: remove it!
    p.advanced = EmptyClass()

    p.istop = 0
    p.iter = 0
    p.graphics.nPointsPlotted = 0
    p.finalIterFcnFinished = False
    #for fn in p.nEvals.keys(): p.nEvals[fn] = 0 # NB! f num is used in LP/QP/MILP/etc stop criteria check

    p.msg = ''
    if not type(p.callback) in (list,  tuple): 
        p.callback = [p.callback]
    if hasattr(p, 'xlabel'): 
        p.graphics.xlabel = p.xlabel
    if p.graphics.xlabel == 'nf': 
        p.iterValues.nf = [] # iter ObjFunc evaluation number

