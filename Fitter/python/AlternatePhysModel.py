from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
import ROOT, os, re
import json
import numpy as np

class AlternatePhysicsModel(PhysicsModel):
  
    def __init__(self):
        PhysicsModel.__init__(self)
        self.poiNames=[]
        self.numOperators = {}
        self.Operators = {}
        self.verbose=False

        self.top22006_sgnl_known = ['ttH','tllq','ttll','ttlnu','tHq','tttt']
        
        # regular expressions for process names. relevant for top-22-006:
        self.sm_re    = re.compile('(?P<proc>.*)_sm')
        self.lin_re   = re.compile('(?P<proc>.*)_lin_(?P<c1>.*)')
        self.quad_re  = re.compile('(?P<proc>.*)_quad_(?P<c1>.*)')
        self.mixed_re = re.compile('(?P<proc>.*)_quad_mixed_(?P<c1>.*)_(?P<c2>.*)') # should go before quad
        
        #the following three options are relevant for top-21-003
        self.linear_only = False
        self.fits = None  
        self.procbins = []
        
        self.wcs = [ # Hardcoded currently...
            'ctW',
            'ctp',
            'cpQM',
            'ctZ',
            #'ctG',
            'cbW',
            'cpQ3',
            'cptb',
            'cpt',
        ]

    def loadOperators(self,fpath):
        print("Loading operators from {fpath}".format(fpath=fpath))
        jsn = open(fpath,'r').read()
        operators = json.loads(jsn)
        self.alloperators = []

        for sig in self.top22006_sgnl_known:
            self.Operators[sig] = operators[sig]
            self.numOperators[sig] = len(operators[sig])
            self.alloperators.extend(operators[sig])
        self.alloperators = list(set(self.alloperators))
        print("Operators: {ops}".format(ops=self.Operators))

    def setPhysicsOptions(self, physOptions):
        for po in physOptions:
            if po.startswith("selectedWCs="):
                selected_wcs_fpath = po.replace("selectedWCs=","").strip()
                self.loadOperators(selected_wcs_fpath)

            if po.startswith("fits="):
                self.fits = po.replace("fits=","").strip() 
                fits = np.load(self.fits)
                self.procbins.extend(fits.item().keys())            #only relevant when using fits=EFTParam_v8.npy
           
    def doParametersOfInterest(self):
        self.modelBuilder.doVar("r[1,-10,10]")
        self.poiNames="r"
        for operator in self.alloperators:
            self.modelBuilder.doVar(operator + "[0,-200,200]")
            self.poiNames += "," + operator

        self.modelBuilder.doSet("POI", self.poiNames)
        
        #top22006 part
        self.quadFactors = []
        for sig in self.top22006_sgnl_known:
            sgnl_ops = self.Operators[sig]
            if self.numOperators[sig] != 1:
                for op in sgnl_ops:
                    op_idx = sgnl_ops.index(op)
                    oplist = sgnl_ops[op_idx:]
                    terms = " ".join(['-@0*@%d'%(i+1) for i,x in enumerate(oplist[1:])])
                    formula = "@0{}".format(terms)
                    expression = "expr::func_{sig}_sm_{op}(\"{formula}\", {oplist})".format(
                        sig=sig,
                        op=op,
                        formula=formula,
                        oplist=','.join(oplist)
                    )
                    if self.verbose:
                        print("SM sub-expr: {}".format(expression))
                    self.modelBuilder.factory_(expression)
 
                terms = "-@".join(str(i+1) for i,x in enumerate(sgnl_ops))
                formula = "@0*(1-@{})".format(terms)
                functions = ", ".join('func_{sig}_sm_{op}'.format(sig=sig,op=op) for op in self.Operators[sig])
                expression = "expr::func_{sig}_sm(\"{formula}\", 1, {funcs})".format(
                    sig=sig,
                    formula=formula,
                    funcs=functions
                )
                if self.verbose:
                    print("SM expr: {}".format(expression))
                self.modelBuilder.factory_(expression)
            else:
                formula = "@0*(1-(@1))"
                expression = "expr::func_{sig}_sm(\"{formula}\", 1, {op})".format(
                    sig=sig,
                    formula=formula,
                    op=sgnl_ops[0]
                )
                if self.verbose:
                    print("SM expr: {}".format(expression))
                self.modelBuilder.factory_(expression)
        # This is the coefficient of "SM + Lin_i + Quad_i"
        for sig in self.top22006_sgnl_known:
            sgnl_ops = self.Operators[sig]
            if self.numOperators[sig] != 1:
                for i in range(0,self.numOperators[sig]):
                    op = sgnl_ops[i]
                    oplist = [str(sgnl_ops[j]) for j in range(len(sgnl_ops)) if i != j]
                    terms = "+@".join([str(j+2) for j in range(len(sgnl_ops) - 1)])
                    formula = "@0*(@1 * (1-(@{terms})))".format(terms=terms)
                    expression = "expr::func_{sig}_sm_linear_quadratic_{op}(\"{formula}\", 1,{op}, {oplist})".format(
                        sig=sig,
                        op=op,
                        formula=formula,
                        oplist=", ".join(oplist)
                    )
                    if self.verbose:
                        print("SM+L+Q expr: {}".format(expression))
                    self.modelBuilder.factory_(expression)

            else:
                formula = "@0*(@1)"
                expression = "expr::func_{sig}_sm_linear_quadratic_(\"{formula}\", 1, {op})".format(
                    sig=sig,
                    formula=formula,
                    op=sgnl_ops[0]
                )
                if self.verbose:
                    print("SM+L+Q expr: {}".format(expression))
                self.modelBuilder.factory_(expression)

        # Quadratic term in each Wilson coefficient
        # e.g. expr::func_sm_linear_quadratic_cH("@0*(@1 * (1-2*(@2+@3) ))", 1,cH, cG, cGtil)
        for sig in self.top22006_sgnl_known:
            sgnl_ops = self.Operators[sig]
            # This is the coefficient of "Quad_i"
            for i in range(0,self.numOperators[sig]):
                formula = "@0*(@1*@1-@1)"
                expression = "expr::func_{sig}_quadratic_{op}(\"{formula}\", 1, {op})".format(
                    sig=sig,
                    op=sgnl_ops[i],
                    formula=formula
                )
                if self.verbose:
                    print("Quad expr: {}".format(expression))
                self.modelBuilder.factory_(expression)
            # (SM + linear) + quadratic + interference between pairs of Wilson coefficients
            if self.numOperators[sig] != 1:
                for i in range(0, self.numOperators[sig]):
                    for j in range(i+1,self.numOperators[sig]):
                        # This is the coefficient of "SM + Lin_i + Lin_j + Quad_i + Quad_j + 2 * M_ij"
                        func_name = "func_{sig}_sm_linear_quadratic_mixed_{op1}_{op2}".format(
                            sig=sig,
                            op1=sgnl_ops[j],    # Note: This is the order as it was before the code cleanup
                            op2=sgnl_ops[i]
                        )
                        self.quadFactors.append(func_name)
                        formula = "@0*@1*@2"
                        expression = "expr::{name_top22006}(\"{formula}\", 1,{op1},{op2})".format(
                            name_top22006=func_name,
                            formula=formula,
                            op1=sgnl_ops[i],
                            op2=sgnl_ops[j]
                        )
                        if self.verbose:
                            print("Mixed expr: {}".format(expression))
                        self.modelBuilder.factory_(expression)
        print(" parameters of interest = {}".format(self.poiNames))
        print(" self.numOperators = {}".format(self.numOperators))


        ##top-21-003 part
        fits = np.load(self.fits)                                          #for EFTParam_v8.npy file, this is a numpy array. Next line converts this to a dictionary.
        fits = dict(zip(fits.item().keys(),fits.item().values()))          #Needed when using EFTParam_v8.npy file. 
        for i,procbin in enumerate(sorted(self.procbins)):
            name = 'r_{proc}_{cat}'.format(proc=procbin[0],cat=procbin[1])
            procbin_name = '_'.join(procbin)
            if not self.modelBuilder.out.function(name):
                # Initialize function pieces
                constant = '{}'.format(fits[procbin][('sm','sm')]) # constant term (should be 1)
                lin_name = procbin_name+"_L" # Name of linear function
                lin_term = [] # Linear term
                lin_args = [] # List of wcs in linear term
                quartic_names = [procbin_name+"_Q"+str(idx) for idx,wc in enumerate(self.wcs)] # Names of quadratic functions
                quartic_terms = [[] for wc in self.wcs] # Quartic terms, but split into chunks
                quartic_args = [[] for wc in self.wcs] # List of wcs in quartic terms
                fit_terms = [constant] # List of fit terms

                # Fill function pieces
                for idx,wc1 in enumerate(self.wcs):
                    #if abs(fits[procbin][('sm',wc1)]) >= 0.001:
                    if abs(fits[procbin][('sm',wc1)]) >= 0.0:
                        if fits[procbin][('sm',wc1)] < 0.0:
                            lin_term.append('({0})*{1}'.format(fits[procbin][('sm',wc1)],wc1))
                        else:
                            lin_term.append('{0}*{1}'.format(fits[procbin][('sm',wc1)],wc1))
                        lin_args.append(wc1)
                    for idy,wc2 in enumerate(self.wcs):
                        if self.linear_only: continue
                        key = (wc1,wc2)
                        if (wc1,wc2) not in fits[procbin]: key = (wc2, wc1)
                        #if (idy >= idx) and (abs(fits[procbin][(wc1,wc2]) >= 0.001):
                        if (idy >= idx) and (abs(fits[procbin][key]) >= 0.000):
                            quartic_terms[idx].append('{0}*{1}*{2}'.format(fits[procbin][key],key[0],key[1])) # wc1, wc2
                            quartic_args[idx].extend([key[0],key[1]]) # wc1 , wc2
                # Compile linear function for combine
                if lin_term:
                    lin_expr = "expr::{lin_name}('{lin_term}',{lin_args})".format(lin_name=lin_name,lin_term="+".join(lin_term),lin_args=",".join(lin_args))
                    lin_func = self.modelBuilder.factory_(lin_expr)
                    self.modelBuilder.out._import(lin_func)
                    fit_terms.append(lin_name)
                # Compile quartic functions separately first
                for idx,fn in enumerate(quartic_terms):
                    if not fn: continue # Skip empty quartic functions
                    quartic_expr = "expr::{quartic_names}('{quartic_terms}',{quartic_args})".format(quartic_names=quartic_names[idx], quartic_terms="+".join(fn), quartic_args=",".join(list(set(quartic_args[idx]))))
                    quartic_func = self.modelBuilder.factory_(quartic_expr)
                    self.modelBuilder.out._import(quartic_func)
                    fit_terms.append(quartic_names[idx])
                # Compile the full function
                fit_function = "expr::{name}('{fit_terms}',{fit_args})".format(name=name,fit_terms="+".join(fit_terms),fit_args=",".join(fit_terms[1:]))
                quadratic = self.modelBuilder.factory_(fit_function)

                # Export fit function
                self.modelBuilder.out._import(quadratic)
 

    def getYieldScale(self, bin, process):
        if bin.startswith("top22006"):
            if any( process.startswith(x) for x in self.top22006_sgnl_known):
                if self.sm_re.search(process):
                    match = self.sm_re.search(process)
                    return "func_{p}_sm".format(p=match.group('proc'))
                elif self.lin_re.search(process):
                    match = self.lin_re.search(process)
                    return "func_{p}_sm_linear_quadratic_{c1}".format(p=match.group('proc'),c1=match.group('c1'))
                elif self.mixed_re.search(process):
                    match = self.mixed_re.search(process)
                    c1 = match.group('c1')
                    c2 = match.group('c2')
                    proc = match.group('proc')
                    name = "func_{p}_sm_linear_quadratic_mixed_{c1}_{c2}".format(p=proc,c1=c1,c2=c2)
                    if name in self.quadFactors:
                        return "func_{p}_sm_linear_quadratic_mixed_{c1}_{c2}".format(p=proc,c1=c1,c2=c2)
                    else:
                        return "func_{p}_sm_linear_quadratic_mixed_{c2}_{c1}".format(p=proc,c1=c1,c2=c2)
                elif self.quad_re.search(process):
                    match = self.quad_re.search(process)
                    c1 = match.group('c1')
                    proc = match.group('proc')
                    return "func_{p}_quadratic_{c1}".format(p=proc,c1=c1)
                else:
                    #raise RuntimeError("Undefined process %s"%process)
                    print("Undefined process {}, probably below threshold and ignored".format(process))
            else:
                print("Process {} not a signal".format(process))
                return 1

        elif bin.startswith("top21003"):
            if (process,bin) not in self.procbins:
                return 1
            else:
                name = 'r_{0}_{1}'.format(process,bin)
                return name
alternatePhysicsModel = AlternatePhysicsModel()
