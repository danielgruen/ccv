"""
This file contains a CCV interface object.
It is initialized with the inputs that we
give to CCV and then we can generate models with it.
It also contains functions needed to replace
its internal variables.

This is sort a work in progress since the pathing will
be all screwed up when this rolls into production.

This file assumes that the proper files have been compiled
with the Makefile in ccv/. This needs to be done for
each redshift or more specifically each
set of annuli.
"""
import os

class CCV_interface(object):
    def __init__(self,output_prefix,zlens,annuli_prefix,\
                 Mass,zsource=None,pz_prefix=None,\
                 amplitudes=None):
        if zsource is None and pz_prefix is None:
            raise Exception("Require either a zsource or a pz_prefix.")
        if zsource is not None and pz_prefix is not None:
            raise Exception("Cannot have both a zsource and a pz_prefix.")
        self.outpath = output_prefix+".fits"
        self.zlens = zlens
        self.annuli_prefix = annuli_prefix
        self.Mass = Mass
        if zsource is not None: self.zsource = zsource
        if pz_prefix is not None: self.pz_prefix = pz_prefix
        if amplitudes is None: 
            self.use_defaults = True
        else:
            self.use_defaults = False
            self.cconc = amplitudes['cconc']
            self.ccorr = amplitudes['ccorr']
            self.cell  = amplitudes['cell']
            self.coff  = amplitudes['coff']
            self.clss  = amplitudes['clss']
        return

    def __str__(self):
        string = ''
        for a in dir(self): 
            if not a.startswith('__'): string+=a+"\n"
        return string

    def call_ccv(self):
        """
        Run the getmodel_ds from ccv and 
        """
        if hasattr(self,'zsource'):
            command = "./src/getmodel_ds %s %s %s %s %s"%(self.outpath,self.zlens,\
                                                         self.annuli_prefix,self.zsource,\
                                                         self.Mass)
        else:
            command = "./src/getmodel_ds %s %s %s %s %s"%(self.outpath,self.zlens,\
                                                         self.annuli_prefix,self.pz_prefix,\
                                                         self.Mass)
        if not self.use_defaults:
            command += " %.2f %.2f %.2f %.2f %.2f"%(self.cconc,self.ccorr,self.cell,self.coff,self.clss)
        print "Running the following command:"
        print command
        os.system(command)
        return

#A basic unit test
if __name__=="__main__":
    output_prefix = "test"
    zlens = "0.24533"
    annuli_prefix = "thetas_z0_l0"
    Mass = "5e14"
    zsource = None
    pz_prefix = "pz"
    amplitudes = {}
    amplitudes['cconc'] = 1.0
    amplitudes['ccorr'] = 1.0
    amplitudes['cell']  = 1.0
    amplitudes['coff']  = 1.0
    amplitudes['clss']  = 1.0
    ccv = CCV_interface(output_prefix,zlens,annuli_prefix,Mass,zsource,pz_prefix,amplitudes)
    ccv.call_ccv()
    print "Unit test complete"
