import numpy as np
from neuron import h
from neuron import load_mechanisms
from . import neuronFunctions as nfx
from . import get_src_dir

src_directory = get_src_dir()
load_mechanisms(f'{src_directory}/mod.files/')
h('objref nil')


class L23(object):
    def __init__(self,cellID=0,cutExperiment=0,dendNa=[6,2,200,False],dendK=[0.048,100,False,False],dxSeg=5,fixDiam=None):
        """
        Initializer for L23 object. cellID determines which cell to load (see "loadHoc" function of this class for file names)
        """
        if cellID>7 or cellID<0:
            raise ValueError(f"CellIDs can only be 0-7, cellID={cellID} was provided.")
        self.cellID = cellID
        self.loadHoc(cellID,cutExperiment) # Load Hoc File
        self.props(dendNa=dendNa,dendK=dendK) # Load properties into class
        self._topol(cellID,dxSeg,fixDiam) # Construct topology (and load names of sections into L23 object)
        self._addAxon() # Add axon
        self._biophys() # Add passive properties
        self.init_active(dendNa=dendNa, dendK=dendK)
        self.defineTargetROIs(cellID)
        
    def loadHoc(self,cellID,cutExperiment):
        hocFileName = {0:'ATL200625a',1:'ATL200626b',2:'ATL200910c',3:'ATL201007b',4:'ATL201009b',5:'ATL201012b',6:'ATL201020a',7:'ATL201021c'}
        cutName = ''
        if cutExperiment!=0:
            cutName = f"_dcut{cutExperiment}"
        #print('Creating cell : {fName}{cName}'.format(fName=hocFileName[cellID],cName=cutName))
        file_name = f"{get_src_dir()}/silentDendrite_uncageMappingMorphologies/{hocFileName[cellID]}{cutName}.hoc"    
        h(f'xopen("{file_name}")');
        
    def props(self,dendNa,dendK):
        # Passive properties
        self.CM = 1.0
        self.RM = 7000.0
        self.RA = 100.0
        self.E_PAS = -75
        self.CELSIUS = 35

        # Active properties
        self.Ek = -90
        self.Ena = 60
        self.Eca = 140

        # From Lorincz & Nusser, 2010
        self.singleChannelNA = 17
        
        self.gna_ais = self.singleChannelNA*150 #Ujfalussy: 1000
        self.gkv_ais = 100 # Ujfalussy: 100
        self.gna_axon = self.singleChannelNA*10
        self.gkv_axon = 33

        self.gna_soma = self.singleChannelNA * dendNa[0] # Ujfalussy: 1000
        self.gkv_soma = 100 
        self.gkm_soma = 2.2 
        self.gkca_soma = 3 
        self.gca_soma = 0.5 
        self.git_soma = 0.0003 

        self.gna_dend = self.singleChannelNA*dendNa[0] # 80 - default
        self.gna_Lambda = dendNa[2] # Drop off parameter
        
        self.gkv_dend = 3 # 3 - default
        self.gkm_dend = 1 # dendK[0][1] # 1 - default
        self.gkca_dend = 3 # dendK[0][2] # 3 - default
        self.gk_Lambda = dendK[1] # Drop off parameter
        self.gka_base = dendK[0]
        
        self.gca_dend = 0.5
        self.git_dend = 0.00015 
        self.gh_dend = 0
        
        self.caTau = 20 # ms - use for decay and 1/gamma
        self.caDepth = 0.1 # µm how much to focus calcium in shell
        self.caMinCai = 75e-6 # mM - minimum calcium concentration
    
    def _topol(self,cellID,dxSeg,fixDiam):            
        self.soma = h.soma
        self.soma.nseg = 1
        self.dends = []
        
        # Load each section (which has already been defined by h(xopen..)) into the self.dends list
        secListL23 = returnSecListL23() #return name of every dendritic segment in list (excludes soma and axon)
        for secName in secListL23:
            self.dends.append(eval('h.{secName}'.format(secName=secName)))
        
        for sec in self.dends:
            sec.diam = sec.diam
            cnseg = np.ceil(sec.L/dxSeg).astype(int) # have a little bit more than 1 segment per every dxSeg
            if cnseg % 2 == 0: cnseg += 1 # make sure nseg is odd!
            sec.nseg = cnseg
            if fixDiam:
                sec.diam = fixDiam
                            
    def _addAxon(self):
        # Add initial segment and axon
        self.ais = h.Section(name='ais')
        self.ais.L = 30
        self.ais.diam = 2
        self.ais.nseg = 6
        
        self.axon = h.Section(name='axon')
        self.axon.L = 300
        self.axon.diam = 1
        self.axon.nseg = 7

        self.ais.connect(self.soma,1,0)
        self.axon.connect(self.ais,1,0)
        
    def _biophys(self):
        # Set up passive biophysics for soma
        self.soma.cm = self.CM
        self.soma.insert('pas')
        self.soma.e_pas = self.E_PAS
        self.soma.g_pas = 1.0/self.RM
        self.soma.Ra = self.RA
        
        # Set up passive biophysics for axon
        self.ais.cm = self.CM
        self.ais.insert('pas')
        self.ais.e_pas = self.E_PAS
        self.ais.g_pas = 1.0/self.RM
        self.ais.Ra = self.RA
        
        self.axon.cm = self.CM
        self.axon.insert('pas')
        self.axon.e_pas = self.E_PAS
        self.axon.g_pas = 1.0/self.RM
        self.axon.Ra = self.RA
        
        # Set up passive biophysics for dendrites
        for sec in self.dends:
            sec.cm = self.CM
            sec.insert('pas')
            sec.e_pas = self.E_PAS
            sec.g_pas = 1.0/self.RM
            sec.Ra = self.RA
   
    def init_active(self, dendNa, dendK):
        self.axon.insert('na'); self.axon.gbar_na = self.gna_axon
        self.axon.insert('kv'); self.axon.gbar_kv = self.gkv_axon
        self.axon.ena = self.Ena
        self.axon.ek = self.Ek
        
        self.ais.insert('na'); self.ais.gbar_na = self.gna_ais
        self.ais.insert('kv'); self.ais.gbar_kv = self.gkv_ais
        self.ais.ena = self.Ena
        self.ais.ek = self.Ek

        self.soma.insert('na'); self.soma.gbar_na = self.gna_soma
        self.soma.insert('kv'); self.soma.gbar_kv = self.gkv_soma
        self.soma.insert('km'); self.soma.gbar_km = self.gkm_soma
        self.soma.insert('kca'); self.soma.gbar_kca = self.gkca_soma
        self.soma.insert('ca'); self.soma.gbar_ca = self.gca_soma
        self.soma.insert('CaT'); self.soma.gbar_CaT = self.git_soma
        self.soma.insert('CaDynamics');
        self.soma.gamma_CaDynamics = 1/self.caTau
        self.soma.decay_CaDynamics = self.caTau
        self.soma.depth_CaDynamics = self.caDepth
        self.soma.minCai_CaDynamics = self.caMinCai
        self.soma.ena = self.Ena
        self.soma.ek = self.Ek
        self.soma.eca = self.Eca
        
        ###
        #if axon:
        #    model.axon.insert('na'); model.axon.gbar_na = model.gna_axon
        #    model.axon.insert('kv'); model.axon.gbar_kv = model.gkv_axon
        #    model.axon.ena = model.Ena
        #    model.axon.ek = model.Ek

        #if soma:
        #    model.soma.insert('na'); model.soma.gbar_na = model.gna_soma
        #    model.soma.insert('kv'); model.soma.gbar_kv = model.gkv_soma
        #    model.soma.insert('km'); model.soma.gbar_km = model.gkm_soma
        #    model.soma.insert('kca'); model.soma.gbar_kca = model.gkca_soma
        #    model.soma.insert('ca'); model.soma.gbar_ca = model.gca_soma
        #    model.soma.insert('it'); model.soma.gbar_it = model.git_soma
        #    model.soma.ena = model.Ena
        #    model.soma.ek = model.Ek
        #    model.soma.eca = model.Eca

        #if dend:
        #    for d in model.dends:
        #        d.insert('na'); d.gbar_na = model.gna_dend*dendNa
        #        d.insert('kv'); d.gbar_kv = model.gkv_dend
        #        d.insert('km'); d.gbar_km = model.gkm_dend
        #        d.insert('kca'); d.gbar_kca = model.gkca_dend
        #        d.insert('ca'); d.gbar_ca = model.gca_dend*dendCa
        #        d.insert('it'); d.gbar_it = model.git_dend*dendCa
        #        #d.insert('cad'); 
        #        d.ena = model.Ena
        #        d.ek = model.Ek
        #        d.eca = model.Eca

        #if calcium:
        #    model.soma.insert('cad')
        #    for d in model.dends:
        #        d.insert('cad')
        ###
        
        # If we're adding the A-Type channel, add it
        if dendK[2] is True:
            self.soma.insert('kap');
            self.soma.gkabar_kap = self.gka_base

        for d in self.dends:
            d.insert('na'); 
            d.insert('kv'); d.gbar_kv = self.gkv_dend
            d.insert('km'); d.gbar_km = self.gkm_dend
            d.insert('kca'); d.gbar_kca = self.gkca_dend
            d.insert('ca'); d.gbar_ca = self.gca_dend
            d.insert('CaT'); d.gbar_CaT = self.git_dend
            d.insert('CaDynamics');
            d.gamma_CaDynamics = 1/self.caTau
            d.decay_CaDynamics = self.caTau
            d.depth_CaDynamics = self.caDepth
            d.minCai_CaDynamics = self.caMinCai
        
            #d.insert('cad'); 
            d.ena = self.Ena
            d.ek = self.Ek
            d.eca = self.Eca
            
            if dendNa[3] is True:
                # Decay conductance exponentially
                for seg in range(d.nseg):
                    segSpot = nfx.returnSegment(d.nseg,seg+1)
                    segDistance = h.distance(self.soma(1),d(segSpot))
                    gDecaySodium = np.exp(-segDistance/self.gna_Lambda)
                    d(segSpot).gbar_na = self.singleChannelNA*((dendNa[0]-dendNa[1])*gDecaySodium+dendNa[1]) # take proximal value and decay exponentially
            else:
                d.gbar_na = self.gna_dend

            # Add K-A linearly with distance from soma
            if dendK[2] is True:
                d.insert('kap');
                d.insert('kad');
                for seg in range(d.nseg):
                    segSpot = nfx.returnSegment(d.nseg,seg+1)
                    segDistance = h.distance(self.soma(1),d(segSpot))
                    increaseAType = 1e-3 * segDistance / dendK[1] # increase 1mA/cm2 every 100µm of distance from soma
                    # Assume we're only using proximal version
                    d(segSpot).gkabar_kap = self.gka_base + increaseAType
                    d(segSpot).gkabar_kad = 0
                    # And switch to distal version at 100µm if requested
                    if dendK[3] is True and segDistance>100:
                        d(segSpot).gkabar_kap = 0
                        d(segSpot).gkabar_kad = self.gka_base + increaseAType
                        
    def defineTargetROIs(self,cellID):
        sectionList = []
        segmentList = []
        silentID = []
        if cellID==0:
            sectionList.append(h.D11111)
            segmentList.append(0.6)
            silentID.append(True)
            
            sectionList.append(h.D112111)
            segmentList.append(0.59)
            silentID.append(False)
            
            sectionList.append(h.D112111)
            segmentList.append(0.78)
            silentID.append(False)
            
        elif cellID==1:
            sectionList.append(h.D1111111)
            segmentList.append(0.95)
            silentID.append(True)

            sectionList.append(h.D11121)
            segmentList.append(0.3)
            silentID.append(False)
            
            sectionList.append(h.D11211)
            segmentList.append(0.73)
            silentID.append(False)
            
            sectionList.append(h.D11221)
            segmentList.append(0.84)
            silentID.append(False)
            
        elif cellID==2:
            sectionList.append(h.D11111)
            segmentList.append(0.99)
            silentID.append(True)

            # The following ROI is a lot farther than the two active ROIs below...
            #sectionList.append(h.D111111)
            #segmentList.append(0.6)
            #silentID.append(True)
            
            sectionList.append(h.D121)
            segmentList.append(0.8)
            silentID.append(False)

            sectionList.append(h.D2111)
            segmentList.append(0.86)
            silentID.append(False)
            
        elif cellID==3:
            sectionList.append(h.D1111111)
            segmentList.append(0.45)
            silentID.append(True)

            # The following ROI is a lot farther than the active ROI...
            #sectionList.append(h.D11111112)
            #segmentList.append(0.2)
            #silentID.append(True)
            
            sectionList.append(h.D11122)
            segmentList.append(0.95)
            silentID.append(False)
            
        elif cellID==4:
            sectionList.append(h.D111111)
            segmentList.append(0.5)
            silentID.append(True)

            sectionList.append(h.D1121112)
            segmentList.append(0.58)
            silentID.append(False)
            
        elif cellID==5:
            sectionList.append(h.D11111)
            segmentList.append(0.85)
            silentID.append(True)

            sectionList.append(h.D1121)
            segmentList.append(0.73)
            silentID.append(False)
            
        elif cellID==6:
            sectionList.append(h.D11111)
            segmentList.append(0.26)
            silentID.append(True)

            sectionList.append(h.D111211)
            segmentList.append(0.92)
            silentID.append(False)
            
        elif cellID==7:
            sectionList.append(h.D11111)
            segmentList.append(0.93)
            silentID.append(True)

            sectionList.append(h.D112)
            segmentList.append(0.97)
            silentID.append(False)
            
            
        self.sectionList = sectionList
        self.segmentList = segmentList
        self.silentID = silentID
        
        
def returnSecListL23():
    secList = []
    for sec in h.allsec():
        if sec!=h.soma:
            secList.append(sec)
    return secList
        
        
        
    
    
    
    
    
    
    
    
    
