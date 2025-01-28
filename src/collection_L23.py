import numpy as np
from neuron import h, load_mechanisms
from . import neuronFunctions as nfx
from . import get_src_dir
src_directory = get_src_dir()

load_mechanisms(f'{src_directory}/mod.files/')
h('objref nil')


class L23(object):
    def __init__(self,cellID=0,dendNa=[6,2,200,False],dendK=[0.048,100,False,False],dxSeg=5,fixDiam=None):
        """
        Initializer for L23 object. cellID determines which cell to load (see "loadHoc" function of this class for file names)
        """
        if cellID>5 or cellID<0:
            raise ValueError(f"CellIDs can only be 0-5, cellID={cellID} was provided.")
        self.cellID = cellID
        self.loadHoc(cellID) # Load Hoc File
        self.props(dendNa=dendNa,dendK=dendK) # Load properties into class
        self._topol(cellID,dxSeg,fixDiam) # Construct topology (and load names of sections into L23 object)
        if cellID==0: self._changeLength(cellID) # if Cat Cell, fix lengths
        self._addAxon() # Add axon
        self._biophys() # Add passive properties
        self.init_active(dendNa=dendNa, dendK=dendK)
        
    def loadHoc(self,cellID):
        hocFileName = {0:'L23',1:'29_CDK170205_registered_D2',2:'layer23_1_ATL',3:'layer23_2_ATL',4:'layer23_3_ATL',5:'layer23_4_ATL'}
        print('Creating cell : {fName}'.format(fName=hocFileName[cellID]))
        file_name = f"{get_src_dir()}/KeyNeuronMorphologies/{hocFileName[cellID]}.hoc"
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
        secListL23 = returnSecListL23(cellID) #return name of every dendritic segment in list (excludes soma and axon)
        for secName in secListL23:
            self.dends.append(eval('h.{secName}'.format(secName=secName)))
        
        # Set the diameter of each section correctly 
        for sec in self.dends:
            sec.diam = sec.diam * 0.7
            sec.nseg = np.ceil(sec.L/dxSeg).astype(int) # have a little bit more than 1 segment per every dxSeg(default=10µm)
            if fixDiam:
                sec.diam = fixDiam
                
    def _changeLength(self,cellID):
        if cellID != 0:
            return
            
        # Update length of all sections because L23.hoc based on cat neuron
        self.soma.L = self.soma.L * 0.7
        for sec in self.dends:
            sec.L = sec.L * 0.7
            
    def _addAxon(self):
        
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
                        

def returnElevAxim(cellID):
    if cellID==0:
        return None
    elif cellID==1:
        return None
        
        
def returnSecListL23(cellID):
    secList = []
    if cellID==0:
        # Cat cell from Hausser Lab
        secList.append("dend1_1")
        secList.append("dend1_11")
        secList.append("dend1_111")
        secList.append("dend1_1111")
        secList.append("dend1_1112")
        secList.append("dend1_112")
        secList.append("dend1_1121")
        secList.append("dend1_1122")
        secList.append("dend1_12")
        secList.append("dend1_121")
        secList.append("dend1_1211")
        secList.append("dend1_1212")
        secList.append("dend1_122")
        secList.append("dend1_1221")
        secList.append("dend1_1222")
        secList.append("dend1_12221")            
        secList.append("dend1_12222")            
        secList.append("dend2_1")                
        secList.append("dend2_11")               
        secList.append("dend2_111")              
        secList.append("dend2_1111")             
        secList.append("dend2_1112")             
        secList.append("dend2_112")              
        secList.append("dend2_1121")             
        secList.append("dend2_1122")             
        secList.append("dend2_12")               
        secList.append("dend2_121")              
        secList.append("dend2_1211")             
        secList.append("dend2_12111")            
        secList.append("dend2_121111")           
        secList.append("dend2_121112")           
        secList.append("dend2_12112")            
        secList.append("dend2_121121")           
        secList.append("dend2_1211211")        
        secList.append("dend2_1211212")       
        secList.append("dend2_12112121")       
        secList.append("dend2_12112122")
        secList.append("dend2_121122")  
        secList.append("dend2_1212") 
        secList.append("dend2_12121")
        secList.append("dend2_121211")
        secList.append("dend2_121212") 
        secList.append("dend2_12122")   
        secList.append("dend2_122")  
        secList.append("dend2_1221") 
        secList.append("dend2_12211") 
        secList.append("dend2_12212")  
        secList.append("dend2_1222") 
        secList.append("dend2_12221") 
        secList.append("dend2_12222")
        secList.append("dend3_1")    
        secList.append("dend3_11")   
        secList.append("dend3_111")  
        secList.append("dend3_1111")  
        secList.append("dend3_1112") 
        secList.append("dend3_11121") 
        secList.append("dend3_11122")   
        secList.append("dend3_112")  
        secList.append("dend3_1121")  
        secList.append("dend3_1122")    
        secList.append("dend3_12")   
        secList.append("dend3_121")  
        secList.append("dend3_1211")  
        secList.append("dend3_1212") 
        secList.append("dend3_12121") 
        secList.append("dend3_12122")
        secList.append("dend3_121221")
        secList.append("dend3_121222")
        secList.append("dend3_1212221")
        secList.append("dend3_1212222")
        secList.append("dend3_12122221")
        secList.append("dend3_12122222")
        secList.append("dend3_121222221")
        secList.append("dend3_1212222211")
        secList.append("dend3_12122222111")
        secList.append("dend3_12122222112")
        secList.append("dend3_1212222212")
        secList.append("dend3_121222222")
        secList.append("dend3_1212222221")
        secList.append("dend3_1212222222")
        secList.append("dend3_12122222221")
        secList.append("dend3_12122222222")
        secList.append("dend3_122")
        secList.append("dend3_1221")
        secList.append("dend3_12211")
        secList.append("dend3_12212")
        secList.append("dend3_1222")
        secList.append("dend4_1")
        secList.append("dend4_11")
        secList.append("dend4_111")
        secList.append("dend4_1111")
        secList.append("dend4_1112")
        secList.append("dend4_11121")
        secList.append("dend4_11122")
        secList.append("dend4_112")
        secList.append("dend4_1121")
        secList.append("dend4_1122")
        secList.append("dend4_12")
        secList.append("dend4_121")
        secList.append("dend4_122")
        secList.append("dend4_1221")
        secList.append("dend4_1222")
        secList.append("dend4_12221")
        secList.append("dend4_12222")
        
    elif cellID==1:
        # EggerSchmitt Cell
        secList.append("apical_1_0")
        secList.append("apical_1_0_0")
        secList.append("apical_1_0_0_0")
        secList.append("apical_1_0_0_0_0")
        secList.append("apical_1_0_0_0_1")
        secList.append("apical_1_0_0_1")
        secList.append("apical_1_0_0_1_0")
        secList.append("apical_1_0_0_1_1")
        secList.append("apical_1_0_0_1_1_0")
        secList.append("apical_1_0_0_1_1_1")
        secList.append("apical_1_0_1")
        secList.append("apical_1_0_1_0")
        secList.append("apical_1_0_1_0_0")
        secList.append("apical_1_0_1_0_1")
        secList.append("apical_1_0_1_1")
        secList.append("apical_1_0_1_1_0")
        secList.append("apical_1_0_1_1_0_0")
        secList.append("apical_1_0_1_1_0_0_0")
        secList.append("apical_1_0_1_1_0_0_1")
        secList.append("apical_1_0_1_1_0_0_1_0")
        secList.append("apical_1_0_1_1_0_0_1_1")
        secList.append("apical_1_0_1_1_0_0_1_1_0")
        secList.append("apical_1_0_1_1_0_0_1_1_0_0")
        secList.append("apical_1_0_1_1_0_0_1_1_0_1")
        secList.append("apical_1_0_1_1_0_0_1_1_1")
        secList.append("apical_1_0_1_1_0_0_1_1_1_0")
        secList.append("apical_1_0_1_1_0_0_1_1_1_0_0")
        secList.append("apical_1_0_1_1_0_0_1_1_1_0_0_0")
        secList.append("apical_1_0_1_1_0_0_1_1_1_0_0_1")
        secList.append("apical_1_0_1_1_0_0_1_1_1_0_1")
        secList.append("apical_1_0_1_1_0_0_1_1_1_0_1_0")
        secList.append("apical_1_0_1_1_0_0_1_1_1_0_1_1")
        secList.append("apical_1_0_1_1_0_0_1_1_1_0_1_1_0")
        secList.append("apical_1_0_1_1_0_0_1_1_1_0_1_1_1")
        secList.append("apical_1_0_1_1_0_0_1_1_1_1")
        secList.append("apical_1_0_1_1_0_0_1_1_1_1_0")
        secList.append("apical_1_0_1_1_0_0_1_1_1_1_0_0")
        secList.append("apical_1_0_1_1_0_0_1_1_1_1_0_1")
        secList.append("apical_1_0_1_1_0_0_1_1_1_1_1")
        secList.append("apical_1_0_1_1_0_1")
        secList.append("apical_1_0_1_1_0_1_0")
        secList.append("apical_1_0_1_1_0_1_1")
        secList.append("apical_1_0_1_1_1")
        secList.append("apical_1_0_1_1_1_0")
        secList.append("apical_1_0_1_1_1_1")
        secList.append("dend_1_0")
        secList.append("dend_1_0_0")
        secList.append("dend_1_0_1")
        secList.append("dend_1_0_1_0")
        secList.append("dend_1_0_1_0_0")
        secList.append("dend_1_0_1_0_1")
        secList.append("dend_1_0_1_1")
        secList.append("dend_1_0_1_1_0")
        secList.append("dend_1_0_1_1_1")
        secList.append("dend_2_0")
        secList.append("dend_2_0_0")
        secList.append("dend_2_0_0_0")
        secList.append("dend_2_0_0_1")
        secList.append("dend_2_0_1")
        secList.append("dend_3_0")
        secList.append("dend_4_0")
        secList.append("dend_4_0_0")
        secList.append("dend_4_0_0_0")
        secList.append("dend_4_0_0_0_0")
        secList.append("dend_4_0_0_0_0_0")
        secList.append("dend_4_0_0_0_0_0_0")
        secList.append("dend_4_0_0_0_0_0_0_0")
        secList.append("dend_4_0_0_0_0_0_0_1")
        secList.append("dend_4_0_0_0_0_0_1")
        secList.append("dend_4_0_0_0_0_1")
        secList.append("dend_4_0_0_0_0_1_0")
        secList.append("dend_4_0_0_0_0_1_1")
        secList.append("dend_4_0_0_0_1")
        secList.append("dend_4_0_0_0_1_0")
        secList.append("dend_4_0_0_0_1_1")
        secList.append("dend_4_0_0_1")
        secList.append("dend_4_0_0_1_0")
        secList.append("dend_4_0_0_1_1")
        secList.append("dend_4_0_1")
        secList.append("dend_4_0_1_0")
        secList.append("dend_4_0_1_0_0")
        secList.append("dend_4_0_1_0_1")
        secList.append("dend_4_0_1_0_1_0")
        secList.append("dend_4_0_1_0_1_0_0")
        secList.append("dend_4_0_1_0_1_0_0_0")
        secList.append("dend_4_0_1_0_1_0_0_1")
        secList.append("dend_4_0_1_0_1_0_1")
        secList.append("dend_4_0_1_0_1_1")
        secList.append("dend_4_0_1_1")
        secList.append("dend_4_0_1_1_0")
        secList.append("dend_4_0_1_1_0_0")
        secList.append("dend_4_0_1_1_0_1")
        secList.append("dend_4_0_1_1_0_1_0")
        secList.append("dend_4_0_1_1_0_1_0_0")
        secList.append("dend_4_0_1_1_0_1_0_1")
        secList.append("dend_4_0_1_1_0_1_1")
        secList.append("dend_4_0_1_1_1")
        secList.append("dend_4_0_1_1_1_0")
        secList.append("dend_4_0_1_1_1_1")
        secList.append("dend_5_0")
        secList.append("dend_5_0_0")
        secList.append("dend_5_0_0_0")
        secList.append("dend_5_0_0_0_0")
        secList.append("dend_5_0_0_0_1")
        secList.append("dend_5_0_0_1")
        secList.append("dend_5_0_1")
        secList.append("dend_5_0_1_0")
        secList.append("dend_5_0_1_1")
        secList.append("dend_5_0_1_1_0")
        secList.append("dend_5_0_1_1_0_0")
        secList.append("dend_5_0_1_1_0_1")
        secList.append("dend_5_0_1_1_1")
        secList.append("dend_5_0_1_1_1_0")
        secList.append("dend_5_0_1_1_1_1")
        secList.append("dend_6_0")
        secList.append("dend_6_0_0")
        secList.append("dend_6_0_0_0")
        secList.append("dend_6_0_0_1")
        secList.append("dend_6_0_1")
        secList.append("dend_6_0_1_0")
        secList.append("dend_6_0_1_1")
        secList.append("dend_7_0")
        secList.append("dend_7_0_0")
        secList.append("dend_7_0_1")
        secList.append("dend_8_0")
        secList.append("dend_8_0_0")
        secList.append("dend_8_0_0_0")
        secList.append("dend_8_0_0_1")
        secList.append("dend_8_0_1")
        secList.append("dend_9_0")
        secList.append("dend_9_0_0")
        secList.append("dend_9_0_1")
        
    elif cellID==2:
        # Poleg-Polsky L23_1
        secList.append("apic[0]")
        secList.append("apic[1]")
        secList.append("apic[2]")
        secList.append("apic[3]")
        secList.append("apic[4]")
        secList.append("apic[5]")
        secList.append("apic[6]")
        secList.append("apic[7]")
        secList.append("apic[8]")
        secList.append("apic[9]")
        secList.append("apic[10]")
        secList.append("apic[11]")
        secList.append("apic[12]")
        secList.append("apic[13]")
        secList.append("apic[14]")
        secList.append("apic[15]")
        secList.append("apic[16]")
        secList.append("apic[17]")
        secList.append("apic[18]")
        secList.append("apic[19]")
        secList.append("apic[20]")
        secList.append("apic[21]")
        secList.append("apic[22]")
        secList.append("apic[23]")
        secList.append("apic[24]")
        secList.append("apic[25]")
        secList.append("apic[26]")
        secList.append("apic[27]")
        secList.append("apic[28]")
        secList.append("apic[29]")
        secList.append("apic[30]")
        secList.append("apic[31]")
        secList.append("apic[32]")
        secList.append("apic[33]")
        secList.append("apic[34]")
        secList.append("apic[35]")
        secList.append("apic[36]")
        secList.append("apic[37]")
        secList.append("apic[38]")
        secList.append("dend[0]")
        secList.append("dend[1]")
        secList.append("dend[2]")
        secList.append("dend[3]")
        secList.append("dend[4]")
        secList.append("dend[5]")
        secList.append("dend[6]")
        secList.append("dend[7]")
        secList.append("dend[8]")
        secList.append("dend[9]")
        secList.append("dend[10]")
        secList.append("dend[11]")
        secList.append("dend[12]")
        secList.append("dend[13]")
        secList.append("dend[14]")
        secList.append("dend[15]")
        secList.append("dend[16]")
        secList.append("dend[17]")
        secList.append("dend[18]")
        secList.append("dend[19]")
        secList.append("dend[20]")
        secList.append("dend[21]")
        secList.append("dend[22]")
        secList.append("dend[23]")
        secList.append("dend[24]")
        secList.append("dend[25]")
        secList.append("dend[26]")
        secList.append("dend[27]")
        secList.append("dend[28]")
    elif cellID==3:
        # Poleg-Polsky L23_2
        secList.append("apic[0]")
        secList.append("apic[1]")
        secList.append("apic[2]")
        secList.append("apic[3]")
        secList.append("apic[4]")
        secList.append("apic[5]")
        secList.append("apic[6]")
        secList.append("apic[7]")
        secList.append("apic[8]")
        secList.append("apic[9]")
        secList.append("apic[10]")
        secList.append("apic[11]")
        secList.append("apic[12]")
        secList.append("apic[13]")
        secList.append("apic[14]")
        secList.append("apic[15]")
        secList.append("apic[16]")
        secList.append("apic[17]")
        secList.append("apic[18]")
        secList.append("apic[19]")
        secList.append("apic[20]")
        secList.append("apic[21]")
        secList.append("apic[22]")
        secList.append("apic[23]")
        secList.append("apic[24]")
        secList.append("apic[25]")
        secList.append("apic[26]")
        secList.append("apic[27]")
        secList.append("apic[28]")
        secList.append("dend[0]")
        secList.append("dend[1]")
        secList.append("dend[2]")
        secList.append("dend[3]")
        secList.append("dend[4]")
        secList.append("dend[5]")
        secList.append("dend[6]")
        secList.append("dend[7]")
        secList.append("dend[8]")
        secList.append("dend[9]")
        secList.append("dend[10]")
        secList.append("dend[11]")
        secList.append("dend[12]")
        secList.append("dend[13]")
        secList.append("dend[14]")
        
    elif cellID==4:
        # Poleg-Polsky L23_3
        secList.append("dend[0]")
        secList.append("dend[1]")
        secList.append("dend[2]")
        secList.append("dend[3]")
        secList.append("dend[4]")
        secList.append("dend[5]")
        secList.append("dend[6]")
        secList.append("dend[7]")
        secList.append("dend[8]")
        secList.append("dend[9]")
        secList.append("dend[10]")
        secList.append("dend[11]")
        secList.append("dend[12]")
        secList.append("dend[13]")
        secList.append("dend[14]")
        secList.append("dend[15]")
        secList.append("dend[16]")
        secList.append("dend[17]")
        secList.append("dend[18]")
        secList.append("dend[19]")
        secList.append("dend[20]")
        secList.append("dend[21]")
        secList.append("dend[22]")
        secList.append("dend[23]")
        secList.append("dend[24]")
        secList.append("dend[25]")
        secList.append("dend[26]")
        secList.append("dend[27]")
        secList.append("dend[28]")
        secList.append("dend[29]")
        secList.append("dend[30]")
        secList.append("dend[31]")
        secList.append("dend[32]")
        secList.append("apic[0]")
        secList.append("apic[1]")
        secList.append("apic[2]")
        secList.append("apic[3]")
        secList.append("apic[4]")
        secList.append("apic[5]")
        secList.append("apic[6]")
        secList.append("apic[7]")
        secList.append("apic[8]")
        secList.append("apic[9]")
        secList.append("apic[10]")
        secList.append("apic[11]")
        secList.append("apic[12]")
        secList.append("apic[13]")
        secList.append("apic[14]")
        secList.append("apic[15]")
        secList.append("apic[16]")
        secList.append("apic[17]")
        secList.append("apic[18]")

    elif cellID==5:
        # Poleg-Polsky L23_4
        secList.append("dend[0]")
        secList.append("dend[1]")
        secList.append("dend[2]")
        secList.append("dend[3]")
        secList.append("dend[4]")
        secList.append("dend[5]")
        secList.append("dend[6]")
        secList.append("dend[7]")
        secList.append("dend[8]")
        secList.append("dend[9]")
        secList.append("dend[10]")
        secList.append("dend[11]")
        secList.append("dend[12]")
        secList.append("dend[13]")
        secList.append("dend[14]")
        secList.append("dend[15]")
        secList.append("dend[16]")
        secList.append("dend[17]")
        secList.append("dend[18]")
        secList.append("dend[19]")
        secList.append("dend[20]")
        secList.append("dend[21]")
        secList.append("dend[22]")
        secList.append("dend[23]")
        secList.append("dend[24]")
        secList.append("dend[25]")
        secList.append("dend[26]")
        secList.append("dend[27]")
        secList.append("dend[28]")
        secList.append("dend[29]")
        secList.append("dend[30]")
        secList.append("dend[31]")
        secList.append("dend[32]")
        secList.append("dend[33]")
        secList.append("dend[34]")
        secList.append("dend[35]")
        secList.append("dend[36]")
        secList.append("dend[37]")
        secList.append("dend[38]")
        secList.append("dend[39]")
        secList.append("dend[40]")
        secList.append("dend[41]")
        secList.append("dend[42]")
        secList.append("dend[43]")
        secList.append("dend[44]")
        secList.append("dend[45]")
        secList.append("dend[46]")
        secList.append("dend[47]")
        secList.append("dend[48]")
        secList.append("dend[49]")
        secList.append("dend[50]")
        secList.append("dend[51]")
        secList.append("dend[52]")
        secList.append("dend[53]")
        secList.append("dend[54]")
        secList.append("dend[55]")
        secList.append("dend[56]")
        secList.append("dend[57]")
        secList.append("dend[58]")
        secList.append("dend[59]")
        secList.append("dend[60]")
        secList.append("dend[61]")
        secList.append("dend[62]")
        secList.append("dend[63]")
        secList.append("dend[64]")
        secList.append("dend[65]")
        secList.append("dend[66]")
        secList.append("dend[67]")
        secList.append("dend[68]")
        secList.append("dend[69]")
        secList.append("dend[70]")
        secList.append("dend[71]")
        secList.append("dend[72]")
        secList.append("dend[73]")
        secList.append("dend[74]")
        secList.append("dend[75]")
        secList.append("dend[76]")
        secList.append("dend[77]")
        secList.append("dend[78]")
        secList.append("dend[79]")
        secList.append("dend[80]")
        secList.append("dend[81]")
        secList.append("dend[82]")
        secList.append("dend[83]")
        secList.append("dend[84]")
        secList.append("dend[85]")
        secList.append("dend[86]")
        secList.append("dend[87]")
        secList.append("dend[88]")
        secList.append("dend[89]")

    return secList
    
    
    
    
    
    
    
    
    
    
