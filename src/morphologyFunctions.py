import numpy as np
from neuron import h
from . import neuronFunctions as nfx

def findSites(soma, dist, method='hoc', dends=None, incDiam=False):
    """
    Take top level section (usually soma) and find all sections connected to it at a certain dist from soma.
    If method = 'hoc', uses soma.subtree(), and finds all hoc sections connected to soma
    If method = 'struct', requires 4th input dends, which is a list of connections to soma (that must be in soma.subtree())
        - just look for sections in dends that are connected at a certain distance
    Return:
        list of section names
        segment value that's the requested distance (i.e. if section is 90-110µm and requested distance is 100µm, return 0.5)
        actual distance (this is based off of how many segments, will be within smallest dx of segment, set elsewhere)
    """
    # Get subtree
    if method=='hoc':
        # Use generic hoc method
        tree = soma.subtree()
    elif method=='struct':
        # Use structure method (requires 4th input)
        tree = dends
    else:
        print('Did not recognize method')
        return None
    
    # Find sections with requested distance from soma, then return segment number for specific distance
    N = len(tree)
    isDistance = np.zeros(N,dtype='bool')
    nSegment = np.zeros(N)
    prevDistance = np.zeros(N) # length of dendrite from previous branch point
    postDistance = np.zeros(N) # length of dendrite after site requested
    distFromBranch = np.zeros(N) # distance of site from previous branch point
    diamAtBranch = np.zeros(N) # diameter at previous branch point
    for n in range(N): 
        proxDist = h.distance(tree[n](0),soma(1))
        L = tree[n].L
        if (proxDist < dist) & (proxDist+L >= dist):
            isDistance[n]=True # record that this one is a valid distance from the soma
            nSegment[n] = (dist-proxDist)/L # record segment at requested distance 
            
            distFromBranch[n] = h.distance(tree[n](nSegment[n]),tree[n](0)) # distance from previous branch point
            diamAtBranch[n] = tree[n].diam
            
            # Measure total dendritic length after site
            currentTree = tree[n].subtree()
            for ct in currentTree:
                if ct == tree[n]:
                    diamAdjustment = 1
                    if incDiam: diamAdjustment = np.pi*ct.diam
                    postDistance[n] += diamAdjustment * (L - h.distance(tree[n](nSegment[n]),tree[n](0))) # add only the distance after the requested site
                else:
                    diamAdjustment = 1
                    if incDiam: diamAdjustment = np.pi*ct.diam
                    postDistance[n] += ct.L * diamAdjustment # add all children distance
                    
            # Measure total dendritic length after previous branch point (only works rn because the parents are always parents of sisters)
            parentTree = tree[n].parentseg().sec.subtree()
            for ct in parentTree:
                if ct!=tree[n].parentseg().sec: # don't include the actual parent section (which is included in the subtree)
                    diamAdjustment = 1
                    if incDiam: diamAdjustment = np.pi*ct.diam
                    prevDistance[n] += ct.L * diamAdjustment

    outSection = [tree[sec] for sec in np.where(isDistance)[0]]
    outSegment = nSegment[isDistance]
    outPost = postDistance[isDistance]
    outPre = prevDistance[isDistance]
    outDistBranch = distFromBranch[isDistance]
    outDiam = diamAtBranch[isDistance]
    return outSection, outSegment, outPost, outPre, outDistBranch, outDiam

def measurePrePostDistance(section,segment):
    N = len(section)
    prevDistance = np.zeros(N)
    postDistance = np.zeros(N)
    for n in range(N):
        cSection = section[n]
        cSegment = segment[n]
        cTree = cSection.subtree()
        # Measure dendritic length after site
        for ct in cTree:
            postDistance[n] += ct.L
            if ct == cSection:
                postDistance[n] -= (ct.L - h.distance(ct(cSegment),ct(0))) # discount everything proximal to ROI
        
        # Measure dendritic length after previous branch point
        parentTree = cSection.parentseg().sec.subtree()
        for ct in parentTree:
            if ct!=cSection.parentseg().sec: # don't include parent section (which is included in the subtree)
                prevDistance[n] += ct.L
    return postDistance, prevDistance

def recordSites(section,segment,recordVariable='_ref_v'):
    """
    Takes a list of section names and segments within each section, and sets up recording vectors for each
    Section & segment must be registered with one another...
    recordVariable determines what to measure from each, (default is membrane voltage)
    """
    tv = h.Vector() # Time stamp vector
    tv.record(h._ref_t)
    
    vsection = [] # list of hoc vectors for each section
    for sec,seg in zip(section, segment):
        vsection.append(h.Vector())
        vsection[-1].record(getattr(sec(seg),recordVariable)) # record voltage... eventually make this a dynamic attribute name
    return vsection,tv

def injectSites(section,segment,stim=None,amplitude=-0.1):
    # Always record time stamps
    tv = h.Vector()
    tv.record(h._ref_t)
    
    # Inject and record voltage in each segment
    N = len(section)
    i = 0
    vrecord = h.Vector()
    vsection = []
    for sec,seg in zip(section,segment):
        i += 1
        #print('Working on section {0}, {1}/{2}'.format(sec,i,N))
        vrecord.record(sec(seg)._ref_v)
        stim = nfx.attachCC(section=sec, delay=50, dur=50, amp=amplitude, loc=seg)
        nfx.simulate(tstop=101,v_init=-76,celsius=37)        
        vsection.append(np.array(vrecord))
        
    return vsection,tv,stim

def injectAlphaSites(section,segment,syn=None,onset=5,tau=2,gmax=0.1,tstop=25):
    # Always record time stamps
    tv = h.Vector()
    tv.record(h._ref_t)
    
    # Inject and record voltage in each segment
    N = len(section)
    i = 0
    vrecord = h.Vector()
    vsomaRec = h.Vector()
    vsection = []
    vsoma = []
    for sec,seg in zip(section,segment):
        i += 1
        #print('Working on section {0}, {1}/{2}'.format(sec,i,N))
        vrecord.record(sec(seg)._ref_v)
        vsomaRec.record(h.soma(0.5)._ref_v)
        syn = nfx.attachAlpha(section=sec, seg=seg, onset=onset, tau=tau, gmax=gmax)
        nfx.simulate(tstop=tstop,v_init=-76,celsius=35)
        vsection.append(np.array(vrecord))
        vsoma.append(np.array(vsomaRec))
        
    return vsection,vsoma,tv,syn

def recordBranchPointDivision(section):
    tv = h.Vector()
    tv.record(h._ref_t)
    
    targetBranch = []
    sisterBranch = []
    for sec in section:
        targetBranch.append([h.Vector(), h.Vector(), 0.0])
        sisterBranch.append([h.Vector(), h.Vector(), 0.0])
        parentRef = h.SectionRef(sec=sec)
        
        # Record at first and second segment of target branch immediately after previous branch point
        targetBranch[-1][0].record(sec(nfx.returnSegment(sec.nseg,1))._ref_v)
        targetBranch[-1][1].record(sec(nfx.returnSegment(sec.nseg,2))._ref_v)
        # Resistance!
        targetAxialResistance = 4*sec.Ra*1e4 / (np.pi * sec.diam**2)
        targetLength = sec.L/sec.nseg
        targetBranch[-1][2] = 1e-6 * targetAxialResistance * targetLength 
        
        # Record at 1st/2nd segment of sister branch
        sref = h.SectionRef(sec=sec.parentseg().sec)
        if sref.nchild()!=2: 
            print('The parent of {0} had more than 2 children! Exiting prematurely.'.format(sec))
            return
        childIdx = 0
        if sref.child[childIdx]==sec: 
            childIdx=1
        sisterSection = sref.child[childIdx]
        sisterBranch[-1][0].record(sisterSection(nfx.returnSegment(sisterSection.nseg,1))._ref_v)
        sisterBranch[-1][1].record(sisterSection(nfx.returnSegment(sisterSection.nseg,2))._ref_v)
        # Resistance!
        sisterAxialResistance = 4*sisterSection.Ra*1e4 / (np.pi * sisterSection.diam**2)
        sisterLength = sisterSection.L / sisterSection.nseg
        sisterBranch[-1][2] = 1e-6 * sisterAxialResistance * sisterLength 
    
    return tv,targetBranch,sisterBranch

def measureConvolvedBranching(section, segment, lengthConstant):
    convolvedLength = []
    for sec, seg in zip(section, segment):
        thisSecOffsets = [] # keep track of distance to branch points
        thisSecLengths = [] # keep track of distance after branch points

        # Start by measuring distance after ROI itself
        thisSecOffsets.append(0.0) # no distance between ROI and itself
        currentTree = sec.subtree()
        currentPost = 0# -h.distance(sec(seg),sec(0)) # start by subtracting distance from previous branch point to ROI (which will be included and offset in following loop)
        for ct in currentTree:
            currentPost += ct.L
            if ct==sec: currentPost -= h.distance(sec(seg),sec(0))
        thisSecLengths.append(currentPost)

        # Next, measure distance after each previous branch point (including soma)
        currentSec = sec
        while True:
            dist2branch = h.distance(sec(seg),currentSec(0)) # measure distance from ROI
            thisSecOffsets.append(dist2branch)
            currentTree = currentSec.subtree()
            currentPost = 0.0
            for ct in currentTree:
                currentPost += ct.L
            thisSecLengths.append(currentPost)
            currentSec = currentSec.parentseg().sec
            if currentSec.parentseg() is None: break
        
        # Now, compute weighted average using exponential decay as convolutional filter
        expPoints = np.exp(-(np.array(thisSecOffsets))/lengthConstant)
        convolvedLength.append(np.dot(thisSecLengths,expPoints)/np.sum(expPoints))
    return convolvedLength
        
      
def measureLocalBranching(section,segment,soma,lengthConstant):
    localBranching = []
    somaTree = soma.subtree()
    for sec,seg in zip(section,segment):
        currentDistance = 0.0
        for ct in somaTree:
            if ct==sec:
                distalLength = h.distance(sec(seg),sec(1))
                multFactor = lengthConstant/distalLength * (np.exp(-0/lengthConstant) - np.exp(-distalLength/lengthConstant))
                currentDistance += distalLength * multFactor
                
                proxLength = h.distance(sec(seg),sec(0))
                multFactor = lengthConstant/proxLength * (np.exp(-0/lengthConstant) - np.exp(-proxLength/lengthConstant))
                currentDistance += proxLength * multFactor
            elif ct!=soma:
                idx1 = h.distance(sec(seg),ct(1))
                idx2 = h.distance(sec(seg),ct(0))
                further = np.max([idx1,idx2])
                closer = np.min([idx1,idx2])
                cLength = further - closer
                multFactor = lengthConstant/(further - closer) * (np.exp(-closer/lengthConstant) - np.exp(-further/lengthConstant))
                currentDistance += cLength * multFactor
            else:
                None
                # Don't add soma...
        localBranching.append(currentDistance)
    return localBranching

    
def measureDiscountedMorphRatio(section,lengthConstant,method='exponential'):
    discountedLength = []
    for sec in section:
        discountedLength.append([0.0,0.0])
        parentSec = sec.parentseg().sec
        parentSecRef = h.SectionRef(sec=parentSec) # Get parent section reference
        if parentSecRef.nchild()!=2: 
            print('The parent of {0} had more than 2 children! Exiting prematurely.'.format(sec)) 
            return
        childIdx = 0 # Try idx 0
        if parentSecRef.child[childIdx]==sec: 
            childIdx=1 # Set idx to sister branch
        sisSec = parentSecRef.child[childIdx]
        
        # Measure total dendritic length, discounted exponentially
        targetTree = sec.subtree()
        for tt in targetTree:
            idx1 = h.distance(parentSec(1),tt(0))
            idx2 = h.distance(parentSec(1),tt(1))
            if method=='exponential':
                # mult factor is integral of exponential between start and end distance of the current section
                multFactor = lengthConstant/(idx2-idx1)*(np.exp(-idx1/lengthConstant) - np.exp(-idx2/lengthConstant)) # average of exponential evaluated between idx1 & idx2
            elif method=='linear':
                multFactor = 1/(idx2-idx1)*1/(2*lengthConstant) * (idx2**2 - idx1**2)
            elif method=='sigmoid':
                if len(lengthConstant)!=3: 
                    print('For sigmoid, must provide 3 terms, exiting now.')
                    return
                mainTerm = lambda x: lengthConstant[0] / (lengthConstant[0] + np.exp(-lengthConstant[1]*x - lengthConstant[2]))
                multFactor = 1/(idx2-idx1) * (1/lengthConstant[1]) * (np.log(mainTerm(idx2)) - np.log(mainTerm(idx1)))
            elif method=='order':
                if len(lengthConstant)!=3: 
                    print('For order, must provide 3 terms, exiting now.')
                    return
                order = 0
                currSec = tt
                while currSec!=sec:
                    order += 1
                    currSec = currSec.parentseg().sec
                multFactor = 1 - 1/(1+lengthConstant[0]*np.exp(lengthConstant[1]*order+lengthConstant[2]))
            discountedLength[-1][0] += tt.L * multFactor
        sisterTree = sisSec.subtree()
        for st in sisterTree:
            idx1 = h.distance(parentSec(1),st(0))
            idx2 = h.distance(parentSec(1),st(1))
            if method=='exponential':
                # mult factor is integral of exponential between start and end distance of the current section
                multFactor = lengthConstant/(idx2-idx1)*(np.exp(-idx1/lengthConstant) - np.exp(-idx2/lengthConstant)) # average of exponential evaluated between idx1 & idx2
            elif method=='linear':
                multFactor = 1/(idx2-idx1)*1/(2*lengthConstant) * (idx2**2 - idx1**2)
            elif method=='sigmoid':
                if len(lengthConstant)!=3: 
                    print('For sigmoid, must provide 3 terms, exiting now.')
                    return
                mainTerm = lambda x: lengthConstant[0] / (lengthConstant[0] + np.exp(-lengthConstant[1]*x - lengthConstant[2]))
                multFactor = 1/(idx2-idx1) * (1/lengthConstant[1]) * (np.log(mainTerm(idx2)) - np.log(mainTerm(idx1)))
            elif method=='order':
                if len(lengthConstant)!=3: 
                    print('For order, must provide 3 terms, exiting now.')
                    return
                order = 0
                currSec = st
                while currSec!=sisSec:
                    order += 1
                    currSec = currSec.parentseg().sec
                multFactor = 1 - 1/(1+lengthConstant[0]*np.exp(lengthConstant[1]*order+lengthConstant[2]))
            discountedLength[-1][1] += st.L * multFactor
    return discountedLength
