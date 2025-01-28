from neuron import h
from neuron.units import ms, mV, µm

h.load_file('stdrun.hoc')

def simulate(tstop=25,v_init=-75,celsius=34):
    # Set initial voltage for all sections
    h.finitialize(v_init * mV)
    h.celsius = celsius
    h.continuerun(tstop * ms)

def record(section,loc,recVar='_ref_v'):
    """Record a section at location. Record membrane voltage unless specified with recVar."""
    tv = h.Vector()
    tv.record(h._ref_t)
    data = h.Vector()
    data.record(getattr(section(loc),recVar))
    return data,tv

def attachCC(section, delay=5, dur=1, amp=.1, loc=0.5):
    stim = h.IClamp(section(loc))
    stim.delay = delay
    stim.dur = dur
    stim.amp = amp
    return stim

def attachAlpha(section, seg, onset=5, tau=1, gmax=.1, eRev=0):
    syn = h.AlphaSynapse(section(seg))
    syn.onset = onset
    syn.tau = tau
    syn.gmax = gmax
    return syn

def returnSegment(nseg,seg):
    dx = 1.0/nseg # it'll be broken up into segments that are this wide
    dxShift = dx/2.0 # shift to center of segment
    if seg==-1: seg=nseg
    return dx*seg-dxShift

def conductanceFromCurrent(i,v,e):
    g = 1e4 * i / (v-e)
    return g

def showNeuronMorphology():
    shape_window = h.PlotShape() # without extra specification, this doesn't show diameters.... so it just looks like a line
    shape_window.exec_menu('Show Diam')

def deleteEverything():
    for sec in h.allsec():
        h.delete_section(sec=sec)
        