
COMMENT
T-type Ca channel 
ca.mod to lead to thalamic ca current inspired by destexhe and huguenrd
Uses fixed eca instead of GHK eqn
changed from (AS Oct0899)
changed for use with Ri18  (B.Kampa 2005)
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX CaT
	USEION ca READ eca WRITE ica
	RANGE m, h, gca, gbar
	RANGE minf, hinf, mtau, htau, inactF, actF
	GLOBAL  vshift,vmin,vmax, v12m, v12h, vwm, vwh, am, ah, vm1, vm2, vh1, vh2, wm1, wm2, wh1, wh2
}

PARAMETER {
	gbar = 0.0008 (mho/cm2)	: 0.12 mho/cm2
	vshift = 0	(mV)		: voltage shift (affects all)

	cao  = 2.5	(mM)	        : external ca concentration
	cai		(mM)
						 
	v 		(mV)
	dt		(ms)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)

	v12m=50         	(mV)
	v12h=78         	(mV)
	vwm =7.4         	(mV)
	vwh=5.0         	(mV)
	am=3         	(mV)
	ah=85         	(mV)
	vm1=25         	(mV)
	vm2=100         	(mV)
	vh1=46         	(mV)
	vh2=405         	(mV)
	wm1=20         	(mV)
	wm2=15         	(mV)
	wh1=4         	(mV)
	wh2=50         	(mV)


}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
	PI	= (pi) (1)
} 

ASSIGNED {
	ica 		(mA/cm2)
	gca		(pS/um2)
	eca		(mV)
	minf 		hinf
	mtau (ms)	htau (ms)
	tadj
}
 

STATE { m h }

INITIAL { 
	trates(v+vshift)
	m = minf
	h = hinf
}

BREAKPOINT {
        SOLVE states
        gca = gbar*m*m*h
	ica = gca * (v - eca)
} 

LOCAL mexp, hexp

PROCEDURE states() {
        trates(v+vshift)      
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
	VERBATIM
	return 0;
	ENDVERBATIM
}


PROCEDURE trates(v) {  
                      
        LOCAL tinc
        TABLE minf, mexp, hinf, hexp
	DEPEND dt	
	FROM vmin TO vmax WITH 199

	rates(v): not consistently executed from here if usetable == 1

        tinc = -dt 

        mexp = 1 - exp(tinc/mtau)
        hexp = 1 - exp(tinc/htau)
}


PROCEDURE rates(v_) {  
        LOCAL  a, b

	minf = 1.0 / ( 1 + exp(-(v_+v12m)/vwm) )
	hinf = 1.0 / ( 1 + exp((v_+v12h)/vwh) )

	mtau = ( am + 1.0 / ( exp((v_+vm1)/wm1) + exp(-(v_+vm2)/wm2) ) ) 
	htau = ( ah + 1.0 / ( exp((v_+vh1)/wh1) + exp(-(v_+vh2)/wh2) ) ) 
}

