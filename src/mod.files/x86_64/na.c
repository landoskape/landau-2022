/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__na
#define _nrn_initial _nrn_initial__na
#define nrn_cur _nrn_cur__na
#define _nrn_current _nrn_current__na
#define nrn_jacob _nrn_jacob__na
#define nrn_state _nrn_state__na
#define _net_receive _net_receive__na 
#define _f_trates _f_trates__na 
#define rates rates__na 
#define states states__na 
#define trates trates__na 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gbar _p[0]
#define gna _p[1]
#define minf _p[2]
#define hinf _p[3]
#define mtau _p[4]
#define htau _p[5]
#define m _p[6]
#define h _p[7]
#define ina _p[8]
#define ena _p[9]
#define Dm _p[10]
#define Dh _p[11]
#define _g _p[12]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_rates(void);
 static void _hoc_trap0(void);
 static void _hoc_trates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_na", _hoc_setdata,
 "rates_na", _hoc_rates,
 "trap0_na", _hoc_trap0,
 "trates_na", _hoc_trates,
 0, 0
};
#define trap0 trap0_na
 extern double trap0( double , double , double , double );
 /* declare global and static user variables */
#define Rg Rg_na
 double Rg = 0.0091;
#define Rd Rd_na
 double Rd = 0.024;
#define Rb Rb_na
 double Rb = 0.124;
#define Ra Ra_na
 double Ra = 0.182;
#define q10 q10_na
 double q10 = 2.3;
#define qinf qinf_na
 double qinf = 6.2;
#define qi qi_na
 double qi = 5;
#define qa qa_na
 double qa = 9;
#define tadj tadj_na
 double tadj = 0;
#define temp temp_na
 double temp = 23;
#define thinf thinf_na
 double thinf = -65;
#define thi2 thi2_na
 double thi2 = -75;
#define thi1 thi1_na
 double thi1 = -50;
#define tha tha_na
 double tha = -35;
#define usetable usetable_na
 double usetable = 1;
#define vshift vshift_na
 double vshift = -10;
#define vmax vmax_na
 double vmax = 100;
#define vmin vmin_na
 double vmin = -120;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_na", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "vshift_na", "mV",
 "tha_na", "mV",
 "qa_na", "mV",
 "Ra_na", "/ms",
 "Rb_na", "/ms",
 "thi1_na", "mV",
 "thi2_na", "mV",
 "qi_na", "mV",
 "thinf_na", "mV",
 "qinf_na", "mV",
 "Rg_na", "/ms",
 "Rd_na", "/ms",
 "temp_na", "degC",
 "vmin_na", "mV",
 "vmax_na", "mV",
 "gbar_na", "pS/um2",
 "gna_na", "pS/um2",
 "mtau_na", "ms",
 "htau_na", "ms",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double m0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "vshift_na", &vshift_na,
 "tha_na", &tha_na,
 "qa_na", &qa_na,
 "Ra_na", &Ra_na,
 "Rb_na", &Rb_na,
 "thi1_na", &thi1_na,
 "thi2_na", &thi2_na,
 "qi_na", &qi_na,
 "thinf_na", &thinf_na,
 "qinf_na", &qinf_na,
 "Rg_na", &Rg_na,
 "Rd_na", &Rd_na,
 "temp_na", &temp_na,
 "q10_na", &q10_na,
 "vmin_na", &vmin_na,
 "vmax_na", &vmax_na,
 "tadj_na", &tadj_na,
 "usetable_na", &usetable_na,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"na",
 "gbar_na",
 0,
 "gna_na",
 "minf_na",
 "hinf_na",
 "mtau_na",
 "htau_na",
 0,
 "m_na",
 "h_na",
 0,
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 13, _prop);
 	/*initialize range parameters*/
 	gbar = 1000;
 	_prop->param = _p;
 	_prop->param_size = 13;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _na_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 13, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 na /Users/landauland/Dropbox/SabatiniLab/neuron-modeling/smithAdaptation/mod.files/x86_64/na.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double _zmexp , _zhexp ;
 static double *_t_minf;
 static double *_t_hinf;
 static double *_t_mtau;
 static double *_t_htau;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_trates(double);
static int rates(double);
static int trates(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_trates(double);
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   trates ( _threadargscomma_ v + vshift ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 trates ( _threadargscomma_ v + vshift ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   trates ( _threadargscomma_ v + vshift ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0 ) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0 ) ) ) / htau ) - h) ;
   }
  return 0;
}
 static double _mfac_trates, _tmin_trates;
 static void _check_trates();
 static void _check_trates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_celsius;
  static double _sav_temp;
  static double _sav_Ra;
  static double _sav_Rb;
  static double _sav_Rd;
  static double _sav_Rg;
  static double _sav_tha;
  static double _sav_thi1;
  static double _sav_thi2;
  static double _sav_qa;
  static double _sav_qi;
  static double _sav_qinf;
  if (!usetable) {return;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_sav_temp != temp) { _maktable = 1;}
  if (_sav_Ra != Ra) { _maktable = 1;}
  if (_sav_Rb != Rb) { _maktable = 1;}
  if (_sav_Rd != Rd) { _maktable = 1;}
  if (_sav_Rg != Rg) { _maktable = 1;}
  if (_sav_tha != tha) { _maktable = 1;}
  if (_sav_thi1 != thi1) { _maktable = 1;}
  if (_sav_thi2 != thi2) { _maktable = 1;}
  if (_sav_qa != qa) { _maktable = 1;}
  if (_sav_qi != qi) { _maktable = 1;}
  if (_sav_qinf != qinf) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_trates =  vmin ;
   _tmax =  vmax ;
   _dx = (_tmax - _tmin_trates)/199.; _mfac_trates = 1./_dx;
   for (_i=0, _x=_tmin_trates; _i < 200; _x += _dx, _i++) {
    _f_trates(_x);
    _t_minf[_i] = minf;
    _t_hinf[_i] = hinf;
    _t_mtau[_i] = mtau;
    _t_htau[_i] = htau;
   }
   _sav_celsius = celsius;
   _sav_temp = temp;
   _sav_Ra = Ra;
   _sav_Rb = Rb;
   _sav_Rd = Rd;
   _sav_Rg = Rg;
   _sav_tha = tha;
   _sav_thi1 = thi1;
   _sav_thi2 = thi2;
   _sav_qa = qa;
   _sav_qi = qi;
   _sav_qinf = qinf;
  }
 }

 static int trates(double _lv){ _check_trates();
 _n_trates(_lv);
 return 0;
 }

 static void _n_trates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_trates(_lv); return; 
}
 _xi = _mfac_trates * (_lv - _tmin_trates);
 if (isnan(_xi)) {
  minf = _xi;
  hinf = _xi;
  mtau = _xi;
  htau = _xi;
  return;
 }
 if (_xi <= 0.) {
 minf = _t_minf[0];
 hinf = _t_hinf[0];
 mtau = _t_mtau[0];
 htau = _t_htau[0];
 return; }
 if (_xi >= 199.) {
 minf = _t_minf[199];
 hinf = _t_hinf[199];
 mtau = _t_mtau[199];
 htau = _t_htau[199];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 minf = _t_minf[_i] + _theta*(_t_minf[_i+1] - _t_minf[_i]);
 hinf = _t_hinf[_i] + _theta*(_t_hinf[_i+1] - _t_hinf[_i]);
 mtau = _t_mtau[_i] + _theta*(_t_mtau[_i+1] - _t_mtau[_i]);
 htau = _t_htau[_i] + _theta*(_t_htau[_i+1] - _t_htau[_i]);
 }

 
static int  _f_trates (  double _lv ) {
   rates ( _threadargscomma_ _lv ) ;
    return 0; }
 
static void _hoc_trates(void) {
  double _r;
    _r = 1.;
 trates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int  rates (  double _lvm ) {
   double _la , _lb ;
 _la = trap0 ( _threadargscomma_ _lvm , tha , Ra , qa ) ;
   _lb = trap0 ( _threadargscomma_ - _lvm , - tha , Rb , qa ) ;
   tadj = pow( q10 , ( ( celsius - temp ) / 10.0 ) ) ;
   mtau = 1.0 / tadj / ( _la + _lb ) ;
   minf = _la / ( _la + _lb ) ;
   _la = trap0 ( _threadargscomma_ _lvm , thi1 , Rd , qi ) ;
   _lb = trap0 ( _threadargscomma_ - _lvm , - thi2 , Rg , qi ) ;
   htau = 1.0 / tadj / ( _la + _lb ) ;
   hinf = 1.0 / ( 1.0 + exp ( ( _lvm - thinf ) / qinf ) ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double trap0 (  double _lv , double _lth , double _la , double _lq ) {
   double _ltrap0;
 if ( fabs ( ( _lv - _lth ) / _lq ) > 1e-6 ) {
     _ltrap0 = _la * ( _lv - _lth ) / ( 1.0 - exp ( - ( _lv - _lth ) / _lq ) ) ;
     }
   else {
     _ltrap0 = _la * _lq ;
     }
   
return _ltrap0;
 }
 
static void _hoc_trap0(void) {
  double _r;
   _r =  trap0 (  *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
 {
   trates ( _threadargscomma_ v + vshift ) ;
   m = minf ;
   h = hinf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gna = tadj * gbar * m * m * m * h ;
   ina = ( 1e-4 ) * gna * ( v - ena ) ;
   }
 _current += ina;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
 _g = _nrn_current(_v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ena = _ion_ena;
 { error =  states();
 if(error){fprintf(stderr,"at line 97 in file na.mod:\n        SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
   _t_minf = makevector(200*sizeof(double));
   _t_hinf = makevector(200*sizeof(double));
   _t_mtau = makevector(200*sizeof(double));
   _t_htau = makevector(200*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/landauland/Dropbox/SabatiniLab/neuron-modeling/smithAdaptation/mod.files/na.mod";
static const char* nmodl_file_text = 
  "\n"
  "COMMENT\n"
  "\n"
  "26 Ago 2002 Modification of original channel to allow variable time\n"
  "step and to correct an initialization error.\n"
  "    Done by Michael Hines(michael.hines@yale.e) and Ruggero\n"
  "Scorcioni(rscorcio@gmu.edu) at EU Advance Course in Computational\n"
  "Neuroscience. Obidos, Portugal\n"
  "11 Jan 2007\n"
  "    Glitch in trap where (v/th) was where (v-th)/q is. (thanks Ronald\n"
  "van Elburg!)\n"
  "\n"
  "na.mod\n"
  "\n"
  "Sodium channel, Hodgkin-Huxley style kinetics.  \n"
  "\n"
  "Kinetics were fit to data from Huguenard et al. (1988) and Hamill et\n"
  "al. (1991)\n"
  "\n"
  "qi is not well constrained by the data, since there are no points\n"
  "between -80 and -55.  So this was fixed at 5 while the thi1,thi2,Rg,Rd\n"
  "were optimized using a simplex least square proc\n"
  "\n"
  "voltage dependencies are shifted approximately from the best\n"
  "fit to give higher threshold\n"
  "\n"
  "Author: Zach Mainen, Salk Institute, 1994, zach@salk.edu\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX na\n"
  "	USEION na READ ena WRITE ina\n"
  "	RANGE m, h, gna, gbar\n"
  "	GLOBAL tha, thi1, thi2, qa, qi, qinf, thinf \n"
  "	RANGE minf, hinf, mtau, htau\n"
  "	GLOBAL Ra, Rb, Rd, Rg\n"
  "	GLOBAL q10, temp, tadj, vmin, vmax, vshift\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	gbar = 1000   	(pS/um2)	: 0.12 mho/cm2\n"
  "	vshift = -10	(mV)		: voltage shift (affects all)\n"
  "								\n"
  "	tha  = -35	(mV)		: v 1/2 for act		(-42)\n"
  "	qa   = 9	(mV)		: act slope		\n"
  "	Ra   = 0.182	(/ms)		: open (v)		\n"
  "	Rb   = 0.124	(/ms)		: close (v)		\n"
  "\n"
  "	thi1  = -50	(mV)		: v 1/2 for inact 	\n"
  "	thi2  = -75	(mV)		: v 1/2 for inact 	\n"
  "	qi   = 5	(mV)	        : inact tau slope\n"
  "	thinf  = -65	(mV)		: inact inf slope	\n"
  "	qinf  = 6.2	(mV)		: inact inf slope\n"
  "	Rg   = 0.0091	(/ms)		: inact (v)	\n"
  "	Rd   = 0.024	(/ms)		: inact recov (v) \n"
  "\n"
  "	temp = 23	(degC)		: original temp \n"
  "	q10  = 2.3			: temperature sensitivity\n"
  "\n"
  "	v 		(mV)\n"
  "	dt		(ms)\n"
  "	celsius		(degC)\n"
  "	vmin = -120	(mV)\n"
  "	vmax = 100	(mV)\n"
  "}\n"
  "\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "	(pS) = (picosiemens)\n"
  "	(um) = (micron)\n"
  "} \n"
  "\n"
  "ASSIGNED {\n"
  "	ina 		(mA/cm2)\n"
  "	gna		(pS/um2)\n"
  "	ena		(mV)\n"
  "	minf 		hinf\n"
  "	mtau (ms)	htau (ms)\n"
  "	tadj\n"
  "}\n"
  " \n"
  "\n"
  "STATE { m h }\n"
  "\n"
  "INITIAL { \n"
  "	trates(v+vshift)\n"
  "	m = minf\n"
  "	h = hinf\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "        SOLVE states METHOD cnexp\n"
  "        gna = tadj*gbar*m*m*m*h\n"
  "	ina = (1e-4) * gna * (v - ena)\n"
  "} \n"
  "\n"
  "LOCAL mexp, hexp \n"
  "\n"
  "DERIVATIVE states {   :Computes state variables m, h, and n \n"
  "        trates(v+vshift)      :             at the current v and dt.\n"
  "        m' =  (minf-m)/mtau\n"
  "        h' =  (hinf-h)/htau\n"
  "}\n"
  "\n"
  "PROCEDURE trates(v) {  \n"
  "                      \n"
  "        \n"
  "        TABLE minf,  hinf, mtau, htau\n"
  "	DEPEND  celsius, temp, Ra, Rb, Rd, Rg, tha, thi1, thi2, qa, qi, qinf\n"
  "	\n"
  "	FROM vmin TO vmax WITH 199\n"
  "\n"
  "	rates(v): not consistently executed from here if usetable == 1\n"
  "\n"
  ":        tinc = -dt * tadj\n"
  "\n"
  ":        mexp = 1 - exp(tinc/mtau)\n"
  ":        hexp = 1 - exp(tinc/htau)\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE rates(vm) {  \n"
  "        LOCAL  a, b\n"
  "\n"
  "	a = trap0(vm,tha,Ra,qa)\n"
  "	b = trap0(-vm,-tha,Rb,qa)\n"
  "\n"
  "        tadj = q10^((celsius - temp)/10)\n"
  "\n"
  "	mtau = 1/tadj/(a+b)\n"
  "	minf = a/(a+b)\n"
  "\n"
  "		:\"h\" inactivation \n"
  "\n"
  "	a = trap0(vm,thi1,Rd,qi)\n"
  "	b = trap0(-vm,-thi2,Rg,qi)\n"
  "	htau = 1/tadj/(a+b)\n"
  "	hinf = 1/(1+exp((vm-thinf)/qinf))\n"
  "}\n"
  "\n"
  "\n"
  "FUNCTION trap0(v,th,a,q) {\n"
  "	if (fabs((v-th)/q) > 1e-6) {\n"
  "	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))\n"
  "	} else {\n"
  "	        trap0 = a * q\n"
  " 	}\n"
  "}	\n"
  ;
#endif
