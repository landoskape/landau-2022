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
 
#define nrn_init _nrn_init__caL3d
#define _nrn_initial _nrn_initial__caL3d
#define nrn_cur _nrn_cur__caL3d
#define _nrn_current _nrn_current__caL3d
#define nrn_jacob _nrn_jacob__caL3d
#define nrn_state _nrn_state__caL3d
#define _net_receive _net_receive__caL3d 
#define _f_rates _f_rates__caL3d 
#define kstates kstates__caL3d 
#define rates rates__caL3d 
 
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
#define a _p[0]
#define b _p[1]
#define C _p[2]
#define O _p[3]
#define ica _p[4]
#define cao _p[5]
#define cai _p[6]
#define DC _p[7]
#define DO _p[8]
#define _g _p[9]
#define _ion_cai	*_ppvar[0]._pval
#define _ion_cao	*_ppvar[1]._pval
#define _ion_ica	*_ppvar[2]._pval
#define _ion_dicadv	*_ppvar[3]._pval
 
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
 static void _hoc_efun(void);
 static void _hoc_ghk(void);
 static void _hoc_rates(void);
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
 "setdata_caL3d", _hoc_setdata,
 "efun_caL3d", _hoc_efun,
 "ghk_caL3d", _hoc_ghk,
 "rates_caL3d", _hoc_rates,
 0, 0
};
#define efun efun_caL3d
#define ghk ghk_caL3d
 extern double efun( double );
 extern double ghk( double , double , double );
 /* declare global and static user variables */
#define Rb Rb_caL3d
 double Rb = 0.2;
#define Ra Ra_caL3d
 double Ra = 1.6;
#define p p_caL3d
 double p = 0.0002;
#define q10 q10_caL3d
 double q10 = 3;
#define q q_caL3d
 double q = 13;
#define tadj tadj_caL3d
 double tadj = 0;
#define temp temp_caL3d
 double temp = 22;
#define th th_caL3d
 double th = 5;
#define usetable usetable_caL3d
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_caL3d", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "p_caL3d", "cm/s",
 "th_caL3d", "mV",
 "q_caL3d", "mV",
 "Ra_caL3d", "/ms",
 "Rb_caL3d", "/ms",
 "temp_caL3d", "degC",
 "a_caL3d", "/ms",
 "b_caL3d", "/ms",
 0,0
};
 static double C0 = 0;
 static double O0 = 0;
 static double delta_t = 1;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "p_caL3d", &p_caL3d,
 "th_caL3d", &th_caL3d,
 "q_caL3d", &q_caL3d,
 "Ra_caL3d", &Ra_caL3d,
 "Rb_caL3d", &Rb_caL3d,
 "temp_caL3d", &temp_caL3d,
 "q10_caL3d", &q10_caL3d,
 "tadj_caL3d", &tadj_caL3d,
 "usetable_caL3d", &usetable_caL3d,
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
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"caL3d",
 0,
 "a_caL3d",
 "b_caL3d",
 0,
 "C_caL3d",
 "O_caL3d",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 10, _prop);
 	/*initialize range parameters*/
 	_prop->param = _p;
 	_prop->param_size = 10;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* cao */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
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

 void _caL3d_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 10, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 caL3d /Users/landauland/Dropbox/SabatiniLab/neuron-modeling/smithAdaptation/mod.files/x86_64/caL3d.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double F = 96485.3;
 static double R = 8.3145;
 static double *_t_a;
 static double *_t_b;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_rates(double);
static int rates(double);
 extern double *_getelm();
 
#define _MATELM1(_row,_col)	*(_getelm(_row + 1, _col + 1))
 
#define _RHS1(_arg) _coef1[_arg + 1]
 static double *_coef1;
 
#define _linmat1  1
 static void* _sparseobj1;
 static void* _cvsparseobj1;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_rates(double);
 static int _slist1[2], _dlist1[2]; static double *_temp1;
 static int kstates();
 
static int kstates ()
 {_reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=1;_i<2;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 /* ~ C <-> O ( a , b )*/
 f_flux =  a * C ;
 b_flux =  b * O ;
 _RHS1( 1) -= (f_flux - b_flux);
 
 _term =  a ;
 _MATELM1( 1 ,1)  += _term;
 _term =  b ;
 _MATELM1( 1 ,0)  -= _term;
 /*REACTION*/
   /* C + O = 1.0 */
 _RHS1(0) =  1.0;
 _MATELM1(0, 0) = 1;
 _RHS1(0) -= O ;
 _MATELM1(0, 1) = 1;
 _RHS1(0) -= C ;
 /*CONSERVATION*/
   } return _reset;
 }
 static double _mfac_rates, _tmin_rates;
 static void _check_rates();
 static void _check_rates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_Ra;
  static double _sav_Rb;
  static double _sav_th;
  static double _sav_celsius;
  static double _sav_temp;
  static double _sav_q10;
  if (!usetable) {return;}
  if (_sav_Ra != Ra) { _maktable = 1;}
  if (_sav_Rb != Rb) { _maktable = 1;}
  if (_sav_th != th) { _maktable = 1;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_sav_temp != temp) { _maktable = 1;}
  if (_sav_q10 != q10) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rates =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_rates)/200.; _mfac_rates = 1./_dx;
   for (_i=0, _x=_tmin_rates; _i < 201; _x += _dx, _i++) {
    _f_rates(_x);
    _t_a[_i] = a;
    _t_b[_i] = b;
   }
   _sav_Ra = Ra;
   _sav_Rb = Rb;
   _sav_th = th;
   _sav_celsius = celsius;
   _sav_temp = temp;
   _sav_q10 = q10;
  }
 }

 static int rates(double _lv){ _check_rates();
 _n_rates(_lv);
 return 0;
 }

 static void _n_rates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rates(_lv); return; 
}
 _xi = _mfac_rates * (_lv - _tmin_rates);
 if (isnan(_xi)) {
  a = _xi;
  b = _xi;
  return;
 }
 if (_xi <= 0.) {
 a = _t_a[0];
 b = _t_b[0];
 return; }
 if (_xi >= 200.) {
 a = _t_a[200];
 b = _t_b[200];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 a = _t_a[_i] + _theta*(_t_a[_i+1] - _t_a[_i]);
 b = _t_b[_i] + _theta*(_t_b[_i+1] - _t_b[_i]);
 }

 
static int  _f_rates (  double _lv ) {
   tadj = pow( q10 , ( ( celsius - temp ) / 10.0 ) ) ;
   a = Ra / ( 1.0 + exp ( - ( _lv - th ) / q ) ) * tadj ;
   b = Rb / ( 1.0 + exp ( ( _lv - th ) / q ) ) * tadj ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
    _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double ghk (  double _lv , double _lci , double _lco ) {
   double _lghk;
 double _lz ;
 _lz = ( 0.001 ) * 2.0 * F * _lv / ( R * ( celsius + 273.15 ) ) ;
   _lghk = ( .001 ) * 2.0 * F * ( _lci * efun ( _threadargscomma_ - _lz ) - _lco * efun ( _threadargscomma_ _lz ) ) ;
   
return _lghk;
 }
 
static void _hoc_ghk(void) {
  double _r;
   _r =  ghk (  *getarg(1) , *getarg(2) , *getarg(3) );
 hoc_retpushx(_r);
}
 
double efun (  double _lz ) {
   double _lefun;
 if ( fabs ( _lz ) < 1e-4 ) {
     _lefun = 1.0 - _lz / 2.0 ;
     }
   else {
     _lefun = _lz / ( exp ( _lz ) - 1.0 ) ;
     }
   
return _lefun;
 }
 
static void _hoc_efun(void) {
  double _r;
   _r =  efun (  *getarg(1) );
 hoc_retpushx(_r);
}
 
/*CVODE ode begin*/
 static int _ode_spec1() {_reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<2;_i++) _p[_dlist1[_i]] = 0.0;}
 /* ~ C <-> O ( a , b )*/
 f_flux =  a * C ;
 b_flux =  b * O ;
 DC -= (f_flux - b_flux);
 DO += (f_flux - b_flux);
 
 /*REACTION*/
   /* C + O = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1() {_reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<2;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 /* ~ C <-> O ( a , b )*/
 _term =  a ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 0 ,1)  -= _term;
 _term =  b ;
 _MATELM1( 1 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
   /* C + O = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE end*/
 
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
  cai = _ion_cai;
  cao = _ion_cao;
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
 _cvode_sparse(&_cvsparseobj1, 2, _dlist1, _p, _ode_matsol1, &_coef1);
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
  cai = _ion_cai;
  cao = _ion_cao;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  C = C0;
  O = O0;
 {
   C = 1.0 ;
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
  cai = _ion_cai;
  cao = _ion_cao;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   rates ( _threadargscomma_ v ) ;
   ica = O * p * ghk ( _threadargscomma_ v , cai , cao ) ;
   }
 _current += ica;

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
  cai = _ion_cai;
  cao = _ion_cao;
 _g = _nrn_current(_v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
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
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
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
  cai = _ion_cai;
  cao = _ion_cao;
 { error = sparse(&_sparseobj1, 2, _slist1, _dlist1, _p, &t, dt, kstates,&_coef1, _linmat1);
 if(error){fprintf(stderr,"at line 98 in file caL3d.mod:\n	SOLVE kstates METHOD sparse\n"); nrn_complain(_p); abort_run(error);}
    if (secondorder) {
    int _i;
    for (_i = 0; _i < 2; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 } }}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_a = makevector(201*sizeof(double));
   _t_b = makevector(201*sizeof(double));
 _slist1[0] = &(O) - _p;  _dlist1[0] = &(DO) - _p;
 _slist1[1] = &(C) - _p;  _dlist1[1] = &(DC) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/landauland/Dropbox/SabatiniLab/neuron-modeling/smithAdaptation/mod.files/caL3d.mod";
static const char* nmodl_file_text = 
  "\n"
  "COMMENT\n"
  "\n"
  "High threshold Ca2+ channel\n"
  "\n"
  "2-state kinetics with sigmoidal voltage-dependence\n"
  "\n"
  "  C<->O\n"
  "\n"
  "Goldman-Hodgkin-Katz equations\n"
  "\n"
  "     # MODEL\n"
  "    |   MODEL AUTHOR  : D.A. McCormick & J. Huguenard\n"
  "    |   MODEL DATE    : 1992\n"
  "    |   MODEL REF     : A model of the electrophysiological properties of \n"
  "thalamocortical relay neurons. J Neurophysiol, 1992 Oct, 68(4):1384-400.\n"
  " \n"
  "    # EXPERIMENT\n"
  "    |   EXP AUTHOR    : Kay AR; Wong RK\n"
  "    |   EXP DATE      : 1987\n"
  "    |   EXP REF       : Journal of Physiology, 1987 Nov, 392:603-16.\n"
  "    |   ANIMAL        : guinea-pig\n"
  "    |   BRAIN REGION  : hippocampus\n"
  "    |   CELL TYPE     : Ca1 pyramidal\n"
  "    |   TECHNIQUE     : slices, whole-cell\n"
  "    |   RECORDING METHOD  : voltage-clamp\n"
  "    |   TEMPERATURE   : 20-22\n"
  " \n"
  "Reference:\n"
  "\n"
  "   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J. Synthesis of models for\n"
  "   excitable membranes, synaptic transmission and neuromodulation using a \n"
  "   common kinetic formalism, Journal of Computational Neuroscience 1: \n"
  "   195-230, 1994.\n"
  "\n"
  "  (electronic copy available at http://cns.iaf.cnrs-gif.fr)\n"
  "\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX caL3d\n"
  "	USEION ca READ cai, cao WRITE ica\n"
  "	RANGE O, C, I\n"
  "	RANGE a,b\n"
  "	GLOBAL Ra, Rb, q, th, p\n"
  "	GLOBAL q10, temp, tadj\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	F = (faraday) (coulomb)\n"
  "	R = (k-mole) (joule/degC)\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "	(pS) = (picosiemens)\n"
  "	(um) = (micron)\n"
  "	(mM) = (milli/liter)\n"
  "} \n"
  "\n"
  "PARAMETER {\n"
  "	p    = 0.2e-3  	(cm/s)		: max permeability\n"
  "	v 		(mV)\n"
  "\n"
  "	th   = 5	(mV)		: v 1/2 for on/off\n"
  "	q   = 13	(mV)		: voltage dependence\n"
  "\n"
  "	: max rates\n"
  "\n"
  "	Ra   = 1.6	(/ms)		: open (v)\n"
  "	Rb   = 0.2	(/ms)		: close (v)\n"
  "\n"
  "	celsius		(degC)\n"
  "	temp = 22	(degC)		: original temp\n"
  "	q10  = 3			: temperature sensitivity\n"
  "} \n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "	ica 		(mA/cm2)\n"
  "	cao		(mM)\n"
  "	cai		(mM)\n"
  "	a (/ms)	b (/ms)\n"
  "	tadj\n"
  "}\n"
  " \n"
  "\n"
  "STATE { C O }\n"
  "\n"
  "INITIAL { \n"
  "	C = 1 \n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "	rates(v)\n"
  "	SOLVE kstates METHOD sparse\n"
  "	ica = O * p * ghk(v,cai,cao)\n"
  "} \n"
  "\n"
  "\n"
  "KINETIC kstates {\n"
  "	~ C <-> O 	(a,b)	\n"
  "	CONSERVE C+O = 1\n"
  "}	\n"
  "	\n"
  "PROCEDURE rates(v(mV)) {\n"
  "	TABLE a, b\n"
  "	DEPEND Ra, Rb, th, celsius, temp, q10\n"
  "	FROM -100 TO 100 WITH 200\n"
  "\n"
  "	tadj = q10 ^ ((celsius - temp)/10 (degC))\n"
  "\n"
  "	a = Ra / (1 + exp(-(v-th)/q)) * tadj\n"
  "	b = Rb / (1 + exp((v-th)/q)) * tadj\n"
  "}\n"
  "\n"
  ": Special gear for calculating the Ca2+ reversal potential\n"
  ": via Goldman-Hodgkin-Katz eqn.\n"
  ": [Ca2+]o \"cao\" and [Ca2+]i \"cai\" are assumed to be set elsewhere\n"
  "\n"
  "\n"
  "FUNCTION ghk(v(mV), ci(mM), co(mM)) (0.001 coul/cm3) {\n"
  "	LOCAL z\n"
  "\n"
  "	z = (0.001)*2*F*v/(R*(celsius+273.15))\n"
  "	ghk = (.001)*2*F*(ci*efun(-z) - co*efun(z))\n"
  "}\n"
  "\n"
  "FUNCTION efun(z) {\n"
  "	if (fabs(z) < 1e-4) {\n"
  "		efun = 1 - z/2\n"
  "	}else{\n"
  "		efun = z/(exp(z) - 1)\n"
  "	}\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "\n"
  ;
#endif
