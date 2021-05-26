#include <RcppArmadillo.h>
#include <array>
#include <vector>
#include <math.h>
using namespace std;

double ifelse(bool condition, const double& result1, const double& result2) {
  if(condition) {
    return result1;
  } else {
    return result2;
  }
}

inline double heaviside(const double& arg) {
  return arg <= 0.0 ? 0.0 : 1.0;
}

inline double dirac(const double& arg) {
  return abs(arg) <= numeric_limits<double>::epsilon()  ? 0 : numeric_limits<double>::infinity();
}

inline double Min(const double& arg1) {
  return arg1;
}

inline double Max(const double& arg1) {
  return arg1;
}

inline double Min(const double& arg1, const double& arg2) {
  return arg1 <= arg2 ? arg1 : arg2;
}

inline double Max(const double& arg1, const double& arg2) {
  return arg1 >= arg2 ? arg1 : arg2;
}

extern "C" {

  array<vector<double>, 2> MiniModelQss(const double& t, const vector<double>& states,
            const vector<double>& params, const vector<double>& forcs) { 


const double sigma2 = params[0];
const double Jmax25 = params[1];
const double DHaJmax = params[2];
const double DHdJmax = params[3];
const double DsJmax = params[4];
const double kD0 = params[5];
const double kDinh = params[6];
const double kf = params[7];
const double kp = params[8];
const double fcyc = params[9];
const double fpseudo = params[10];
const double theta = params[11];
const double gamma1 = params[12];
const double gamma2 = params[13];
const double gamma3 = params[14];
const double PhiqEmax = params[15];
const double Kinh0 = params[16];
const double fprot = params[17];
const double Krep25 = params[18];
const double DHaKrep = params[19];
const double DHdKrep = params[20];
const double DsKrep = params[21];
const double alphar_alpha25 = params[22];
const double DHaAlphar = params[23];
const double Iac = params[24];
const double alpharac = params[25];
const double alpharav = params[26];
const double thetaalphar = params[27];
const double fR0 = params[28];
const double alphafR = params[29];
const double thetafR = params[30];
const double Vrmax = params[31];
const double KmPGA = params[32];
const double RB = params[33];
const double Kc25 = params[34];
const double DHaKc = params[35];
const double Ko25 = params[36];
const double DHaKo = params[37];
const double Kmc25 = params[38];
const double DHaKmc = params[39];
const double Kmo25 = params[40];
const double DHaKmo = params[41];
const double KaRCA = params[42];
const double ac = params[43];
const double bc = params[44];
const double fRBmin = params[45];
const double O2 = params[46];
const double RCA = params[47];
const double DHdRCA = params[48];
const double ToRCA = params[49];
const double DHaRCA = params[50];
const double KmRuBP = params[51];
const double Vch = params[52];
const double KiPGA = params[53];
const double TPU25 = params[54];
const double DHaTPU = params[55];
const double DHdTPU = params[56];
const double DsTPU = params[57];
const double DHaGc = params[58];
const double DHdGc = params[59];
const double DsGc = params[60];
const double DHaGw = params[61];
const double DHdGw = params[62];
const double DsGw = params[63];
const double Rm25 = params[64];
const double DHaRm = params[65];
const double Vref = params[66];
const double Scm = params[67];
const double Sm = params[68];
const double falphaSc = params[69];
const double gcm25 = params[70];
const double gw25 = params[71];
const double D0 = params[72];
const double fI0 = params[73];
const double gswm = params[74];
const double alphafI = params[75];
const double thetafI = params[76];
const double gbw = params[77];
const double volume_chamber = params[78];
const double leaf_surface = params[79];
const double Flow = params[80];
const double alphab = params[81];
const double alphag = params[82];
const double alphared = params[83];
const double alphabp = params[84];
const double alphagp = params[85];
const double alpharp = params[86];
const double PGA = states[0];
const double RuBP = states[1];
const double Cc = states[2];
const double Ccyt = states[3];
const double Ci = states[4];
const double Ca = states[5];
const double H2OS = states[6];
const double Ib = forcs[0];
const double Ig = forcs[1];
const double Ir = forcs[2];
const double Ta = forcs[3];
const double Tl = forcs[4];
const double H2OR = forcs[5];
const double CO2R = forcs[6];
const double zeroKinh = 0.0;
const double zero_pressure = 100000.0;
const double air_pressure = 1.01e8;
const double Flux0 = 9.999999999999999e-23;
const double R = 8310.0;
const double Tzero = 273.15;
const double es0 = 610780.0;
const double es_k = 17.269;
const double es_Tref = 237.3;
const double Tref = 298.15;
const double PARaP = Ib * alphabp + Ig * alphagp + Ir * alpharp;
const double PAR = Ib + Ig + Ir;
const double Photo = (Flow * CO2R - Flow * Ca) / leaf_surface;
const double Trmmol = (Flow * (H2OS - H2OR)) / (leaf_surface * (1 - H2OS));
const double Kmapp_RuBP = KmRuBP * (1 + PGA / (Vch * KiPGA));
const double fRCA = (DHdRCA * exp(((Tl - ToRCA) * DHaRCA) / (ToRCA * R * Tl))) / (DHdRCA - DHaRCA * (1 - exp((DHdRCA * (Tl - ToRCA)) / (ToRCA * R * Tl))));
const double Kc = Kc25 * exp(((Tl - Tref) * DHaKc) / (Tref * R * Tl));
const double Kmc = Kmc25 * exp(((Tl - Tref) * DHaKmc) / (Tref * R * Tl));
const double alphar_alpha = alphar_alpha25 * exp((-((Tl - Tref)) * DHaAlphar) / (Tref * R * Tl));
const double Krep = (Krep25 * exp(((Tl - Tref) * DHaKrep) / (Tref * R * Tl)) * (1 + exp((Tref * DsKrep - DHdKrep) / (Tref * R)))) / (1 + exp((Tl * DsKrep - DHdKrep) / (Tl * R)));
const double Fm_a = kf / (kf + kD0);
const double Jmax = (Jmax25 * exp(((Tl - Tref) * DHaJmax) / (Tref * R * Tl)) * (1 + exp((Tref * DsJmax - DHdJmax) / (Tref * R)))) / (1 + exp((Tl * DsJmax - DHdJmax) / (Tl * R)));
const double TPU = (TPU25 * exp(((Tl - Tref) * DHaTPU) / (Tref * R * Tl)) * (1 + exp((Tref * DsTPU - DHdTPU) / (Tref * R)))) / (1 + exp((Tl * DsTPU - DHdTPU) / (Tl * R)));
const double Rm = Rm25 * exp(((Tl - Tref) * DHaRm) / (Tref * R * Tl));
const double gbc = gbw / 1.37;
const double ea = H2OR * air_pressure;
const double gw = (gw25 * exp(((Tl - Tref) * DHaGw) / (Tref * R * Tl)) * (1 + exp((Tref * DsGw - DHdGw) / (Tref * R)))) / (1 + exp((Tl * DsGw - DHdGw) / (Tl * R)));
const double Mv = (R * Ta) / air_pressure;
const double PARa = Ib * alphab + Ig * alphag + Ir * alphared;
const double fRBss_nr = Min(1.0,ac + bc * Cc);
const double Ko = Ko25 * exp(((Tl - Tref) * DHaKo) / (Tref * R * Tl));
const double Fo_a = kf / (kf + kD0 + kp);
const double fI_a = thetafI;
const double es_leaf = es0 * exp((es_k * (Tl - Tzero)) / (es_Tref + (Tl - Tzero)));
const double Kmo = Kmo25 * exp(((Tl - Tref) * DHaKmo) / (Tref * R * Tl));
const double fRBmax = (RCA * fRCA) / (RCA * fRCA + KaRCA);
const double fI_b = -((1 + fI0 + alphafI * PAR));
const double fRuBP = (1 / ((2 * RB) / Vch)) * ((RB / Vch + Kmapp_RuBP + RuBP / Vch) - sqrt(pow(RB / Vch + Kmapp_RuBP + RuBP / Vch,2) - (4 * RB * RuBP) / pow(Vch,2)));
const double alpharss = Min(1.0 + (Ib / Iac) * alpharac,(1.0 + alpharac) - ((alphar_alpha * (Ib - Iac) + alpharav) - sqrt(pow(alphar_alpha * (Ib - Iac) + alpharav,2) - 4 * alphar_alpha * thetaalphar * alpharav * (Ib - Iac))) / (2 * thetaalphar));
const double Fm = Fm_a;
const double fRss = fR0 + ((alphafR * PARa + (1 - fR0)) - sqrt(pow(alphafR * PARa + (1 - fR0),2) - 4 * alphafR * thetafR * PARa * (1 - fR0))) / (2 * thetafR);
const double VPDleaf = Max(es_leaf - ea,zero_pressure);
const double PhiIId = (Fm_a - Fo_a) / Fm_a;
const double fI_c = fI0 + alphafI * PAR;
const double phi = (Kmc * Ko * O2) / (Kmo * Kc * Cc);
const double Fo = Fo_a;
const double alphar = alpharss;
const double fR = fRss;
const double fvpd = 1 / (1 + VPDleaf / D0);
const double gtw = (Trmmol * (air_pressure - (es_leaf + ea) / 2.0)) / VPDleaf;
const double VrTPU = (3.0 * TPU * (2.0 + 1.5 * phi)) / (1 - 0.5 * phi);
const double fI = (-fI_b - sqrt(pow(fI_b,2) - 4 * fI_a * fI_c)) / (2 * fI_a);
const double PhiIIop = (Fm - Fo) / Fm;
const double PARaP2 = sigma2 * alphar * PARaP;
const double Fmp_d = (alphar * kf) / (kf + kDinh);
const double VrE = (fR * Vrmax * PGA) / (PGA + KmPGA);
const double Sc = Scm * falphaSc * alphar;
const double Cond = 1 / (1 / gtw - 1 / gbw);
const double Fop_d = (alphar * kf) / (kf + kDinh);
const double gss = fI * fvpd * gswm;
const double J2pm = (((Min(VrTPU,VrE) / (1 - fpseudo / (1 - fcyc))) * 2.0) / (2.0 + 1.5 * phi)) * (2.0 + 2.0 * phi);
const double gc = ((Sc / Sm) * gcm25 * exp(((Tl - Tref) * DHaGc) / (Tref * R * Tl)) * (1 + exp((Tref * DsGc - DHdGc) / (Tref * R)))) / (1 + exp((Tl * DsGc - DHdGc) / (Tl * R)));
const double gsw = gss;
const double J2pp = ((PhiIIop * PARaP2 + Jmax) - sqrt(pow(PhiIIop * PARaP2 + Jmax,2) - 4 * PhiIIop * theta * Jmax * PARaP2)) / (2 * theta);
const double gsc = gsw / 1.56;
const double qPp = ifelse(PARaP2 > Flux0,(J2pp / PARaP2) / PhiIIop,1);
const double transpiration = ((VPDleaf / (air_pressure - (es_leaf + ea) / 2.0)) * 1) / (1 / gsw + 1 / gbw);
const double d_Ci_dt = (((Ca - Ci) / (1 / gsc + 1 / gbc) - (Ci - Ccyt) * gw) * Mv) / Vref;
const double A = (Ca - Ci) / (1 / gsc + 1 / gbc);
const double qPm = (J2pm / J2pp) * qPp;
const double d_H2OS_dt = (((-((Flow + leaf_surface * transpiration)) * H2OS + Flow * H2OR + leaf_surface * transpiration) * R * Ta) / volume_chamber) / air_pressure;
const double gm = A / (Ci - Cc);
const double d_Ca_dt = ((((-Flow * Ca + Flow * CO2R) - leaf_surface * A) * R * Ta) / volume_chamber) / air_pressure;
const double qPno_qD = Min(qPm,qPp);
const double fqEss = 1 - qPno_qD;
const double fZ = fqEss;
const double PhiqEss = (fqEss * gamma1 + fqEss * fqEss * gamma2 + fqEss * gamma3) * PhiqEmax;
const double fP = fqEss;
const double Kinh = Max(Kinh0 - fprot * PhiqEss,zeroKinh);
const double PhiIIoss = PhiIIop - PhiqEss;
const double PhiqE = (fP * gamma1 + fP * fZ * gamma2 + fZ * gamma3) * PhiqEmax;
const double PSIId = (PARa * alphar * Kinh) / (PARa * alphar * Kinh + Krep);
const double kD = (kp / (PhiIId - PhiqE) - kf) - kp;
const double PhiIIo = PhiIIop - PhiqE;
const double qI = Fm_a / ((1 - PSIId) * Fm_a + PSIId * Fmp_d) - 1;
const double J2qE = ifelse(PhiIIo < PhiIIoss,(1 - (PhiIIoss - PhiIIo) / PhiIIoss) * J2pp,J2pp);
const double Fmp_a = (alphar * kf) / (kf + kD);
const double Fop_a = (alphar * kf) / (kf + kD + kp);
const double Fmp = (1 - PSIId) * Fmp_a + PSIId * Fmp_d;
const double VrJ = (((J2qE * (1 - fpseudo / (1 - fcyc))) / 2.0) * (2.0 + 1.5 * phi)) / (2.0 + 2.0 * phi);
const double J2 = Min(J2qE,J2pm);
const double fRBss_r = Min(ifelse(PAR > Flux0,(Min(VrTPU,VrJ) / (2.0 + 1.5 * phi)) / ((Kc * RB * Cc) / (Cc + Kmc * (1.0 + O2 / Kmo))),fRBmin),fRBmax);
const double Vr = Min(VrJ,Min(VrTPU,VrE));
const double NPQ = Fm_a / Fmp - 1;
const double qP = (J2 / PARaP2) / PhiIIo;
const double PhiII = ifelse(PARaP2 > Flux0,J2 / PARaP2,PhiIIo);
const double qM = NPQ - ((Fm_a / Fmp) * alphar - 1);
const double reg_limit = ifelse(abs(VrJ - Vr) < Flux0,1.0,ifelse(abs(VrTPU - Vr) < Flux0,2.0,3.0));
const double fRBss = Min(fRBss_nr,fRBss_r);
const double qE = (NPQ - qM) - qI;
const double fRB = fRBss;
const double Vc = (fRB * fRuBP * Kc * RB * Cc) / (Cc + Kmc * (1.0 + O2 / Kmo));
const double Rp = 0.5 * Vc * phi;
const double d_RuBP_dt = ((1.0 + phi) / (2.0 + 1.5 * phi)) * Vr - Vc * (1.0 + phi);
const double d_Cc_dt = (((Ccyt - Cc) * gc - Vc) * Mv) / Vref;
const double d_PGA_dt = (2.0 * Vc + 1.5 * Vc * phi) - Vr;
const double d_Ccyt_dt = ((((Ci - Ccyt) * gw + Rp + Rm) - (Ccyt - Cc) * gc) * Mv) / Vref;

vector<double> derivatives{d_PGA_dt,d_RuBP_dt,d_Cc_dt,d_Ccyt_dt,d_Ci_dt,d_Ca_dt,d_H2OS_dt};
 vector<double> observed{Vr, Vc, fRuBP, fRB, PARaP, PARa, PAR, fqEss, fP, fZ, PhiqE, PhiIIoss, PhiIIo, qP, PhiII, alphar, PSIId, VrJ, NPQ, qI, qE, qM, fR, Rp, A, gsw, gss, VPDleaf, Sc, Photo, transpiration, Trmmol, Cond, gm, reg_limit};
 array<vector<double>,2> output{derivatives, observed};
 return output;
}


};

