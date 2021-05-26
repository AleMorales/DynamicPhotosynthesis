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

  array<vector<double>, 2> MiniModelQssfR(const double& t, const vector<double>& states,
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
const double KiqEp = params[16];
const double KdqEp = params[17];
const double KiqEz = params[18];
const double KdqEz = params[19];
const double Kinh0 = params[20];
const double fprot = params[21];
const double Krep25 = params[22];
const double DHaKrep = params[23];
const double DHdKrep = params[24];
const double DsKrep = params[25];
const double alphar_alpha25 = params[26];
const double DHaAlphar = params[27];
const double DHaKalpha = params[28];
const double DsKalpha = params[29];
const double DHdKalpha = params[30];
const double Iac = params[31];
const double alpharac = params[32];
const double alpharav = params[33];
const double thetaalphar = params[34];
const double Kialpha25 = params[35];
const double Kdalpha25 = params[36];
const double fR0 = params[37];
const double alphafR = params[38];
const double thetafR = params[39];
const double Vrmax = params[40];
const double KmPGA = params[41];
const double RB = params[42];
const double Kc25 = params[43];
const double DHaKc = params[44];
const double Ko25 = params[45];
const double DHaKo = params[46];
const double Kmc25 = params[47];
const double DHaKmc = params[48];
const double Kmo25 = params[49];
const double DHaKmo = params[50];
const double KaRCA = params[51];
const double ac = params[52];
const double bc = params[53];
const double KdRB = params[54];
const double Krca = params[55];
const double fRBmin = params[56];
const double O2 = params[57];
const double RCA = params[58];
const double DHdRCA = params[59];
const double ToRCA = params[60];
const double DHaRCA = params[61];
const double KmRuBP = params[62];
const double Vch = params[63];
const double KiPGA = params[64];
const double TPU25 = params[65];
const double DHaTPU = params[66];
const double DHdTPU = params[67];
const double DsTPU = params[68];
const double DHaGc = params[69];
const double DHdGc = params[70];
const double DsGc = params[71];
const double DHaGw = params[72];
const double DHdGw = params[73];
const double DsGw = params[74];
const double Rm25 = params[75];
const double DHaRm = params[76];
const double kPR = params[77];
const double Vref = params[78];
const double Scm = params[79];
const double Sm = params[80];
const double falphaSc = params[81];
const double gcm25 = params[82];
const double gw25 = params[83];
const double D0 = params[84];
const double fI0 = params[85];
const double gswm = params[86];
const double alphafI = params[87];
const double thetafI = params[88];
const double Kgsi = params[89];
const double Kgsd = params[90];
const double gbw = params[91];
const double volume_chamber = params[92];
const double leaf_surface = params[93];
const double Flow = params[94];
const double alphab = params[95];
const double alphag = params[96];
const double alphared = params[97];
const double alphabp = params[98];
const double alphagp = params[99];
const double alpharp = params[100];
const double PGA = states[0];
const double RuBP = states[1];
const double fRB = states[2];
const double fP = states[3];
const double fZ = states[4];
const double alphar = states[5];
const double PSIId = states[6];
const double PR = states[7];
const double Cc = states[8];
const double Ccyt = states[9];
const double Ci = states[10];
const double Ca = states[11];
const double H2OS = states[12];
const double gsw = states[13];
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
const double PhiqE = (fP * gamma1 + fP * fZ * gamma2 + fZ * gamma3) * PhiqEmax;
const double Rp = 0.5 * PR * kPR;
const double Sc = Scm * falphaSc * alphar;
const double Trmmol = (Flow * (H2OS - H2OR)) / (leaf_surface * (1 - H2OS));
const double Kmapp_RuBP = KmRuBP * (1 + PGA / (Vch * KiPGA));
const double fRCA = (DHdRCA * exp(((Tl - ToRCA) * DHaRCA) / (ToRCA * R * Tl))) / (DHdRCA - DHaRCA * (1 - exp((DHdRCA * (Tl - ToRCA)) / (ToRCA * R * Tl))));
const double Kc = Kc25 * exp(((Tl - Tref) * DHaKc) / (Tref * R * Tl));
const double Kmc = Kmc25 * exp(((Tl - Tref) * DHaKmc) / (Tref * R * Tl));
const double alphar_alpha = alphar_alpha25 * exp((-((Tl - Tref)) * DHaAlphar) / (Tref * R * Tl));
const double Kdalpha = (Kdalpha25 * exp(((Tl - Tref) * DHaKalpha) / (Tref * R * Tl)) * (1 + exp((Tref * DsKalpha - DHdKalpha) / (Tref * R)))) / (1 + exp((Tl * DsKalpha - DHdKalpha) / (Tl * R)));
const double Krep = (Krep25 * exp(((Tl - Tref) * DHaKrep) / (Tref * R * Tl)) * (1 + exp((Tref * DsKrep - DHdKrep) / (Tref * R)))) / (1 + exp((Tl * DsKrep - DHdKrep) / (Tl * R)));
const double Fm_d = kf / (kf + kDinh);
const double Fm_a = kf / (kf + kD0);
const double Fmp_d = (alphar * kf) / (kf + kDinh);
const double Jmax = (Jmax25 * exp(((Tl - Tref) * DHaJmax) / (Tref * R * Tl)) * (1 + exp((Tref * DsJmax - DHdJmax) / (Tref * R)))) / (1 + exp((Tl * DsJmax - DHdJmax) / (Tl * R)));
const double TPU = (TPU25 * exp(((Tl - Tref) * DHaTPU) / (Tref * R * Tl)) * (1 + exp((Tref * DsTPU - DHdTPU) / (Tref * R)))) / (1 + exp((Tl * DsTPU - DHdTPU) / (Tl * R)));
const double Rm = Rm25 * exp(((Tl - Tref) * DHaRm) / (Tref * R * Tl));
const double gbc = gbw / 1.37;
const double ea = H2OR * air_pressure;
const double Mv = (R * Ta) / air_pressure;
const double Photo = (Flow * CO2R - Flow * Ca) / leaf_surface;
const double fRBss_nr = Min(1.0,ac + bc * Cc);
const double Ko = Ko25 * exp(((Tl - Tref) * DHaKo) / (Tref * R * Tl));
const double Fo_d = kf / (kf + kDinh);
const double Fo_a = kf / (kf + kD0 + kp);
const double Fop_d = (alphar * kf) / (kf + kDinh);
const double gsc = gsw / 1.56;
const double es_leaf = es0 * exp((es_k * (Tl - Tzero)) / (es_Tref + (Tl - Tzero)));
const double PARa = Ib * alphab + Ig * alphag + Ir * alphared;
const double Kmo = Kmo25 * exp(((Tl - Tref) * DHaKmo) / (Tref * R * Tl));
const double Kialpha = (Kialpha25 * exp(((Tl - Tref) * DHaKalpha) / (Tref * R * Tl)) * (1 + exp((Tref * DsKalpha - DHdKalpha) / (Tref * R)))) / (1 + exp((Tl * DsKalpha - DHdKalpha) / (Tl * R)));
const double gw = (gw25 * exp(((Tl - Tref) * DHaGw) / (Tref * R * Tl)) * (1 + exp((Tref * DsGw - DHdGw) / (Tref * R)))) / (1 + exp((Tl * DsGw - DHdGw) / (Tl * R)));
const double fI_a = thetafI;
const double fRBmax = (RCA * fRCA) / (RCA * fRCA + KaRCA);
const double PARaP2 = sigma2 * alphar * PARaP;
const double fI_b = -((1 + fI0 + alphafI * PAR));
const double gc = ((Sc / Sm) * gcm25 * exp(((Tl - Tref) * DHaGc) / (Tref * R * Tl)) * (1 + exp((Tref * DsGc - DHdGc) / (Tref * R)))) / (1 + exp((Tl * DsGc - DHdGc) / (Tl * R)));
const double fRuBP = (1 / ((2 * RB) / Vch)) * ((RB / Vch + Kmapp_RuBP + RuBP / Vch) - sqrt(pow(RB / Vch + Kmapp_RuBP + RuBP / Vch,2) - (4 * RB * RuBP) / pow(Vch,2)));
const double qI = Fm_a / ((1 - PSIId) * Fm_a + PSIId * Fmp_d) - 1;
const double alpharss = Min(1.0 + (Ib / Iac) * alpharac,(1.0 + alpharac) - ((alphar_alpha * (Ib - Iac) + alpharav) - sqrt(pow(alphar_alpha * (Ib - Iac) + alpharav,2) - 4 * alphar_alpha * thetaalphar * alpharav * (Ib - Iac))) / (2 * thetaalphar));
const double Kinh = Max(Kinh0 - fprot * PhiqE,zeroKinh);
const double fI_c = fI0 + alphafI * PAR;
const double A = (Ca - Ci) / (1 / gsc + 1 / gbc);
const double VPDleaf = Max(es_leaf - ea,zero_pressure);
const double PhiIId = (Fm_a - Fo_a) / Fm_a;
const double Fo = (1 - PSIId) * Fo_a + PSIId * Fo_d;
const double fRss = fR0 + ((alphafR * PARa + (1 - fR0)) - sqrt(pow(alphafR * PARa + (1 - fR0),2) - 4 * alphafR * thetafR * PARa * (1 - fR0))) / (2 * thetafR);
const double phi = (Kmc * Ko * O2) / (Kmo * Kc * Cc);
const double Fm = (1 - PSIId) * Fm_a + PSIId * Fm_d;
const double d_Ci_dt = (((Ca - Ci) / (1 / gsc + 1 / gbc) - (Ci - Ccyt) * gw) * Mv) / Vref;
const double gm = A / (Ci - Cc);
const double kD = (kp / (PhiIId - PhiqE) - kf) - kp;
const double fvpd = 1 / (1 + VPDleaf / D0);
const double gtw = (Trmmol * (air_pressure - (es_leaf + ea) / 2.0)) / VPDleaf;
const double d_alphar_dt = ifelse(alpharss > alphar,(alpharss - alphar) * Kialpha,(alpharss - alphar) * Kdalpha);
const double d_Ca_dt = ((((-Flow * Ca + Flow * CO2R) - leaf_surface * A) * R * Ta) / volume_chamber) / air_pressure;
const double Vc = (fRB * fRuBP * Kc * RB * Cc) / (Cc + Kmc * (1.0 + O2 / Kmo));
const double fR = fRss;
const double transpiration = ((VPDleaf / (air_pressure - (es_leaf + ea) / 2.0)) * 1) / (1 / gsw + 1 / gbw);
const double d_PSIId_dt = (1 - PSIId) * PARa * alphar * Kinh - PSIId * Krep;
const double PhiIIop = (Fm - Fo) / Fm;
const double VrTPU = (3.0 * TPU * (2.0 + 1.5 * phi)) / (1 - 0.5 * phi);
const double d_Ccyt_dt = ((((Ci - Ccyt) * gw + Rp + Rm) - (Ccyt - Cc) * gc) * Mv) / Vref;
const double fI = (-fI_b - sqrt(pow(fI_b,2) - 4 * fI_a * fI_c)) / (2 * fI_a);
const double Fmp_a = (alphar * kf) / (kf + kD);
const double VrE = (fR * Vrmax * PGA) / (PGA + KmPGA);
const double d_Cc_dt = (((Ccyt - Cc) * gc - Vc) * Mv) / Vref;
const double d_H2OS_dt = (((-((Flow + leaf_surface * transpiration)) * H2OS + Flow * H2OR + leaf_surface * transpiration) * R * Ta) / volume_chamber) / air_pressure;
const double Cond = 1 / (1 / gtw - 1 / gbw);
const double PhiIIo = PhiIIop - PhiqE;
const double J2pp = ((PhiIIop * PARaP2 + Jmax) - sqrt(pow(PhiIIop * PARaP2 + Jmax,2) - 4 * PhiIIop * theta * Jmax * PARaP2)) / (2 * theta);
const double Fop_a = (alphar * kf) / (kf + kD + kp);
const double d_PR_dt = Vc * phi - PR * kPR;
const double gss = fI * fvpd * gswm;
const double Fmp = (1 - PSIId) * Fmp_a + PSIId * Fmp_d;
const double J2pm = (((Min(VrTPU,VrE) / (1 - fpseudo / (1 - fcyc))) * 2.0) / (2.0 + 1.5 * phi)) * (2.0 + 2.0 * phi);
const double qPp = ifelse(PARaP2 > Flux0,(J2pp / PARaP2) / PhiIIop,1);
const double d_gsw_dt = ifelse(gss > gsw,(gss - gsw) * Kgsi,(gss - gsw) * Kgsd);
const double NPQ = Fm_a / Fmp - 1;
const double qPm = (J2pm / J2pp) * qPp;
const double qM = NPQ - ((Fm_a / Fmp) * alphar - 1);
const double qPno_qD = Min(qPm,qPp);
const double qE = (NPQ - qM) - qI;
const double fqEss = 1 - qPno_qD;
const double PhiqEss = (fqEss * gamma1 + fqEss * fqEss * gamma2 + fqEss * gamma3) * PhiqEmax;
const double d_fP_dt = ifelse(fqEss > fP,(fqEss - fP) * KiqEp,(fqEss - fP) * KdqEp);
const double d_fZ_dt = ifelse(fqEss > fZ,(fqEss - fZ) * KiqEz,(fqEss - fZ) * KdqEz);
const double PhiIIoss = PhiIIop - PhiqEss;
const double J2qE = ifelse(PhiIIo < PhiIIoss,(1 - (PhiIIoss - PhiIIo) / PhiIIoss) * J2pp,J2pp);
const double VrJ = (((J2qE * (1 - fpseudo / (1 - fcyc))) / 2.0) * (2.0 + 1.5 * phi)) / (2.0 + 2.0 * phi);
const double J2 = Min(J2qE,J2pm);
const double fRBss_r = Min(ifelse(PAR > Flux0,(Min(VrTPU,VrJ) / (2.0 + 1.5 * phi)) / (Vc / (fRB * fRuBP)),fRBmin),fRBmax);
const double Vr = Min(VrJ,Min(VrTPU,VrE));
const double qP = (J2 / PARaP2) / PhiIIo;
const double PhiII = ifelse(PARaP2 > Flux0,J2 / PARaP2,PhiIIo);
const double reg_limit = ifelse(abs(VrJ - Vr) < Flux0,1.0,ifelse(abs(VrTPU - Vr) < Flux0,2.0,3.0));
const double d_RuBP_dt = ((1.0 + phi) / (2.0 + 1.5 * phi)) * Vr - Vc * (1.0 + phi);
const double fRBss = Min(fRBss_nr,fRBss_r);
const double d_PGA_dt = (2.0 * Vc + 1.5 * Vc * phi) - Vr;
const double d_fRB_dt = ifelse(fRBss > fRB,(fRBss - fRB) * Krca * RCA * fRCA,(fRBss - fRB) * KdRB);

vector<double> derivatives{d_PGA_dt,d_RuBP_dt,d_fRB_dt,d_fP_dt,d_fZ_dt,d_alphar_dt,d_PSIId_dt,d_PR_dt,d_Cc_dt,d_Ccyt_dt,d_Ci_dt,d_Ca_dt,d_H2OS_dt,d_gsw_dt};
 vector<double> observed{Vr, Vc, fRuBP, PARaP, PARa, PAR, fqEss, PhiqE, PhiIIoss, PhiIIo, qP, PhiII, VrJ, NPQ, qI, qE, qM, fR, Rp, A, gss, VPDleaf, Sc, Photo, transpiration, Trmmol, Cond, gm, reg_limit};
 array<vector<double>,2> output{derivatives, observed};
 return output;
}


};

