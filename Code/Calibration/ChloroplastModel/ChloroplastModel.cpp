#include <array>
#include <vector>
#include <math.h>

using namespace std;

extern "C" {
array<vector<double>, 2> chloroplast(const double& t, const vector<double>& states, const vector<double>& params, const vector<double>& forcs) { 
  // Parameters
  const double alphar_alpha25 = params[0];
  const double DHaAlphar = params[1];
  const double Iac = params[2];
  const double alpharac = params[3];
  const double alpharav = params[4];
  const double thetaalphar = params[5];
  const double alphar_K_i = params[6];
  const double alphar_K_d = params[7];
  const double HaKalpha = params[8];
  const double DsKalpha = params[9];  
  const double HdKalpha = params[10];
//const double falpharav = params[11];    
  const double T = params[12];
  // States
  const double alphar = states[0];
  // Forcings
  const double Ib = forcs[0];
  
  // Model
  constexpr double R = 8.31;
  constexpr double Tref = 298.15;
  const double alphar_alpha = alphar_alpha25*exp(-(T - Tref)*DHaAlphar/(Tref*R*T));
  const double alpharss = min(1.0 + Ib/Iac*alpharac, 1.0 + alpharac - (alphar_alpha*(Ib - Iac) + alpharav - 
                              sqrt(pow(alphar_alpha*(Ib - Iac) + alpharav, 2.0) - 4.0*alphar_alpha*thetaalphar*alpharav*(Ib - Iac)))/(2*thetaalphar));
  const double fkT = exp((T - Tref)*HaKalpha/(Tref*R*T))*(1.0 + exp((Tref*DsKalpha - HdKalpha)/(Tref*R)))/(1 + exp((T*DsKalpha - HdKalpha)/(T*R)));
  const double d_alphar_dt = (alpharss > alphar) ? ((alpharss - alphar)*alphar_K_i*fkT) : ((alpharss - alphar)*alphar_K_d*fkT);
    
  // Return derivatives
  vector<double> derivatives{d_alphar_dt};
  vector<double> observed{};
  array<vector<double>,2> output{derivatives, observed};
  return output;
}
}