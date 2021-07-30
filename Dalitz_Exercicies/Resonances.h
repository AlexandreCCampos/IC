#ifndef RESONANCES_H
#define RESONANCES_H
#endif

double momentum(int spin);
double cos12(double m0, double m1, double m2, double m3, double s12, double s13);
double Form_Blatt_Weisskopf(double p, double r, int J);
TComplex BW(double m0, double s12, double gamma);
double lambda(double x, double y, double z);
double s_12_lim(double M, double m_1, double m_2, double m_3, double s23, bool c);

// Ressonance
class Resonance {
public:
    double mass;
    int spin;
    int pair; //need to know which pair the resonance represents in the decay for calculate the cosine - not sure if its needed. Must receive 12, 13 or 23.
    double width;
    double amplitude;
    double phase;

    Resonance(double m, double J, double w, double A, double phi, int p){
        mass = m;
        spin = J;
        amplitude = A;
        phase = phi;
        width = w;
        pair = p;
    }

    // Returns the coef for this resonance to calculate the pdf
    TComplex getCoeficient(double m0, double m1, double m2, double m3, double s12, double s13){
        double F_D = 1., F_R = 1.; //for now

        //taking off legendre pol because its not working and for spin 1, legpoly is x.
        //double a = F_D * F_R * momentum(spin) * std::legendre(spin, cos13(m0, m1, m2, m3, s12, s13)) * BW(m0, s12, width);
        
        TComplex a = TComplex(F_D * F_R * cos12(m0, m1, m2, m3, s12, s13), 0) * BW(mass, s12, width);
        return a;
    }

};