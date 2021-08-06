#ifndef RESONANCES_H
#define RESONANCES_H
#endif

double momentum(int spin);
double cos12(double m0, double m1, double m2, double m3, double s12, double s13);
double Form_Blatt_Weisskopf(double p, double r, int J);
TComplex BW(double m0, double s12, double gamma);
double lambda(double x, double y, double z);
double s_12_lim(double M, double m_1, double m_2, double m_3, double s23, bool c);
double cos13(double m0, double m1, double m2, double m3, double s12, double s13);

// Ressonance
class Resonance {
public:
    double mass;
    int spin;
    int pair; //need to know which pair the resonance represents in the decay for calculate the cosine - not sure if its needed. Must receive 12, 13 or 23.
    double width;
    double magnitude; //magnitude
    double phase;

    Resonance(double m, double J, double w, double A, double phi, int p){
        mass = m;
        spin = J;
        magnitude = A;
        phase = phi;
        width = w;
        pair = p;
    }

    // Returns the coef for this resonance to calculate the pdf
    TComplex getAmplitude(double m0, double m1, double m2, double m3, double s12, double s13){
        double F_D = 1., F_R = 1.; //for now
        TComplex a;
        if(pair == 12){
            //a = TComplex(F_D * F_R) * BW(mass, s12, width);
            a = TComplex(F_D * F_R * ROOT::Math::legendre(spin, cos13(m0, m1, m2, m3, s12, s13)), 0) * BW(mass, s12, width);
            //cout << "spin = " << spin << endl;
            //cout << "cos12(m0, m1, m2, m3, s12, s13) = " << cos12(m0, m1, m2, m3, s12, s13) << endl;
            //cout << "ROOT::Math::legendre(spin, cos12(m0, m1, m2, m3, s12, s13)), 0)  = " << ROOT::Math::legendre(spin, cos12(m0, m1, m2, m3, s12, s13))  << endl;
        } else if (pair == 13){
            a = TComplex(F_D * F_R * ROOT::Math::legendre(spin, cos12(m0, m1, m2, m3, s12, s13)), 0) * BW(mass, s13, width);
            //cout << "F_D * F_R * ROOT::Math::legendre(spin, cos12(m0, m1, m2, m3, s12, s13)) = " << F_D * F_R * ROOT::Math::legendre(spin, cos12(m0, m1, m2, m3, s12, s13)) << endl;
            //cout << "BW = " << BW(mass, s12, width) << endl;
            //cout << "a = " << a << endl;
        }

        return a*TComplex(magnitude, 0);
    }

};