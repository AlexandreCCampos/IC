/*#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TRandom3.h"
#include <complex> //complex lib of cpp. Root has a lib too.
*/
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TROOT.h>
#include <TNtuple.h>
#include <TComplex.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TH2Poly.h>
#include "Math/IFunction.h"
#include <cmath>
#include <iostream>
#include "Resonances.h"

//using namespace Resonances;
//using namespace ROOT;

#define POW2(x)(pow(x, 2))
/*
Resonances for DKKP - for reference
      ks_892
      ks_1430  
      phi_1020
      phi_1680
      a0_1450
      kappa
      k2s_1430
*//*
Data from PDG
      ks_892   - mass = 891.67   +- 0.26;    width = 51.4  +- 0.8
      phi_1020 - mass = 1019.461 +- 0.016;   width = 4.249 +- 0.013
      phi_1680 - mass = 1680     +- 20;      width = 150   +- 50 (???)
*/

//momentum term on amplitude's coefficient. not returning the right momentum yet.
//(-2*p1*p3)**J
double momentum(int spin){

   return pow(-2, spin);
}

// cosine between particle
//function on Kajanti - equivalent
/*double cos12(double m0, double m1, double m2, double m3, double s12, double s13){
   double s = POW2(m0);
   double s23 = s + POW2(m1) + POW2(m2) + POW2(m3) - s12 - s13;
   double a = 2*s*(POW2(m1) + POW2(m2) - s12);
   double b = (s - s23 + POW2(m1)) * (s - s13 + POW2(m2));
   double c = lambda(s, POW2(m1), s23);
   double d = lambda(s, POW2(m2), s13);
   return (a+b)/pow(c * d, 1./2);
}*/
  //cos pedro
double cos12(double m0, double m1, double m2, double m3, double s12, double s13){
   double s = POW2(m0);
   double s23 = s + POW2(m1) + POW2(m2) + POW2(m3) - s12 - s13;
   double a = 2*s13*(POW2(m1) + POW2(m2) - s12);
   double b = (s13 - POW2(m3) + POW2(m1)) * (s - s13 - POW2(m2));
   double c = lambda(s13, POW2(m1), POW2(m3));
   double d = lambda(s, s13, POW2(m2));
   return (a+b)/pow(c * d, 1./2);
}

double cos13(double m0, double m1, double m2, double m3, double s12, double s13){
   double s = POW2(m0);
   double s23 = s + POW2(m1) + POW2(m2) + POW2(m3) - s12 - s13;
   double a = 2*s12*(POW2(m1) + POW2(m3) - s13);
   double b = (s12 - POW2(m2) + POW2(m1)) * (s - s12 - POW2(m3));
   double c = lambda(s12, POW2(m2), POW2(m1));
   double d = lambda(s, s12, POW2(m3));
   return (a+b)/pow(c * d, 1./2);
}

// not using
// p and r ??
double Form_Blatt_Weisskopf(double p, double r, int J){
   if(J==0){
      return 1;
   } else if(J==1){
      return 1/pow(1+pow(r, 2)*pow(p, 2), 1./2);
   } else if(J==2){
      return 1/pow(9+3*pow(r, 2)*pow(p, 2)+pow(r, 4)*pow(p, 4), 1./2);
   }

   return 1; //just in case of error
}

/*
// not using
// Breit-Wigner's Gamma function 
double G_BW(double m0, double m12){
   double dummy = Form_Blatt_Weisskopf(p, r, J)/Form_Blatt_Weisskopf(p, r, J);//not right
   double dummy2 = pow(p/p0, 2*J+1);
   double dummy3 = m0*Form_Blatt_Weisskopf(p, r, J)/m12
}*/

// Breit-Wigner's function
//use gamma from pdg
TComplex BW(double m0, double s12, double gamma){
   double a = POW2(m0) - s12, b = m0 * gamma;
   double realBW = a/(POW2(a)+POW2(b));
   double imBW = b/(POW2(a)+POW2(b));
   return TComplex(realBW, imBW);
}

// pedro's BW
/*TComplex BW(double m0, double s12, double gamma){
   s12 = pow(s12, 1./2);
   double a = 1/(M_PI*gamma), b = POW2(gamma)/POW2(s12-m0);
   double realBW = a*b+POW2(gamma);
   double imBW = 0;
   return TComplex(realBW, imBW);
}*/

// Kallen's function 
double lambda(double x, double y, double z) {
   return pow(x-y-z,2)-4*y*z;
}

//abs val of "A" cursive sqrd on isobaric model. Expressed as a function of s_ij and s_ik.
TComplex pdf(Resonance res[], double m0, double m1, double m2, double m3, double s12, double s13, int dim_array){//m0 = mother mass
   TComplex A_total = TComplex(0, 0);
   for (int j = 0; j<dim_array; j++){
      A_total += res[j].getAmplitude(m0, m1, m2, m3, s12, s13);
   }
   return A_total;
}

//cinematics limits
//implement sqr mass as variable
//c is 1 for upper and 0 for lower
double s_12_lim(double M, double m_1, double m_2, double m_3, double s23, bool c){
   double dummy, s=pow(M,2);
   
   if(c){
      dummy = ((s23-s+pow(m_1,2))*(s23+pow(m_2,2)-pow(m_3,2))-pow(lambda(s23,s,pow(m_1,2)),1./2)*pow(lambda(s23,pow(m_2,2),pow(m_3,2)),1./2));
   }else{
      dummy = ((s23-s+pow(m_1,2))*(s23+pow(m_2,2)-pow(m_3,2))+pow(lambda(s23,s,pow(m_1,2)),1./2)*pow(lambda(s23,pow(m_2,2),pow(m_3,2)),1./2));
   }/* else{
      cout << "Not a valid cinematic limit: c must be 0 or 1." << endl;
   }*/

   return pow(m_1, 2)+pow(m_2,2)-dummy/(2*s23);
}

//dalitz_isobaric
void dalitz(int nEntries=10){
   //PDG
   double m_K493=493.67, m_pi139=139.57, m_Dplus=1869.66;
   //number of terms on pdf
   
   //implement mem alloc
   int nRes = 10, nBins =100, offset12, offset13;

   //Resonance res[] = {Resonance(891.67e-3, 1, 51.4, 1, 0, 13), Resonance(1019.461e-3, 1, 4.249, 1, 0, 12)};
   //Resonance res[] = {Resonance(891.67, 1, 51.4, 1, 0, 13)};
   Resonance res[] = {Resonance(1019.461, 1, 4.249, 1, 0, 12)};

   // KKP for now
   double s13_min = pow(m_K493+m_pi139, 2);
   double s13_max = pow(m_Dplus-m_K493, 2);
   double s12_min = pow(m_K493+m_K493, 2);
   double s12_max = pow(m_Dplus-m_pi139, 2);

   // cosmetic
   offset12 = 0.1*s12_max;
   offset13 = 0.1*s13_max;

   TH1F h_s12("h_s12", "", nBins, s12_min-offset12, s12_max+offset12);
   TH1F h_s13("h_s13", "", nBins, s13_min-offset13, s13_max+offset13);
   TH2F h_dalitz("h_dalitz", "; s13; s12", nBins, s13_min-offset13, s13_max+offset13, nBins, s12_min-offset12, s12_max+offset12);

   TComplex A_total = TComplex(1,0);

   // making a loop to get the maximum value of the pdf.
   double max_pdf[3] = {0};
   double i_s12 = 0, i_s13 = 0, ds13 = (s13_max - s13_min)/nBins, ds12 = (s12_max - s12_min)/nBins;
   for(int k = 0; k<nBins; k++){
      for(int l = 0; l<nBins; l++){
         i_s12 = s12_min + k*ds12 + ds12/2;
         i_s13 = s13_min + l*ds13 + ds13/2;
         if(i_s12<=s_12_lim(m_Dplus, m_K493, m_K493, m_pi139, i_s13, 1) && i_s12>=s_12_lim(m_Dplus, m_K493, m_K493, m_pi139, i_s13, 0)){//force not to take a bin outside the dalitz. this would return a cosine greater ou lesser than 1.
            A_total = pdf(res, m_Dplus, m_K493, m_K493, m_pi139, i_s12, i_s13, int(sizeof(res)/sizeof(*res)));
            double A_total_sqr = A_total.Rho2();

            if(A_total_sqr > max_pdf[0]){
               max_pdf[0] = A_total_sqr;
               max_pdf[1] = i_s12; //bin s12
               max_pdf[2] = i_s13; //bin s13
            }
         }
      }
   }

   TRandom3 r(1);
   cout << "Máximo valor da pdf: " << max_pdf[0] << endl;

   // Generating a toy under the cinematics limits
   int i=0;
   while(i<nEntries){

      double s12 = r.Uniform(s12_min, s12_max);
      double s13 = r.Uniform(s13_min, s13_max);

      if(s12<=s_12_lim(m_Dplus, m_K493, m_K493, m_pi139, s13, 1) && s12>=s_12_lim(m_Dplus, m_K493, m_K493, m_pi139, s13, 0)){
         
         //calculating pdf
         A_total = pdf(res, m_Dplus, m_K493, m_K493, m_pi139, s12, s13, int(sizeof(res)/sizeof(*res)));

         double A_total_sqr = A_total.Rho2();

         if(r.Uniform(0, max_pdf[0])<=A_total_sqr){ //could generate a rnd number uppon the probability

            if (i%50000 == 0) cout << "Processing event " << i << endl;

            //fill histograms
            i++;
            h_s12.Fill(s12);
            h_s13.Fill(s13);
            h_dalitz.Fill(s13,s12);
         }
      }
   }

   // Let's open a TFile
   TFile out_file("dalitz.root","RECREATE");
   
   // Write the histogram in the file
   h_s12.Write();
   h_s13.Write();
   h_dalitz.Write();

   // Close the file
   out_file.Close();

   //TCanvas c = new TCanvas();
   //h_dalitz.Draw("COLZ");

}