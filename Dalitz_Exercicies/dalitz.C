#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TRandom3.h"


/*
Resonances for DKKP - for reference
      ks_892
      ks_1430  
      phi_1020
      phi_1680
      a0_1450
      kappa
      k2s_1430
*/
/*
Make a .h file with these function

*/

//p and r ??
double Form_Blatt_Wieskopf(double p, double r, int J){
   if(J==0){
      return 1;
   } else if(J==1){
      return 1/pow(1+pow(r, 2)*pow(p, 2), 1./2);
   } else if(J==2){
      return 1/pow(9+3*pow(r, 2)*pow(p, 2)+pow(r, 4)*pow(p, 4), 1./2);
   }

}

// Brit-Wigner's Gamma function 
double G_BW(double m0, double m12){

}

// Briet-Wigner's function
double BW(){
   1/(pow(m0, 2)-pow(m12, 2)-IMAGINARY * m0 * G_BW(m0, m12));
}
//ROOT::Math::legendre([0],x)

//Kallen's function 
double lambda(double x, double y, double z) {
   return pow(x-y-z,2)-4*y*z;
}


/* // not using
//invariant mass sqr
double s_ij(double p1[], double p2[]){
   double m_12 = p1+p2;
   return pow(m12[0],2)-pow(m12[1],2)-pow(m12[2],2)-pow(m12[3],2);
}*/

//"A" cursive on isobaric model. Expressed as a function of s_ij and s_ik.
/*double pdf(){
   return ;
}

double pdf_coef(){
   return ; 
}
*/


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

void dalitz_isobaric(int nEntries=1000){
   //PDG
   double m_K493=493.67, m_pi139=139.57, m_Dplus=1869.66;

   //KKP for now
   double s23_min = pow(m_K493+m_pi139, 2);
   double s23_max = pow(m_Dplus-m_K493, 2);
   double s12_min = pow(m_K493+m_K493, 2);
   double s12_max = pow(m_Dplus-m_pi139, 2);

   TH1F h_s12("h_s12","",100, s12_min, s12_max);
   TH1F h_s23("h_s23","",100, s23_min, s23_max);
   TH2F h_dalitz("h_dalitz","",100, s12_min, s12_max, 100, s23_min, s23_max);

   TRandom3 r(1);
   
   // Generating a toy under the cinematics limits
   for(int i=0; i<nEntries; i++){

      double s12 = r.Uniform(s12_min, s12_max);
      double s23 = r.Uniform(s23_min, s23_max);

      if(s12<=s_12_lim(m_Dplus, m_K493, m_K493, m_pi139, s23, 1) && s12>=s_12_lim(m_Dplus, m_K493, m_K493, m_pi139, s23, 0)){
         //fill histograms
         h_s12.Fill(s12);
         h_s23.Fill(s23);
         h_dalitz.Fill(s12,s23);
      }
   }

   //not working - look for syntax
   //h_dalitz.Fill(h_s12,h_s23);

   // Let's open a TFile
   TFile out_file("dalitz.root","RECREATE");
   
   // Write the histogram in the file
   h_s12.Write();
   h_s23.Write();
   h_dalitz.Write();

   // Close the file
   out_file.Close();

   //TCanvas c = new TCanvas();
   //h_dalitz.Draw("COLZ");

}

void dalitz(int nEntries=1000){
   // PDG
   double m_K493=493.67, m_pi139=139.57, m_Dplus=1869.66;

   // KKP for now
   double s23_min = pow(m_K493+m_pi139, 2);
   double s23_max = pow(m_Dplus-m_K493, 2);
   double s12_min = pow(m_K493+m_K493, 2);
   double s12_max = pow(m_Dplus-m_pi139, 2);

   TH1F h_s12("h_s12","",100, s12_min, s12_max);
   TH1F h_s23("h_s23","",100, s23_min, s23_max);
   TH2F h_dalitz("h_dalitz","",100, s12_min, s12_max, 100, s23_min, s23_max);

   TRandom3 r(1);
   
   // Generating a toy under the cinematics limits
   for(int i=0; i<nEntries; i++){

      double s12 = r.Uniform(s12_min, s12_max);
      double s23 = r.Uniform(s23_min, s23_max);

      if(s12<=s_12_lim(m_Dplus, m_K493, m_K493, m_pi139, s23, 1) && s12>=s_12_lim(m_Dplus, m_K493, m_K493, m_pi139, s23, 0)){
         //fill histograms
         h_s12.Fill(s12);
         h_s23.Fill(s23);
         h_dalitz.Fill(s12,s23);
      }
   }

   //not working - look for syntax
   //h_dalitz.Fill(h_s12,h_s23);

   // Let's open a TFile
   TFile out_file("dalitz.root","RECREATE");
   
   // Write the histogram in the file
   h_s12.Write();
   h_s23.Write();
   h_dalitz.Write();

   // Close the file
   out_file.Close();

   //TCanvas c = new TCanvas();
   //h_dalitz.Draw("COLZ");

}


/*
MC - Dalitz plot

pt1 - Gerando sem as ressonâncias.

Definir o limites do dalitz de acordo com as massas do decaimento.

Criar os histogramas que serão preechidos.

Começo de um loop.

   Sortear as variaveis s12 e s23, que são as massas invariante ao quadrado do pares 12 e 23, dentro dos limites de maneira uniforme.

   Calcular de acordo com o limite cinemático se o ponto no dalitz está dentro o fora da área. O limite das bordas do dalitz é dado por:

   Se está dentro, preecher os histogramas dos pares s12, s23 e do dalitz.

Fim do loop

Salvar em um arquivo .root.

pt2 - Adicionando as ressonancias.

*/