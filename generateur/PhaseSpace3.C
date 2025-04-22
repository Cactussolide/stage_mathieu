#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include <array>
#include <iostream>
#include <vector>
#include <tuple>
#include <TString.h>
#include <fstream>
#include <cmath>
#include <iomanip>
using namespace std;



// Distribution E_gamma

double N_EPA(double Eb, double Eg, double Q2_max) {
   const double alpha = 1.0 / 137.0;
   const double PI = 3.141592653589793;


   double x = Eg / Eb;
   double me = 0.00051;
   double Q2_min = me * me * x * x / (1 - x);


   return (1 / Eb) * alpha / (PI * x) * ((1 - x + x * x / 2) * log(Q2_max / Q2_min) - (1 - x));
}

double N_Brem(double Eg, double Eb, double d, double X0) {
   return (0.5 * d / X0) * (1 / Eg) * ((4.0 / 3.0) - (4.0 / 3.0) * (Eg / Eb) + Eg * Eg / (Eb * Eb));
}

// Classe section efficace

class Physics {


   public:


      double cross_section_holo_2(double t, double Epho, double D_0, double m_A, double m_D, int expo = 3, double A_0 = 0.430, double N = 2.032, double Mn = 0.935, double Mv = 1.019) {
         // Paramètre de normalisation
         double N_e = N;  // En nb GeV^(-2)
   
         // Calcul de la variable de Mandelstam s
         double s = pow(Mn, 2) + 2 * Mn * Epho;
   
         // Facteur tilde F
         double tilde_F = pow((s - pow(Mn, 2)) / 2.0, 4);
   
         // Paramètre eta
         double eta = 0;
   
         // Facteurs A(t) et D(t)
         double A_t = A_0 / pow(1 - (t / pow(m_A, 2)), expo);
         double D_t = D_0 / pow(1 - (t / pow(m_D, 2)), expo);
   
         // Calcul de la section efficace
         double sigma = pow(N_e, 2) * (1.0 / (64 * M_PI * pow(s - pow(Mn, 2), 2))) *
                        pow((A_t + pow(eta, 2) * D_t) / A_0, 2) * tilde_F * 8;
   
         return sigma;
      }      


      double Lambda(double x, double y, double z) {
         return (x - y - z)*(x - y - z) - 4 * y*z;
      }


       //From Byukling Kayanti Formula (5.14) Page 86

      double T_min(double ma_2, double mb_2, double m1_2, double m2_2, double s) // arguments are squares of masses of particles in the reaction a+b->1+2, and s is the square of the total c.m. energy i.e. (a+b)^2
      {
         return ma_2 + m1_2 - (1 / (2 * s))*((s + ma_2 - mb_2)*(s + m1_2 - m2_2) - sqrt(Lambda(s, ma_2, mb_2) * Lambda(s, m1_2, m2_2)));
      }

      //From Byukling Kayanti Formula (5.14) page 86

      double T_max(double ma_2, double mb_2, double m1_2, double m2_2, double s) {
         return ma_2 + m1_2 - (1 / (2 * s))*((s + ma_2 - mb_2)*(s + m1_2 - m2_2) + sqrt(Lambda(s, ma_2, mb_2) * Lambda(s, m1_2, m2_2)));
      }

   };



// Generateur 
void PhaseSpace3(double nb_fichier = 1000, double nb_event = 1000, double Eb = 10.6) {


   //paramètres distrib E_gamma EPA et BREM
  
   double Q2_max = 1; 
   double d = 5.0;     
   double X0 = 929.0;

   // plage de calcul de la distrib E_gamma EPA et BREM

   double Eg_min = 1.58;
   double Eg_max = Eb;

   // paramettres section efficace HOLO

   double D_0 = -1.52;
   double m_A = 1.612;
   double m_D = 1.206;
   double A_O = 0.430;
   double N = 2.032;

   // particules
   double Mp = 0.938;
   double Mv = 1.019;
   double me = 0.00051;

   // variables utiles
   Double_t vz = 0; // indice position vertex
   double BR = 2.9*1e-4; // branching ratio

   //variable plot
   double echelle_plot_couleur = 0.00001;

   // NaN
   double counter_nan = 0;


   Physics phys;
   

   // Définition du proton cible au repos
   TLorentzVector target(0.0, 0.0, 0.0, Mp);  
  
   Double_t masses[2] = {Mp, Mv}; // masse proton et phi

   TGenPhaseSpace event;

   // initialisation des histo

   TH1F *hMphi = new TH1F("hMphi", "mass of phi", 100, 0, 1.05);
   TH1F *hMee = new TH1F("hMee", "invariant mass of e+ e-", 100, 0.9, 1.1);
   TH1F *hE = new TH1F("hE", "Histogram E_gamma", 100, 1.50, 12);


   TH2F *h2D = new TH2F("h2D", "Histogram |t| and E gamma", 600, 0, 12, 1000, -0.1, 20);  
   TH2F *h2D2 = new TH2F("h2D2", "Histogram p and theta for phi", 1000, 0, 90, 500, -0.1, 11); 
   TH2F *h2D3 = new TH2F("h2D3", "Histogram p and theta for e-", 1000, 0, 180, 500, -0.1, 11); 
   TH2F *h2D4 = new TH2F("h2D4", "Histogram p and theta for e+", 1000, 0, 180, 500, -0.1, 11);  
   TH2F *h2D5 = new TH2F("h2D5", "Histogram p and theta for proton", 1000, 0, 90, 500, -0.1, 11);  


   std::vector<std::ofstream> files(nb_fichier);

   
   for (int i = 0; i < nb_fichier; ++i) {

       // Création et ouverture des fichiers

       std::string filename = "PhiGen_" + std::to_string(i) + ".txt";
       files[i].open(filename);

       if (!files[i]) {
           std::cerr << "Erreur lors de l'ouverture du fichier : " << filename << std::endl;
       }


       for (Int_t n = 0; n < nb_event; n++) {

         // Remplir un fichier avec nb_event

         double indice_evenement = n;

         Double_t E_gamma = gRandom->Uniform(Eg_min, Eg_max);
         Double_t Poids_gamma = (N_EPA(Eb, E_gamma, Q2_max) + N_Brem(E_gamma, Eb, d, X0));
         TLorentzVector beam(0.0, 0.0, E_gamma, E_gamma);
         TLorentzVector W = beam + target;   
         event.SetDecay(W, 2, masses);
         Double_t weight = event.Generate();
  
         TLorentzVector *pProton = event.GetDecay(0); // Proton final
         TLorentzVector *pPhi    = event.GetDecay(1); //  phi
  
         Double_t t_3 = Mv*Mv - 2*(E_gamma*pPhi->E() - E_gamma*pPhi->Pz()); // calcul de t 
  
         double weight_crosssection = phys.cross_section_holo_2(t_3, E_gamma, D_0, m_A, m_D); //dsigma/dt
  
         double s = Mp*Mp + 2*Mp*E_gamma;
  
         double t_min = phys.T_min(0., Mp*Mp, Mv*Mv, Mp*Mp, s);
         double t_max = phys.T_max(0., Mp*Mp, Mv*Mv, Mp*Mp, s);
  
         double correction_weight = abs(t_min - t_max); // corection weight

         double weight_range_photon = Eg_max - Eg_min; // correction tirage range photon
  
         double tot_weight = Poids_gamma*correction_weight*weight_crosssection*BR*weight_range_photon; // poid total
  
         h2D->Fill(E_gamma, abs(t_3), tot_weight);
  
         hE->Fill(E_gamma, Poids_gamma);
         hMphi->Fill(pPhi->M()); 
  
  
         // deuxieme decay : désintégration du phi en e+ e-
  
         Double_t masses_elec[2] = { me, me };
         TGenPhaseSpace decay;
         decay.SetDecay(*pPhi, 2, masses_elec);
         Double_t weight2 = decay.Generate();
  
         TLorentzVector *pElectron = decay.GetDecay(0);
         TLorentzVector *pPositron = decay.GetDecay(1);
         
  
         // masse invariante e+e-
  
         TLorentzVector Mee = *pElectron + *pPositron;
         hMee->Fill(Mee.M());
  
         // histo p et thetha
  
         Double_t theta_phi = pPhi->Theta();
         Double_t theta_elec = pElectron->Theta();
         Double_t theta_posi = pPositron->Theta();
         Double_t theta_proton = pProton->Theta();
  
  
         Double_t P_phi = pPhi->P();
         Double_t P_elec = pElectron->P();
         Double_t P_posi = pPositron->P();
         Double_t P_proton = pProton->P();
  
         // fill les autres histo
  
         h2D2->Fill(theta_phi*(180/3.14), P_phi, tot_weight);
         h2D3->Fill(theta_elec*(180/3.14), P_elec, tot_weight);
         h2D4->Fill(theta_posi*(180/3.14), P_posi, tot_weight);
         h2D5->Fill(theta_proton*(180/3.14), P_proton, tot_weight);
  
         if (n % 1000 == 0) {
  
           std::cout << "Événement " << n << " : E_gamma = " << E_gamma << std::endl;
           std::cout << "Événement " << n << " : t3 = " << t_3 << std::endl;
           cout << "Section efficace : " << weight_crosssection << " GeV^-2" << endl;
        
         }

         // il y a t il des NaN ? 

         if (std::isnan(pPhi->Px()) || std::isnan(pPhi->Py()) || std::isnan(pPhi->Pz()) 
            || std::isnan(pProton->Px()) || std::isnan(pProton->Py()) || std::isnan(pProton->Pz())
            || std::isnan(pElectron->Px()) || std::isnan(pElectron->Py()) || std::isnan(pElectron->Pz())
            || std::isnan(pPositron->Px()) || std::isnan(pPositron->Py()) || std::isnan(pPositron->Pz())) {


            std::cerr << "Une des composantes est NaN" << std::endl;
            std::cout << "Événement " << n << " concerne avec E_gamma = " << E_gamma << std::endl;
            counter_nan += 1;

         }

         
         // Ecriture fichier

         //header
         files[i] << 3 << setw(15) << 1 << setw(5) << 1 << setw(15) << 0 << setw(15) << 0 << setw(15) << 22 << setw(15) << E_gamma << setw(15) << 0 << setw(15) << indice_evenement << setw(15) << tot_weight << endl;
  
         // proton
         files[i] << 1 << setw(5) << 1 << setw(5) << 1 << setw(7) << 2212 << setw(5) << 0 << setw(5) << 0 << setw(15) << pProton->Px() << setw(15) << pProton->Py() << setw(15) << pProton->Pz();
         files[i] << setw(15) << pProton->E() << setw(15) << Mp << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << endl;
  
         // electron
         files[i] << 2 << setw(5) << -1 << setw(5) << 1 << setw(7) << 11 << setw(5) << 0 << setw(5) << 0 << setw(15) << pElectron->Px() << setw(15) << pElectron->Py() << setw(15) << pElectron->Pz();
         files[i] << setw(15) << pElectron->E() << setw(15) << me << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << endl;
  
         // positron
         files[i] << 3 << setw(5) << 1 << setw(5) << 1 << setw(7) << -11 << setw(5) << 0 << setw(5) << 0 << setw(15) << pPositron->Px() << setw(15) << pPositron->Py() << setw(15) << pPositron->Pz();
         files[i] << setw(15) << pPositron->E() << setw(15) << me << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << endl;
  
  
     }

     files[i].close();

   }

   // ajustement graphique des histos


   hE->GetXaxis()->SetTitle("E_gamma [GeV]");  
   hE->GetYaxis()->SetTitle("weighted distribution");

   h2D->GetXaxis()->SetTitle("E_gamma [GeV]");  
   h2D->GetYaxis()->SetTitle("|t| [GeV^{2}]");      

   h2D2->GetXaxis()->SetTitle("#theta [degree]"); 
   h2D2->GetYaxis()->SetTitle("p [GeV]");        

   h2D3->GetXaxis()->SetTitle("#theta [degree]"); 
   h2D3->GetYaxis()->SetTitle("p [GeV]");        

   h2D4->GetXaxis()->SetTitle("#theta [degree]"); 
   h2D4->GetYaxis()->SetTitle("p [GeV]");   

   h2D5->GetXaxis()->SetTitle("#theta [degree]");  
   h2D5->GetYaxis()->SetTitle("p [GeV]");   

   TCanvas *c1 = new TCanvas("c1", "Masse du Phi", 800, 600);
   hMphi->Draw();

   TCanvas *c2 = new TCanvas("c2", "Masse e+e-", 800, 600);
   hMee->Draw();

   TCanvas *c3 = new TCanvas("c3", "Histogramme |t| et E gamma", 800, 600);
   h2D->Draw("COLZ");  // "COLZ" palette de couleurs

   h2D->SetMinimum(0);  // échelle de couleur
   h2D->SetMaximum(echelle_plot_couleur);

   TCanvas *c4 = new TCanvas("c4", "Histogramme E_gamma", 800, 600);
   hE->Draw();  

   TCanvas *c5 = new TCanvas("c5", "Histogramme p et theta pour le phi", 800, 600);
   h2D2->Draw("COLZ");  

   h2D2->SetMinimum(0);  
   h2D2->SetMaximum(echelle_plot_couleur);


   TCanvas *c6 = new TCanvas("c6", "Histogramme p et theta pour e-", 800, 600);
   h2D3->Draw("COLZ");  

   h2D3->SetMinimum(0);  
   h2D3->SetMaximum(echelle_plot_couleur);


   TCanvas *c7 = new TCanvas("c7", "Histogramme p et theta pour e+", 800, 600);
   h2D4->Draw("COLZ");  

   h2D4->SetMinimum(0);  
   h2D4->SetMaximum(echelle_plot_couleur);


   TCanvas *c8 = new TCanvas("c8", "Histogramme p et theta pour le proton", 800, 600);
   h2D5->Draw("COLZ"); 
   
   h2D5->SetMinimum(0);  
   h2D5->SetMaximum(echelle_plot_couleur);

   std::cout << nb_fichier << " fichiers de " << nb_event << " events generés" << std::endl;
   std::cout << counter_nan << " evenements corrompus" << std::endl;
}
