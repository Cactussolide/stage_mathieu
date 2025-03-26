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
#include <iostream>
#include <fstream>
#include <cmath>
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


       double cross_section_holo_2(double t, double Epho, double D_0, double m_A, double m_D, int expo = 3, double A_0 = 0.414, double N = 2.032, double Mn = 0.935, double Mv = 0) {
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
   };





// Generateur 
void PhaseSpace2() {


   //paramètres distrib E_gamma 
   double Eb = 10.6;   
   double Q2_max = 1; 
   double d = 5.0;     
   double X0 = 929.0;    

   // plage de calcul de la distrib
   double Eg_min = 0.1;
   double Eg_max = 10.6;
   int points = 1000;

   // fichier de sortie pour la distrib
   std::ofstream file("data.dat");

   // Générer les données pour N_EPA, N_Brem et leur somme dans un txt (pour tracer sur python apres)
   for (int i = 0; i <= points; ++i) {

       double Eg = Eg_min + i * (Eg_max - Eg_min) / points; // Calcul de Eg a chaque iteration

       // Calcul des valeurs pour N_EPA, N_Brem et leur somme
       double n_epa = N_EPA(Eb, Eg, Q2_max); // poids  EPA
       double n_brem = N_Brem(Eg, Eb, d, X0); // poids  Brem
       double n_sum = n_epa + n_brem; // somme poids 

       file << Eg << "\t" << n_epa << "\t" << n_brem << "\t" << n_sum << "\n";
   }

   file.close();

   std::cout << "Les données ont été enregistrées dans data.dat" << std::endl;
   


   // Définition du proton cible au repos
   TLorentzVector target(0.0, 0.0, 0.0, 0.938);  
  
   Double_t masses[2] = {0.938, 1.019}; // masse proton et phi

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

   Physics phys;

   // paramettres section efficace

   double D_0 = -1.52;
   double m_A = 1.612;
   double m_D = 1.206;
   double A_O = 0.430;
   double N = 2.032;

   for (Int_t n = 0; n < 100000; n++) {



       Double_t E_gamma = gRandom->Uniform(1.57,10.6);
       Double_t Poids = (N_EPA(Eb, E_gamma, Q2_max) + N_Brem(E_gamma, Eb, d, X0));
       TLorentzVector beam(0.0, 0.0, E_gamma, E_gamma);
       TLorentzVector W = beam + target;   
       event.SetDecay(W, 2, masses); 
       Double_t weight = event.Generate();

       TLorentzVector *pProton = event.GetDecay(0); // Proton final
       TLorentzVector *pPhi    = event.GetDecay(1); //  phi

       Double_t t = (beam - *pPhi).M()*(beam - *pPhi).M();  //calcul de t d'une maniere (attention renvoi |t| et pas t)
       Double_t t_3 = 1.019*1.019 - 2*(E_gamma*pPhi->E() - E_gamma*pPhi->Pz()); // calcul de t d'une autre maniere (renvoi bien t)

       double weight_crosssection = phys.cross_section_holo_2(t_3, E_gamma, D_0, m_A, m_D);


       Double_t E = E_gamma - pPhi->E();
       Double_t px = 0 - pPhi->Px();
       Double_t py = 0 - pPhi->Py();
       Double_t pz = E_gamma - pPhi->Pz();

       Double_t t_2 = E*E - px*px -py*py - pz*pz; //calcul de t d'une 3eme maniere 

       // verification de la conservation de l'énergie
       Double_t E_ini = E_gamma + 0.938;
       Double_t E_fin = pPhi->E() + pProton->E();

       // verification que somme impulsion ini = somme impulsion final
       // normalement tout les P_diff valent 0

       Double_t Px_diff = 0 + 0 - pProton->Px() - pPhi->Px();
       Double_t Py_diff = 0 + 0 - pProton->Py() - pPhi->Py();
       Double_t Pz_diff = E_gamma + 0 - pProton->Pz() - pPhi->Pz();

       // print plein de choses pour vérif

       if (n % 1000 == 0) {
        std::cout << "Événement " << n << " : E_gamma = " << E_gamma << std::endl;
        std::cout << "Événement " << n << " : m = " << pPhi->M() << std::endl;
        std::cout << "Calcul de t " << std::endl;
        std::cout << "Événement " << n << " : t = " << t << std::endl;
        std::cout << "Événement " << n << " : t2 = " << t_2 << std::endl;
        std::cout << "Événement " << n << " : t3 = " << t_3 << std::endl;
        std::cout << "Verification conservation de l'énergie " << std::endl;
        std::cout << "Événement " << n << " : E ini = " << E_ini << std::endl;
        std::cout << "Événement " << n << " : E fini = " << E_fin << std::endl;
        std::cout << "Verification somme impulsion ini = somme impulsion fini " << std::endl;
        std::cout << "Événement " << n << " : Px_diff = " << Px_diff << std::endl;
        std::cout << "Événement " << n << " : Py_diff = " << Py_diff << std::endl;
        std::cout << "Événement " << n << " : Pz_diff = " << Pz_diff << std::endl;
        cout << "Section efficace : " << weight_crosssection << " GeV^-2" << endl;
        std::cout << "pPhi: (E=" << pPhi->E() 
        
          << ", px=" << pPhi->Px() 
          << ", py=" << pPhi->Py() 
          << ", pz=" << pPhi->Pz() 
          << ", M=" << pPhi->M() << ")" << std::endl;

       }

       // histo fill en pondérant par Poids (celui des E_gamma) et weigh_crosssection (qui vient du modèle de section efficace)

       // j'ai pas mis les poids "weight" car pour l'instant ils valent 1.

       h2D->Fill(E_gamma, t, Poids*weight_crosssection);

       hE->Fill(E_gamma, Poids);


       hMphi->Fill(pPhi->M()); 


       // deuxieme decay : désintégration du phi en e+ e-

       Double_t masses_elec[2] = { 0.000511, 0.000511 };
       TGenPhaseSpace decay;
       decay.SetDecay(*pPhi, 2, masses_elec);

       Double_t weight2 = decay.Generate();

       TLorentzVector *pElectron = decay.GetDecay(0);
       TLorentzVector *pPositron = decay.GetDecay(1);
       

       Double_t poid_total = weight*weight2;

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

       h2D2->Fill(theta_phi*(180/3.14), P_phi, Poids*weight_crosssection);
       h2D3->Fill(theta_elec*(180/3.14), P_elec, Poids*weight_crosssection);
       h2D4->Fill(theta_posi*(180/3.14), P_posi, Poids*weight_crosssection);
       h2D5->Fill(theta_proton*(180/3.14), P_proton, Poids*weight_crosssection);


       // bizarre les poids weight et weight2 valent tout le temps 1

       if (n % 1000 == 0) {

         std::cout << "Événement " << n << " : poids = " << weight << std::endl;

         std::cout << "Événement " << n << " : poids deuxieme decay = " << weight2 << std::endl;
      
       }


   }


   // ajustement graphique des histos

   h2D->GetXaxis()->SetTitle("E_gamma (GeV)");  
   h2D->GetYaxis()->SetTitle("|t|");      

   h2D2->GetXaxis()->SetTitle("Theta (#degree)"); 
   h2D2->GetYaxis()->SetTitle("p");        

   h2D3->GetXaxis()->SetTitle("Theta (#degree)"); 
   h2D3->GetYaxis()->SetTitle("p");        

   h2D4->GetXaxis()->SetTitle("Theta (#degree)"); 
   h2D4->GetYaxis()->SetTitle("p");   

   h2D5->GetXaxis()->SetTitle("Theta (#degree)");  
   h2D5->GetYaxis()->SetTitle("p");        


   TCanvas *c1 = new TCanvas("c1", "Masse du Phi", 800, 600);
   hMphi->Draw();

   TCanvas *c2 = new TCanvas("c2", "Masse e+e-", 800, 600);
   hMee->Draw();

   TCanvas *c3 = new TCanvas("c3", "Histogramme |t| et E gamma", 800, 600);
   h2D->Draw("COLZ");  // "COLZ" palette de couleurs

   h2D->SetMinimum(0);  // échelle de couleur
   h2D->SetMaximum(0.001);

   TCanvas *c4 = new TCanvas("c4", "Histogramme E_gamma", 800, 600);
   hE->Draw();  

   TCanvas *c5 = new TCanvas("c5", "Histogramme p et theta pour le phi", 800, 600);
   h2D2->Draw("COLZ");  

   h2D2->SetMinimum(0);  
   h2D2->SetMaximum(0.001);


   TCanvas *c6 = new TCanvas("c6", "Histogramme p et theta pour e-", 800, 600);
   h2D3->Draw("COLZ");  

   h2D3->SetMinimum(0);  
   h2D3->SetMaximum(0.001);


   TCanvas *c7 = new TCanvas("c7", "Histogramme p et theta pour e+", 800, 600);
   h2D4->Draw("COLZ");  

   h2D4->SetMinimum(0);  
   h2D4->SetMaximum(0.001);


   TCanvas *c8 = new TCanvas("c8", "Histogramme p et theta pour le proton", 800, 600);
   h2D5->Draw("COLZ"); 
   
   h2D5->SetMinimum(0);  
   h2D5->SetMaximum(0.001);

}
