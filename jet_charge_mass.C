#include "TFile.h"
#include "TSystem.h"
#include "TPythia8.h"
#include "Pythia8/Pythia.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TTree.h"
#include "TError.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "THStack.h"
#include <TLorentzVector.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#ifndef __CINT__
#include "TCanvas.h"
#include "TStyle.h"
#include "fastjet/ClusterSequence.hh"
#include "TLegend.h"

using namespace std;
#endif
using namespace Pythia8;
using namespace fastjet;

double calculate_energy(double px,double py,double pz,double m = 0.139570)
{

  double e = sqrt(pow(m,2)+pow(px,2)+pow(py,2)+pow(pz,2));

  return e;
}

void jet_charge_mass()
{
  int xBins = 6;
  int yBins = 10;

  TH2* jet_c_m = new TH2D("Statistics","PYTHIA pp 200 GeV", 12, 0., 6., 8, -2., 2.);

  TH1* jetcharge = new TH1D("","PYTHIA pp 200 GeV",yBins, -2.5,2.5);
  TH1* jetmass = new TH1D("","PYTHIA pp 200 GeV", xBins, 0., 6.);

  Pythia8::Pythia pythia;

  Settings& settings = pythia.settings;
  Pythia8::Info& info = pythia.info;


  pythia.readString("HardQCD:all = on");
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 200.");
  pythia.readString("PhaseSpace:pTHatMin = 20.");

  pythia.init();

//  vector<PseduoJet> injets;

  for (int iEv = 0; iEv <10000; ++iEv)
  {
    if (!pythia.next()) continue;
    vector<PseudoJet> fjInputs;

    for (int i = 0; i < pythia.event.size(); ++i)
    {
       if (pythia.event[i].isFinal() && abs(pythia.event[i].eta())<1 && pythia.event[i].pT()>1 && pythia.event[i].isCharged())
       {
         double energy = calculate_energy(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz());
         //double energy = pythia.event[i].e();

         Int_t charge = pythia.event[i].charge();
         fastjet::PseudoJet particle_Track(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), energy);
         particle_Track.set_user_index(charge);
         fjInputs.push_back(particle_Track);
       }
     }

    double R = 0.4;

    if (fjInputs.size()==0) continue;

  //  injets.clear();

    JetDefinition jet_def(antikt_algorithm, R);
    ClusterSequence cluseq(fjInputs, jet_def);

    vector<PseudoJet> jets = sorted_by_pt(cluseq.inclusive_jets());

    Int_t N = jets.size();


    for (unsigned i = 0; i < jets.size(); i++)
    {

      double jet_energy = jets[i].E();
      double jet_pt = jets[i].pt();
      double jet_pz = jets[i].pz();

      double jet_mass = jets[i].m();



      if (jet_pt>10. && jet_pt<80. && abs(jets[i].eta())<0.6 && jet_mass>0.)
      {
        vector <PseudoJet> constituents = jets[i].constituents();
        vector <PseudoJet> sortedconstituents = sorted_by_pt(constituents);

        jetmass->Fill(jet_mass);
        double chargesum[] = {0.0};

        for (unsigned j = 0; j < sortedconstituents.size(); j++)
        {
          if (sortedconstituents.size() == 1) continue;

          Int_t track_charge = sortedconstituents[j].user_index();

            Float_t z = (sortedconstituents[j].pt()) / jet_pt;
            Float_t jet_charge_5 = (pow(z,0.5) * track_charge);
            chargesum[i] += jet_charge_5;

          }
          for (unsigned k = 0; k < jets.size(); k++)
          {
            jet_c_m->Fill(jet_mass,chargesum[k]);
            jetcharge->Fill(chargesum[k]);
          }
        }
      }

    }

TCanvas* c1 = new TCanvas("c1","inclusive jet charge vs mass");
c1->cd();
jet_c_m->Draw("SURF3");
jet_c_m->SetStats(0);


TCanvas* c2 = new TCanvas("c2","inclusive jet mass");
c2->cd();
jetmass->Draw();
jetmass->SetStats(0);
jetmass->SetFillColor(kGray);
jetmass->SetYTitle("(1/N_{jets}) dN/dM_{j}");
jetmass->SetXTitle("Charged Jet Mass (GeV/c^{2})");
jetmass->Scale(1.0/((jetmass->GetEntries())*(jetmass->GetBinWidth(xBins))));


TCanvas* c3 = new TCanvas("c3","inclusive jet charge");
c3->cd();
jetcharge->Draw();
jetcharge->SetStats(0);
jetcharge->SetFillColor(kGray);
jetcharge->SetYTitle("Normalized");
jetcharge->SetXTitle("Jet Charge #kappa = 0.5");
jetcharge->Scale(1./(jetcharge->GetEntries()));


}
