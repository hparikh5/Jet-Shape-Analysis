// includes
#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"

#ifndef __CINT__
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "Riostream.h"
#include <cstdlib>
#include "TH3F.h"
#include "TH2F.h"
#include "TMath.h"
//#include "TMCParticle.h"
#include "THnSparse.h"
#include <iostream>
using namespace std;
#endif

// histogram classes
class TH3F;
class TH1F;

//fastjet interface

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/Selector.hh"

// set namespace as Pythia8
using namespace Pythia8;

//function for annuli bin

Int_t GetAnnuliBin(Float_t deltaR)
{
  //initialize annuli bin
  Int_t annuliBin = -99;

  // get annuli bin number
  if(deltaR >= 0.00 && deltaR <= 0.05)     { annuliBin = 0; }
  else if(deltaR > 0.05 && deltaR <= 0.10) { annuliBin = 1; }
  else if(deltaR > 0.10 && deltaR <= 0.15) { annuliBin = 2; }
  else if(deltaR > 0.15 && deltaR <= 0.20) { annuliBin = 3; }
  else if(deltaR > 0.20 && deltaR <= 0.25) { annuliBin = 4; }
  else if(deltaR > 0.25 && deltaR <= 0.30) { annuliBin = 5; }
  else if(deltaR > 0.30 && deltaR <= 0.35) { annuliBin = 6; }
  else if(deltaR > 0.35 && deltaR <= 0.40) { annuliBin = 7; }

  return annuliBin;

}

//function to get jet pt bin

Int_t GetJetPtBin(Double_t jetpt)
{
  // initialize jet pt bin
  Int_t jetPtBin = -99;

  // get jet pt bin number
  if(jetpt >= 15.0 /* && jetpt < 15.0 */)      { jetPtBin = 0; }
  //else if(jetpt >= 15.0 && jetpt < 20.0) { jetPtBin = 1; }
  //else if(jetpt >= 20.0 && jetpt < 40.0) { jetPtBin = 2; }
  //else if(jetpt >= 40.0 && jetpt < 60.0) { jetPtBin = 3; }
  //else if(jetpt >= 60.0)                 { jetPtBin = 4; }

  return jetPtBin;
}

//function to get track pt bin

Int_t GetTrackPtBin(Double_t pt)
{
  // initialize track pt bin
  Int_t trackPtBin = -99;

  // additional pt selection when doing pt associated bin method
  if((pt >= 0.20) && (pt < 0.5)) trackPtBin = 0;  // 0.20 - 0.5 GeV assoc bin used for jet shapes
  if((pt >= 0.50) && (pt < 1.0)) trackPtBin = 1;  // 0.50 - 1.0 GeV assoc bin used for jet shapes
  if((pt >= 1.00) && (pt < 1.5)) trackPtBin = 2;  // 1.00 - 1.5 GeV assoc bin used for jet shapes
  if((pt >= 1.50) && (pt < 2.0)) trackPtBin = 3;  // 1.50 - 2.0 GeV assoc bin used for jet shapes
  if((pt >= 2.00) && (pt < 3.0)) trackPtBin = 4;  // 2.00 - 3.0 GeV assoc bin used for jet shapes
  if((pt >= 3.00) && (pt < 4.0)) trackPtBin = 5;  // 3.00 - 4.0 GeV assoc bin used for jet shapes
  if((pt >= 4.00) && (pt < 6.0)) trackPtBin = 6;  // 4.00 - 6.0 GeV assoc bin used for jet shapes
  if(pt >=  6.0)                 trackPtBin = 7;  //       6.0+ GeV assoc bin used for jet shapes

  return trackPtBin;
}

//function to get jet mass bin
/*Int_t GetJetMassBin(Double_t jetmass)
{
  //initialize jet mass bin
  Int_t jetMassBin = -99;

  //get jet mass bin number
  if (jetmass >= 0. && jetmass < 1.0)       { jetMassBin = 0; }
  else if (jetmass >= 1.0 && jetmass < 2.0) { jetMassBin = 1; }
  else if (jetmass >= 2.0 && jetmass < 3.0) { jetMassBin = 2; }
  else if (jetmass >= 3.0 && jetmass < 4.0) { jetMassBin = 3; }
  else if (jetmass >= 4.0 && jetmass < 5.0) { jetMassBin = 4; }
  else if (jetmass >= 5.0 && jetmass < 6.0) { jetMassBin = 5; }
  else if (jetmass >= 6.0)                  { jetMassBin = 6; }

  return jetMassBin;

}*/

//function for deltaR

Double_t GetDeltaR(fastjet::PseudoJet jet, fastjet::PseudoJet trk) {

  // constants
  Float_t deltaR = -99.;
  Float_t pi = 1.0*TMath::Pi();

  // track variables
  Float_t tphi = trk.phi();
  if(tphi > 2.0*pi) tphi -= 2.0*pi;
  if(tphi < 0.0   ) tphi += 2.0*pi;
  Float_t teta = trk.eta();
  //Float_t teta = trk.rap(); // this is what FastJet uses

  // jet variables
  Float_t jphi = jet.phi();
  Float_t jeta = jet.eta();
  //Float_t jeta = jet.rap(); // this is what FastJet uses

  // calculate radial distance
  Float_t deltaEta = 1.0*TMath::Abs(jeta - teta);
  Float_t deltaPhi = 1.0*TMath::Abs(jphi - tphi);
  if(deltaPhi > 1.0*pi) deltaPhi = 2.0*pi - deltaPhi;
  deltaR = 1.0*TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

  return deltaR;
}

double calculate_energy(double px,double py,double pz,double m = 0.139570)
{

  double e = sqrt(pow(m,2)+pow(px,2)+pow(py,2)+pow(pz,2));

  return e;
}

void jet_shape_final()
{

  Int_t jobID = 0;
  Int_t tune = 5;
  Int_t power = -1;
  Float_t jetResolutionR = 0.40;


  //initialize pythia
  Pythia8::Pythia pythia;

  Settings& settings = pythia.settings;
  Pythia8::Info& info = pythia.info;


  //configure
  pythia.readString("HardQCD:all = on");
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 200.");
  pythia.readString("PhaseSpace:pTHatMin = 20.");

  pythia.init();

  //histogram time!

  //basic histos for essential distributions
  TH1F* pt_jet = new TH1F("", "PYTHIA pp 200 GeV", 84, 8.,50.);
  TH1F* pt_part = new TH1F("","PYTHIA pp 200 GeV", 100, 0., 50.);

  //jet shape histos
  TH1F *hFragFunc = new TH1F("hFragFunc", "Z = p_{T} / E^{jet}_{T}", 100, 0, 1);
  //arrays are used to store info according to pT bins and annuli bins
  TH1F* hJetShape[1][9] = {0x0};
  //TH1F* hJetMass[5][9] = {0x0};
  TH1F* hJetShapeALL = new TH1F("hJetShapeALL", "Jet shape #rho(r)", 8, 0.0, 0.40);

  for(Int_t i = 0; i < 1; i++)
  {
    for(Int_t j = 0; j < 9; j++)
    {
      if (j < 8)
      {
      const char *name = Form("JetShape_%i_%i", i, j);
      const char *title = Form("Jet Shape #rho(r) - Jet p_{T} bin %i, Associated p_{T} bin %i", i, j);
      hJetShape[i][j]  = new TH1F(name, title, 8, 0.0, 0.40);
    }
     if (j==8)
     {
       const char *name1 = Form("Jet Shape Stats");
       const char *title1 = Form("PYTHIA pp 200 GeV");
       hJetShape[i][j]  = new TH1F(name1, title1, 8, 0.0, 0.40);
     }
    }
  }



  fastjet::Strategy              strategy = fastjet::Best;
  fastjet::RecombinationScheme   recombScheme = fastjet::BIpt_scheme;
  fastjet::JetDefinition        *jetDef = NULL;
  fastjet::JetAlgorithm          algorithm;

  if (power == -1)      algorithm = fastjet::antikt_algorithm;
  if (power ==  0)      algorithm = fastjet::cambridge_algorithm;
  if (power ==  1)      algorithm = fastjet::kt_algorithm; // was set as the default - use for background

  //set jet definition
  jetDef = new fastjet::JetDefinition(algorithm, jetResolutionR, recombScheme, strategy);


  //make the pseudojet vector
  std::vector <fastjet::PseudoJet> fjInputs;

  //event time! generate pythia events; increase statistics for final run
  Int_t nev = 10000;

  int nJets = 0;
  for (Int_t event = 0; event < nev; event++)
  {
    if (!pythia.next()) continue;

    // Reset Fastjet input
    fjInputs.resize(0);

    // loop over particles
    for (Int_t part=0; part<pythia.event.size(); ++part)
    {
      Double_t pt= pythia.event[part].pT();
      Double_t phi = pythia.event[part].phi();
      Double_t eta = pythia.event[part].eta();

      Int_t charge = pythia.event[part].charge();

      //apply constituent cuts
      if (pythia.event[part].isFinal() && abs(pythia.event[part].eta())<1 && pythia.event[part].pT()>1.0 && pythia.event[part].isCharged())
      {
        double energy = calculate_energy(pythia.event[part].px(), pythia.event[part].py(), pythia.event[part].pz());
      // Store as input to Fastjet
      fastjet::PseudoJet particle_Track(pythia.event[part].px(), pythia.event[part].py(), pythia.event[part].pz(), energy);
      particle_Track.set_user_index(charge);
      fjInputs.push_back(particle_Track);

      pt_part->Fill(pythia.event[part].pT());

      }

    }//end particle loop

    if (fjInputs.size() == 0) continue;

    // Run Fastjet algorithm
    vector <fastjet::PseudoJet> inclusiveJets, sortedJets, selJets;
    fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);

    Float_t minJetPt = 10.0;
    inclusiveJets = clustSeq.inclusive_jets(minJetPt);

    //apply jet cuts
    fastjet::Selector select_eta = fastjet::SelectorEtaRange(-0.6, 0.6);
    fastjet::Selector select_phi = fastjet::SelectorPhiRange(0.0, 2.0*TMath::Pi());
    fastjet::Selector select_pt = fastjet::SelectorPtRange(10, 1000); // selecta 10 < pt < 1000
    fastjet::Selector select_Nhard = fastjet::SelectorNHardest(2);
    fastjet::Selector select_eta_phi_2hard = select_eta && select_phi;

    // subselection of jets
    selJets = sorted_by_pt(inclusiveJets);
    sortedJets = select_eta_phi_2hard(selJets);

    // constants
    Float_t pi = 1.0*TMath::Pi();
    Float_t rbinSize = 0.05;

    //annuli jet mass distribution
  /*  Float_t rJetMassSum0GeV[10] = {0.0};
    Float_t rJetMassSum[7][10] = {0.0};*/

    //loop over reconstructed jets
    for(unsigned i = 0; i < sortedJets.size(); i++)
    {
      //jet properties
      Float_t jetpt = sortedJets[i].pt();
      Float_t jetet = sortedJets[i].Et();
      Float_t jeteta = sortedJets[i].eta();
      Float_t jetphi = sortedJets[i].phi();

      Double_t jetmass = sortedJets[i].m();

      pt_jet->Fill(jetpt);

      //get jet constituent vectors
      vector <fastjet::PseudoJet> constituents = sortedJets[i].constituents();
      vector <fastjet::PseudoJet> sortedconstituents = sorted_by_pt(constituents);

      //check for jets with only 1 constituent
      int nJetConstituents = sortedconstituents.size();
      if(nJetConstituents == 1) continue;

      //get jet pt bin
      Int_t jetPtBin = GetJetPtBin(jetpt);
      if(jetPtBin < 0) {continue;}

      //get jet mass bin
    /*  Int_t jetMassBin = GetJetMassBin(jetmass);
      if(jetMassBin < 0) {continue;} */


      //annuli pt sums - initialize
      Float_t rsumALL[10] = {0.0};
      Float_t rsum1GeV[10] = {0.0};
      Float_t rsum[8][10] = {0.0};

      nJets += 1;

      //loop over jet constituents
      for (unsigned j = 0; j < sortedconstituents.size(); j++)
      {
        // jet constituents info
        Float_t pt = sortedconstituents[j].pt();
        Float_t eta = sortedconstituents[j].eta();
        Float_t phi = sortedconstituents[j].phi();

        // fill fragmentation function plot
        Double_t z = pt / jetet;
        hFragFunc->Fill(z);

        // get track pt bin
        Int_t trackPtBin = GetTrackPtBin(pt);
        if(trackPtBin < 0) {continue;}

        //get radial distance between track and jet axis
        Float_t deltaR = GetDeltaR(sortedJets[i], sortedconstituents[j]);

        // get annuli bin
        Int_t annuliBin = GetAnnuliBin(deltaR);
        if(annuliBin < 0) {continue;}

      /*  rJetMassSum[trackPtBin][annuliBin] = jetmass;
        if (jetmass > 0) rJetMassSum0GeV[annuliBin] = jetmass; */

        // calculate radial pt sum
        rsum[trackPtBin][annuliBin] += pt;
        rsumALL[annuliBin] += pt;
        if(pt >=  1.0) rsum1GeV[annuliBin] += pt;

      } //end constituent loop


      //fill histos for jet shapes - loop over annuli bins while filling
      for(Int_t k = 0; k < 9; k++)
      { // associated pt loop
         for(Int_t i = 0; i < 10; i++)
         { // annuli bin loop
           if(k <9) hJetShape[jetPtBin][k]->Fill(i*rbinSize + 1e-3, 1.0*rsum[k][i]/(jetpt));
           if(k==9) hJetShape[jetPtBin][9]->Fill(i*rbinSize + 1e-3, 1.0*rsum1GeV[i]/jetpt); // 1GeV+ bin
           //Int_t ent = hJetShape[jetPtBin][9]->GetEntries();
           //cout<<"entries: "<<ent<<endl;}
         }
       }//end jet shape histo loop

       //fill histos for jet mass - loop over annuli bins while filling
    /*   for(Int_t y = 0; y < 9; y++)
       { //associated pt loop
         for(Int_t h = 0; h < 10; h++)
         {
           if(y <9) hJetMass[jetMassBin][y]->Fill(1.0*rJetMassSum[y][h]);
           if(y==9) hJetMass[jetMassBin][9]->Fill(rJetMassSum0GeV[h]);
         }
       } */

       for(Int_t i = 0; i < 10; i++)
       { // annuli bin loop
         hJetShapeALL->Fill(i*rbinSize + 1e-3, 1.0*rsumALL[i]/(jetpt));
       }//end jet shape - ALL histo

    }//end sorted jets loop

  }//end event loop


    //write objects to root file

    //create root file
    TFile *MyFile = new TFile("Jet_Charge_Shape.root","RECREATE");
    MyFile->cd();

    pt_part->Scale(1.0/(pt_part->GetEntries()));
    pt_part->GetYaxis()->SetTitle("Normalized counts");
    pt_part->GetXaxis()->SetTitle("p_{T} (GeV/c) - Constituents");
    pt_part->SetStats(0);
    pt_part->Write();
    pt_jet->Scale(1.0/((pt_jet->GetEntries())*(pt_jet->GetBinWidth(84))));
    pt_jet->GetYaxis()->SetTitle("1/N_{jets} dN/dp_{T,jets}");
    pt_jet->GetXaxis()->SetTitle("p_{T} (GeV/c) - Jets");
    pt_jet->SetStats(0);
    pt_jet->Write();

    for(Int_t j = 0; j < 1; j++)
    {
      for(Int_t i = 0; i < 9; i++)
      {
        // hJetShape[j][i]->Write();
        if (i==8)
        {
          Int_t entries = hJetShape[j][8]->GetEntries();
          hJetShape[j][8]->Scale(1.0/(0.05*entries));
          hJetShape[j][8]->GetYaxis()->SetTitle("#rho(r)");
          hJetShape[j][8]->GetXaxis()->SetTitle("r");
          hJetShape[j][8]->SetStats(0);
          hJetShape[j][8]->Write();
          //hJetShape[j][8]->SetLogy();
        }
      }
    }//end write file loop


  /*  for(Int_t j = 0; j < 7; j++)
    {
      for(Int_t i = 0; i < 9; i++)
      {
        // hJetShape[j][i]->Write();
        if (i==8)
        {
          Int_t entries = hJetMass[j][8]->GetEntries();
          hJetMass[j][8]->Scale(1.0/(0.05*entries));
          hJetMass[j][8]->GetYaxis()->SetTitle("#rho(r)");
          hJetMass[j][8]->GetXaxis()->SetTitle("r");
          hJetMass[j][8]->Write();
          //hJetShape[j][8]->SetLogy();
        }
      }
    }//end write file loop */

}//THE END
