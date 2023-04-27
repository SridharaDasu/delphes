#include "modules/TauReconstructor.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "TLorentzVector.h"

TauReconstructor::TauReconstructor() : fMinTauSeedPT(5.), fMaxTauSeedEta(2.5), fMaxTauIsolDeltaR(0.5), fMaxTauCoreDeltaR(0.3) {
}

TauReconstructor::~TauReconstructor() {return;}

void TauReconstructor::Init() {
  fMinTauSeedPT = GetDouble("MinTauSeedPT", 5.0);
  fMaxTauSeedEta = GetDouble("MaxTauSeedEta", 2.5);
  fMaxTauIsolDeltaR = GetDouble("MaxTauIsolDeltaR", 0.5);
  fMaxTauCoreDeltaR = GetDouble("MaxTauCoreDeltaR", 0.3);
  fInputArray = ImportArray(GetString("InputArray", "TrackSmearing/tracks"));
  fItInputArray = fInputArray->MakeIterator();
  fOutputArray = ExportArray(GetString("OutputArray", "taus"));
}

void TauReconstructor::Process() {
  std::vector <Candidate *> trackSeeds; 
  fItInputArray->Reset();
  while(Candidate *track = static_cast<Candidate *>(fItInputArray->Next())) {
    if(track->Momentum.Pt() < fMinTauSeedPT || fabs(track->Momentum.Eta()) > fMaxTauSeedEta)
      continue;
    else
      trackSeeds.push_back(track);
  }
  std::vector <Tau> recTaus;
  TLorentzVector recoTau;
  for (UInt_t i = 0; i < trackSeeds.size(); i++) {
    Candidate *track = trackSeeds[i];
    fItInputArray->Reset();
    bool selected = true;
    recoTau.SetPtEtaPhiM(track->Momentum.Pt(), track->Momentum.Eta(), track->Momentum.Phi(), track->Momentum.M());
    double isolation = 0;
    int charge = track->Charge;
    int nProngs = 1;
    while(Candidate* fragment = static_cast<Candidate *>(fItInputArray->Next())) {
      // if (entry < 10) cout << "Max: EFTrack[" << t << "] = (" << track->PT << ", " << track->Eta << ", " << track->Phi << ", " << track->Mass << ")";
      if (track != fragment) {
	double deltaR = recoTau.DeltaR(fragment->Momentum);
	if(deltaR < fMaxTauIsolDeltaR) {
	  if (fragment->Momentum.Pt() > track->Momentum.Pt()) {
	    // if (entry < 10) cout << " - Deselected" << endl;
	    selected = false;
	    break;
	  }
	  else {
	    if (deltaR < fMaxTauCoreDeltaR) {
	      recoTau += fragment->Momentum;
	      // if (entry < 10) cout << endl << "Frg: EFTrack[" << a << "] = (" << fragment->PT << ", " << fragment->Eta << ", " << fragment->Phi << ", " << fragment->Mass << ")";
	      charge += fragment->Charge;
	      nProngs += 1;
	    }
	    else {
	      isolation += fragment->Momentum.Pt();
	    }
	  }
	}
      }
    }
    // Add in photons (mostly from pizero decays) if within 0.3, and for isolation if between 0.3-0.5
    int nPhotons = 0;
    /*
    for (Int_t a = t + 1; a < branchEFPhotons->GetEntries(); ++a) {
      Tower *fragment = (Tower*) branchEFPhotons->At(a);
      fragmentP4.SetPtEtaPhiM(fragment->ET, fragment->Eta, fragment->Phi, 0.);
      double deltaR = recoTau.DeltaR(fragmentP4);
      if(deltaR < fMaxTauIsolDeltaR) {
	if (deltaR < fMaxTauCoreDeltaR) {
	  recoTau += fragmentP4;
	  nPhotons++;
	}
	else {
	  isolation += fragment->ET;
	}
      }
    }
    */
    // Consider non-pizero remnant neutral hadrons only for isolation
    int nNHadrons = 0;
    /*
    for (Int_t a = t + 1; a < branchEFNHadrons->GetEntries(); ++a) {
      Tower *fragment = (Tower*) branchEFNHadrons->At(a);
      fragmentP4.SetPtEtaPhiM(fragment->ET, fragment->Eta, fragment->Phi, 0.);
      double deltaR = recoTau.DeltaR(fragmentP4);
      if(deltaR < fMaxTauIsolDeltaR) {
	if (deltaR < fMaxTauCoreDeltaR) {
	  recoTau += fragmentP4;
	  nNHadrons++;
	}
	else {
	  isolation += fragment->ET;
	}
      }
    }
    */
    if (selected) {
      // if (entry < 10) cout << endl << "SelRecoTau = (" << recoTau.Pt() << ", " << recoTau.Eta() << ", " << recoTau.Phi() << ", "<< recoTau.M() << ")" << endl;
      if (nProngs <= 5) {
	Tau tau;
	tau.PT = recoTau.Pt();
	tau.Eta = recoTau.Eta();
	tau.Phi = recoTau.Phi();
	tau.Mass = recoTau.M();
	// tau.T = track->T;
	tau.Charge = charge;
	tau.nProngs = nProngs;
	tau.nPhotons = nPhotons;
	tau.nNHadrons = nNHadrons;
	tau.isolation = isolation;
	recTaus.push_back(tau);
      }
    } // Selected objects
  } // All seeds above 5 GeV
  
}

void TauReconstructor::Finish() {
}
