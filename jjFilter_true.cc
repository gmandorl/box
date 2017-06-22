// -*- C++ -*-
//
// Package:    filter/jjFilter
// Class:      jjFilter
// 
/**\class jjFilter jjFilter.cc filter/jjFilter/plugins/jjFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giulio Mandorli
//         Created:  Wed, 21 Jun 2017 13:13:12 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TLorentzVector.h>

#include <DataFormats/JetReco/interface/GenJet.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>

//
// class declaration
//

class jjFilter : public edm::EDFilter {
   public:
      explicit jjFilter(const edm::ParameterSet&);
      ~jjFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
jjFilter::jjFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


jjFilter::~jjFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
jjFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

        return true;

        using namespace edm;
        using namespace std;

        int numberOfMuAndTauHard = 0;
        int numberOfJet = 0;



        Handle<std::vector<reco::GenParticle>> particelesCollection;
        iEvent.getByLabel("genParticles", particelesCollection);
        std::vector<reco::GenParticle>  particles = *particelesCollection.product();
//         std::cout <<  "Partile size:  " << particles.size() << std::cout;

        Handle<std::vector<reco::GenJet>> jetsCollection;
        iEvent.getByLabel("ak4GenJets", jetsCollection);
        std::vector<reco::GenJet>  jets = *jetsCollection.product();




        std::vector<TLorentzVector> mu_and_tau_vector;
        for(std::vector<reco::GenParticle>::const_iterator pit = particles.begin() ; pit != particles.end() ; ++pit) {
            if (((abs(pit->pdgId()) == 11) || (abs(pit->pdgId()) == 13) || (abs(pit->pdgId()) == 15)) && pit->isHardProcess()) {
                 TLorentzVector tmpVector;
                 tmpVector.SetPtEtaPhiM(pit->p4().pt(), pit->p4().Eta(), pit->p4().Phi(), pit->p4().M());
                 mu_and_tau_vector.push_back(tmpVector);
            }
        }
        numberOfMuAndTauHard =  mu_and_tau_vector.size();


        if (numberOfMuAndTauHard < 2) return false; // == 2 ?


        std::vector<TLorentzVector> jets_without_lep;
        for(std::vector<reco::GenJet>::const_iterator jit = jets.begin() ; jit != jets.end() ; ++jit) {
             TLorentzVector tmpVector;
             tmpVector.SetPtEtaPhiM(jit->pt(),jit->eta(),jit->phi(),jit->mass());
             if (mu_and_tau_vector[0].DeltaR(tmpVector) > 0.3  && mu_and_tau_vector[1].DeltaR(tmpVector) > 0.3) jets_without_lep.push_back(tmpVector);
        }
        numberOfJet = jets_without_lep.size();


        if (numberOfJet < 2) return false;


        float Mll = (mu_and_tau_vector[0] + mu_and_tau_vector[1]).M();
        float Mjj = (jets_without_lep[0] + jets_without_lep[1]).M();

//        std::cout << "jets size: "  << jets.size() << " \t jets_without_lep size: "  << jets_without_lep.size() << "difference: "  << jets.size() - jets_without_lep.size() << std::endl;

        if (Mll > 100 && Mjj > 200) {
            std::cout << "jets size: "  << jets.size() << " \t jets_without_lep size: "  << jets_without_lep.size() << " \t difference: "  << jets.size() - jets_without_lep.size() << std::endl;
            std::cout << "Mjj: "  << Mjj << " \t Mll: "  << Mll << std::endl;
            
            std::cout << "mu_and_tau_vector[0].pt: "  << mu_and_tau_vector[0].Pt() << " \t mu_and_tau_vector[0].eta: "  << mu_and_tau_vector[0].Eta ()  << std::endl;
            return true;
        }
        else return false;

}



// ------------ method called once each job just before starting event loop  ------------
void 
jjFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
jjFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
jjFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
jjFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
jjFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
jjFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
jjFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(jjFilter);
