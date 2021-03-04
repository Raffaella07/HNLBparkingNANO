///////////////////////// Code to produce K* candidates ////////////////////////


#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"




template<typename Lepton>
class HNLLeptonBuilder : public edm::global::EDProducer<> {

  
public:

  typedef std::vector<Lepton> LeptonCollection;
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  
  explicit HNLLeptonBuilder(const edm::ParameterSet &cfg):
    lepton_selection_{cfg.getParameter<std::string>("leptonSelection")},
    trk_selection_{cfg.getParameter<std::string>("trkSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    pions_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("pions") )},
    ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("pionsTransientTracks") )}, 
    src_{consumes<LeptonCollection>( cfg.getParameter<edm::InputTag>("src") )},
    leptons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("leptonsTransientTracks") )},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} { 
      //output
       produces<pat::CompositeCandidateCollection>();

    }

  ~HNLLeptonBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<Lepton> lepton_selection_; // cuts on leading cand
  const StringCutObjectSelector<pat::CompositeCandidate> trk_selection_; // sub-leading cand
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> pions_; //input PF cands this is sorted in pT in previous step
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_; //input TTracks of PF cands
  const edm::EDGetTokenT<LeptonCollection> src_;
  const edm::EDGetTokenT<TransientTrackCollection> leptons_ttracks_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};


template<typename Lepton>
void HNLLeptonBuilder<Lepton>::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //inputs  
  edm::Handle<LeptonCollection> leptons;
  evt.getByToken(src_, leptons);
 
  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(leptons_ttracks_, leptons_ttracks);
 
  edm::Handle<pat::CompositeCandidateCollection> pions;
  evt.getByToken(pions_, pions);
  
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_, ttracks);
 
  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  
  
  bool debug = false; 

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> hnl_out(new pat::CompositeCandidateCollection());

 if (debug)std::cout << "___________________________HNL MAIN LOOP"<< std::endl;

  // main loop
  for(size_t lep_idx = 0; lep_idx < leptons->size(); ++lep_idx) {


   edm::Ptr<Lepton> l_ptr(leptons, lep_idx);
   if(!lepton_selection_(*l_ptr)) continue;
 if (debug)std::cout << "___________________________DEBUG LEPTONCYCLE " << l_ptr->pt() <<" "<< l_ptr->isPF()<< std::endl;
   for(size_t trk_idx = 0; trk_idx < pions->size(); ++trk_idx ){

     edm::Ptr<pat::CompositeCandidate> trk_ptr( pions, trk_idx );
     if(!trk_selection_(*trk_ptr)) continue; 
     if (trk_ptr->charge()+l_ptr->charge()!=0 ) continue; 
          
// std::cout << "___________________________DEBUG trackCYCLE"<< std::endl;
     // create a HNL* candidate; add first quantities that can be used for pre fit selection
     pat::CompositeCandidate hnl_cand;
     auto trk_p4=trk_ptr->polarP4();
     auto lep_p4=l_ptr->polarP4();
     trk_p4.SetM(PI_MASS);
     lep_p4.SetM(l_ptr->mass());

// std::cout << "___________________________DEBUG INDEXES"<< std::endl;
     //adding stuff for pre fit selection
     hnl_cand.setP4(trk_p4 + lep_p4);
     hnl_cand.setCharge(trk_ptr->charge() + l_ptr->charge());
     hnl_cand.addUserFloat("trk_deltaR", reco::deltaR(*trk_ptr, *l_ptr));

     // save indices
     hnl_cand.addUserInt("trk_idx", trk_idx );
     hnl_cand.addUserInt("lep_idx", lep_idx );

     // save cands      
     hnl_cand.addUserCand("trk", trk_ptr );
     hnl_cand.addUserCand("lep", l_ptr );

     //second mass hypothesis
  //   trk1_p4.SetM(PI_MASS);
  //   trk2_p4.SetM(K_MASS);
   //  hnl_cand.addUserFloat("barMass", (trk1_p4 + trk2_p4).M() );
     
     // selection before fit
     if( !pre_vtx_selection_(hnl_cand) ) continue;
           
 if (debug)std::cout << "___________________________DEBUG prefit " << leptons_ttracks->at(lep_idx).isValid() << " " <<  ttracks->at(trk_idx).isValid() << std::endl;

     KinVtxFitter fitter(
       {leptons_ttracks->at(lep_idx), ttracks->at(trk_idx)},
       { l_ptr->mass(), PI_MASS },
       {LEP_SIGMA, PI_SIGMA} //K and PI sigma equal...
        );
      if ( !fitter.success() ) continue;           

      // save quantities after fit
      hnl_cand.setVertex( 
          reco::Candidate::Point( 
            fitter.fitted_vtx().x(),
            fitter.fitted_vtx().y(),
            fitter.fitted_vtx().z()
          )  
        );

 if (debug)std::cout << "___________________________DEBUG postfit"<< std::endl;
     auto fit_p4 = fitter.fitted_p4();
     auto lxy    = l_xy(fitter, *beamspot);

     hnl_cand.addUserInt  ("hnl_vtx_OK"             , fitter.success()                                                        );
     hnl_cand.addUserFloat("hnl_vtx_chi2"           , fitter.chi2()                                                           );
     hnl_cand.addUserFloat("hnl_vtx_ndof"           , fitter.dof()                                                            ); // float??
     hnl_cand.addUserFloat("hnl_vtx_prob"           , fitter.prob()                                                           );
     hnl_cand.addUserFloat("hnl_fitted_pt"          , fit_p4.pt()                                                             ); 
     hnl_cand.addUserFloat("hnl_fitted_eta"         , fit_p4.eta()                                                            );
     hnl_cand.addUserFloat("hnl_fitted_phi"         , fit_p4.phi()                                                            );
     hnl_cand.addUserFloat("hnl_fitted_mass"        , fitter.fitted_candidate().mass()                                        );      
     hnl_cand.addUserFloat("hnl_fitted_massErr"     , sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)));      
     hnl_cand.addUserFloat("hnl_cos_theta_2D"       , cos_theta_2D(fitter, *beamspot, hnl_cand.p4())                          );
     hnl_cand.addUserFloat("hnl_fitted_cos_theta_2D", cos_theta_2D(fitter, *beamspot, fit_p4)                                 );
     hnl_cand.addUserFloat("hnl_l_xy"               , lxy.value()                                                             );
     hnl_cand.addUserFloat("hnl_l_xy_unc"           , lxy.error()                                                             );
     hnl_cand.addUserFloat("hnl_ls_xy"              , lxy.value()/lxy.error()                                                 );
     hnl_cand.addUserFloat("hnl_vtx_x"              , hnl_cand.vx()                                                           );
     hnl_cand.addUserFloat("hnl_vtx_y"              , hnl_cand.vy()                                                           );
     hnl_cand.addUserFloat("hnl_vtx_z"              , hnl_cand.vz()                                                           );
     hnl_cand.addUserFloat("hnl_vtx_ex"             , sqrt(fitter.fitted_vtx_uncertainty().cxx())                             );
     hnl_cand.addUserFloat("hnl_vtx_ey"             , sqrt(fitter.fitted_vtx_uncertainty().cyy())                             );
     hnl_cand.addUserFloat("hnl_vtx_ez"             , sqrt(fitter.fitted_vtx_uncertainty().czz())                             );
     hnl_cand.addUserFloat("hnl_fitted_l_pt"       , fitter.daughter_p4(0).pt()                                              ); 
     hnl_cand.addUserFloat("hnl_fitted_l_eta"      , fitter.daughter_p4(0).eta()                                             );
     hnl_cand.addUserFloat("hnl_fitted_l_phi"      , fitter.daughter_p4(0).phi()                                             );
     hnl_cand.addUserFloat("hnl_fitted_pi_pt"       , fitter.daughter_p4(1).pt()                                              ); 
     hnl_cand.addUserFloat("hnl_fitted_pi_eta"      , fitter.daughter_p4(1).eta()                                             );
     hnl_cand.addUserFloat("hnl_fitted_pi_phi"      , fitter.daughter_p4(1).phi()                                             ); 

      // second mass hypothesis
   //   auto fitted_trk1= fitter.daughter_p4(0);
   //   auto fitted_trk2= fitter.daughter_p4(1);
     // fitted_trk1.SetM(PI_MASS);
    //  fitted_trk2.SetM(K_MASS);
   //   hnl_cand.addUserFloat("fitted_barMass", (fitted_trk1+fitted_trk2).M() );
                    
 if (debug)std::cout << "___________________________DEBUG variabless"<< std::endl;
      // after fit selection
      if( !post_vtx_selection_(hnl_cand) ) continue;
      hnl_out->emplace_back(hnl_cand);
      }
  }
  
  evt.put(std::move(hnl_out));
}

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
//typedef HNLLeptonBuilder<pat::Muon> HNLMuonBuilder;
typedef HNLLeptonBuilder<pat::Electron> HNLElectronBuilder;


#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(HNLMuonBuilder);
DEFINE_FWK_MODULE(HNLElectronBuilder);

