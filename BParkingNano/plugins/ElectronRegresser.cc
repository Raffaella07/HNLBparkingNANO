// Code to apply energy regression

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "CommonTools/CandAlgos/interface/ModifyObjectValueBase.h"
#include "helper.h"
#include "TVector3.h"

class ElectronRegresser : public edm::global::EDProducer<> {

public:
  bool debug=true; 

  explicit ElectronRegresser(const edm::ParameterSet &cfg):
    lowpt_src_{ consumes<pat::ElectronCollection>( cfg.getParameter<edm::InputTag>("lowptSrc") )},
    pf_src_{ consumes<pat::ElectronCollection>( cfg.getParameter<edm::InputTag>("pfSrc") )}
    {

      // LPT regression stuff                                                                                                         
      if( cfg.existsAs<edm::ParameterSet>("lowPtRegressionConfig") ) {
	const edm::ParameterSet& iconf = cfg.getParameterSet("lowPtRegressionConfig");
	const std::string& mname = iconf.getParameter<std::string>("modifierName");
	ModifyObjectValueBase* plugin =
	  ModifyObjectValueFactory::get()->create(mname,iconf);
	regression_.reset(plugin);
	edm::ConsumesCollector sumes = consumesCollector();
	regression_->setConsumes(sumes);
      } else {
	regression_.reset(nullptr);
      }

      // PF regression                                                                                                           
      if( cfg.existsAs<edm::ParameterSet>("gsfRegressionConfig") ) {
	const edm::ParameterSet& iconf = cfg.getParameterSet("gsfRegressionConfig");
	const std::string& mname = iconf.getParameter<std::string>("modifierName");
      ModifyObjectValueBase* plugin =
        ModifyObjectValueFactory::get()->create(mname,iconf);
      regressionGsf_.reset(plugin);
      edm::ConsumesCollector sumes = consumesCollector();
      regressionGsf_->setConsumes(sumes);
      } else {
	regressionGsf_.reset(nullptr);
      }

      produces<pat::ElectronCollection>("regressedElectrons");
      produces<pat::ElectronCollection>("regressedLowPtElectrons");
    }

  ~ElectronRegresser() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const edm::EDGetTokenT<pat::ElectronCollection> lowpt_src_;
  const edm::EDGetTokenT<pat::ElectronCollection> pf_src_;

  // regression stuff                                                                                                                                         
  std::unique_ptr<ModifyObjectValueBase> regression_; // Low pt                                                                                               
  std::unique_ptr<ModifyObjectValueBase> regressionGsf_; // Gsf                                                                                               
};

void ElectronRegresser::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const & iSetup) const {

  //input
  edm::Handle<pat::ElectronCollection> lowpt;
  evt.getByToken(lowpt_src_, lowpt);
  edm::Handle<pat::ElectronCollection> pf;
  evt.getByToken(pf_src_, pf);

  // regression stuff
  regression_->setEvent(evt);
  regression_->setEventContent(iSetup);
  regressionGsf_->setEvent(evt);
  regressionGsf_->setEventContent(iSetup);

  // output
  std::unique_ptr<pat::ElectronCollection>  ele_out_pf      (new pat::ElectronCollection );
  std::unique_ptr<pat::ElectronCollection>  ele_out_lpt      (new pat::ElectronCollection );

  // PF regression
  size_t ipfele=-1;
  for(auto ele : *pf) {
   ipfele++;
   if(debug) {
     std::cout << "ElectronRegresser, Event " << (evt.id()).event() 
	       << " => Pre regression, PF: ele.superCluster()->rawEnergy() = " << ele.superCluster()->rawEnergy()
	       << ", ele.correctedEcalEnergy() = " << ele.correctedEcalEnergy()
	       << ", ele.trackMomentumAtVxt() = " << ele.eSuperClusterOverP()
	       << ", ele.eEleClusterOverP() = " << ele.eEleClusterOverPout()
	       << ", ele gsf track chi2 = " << ele.core()->gsfTrack()->normalizedChi2()
	       << ", ele.p = " << ele.p() << std::endl;
   }
  /* float temp = ele.eSuperClusterOverP(); 
   struct reco::GsfElectron::TrackExtrapolations TrackExtrCopy;
   struct reco::GsfElectron::TrackClusterMatching TrackClusterMatchingCopy;
   struct reco::GsfElectron::FiducialFlags FiducialFlagsCopy;
   struct reco::GsfElectron::ShowerShape ShowerShapeCopy;
   struct reco::GsfElectron::ConversionRejection ConversionRejectionCopy;
   struct reco::GsfElectron::SaturationInfo SaturationInfoCopy;
   reco::GsfElectronCoreRef coreCopy;
   TVector3 computedPAtVtx(0,0,ele.superCluster()->energy()/ele.eSuperClusterOverP());
    
   TrackClusterMatchingCopy.electronCluster = ele.electronCluster() ;  // basic cluster best matching gsf track
   TrackClusterMatchingCopy.eSeedClusterOverP = ele.eSeedClusterOverP() ;         // the seed cluster energy / track momentum at the PCA to the beam spot
   TrackClusterMatchingCopy.eSeedClusterOverPout = ele.eSeedClusterOverPout() ;      // the seed cluster energy / track momentum at calo extrapolated from the outermost track state
   TrackClusterMatchingCopy.eEleClusterOverPout = ele.eEleClusterOverPout();       // the electron cluster energy / track momentum at calo extrapolated from the outermost track state
   TrackClusterMatchingCopy.deltaEtaSuperClusterAtVtx = ele.deltaEtaSuperClusterTrackAtVtx() ; // the supercluster eta - track eta position at calo extrapolated from innermost track state
   TrackClusterMatchingCopy.deltaEtaSeedClusterAtCalo = ele.deltaEtaSeedClusterTrackAtCalo() ; // the seed cluster eta - track eta position at calo extrapolated from the outermost track state
   TrackClusterMatchingCopy.deltaEtaEleClusterAtCalo = ele.deltaEtaEleClusterTrackAtCalo() ;  // the electron cluster eta - track eta position at calo extrapolated from the outermost state
   TrackClusterMatchingCopy.deltaPhiEleClusterAtCalo = ele.deltaPhiEleClusterTrackAtCalo();  // the electron cluster phi - track phi position at calo extrapolated from the outermost track state
   TrackClusterMatchingCopy.deltaPhiSuperClusterAtVtx = ele.deltaPhiSuperClusterTrackAtVtx() ; // the supercluster phi - track phi position at calo extrapolated from the innermost track state
   TrackClusterMatchingCopy.deltaPhiSeedClusterAtCalo = ele.deltaPhiSeedClusterTrackAtCalo(); // the seed cluster phi - track phi position at calo extrapolated from the outermost track state

   TrackExtrCopy.positionAtVtx = ele.trackPositionAtVtx() ;     // the track PCA to the beam spot
   TrackExtrCopy.positionAtCalo = ele.trackPositionAtCalo() ;    // the track PCA to the supercluster position
   TrackExtrCopy.momentumAtVtx = computedPAtVtx;     // the track momentum at the PCA to the beam spot
   TrackExtrCopy.momentumAtCalo =ele.trackMomentumAtCalo();    // the track momentum extrapolated at the supercluster position from the innermost track state
   TrackExtrCopy.momentumOut = ele.trackMomentumOut() ;       // the track momentum extrapolated at the seed cluster position from the outermost track state
   TrackExtrCopy.momentumAtEleClus = ele.trackMomentumAtEleClus() ; // the track momentum extrapolated at the ele cluster position from the outermost track state
   TrackExtrCopy.momentumAtVtxWithConstraint =  ele.trackMomentumAtVtxWithConstraint();     // the track momentum at the PCA to the beam spot using bs constraint
  
   coreCopy= ele.core(); 

    SaturationInfoCopy.nSaturatedXtals = ele.nSaturatedXtals() ;
    SaturationInfoCopy.isSeedSaturated = ele.isSeedSaturated();

    ConversionRejectionCopy.flags = ele.convFlags() ;  // -max:not-computed, other: as computed by Puneeth conversion code
    ConversionRejectionCopy.dist = ele.convDist() ; // distance to the conversion partner
    ConversionRejectionCopy.dcot = ele.convDcot() ; // difference of cot(angle) with the conversion partner track
    ConversionRejectionCopy.radius = ele.convRadius() ; // signed conversion radius
    //ConversionRejectionCopy.vtxFitProb  = ele.convVtxFitProb(); //fit probablity (chi2/ndof) of the matched conversion vtx
   
   
    ShowerShapeCopy.sigmaEtaEta = ele.sigmaEtaEta()  ;        // weighted cluster rms along eta and inside 5x5 (absolute eta)
    ShowerShapeCopy.sigmaIetaIeta = ele.sigmaIetaIeta() ;      // weighted cluster rms along eta and inside 5x5 (Xtal eta)
    ShowerShapeCopy.sigmaIphiIphi = ele.sigmaIphiIphi() ;      // weighted cluster rms along phi and inside 5x5 (Xtal phi)
    ShowerShapeCopy.e1x5 = ele.e1x5() ;               // energy inside 1x5 in etaxphi around the seed Xtal
    ShowerShapeCopy.e2x5Max = ele.e2x5Max() ;            // energy inside 2x5 in etaxphi around the seed Xtal (max bwt the 2 possible sums)
    ShowerShapeCopy.e5x5 = ele.e5x5() ;               // energy inside 5x5 in etaxphi around the seed Xtal
    ShowerShapeCopy.r9 = ele.r9();                 // ratio of the 3x3 energy and supercluster energy
    ShowerShapeCopy.hcalDepth1OverEcal = ele.hcalDepth1OverEcal() ; // hcal over ecal seed cluster energy using 1st hcal depth (using hcal towers within a cone)
    ShowerShapeCopy.hcalDepth2OverEcal = ele.hcalDepth2OverEcal() ; // hcal over ecal seed cluster energy using 2nd hcal depth (using hcal towers within a cone)
    ShowerShapeCopy.hcalDepth1OverEcalBc =  ele.hcalDepth1OverEcalBc(); // hcal over ecal seed cluster energy using 1st hcal depth (using hcal towers behind clusters)
    ShowerShapeCopy.hcalDepth2OverEcalBc =  ele.hcalDepth2OverEcalBc(); // hcal over ecal seed cluster energy using 2nd hcal depth (using hcal towers behind clusters)
    //ShowerShapeCopy.invalidHcal = ele.invalidHcal();          // set to true if the hcal energy estimate is not valid (e.g. the corresponding tower was off or masked)
    ShowerShapeCopy.sigmaIetaIphi = ele.full5x5_sigmaIetaIphi();
    ShowerShapeCopy.e2x5Max = ele.full5x5_e2x5Max();
    ShowerShapeCopy.eTop = ele.full5x5_eTop();
    ShowerShapeCopy.eLeft = ele.full5x5_eLeft();
    ShowerShapeCopy.eRight = ele.full5x5_eRight();
    ShowerShapeCopy.eBottom = ele.full5x5_eBottom(); 
    ShowerShapeCopy.e2x5Top = ele.full5x5_e2x5Top();
    ShowerShapeCopy.e2x5Left = ele.full5x5_e2x5Left();
    ShowerShapeCopy.e2x5Right = ele.full5x5_e2x5Right();
    ShowerShapeCopy.e2x5Bottom = ele.full5x5_e2x5Bottom(); 
   





   FiducialFlagsCopy.isEB = ele.isEB() ;        // true if particle is in ECAL Barrel
   FiducialFlagsCopy.isEE = ele.isEE() ;        // true if particle is in ECAL Endcaps
   FiducialFlagsCopy.isEBEEGap = ele.isEBEEGap() ;   // true if particle is in the crack between EB and EE
   FiducialFlagsCopy.isEBEtaGap = ele.isEBEtaGap() ;  // true if particle is in EB, and inside the eta gaps between modules
   FiducialFlagsCopy.isEBPhiGap = ele.isEBPhiGap() ;  // true if particle is in EB, and inside the phi gaps between modules
   FiducialFlagsCopy.isEEDeeGap = ele.isEEDeeGap();  // true if particle is in EE, and inside the gaps between dees
   FiducialFlagsCopy.isEERingGap = ele.isEERingGap(); // true if particle is in EE, and inside the gaps between rings
*/
   regressionGsf_->modifyObject(ele);
/*
   TrackClusterMatchingCopy.eSuperClusterOverP = ele.eSuperClusterOverP()*ele.ecalEnergy()/ele.superCluster()->energy() ;
  // ele.trackClusterMatcing = 	
        // the supercluster energy / track momentum at the PCA to the beam spot
   ele.setTrackExtrapolations(TrackExtrCopy);


   reco::GsfElectron gsfNewEle(
      ele.charge(),
      ele.chargeInfo() ,
      coreCopy ,
      TrackClusterMatchingCopy ,
      TrackExtrCopy ,
      ele.closestCtfTrack() ,
      FiducialFlagsCopy ,
      ShowerShapeCopy ,
      ele.showerShape() ,
      ConversionRejectionCopy ,
      SaturationInfoCopy
     ) ;
 
   pat::Electron NewEle(gsfNewEle);*/
   if(debug) { 
     std::cout << "ElectronRegresser, Event " << (evt.id()).event() 
	       << " => Post regression, PF: ele.superCluster()->rawEnergy() = " << ele.superCluster()->rawEnergy()
	       << ", ele.correctedEcalEnergy() = " << ele.correctedEcalEnergy()
	       << ", ele.trackMomentumAtVtx() = " << ele.eSuperClusterOverP()
	       << ", ele.eEleClusterOverP() = " << ele.eEleClusterOverPout()
	       << ", ele gsf track chi2 = " << ele.core()->gsfTrack()->normalizedChi2()
	       << ", ele.p = " << ele.p() << std::endl;
   }

  // std::cout << "check on temp _________________" << temp << std::endl;
   ele_out_pf -> emplace_back(ele);
  }

  // LowPt regression
  size_t iele=-1;
  for(auto ele : *lowpt) {
    iele++;
  /* if(debug){ 
     std::cout << "ElectronRegresser, Event " << (evt.id()).event() 
	       << " => Pre regression, LPT: ele.superCluster()->rawEnergy() = " << ele.superCluster()->rawEnergy()
	       << ", ele.correctedEcalEnergy() = " << ele.correctedEcalEnergy()
	       << ", ele gsf track chi2 = " << ele.core()->gsfTrack()->normalizedChi2()
	       << ", ele.p = " << ele.p() << std::endl;
   }
*/
   regression_->modifyObject(ele);
/*
   if(debug) {
     std::cout << "ElectronRegresser, Event " << (evt.id()).event() 
	       << " => Post regression, LPT: ele.superCluster()->rawEnergy() = " << ele.superCluster()->rawEnergy()
	       << ", ele.correctedEcalEnergy() = " << ele.correctedEcalEnergy()
	       << ", ele gsf track chi2 = " << ele.core()->gsfTrack()->normalizedChi2()
	       << ", ele.p = " << ele.p() << std::endl;
   }
*/
   ele_out_lpt -> emplace_back(ele);
  }
   
  // put collections in the event
  evt.put(std::move(ele_out_lpt),  "regressedLowPtElectrons");
  evt.put(std::move(ele_out_pf),  "regressedElectrons");
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronRegresser);
