import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *


########## inputs preparation ################

HNLToPiMu = cms.EDProducer(
       'HNLMuonBuilder',
        src= cms.InputTag('muonsForAnalysis', 'SelectedMuons'),
        leptonsTransientTracks= cms.InputTag('muonsForAnalysis', 'SelectedTransientMuons'),
        pions = cms.InputTag('tracksBPark', 'SelectedTracks'),
        pionsTransientTracks= cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
        leptonSelection = cms.string('pt > 1.5 && abs(eta)<2.4'), #need optimization   
        trkSelection = cms.string('pt > 0.5 && abs(eta)<2.4'), #need optimization
    	beamSpot = cms.InputTag("offlineBeamSpot"),
        preVtxSelection = cms.string('abs(userCand("trk1").vz - userCand("trk2").vz)<1.0' 
        ' &&  pt()>2.0 && ( (mass() < 1.042 && mass() > 0.742)'
        ' || (userFloat("barMass") < 1.042 && userFloat("barMass") > 0.742) ) '
        ),
        postVtxSelection = cms.string('userFloat("sv_prob") > 1.e-5'
        ' && (  (userFloat("fitted_mass")<1.042 && userFloat("fitted_mass")>0.742)'
        ' || (userFloat("fitted_barMass")<1.042 && userFloat("fitted_barMass")>0.742)  )'
)
)

HNLToPiE = cms.EDProducer(
       'HNLElectronBuilder',
        src= cms.InputTag('electronsForAnalysis', 'SelectedElectrons'),
        leptonsTransientTracks= cms.InputTag('electronsForAnalysis', 'SelectedTransientElectrons'),
        pions = cms.InputTag('tracksBPark', 'SelectedTracks'),
        pionsTransientTracks= cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
        leptonSelection = cms.string('pt > 2 && abs(eta)<2.4'), #need optimization   
        trkSelection = cms.string('pt > 0.5 && abs(eta)<2.4'), #need optimization
    	beamSpot = cms.InputTag("offlineBeamSpot"),
        preVtxSelection = cms.string( ' mass > 0.2 && mass < 7.0'#'abs(userCand("trk").vz - userCand("lep").vz)<1.0' 
      #  ' &&  pt()>2.0 && ( (mass() < 1.042 && mass() > 0.742)'
    #    ' || (userFloat("barMass") < 1.042 && userFloat("barMass") > 0.742) ) '
        ),
        postVtxSelection = cms.string( 'userInt("hnl_vtx_OK") == 1 && userFloat("hnl_vtx_prob") > 0.0001 && userFloat("hnl_fitted_cos_theta_2D") >= 0.5 && userFloat("hnl_fitted_mass") > 0.5 && userFloat("hnl_fitted_mass") < 6.5'#'userFloat("sv_prob") > 1.e-5'
#        ' && (  (userFloat("fitted_mass")<1.042 && userFloat("fitted_mass")>0.742)'
 #       ' || (userFloat("fitted_barMass")<1.042 && userFloat("fitted_barMass")>0.742)  )'
)
)
DToKPi = cms.EDProducer(
       'DBuilder',
        pfcands= cms.InputTag('tracksBPark', 'SelectedTracks'),
        transientTracks= cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
        trk1Selection = cms.string('pt > 1.5 && abs(eta)<2.4'), #need optimization   
        trk2Selection = cms.string('pt > 1.0 && abs(eta)<2.4'), #need optimization
        preVtxSelection = cms.string('abs(userCand("trk1").vz - userCand("trk2").vz)<1.0' 
        ' &&  pt()>2.0 && ( (mass() < 1.042 && mass() > 0.742)'
        ' || (userFloat("barMass") < 1.042 && userFloat("barMass") > 0.742) ) '
        ),
        postVtxSelection = cms.string('userFloat("sv_prob") > 1.e-5'
        ' && (  (userFloat("fitted_mass")<1.042 && userFloat("fitted_mass")>0.742)'
        ' || (userFloat("fitted_barMass")<1.042 && userFloat("fitted_barMass")>0.742)  )'
)
)



########################### B-> Dl lPi ##########################
BToDPiMuMu = cms.EDProducer(
    'BToDPiLLBuilder',
    leptons = cms.InputTag('muonsForAnalysis', 'SelectedMuons'),
    leptonsTransientTracks = HNLToPiMu.leptonsTransientTracks,
    dmesons = cms.InputTag('DToKPi'),
    DTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    HNLs = cms.InputTag('HNLToPiMu'),
    HNLLepTransientTracks = HNLToPiMu.leptonsTransientTracks,
    HNLTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    lostTracks = cms.InputTag("lostTracks"),
    isoTracksSelection = cms.string('pt > 0.7 && abs(eta)<2.5'),
    
    beamSpot = cms.InputTag("offlineBeamSpot"),
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
        '&& ( (mass < 7. && mass > 4.) '
        '|| (userFloat("barMass")<7. && userFloat("barMass")>4.) )'
        ),
    postVtxSelection = cms.string(
        'userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& ( (userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.)'
        '|| (userFloat("fitted_barMass") > 4.5 && userFloat("fitted_barMass") < 6.)  )'
    )
)

BToDPiEE = cms.EDProducer(
    'BToDPiLLBuilder',
    leptons = cms.InputTag('electronsForAnalysis', 'SelectedElectrons'),
    leptonsTransientTracks = HNLToPiE.leptonsTransientTracks,
    dmesons = cms.InputTag('DToKPi'),
    DTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    HNLs = cms.InputTag('HNLToPiE'),
    HNLLepTransientTracks = HNLToPiE.leptonsTransientTracks,
    HNLTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    lostTracks = cms.InputTag("lostTracks"),
    isoTracksSelection = cms.string('pt > 0.7 && abs(eta)<2.5'),
    
    beamSpot = cms.InputTag("offlineBeamSpot"),
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
        '&& ( (mass < 7. && mass > 4.) '
        '|| (userFloat("barMass")<7. && userFloat("barMass")>4.) )'
        ),
    postVtxSelection = cms.string(
        'userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& ( (userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.)'
        '|| (userFloat("fitted_barMass") > 4.5 && userFloat("fitted_barMass") < 6.)  )'
    )
)
BToDPiMuE = cms.EDProducer(
    'BToDPiLLBuilder',
    leptons = cms.InputTag('muonsForAnalysis', 'SelectedMuons'),
    leptonsTransientTracks = HNLToPiMu.leptonsTransientTracks,
    dmesons = cms.InputTag('DToKPi'),
    dSelection = cms.string(''),
    piSelection = cms.string(''),
    DTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    HNLs = cms.InputTag('HNLToPiE'),
    HNLLepTransientTracks = HNLToPiE.leptonsTransientTracks,
    HNLTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    isoTracksSelection = cms.string('pt > 0.7 && abs(eta)<2.5'),
    
    beamSpot = cms.InputTag("offlineBeamSpot"),
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
        '&& ( (mass < 7. && mass > 4.) '
        '|| (userFloat("barMass")<7. && userFloat("barMass")>4.) )'
        ),
    postVtxSelection = cms.string(
        'userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& ( (userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.)'
        '|| (userFloat("fitted_barMass") > 4.5 && userFloat("fitted_barMass") < 6.)  )'
    )
)
BToDPiEMu = cms.EDProducer(
    'BToDPiLLBuilder',
    leptons = cms.InputTag('electronForAnalysis', 'SelectedElectrons'),
    leptonsTransientTracks = HNLToPiE.leptonsTransientTracks,
    dmesons = cms.InputTag('DToKPi'),
    DTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    HNLs = cms.InputTag('HNLToPiMu'),
    HNLLepTransientTracks = HNLToPiMu.leptonsTransientTracks,
    HNLTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    isoTracksSelection = cms.string('pt > 0.7 && abs(eta)<2.5'),
    
    beamSpot = cms.InputTag("offlineBeamSpot"),
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
        '&& ( (mass < 7. && mass > 4.) '
        '|| (userFloat("barMass")<7. && userFloat("barMass")>4.) )'
        ),
    postVtxSelection = cms.string(
        'userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& ( (userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.)'
        '|| (userFloat("fitted_barMass") > 4.5 && userFloat("fitted_barMass") < 6.)  )'
    )
)

################################### Tables #####################################

DToKPiTable = cms.EDProducer(
   'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("DToKPi"),
    cut = cms.string(""),
    name = cms.string("D"),
    doc = cms.string("D Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
      CandVars,
      barMass = ufloat('barMass'),
      fitted_mass = ufloat('fitted_mass'),
      fitted_barMass = ufloat('fitted_barMass'),
      fitted_pt = ufloat('fitted_pt'),
      fitted_eta = ufloat('fitted_eta'),
      fitted_phi = ufloat('fitted_phi'),
      svprob = ufloat('sv_prob'),         
      trk_deltaR = ufloat('trk_deltaR'),
      trk1_idx = uint('trk1_idx'),
      trk2_idx = uint('trk2_idx')
    )
)

HNLToPiMuTable = cms.EDProducer(
   'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("HNLToPiMu"),
    cut = cms.string(""),
    name = cms.string("HNLToPiMu"),
    doc = cms.string("HNLToPiMu Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
      CandVars,
      ltrk_deltaR          =  ufloat('trk_deltaR'),
      lep_idx               =  uint('lep_idx'),
      trk_idx             =  uint('trk_idx'),
      vtx_OK              =  uint('hnl_vtx_OK'),  
      vtx_chi2            =  ufloat('hnl_vtx_chi2'),
      vtx_ndof            =  ufloat('hnl_vtx_ndof'),
      vtx_prob            =  ufloat('hnl_vtx_prob'),
      fitted_pt           =  ufloat('hnl_fitted_pt'),
      fitted_eta          =  ufloat('hnl_fitted_eta'),
      fitted_phi          =  ufloat('hnl_fitted_phi'),
      fitted_mass         =  ufloat('hnl_fitted_mass'),
      fitted_massErr      =  ufloat('hnl_fitted_massErr'),
      cos_theta_2D        =  ufloat('hnl_cos_theta_2D'),
      fitted_cos_theta    =  ufloat('hnl_fitted_cos_theta_2D'),
      l_xy                =  ufloat('hnl_l_xy'),
      l_xy_unc            =  ufloat('hnl_l_xy_unc'),
      ls_xy               =  ufloat('hnl_ls_xy'),
      vtx_x               =  ufloat('hnl_vtx_x'),
      vtx_y               =  ufloat('hnl_vtx_y'),
      vtx_z               =  ufloat('hnl_vtx_z'),
      vtx_ex              =  ufloat('hnl_vtx_ex'),
      vtx_ey              =  ufloat('hnl_vtx_ey'),
      vtx_ez              =  ufloat('hnl_vtx_ez'),
      fitted_l_pt        =  ufloat('hnl_fitted_l_pt'),
      fitted_l_eta       =  ufloat('hnl_fitted_l_eta'),
      fitted_l_phi       =  ufloat('hnl_fitted_l_phi'),
      fitted_pi_pt        =  ufloat('hnl_fitted_pi_pt'),
      fitted_pi_eta       =  ufloat('hnl_fitted_pi_eta'),
      fitted_pi_phi       =  ufloat('hnl_fitted_pi_phi') 
    )
)

HNLToPiETable = HNLToPiMuTable.clone(
    src = cms.InputTag("HNLToPiE"),
    name = cms.string("HNLToPiE"),
    doc = cms.string("HNLToPiE Variables")
)





BToDPiMuMuTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToDPiMuMu"),
    cut = cms.string(""),
    name = cms.string("BToDPiMuMu"),
    doc = cms.string("BToDPiMuMu Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        l_idx = uint('l_idx'),
        d_trk1_idx = uint('trk1_idx'),
        d_trk2_idx = uint('trk2_idx'),
        d_idx = uint('kstar_idx'),
        hnl_idx = uint('hnl_idx'),
        hnl_l_idx = uint('l_idx'),
        hnl_trk_idx = uint('trk_idx'),
        min_dr = ufloat('min_dr'),
        max_dr = ufloat('max_dr'),
        # fit and vtx info
        # additional mass hypothesis
        Mass = ufloat ('Mass'),
        barMass = ufloat ('barMass'),
        fit_barMass = ufloat('fitted_barMass'),
        # post-fit tracks/leptons
        #l1
      # fit_l1_pt  = ufloat('fitted_l1_pt'),
      # fit_l1_eta = ufloat('fitted_l1_eta'),
      # fit_l1_phi = ufloat('fitted_l1_phi'),
      # #l2
      # fit_l2_pt  = ufloat('fitted_l2_pt'),
      # fit_l2_eta = ufloat('fitted_l2_eta'),
      # fit_l2_phi = ufloat('fitted_l2_phi'),
      # #trk1
      # fit_trk1_pt  = ufloat('fitted_trk1_pt'),
      # fit_trk1_eta = ufloat('fitted_trk1_eta'),
      # fit_trk1_phi = ufloat('fitted_trk1_phi'),
      # #trk2
      # fit_trk2_pt  = ufloat('fitted_trk2_pt'),
      # fit_trk2_eta = ufloat('fitted_trk2_eta'),
      # fit_trk2_phi = ufloat('fitted_trk2_phi'),
        # isolation 
        l1_iso03 = ufloat('l1_iso03'),
        l1_iso04 = ufloat('l1_iso04'),
        l2_iso03 = ufloat('l2_iso03'),
        l2_iso04 = ufloat('l2_iso04'),
        tk1_iso03 = ufloat('tk1_iso03'),
        tk1_iso04 = ufloat('tk1_iso04'),
        tk2_iso03 = ufloat('tk2_iso03'),
        tk2_iso04 = ufloat('tk2_iso04'),
        b_iso03  = ufloat('b_iso03'),
        b_iso04  = ufloat('b_iso04'),
    )
)

BToDPiEETable = BToDPiMuMuTable.clone(
    src = cms.InputTag("BToDPiEE"),
    name = cms.string("BToDPiEE"),
    doc = cms.string("BToDPiEE Variables")
)

BToDPiEMuTable = BToDPiMuMuTable.clone(
    src = cms.InputTag("BToDPiEMu"),
    name = cms.string("BToDPiEMu"),
    doc = cms.string("BToDPiEMu Variables")
)
BToDPiMuETable = BToDPiMuMuTable.clone(
    src = cms.InputTag("BToDPiMuE"),
    name = cms.string("BToDPiMuE"),
    doc = cms.string("BToDPiMuE Variables")
)



########################### Sequencies  ############################

HNLToPiMuSequence = cms.Sequence(  HNLToPiMu )
HNLToPiESequence = cms.Sequence(  HNLToPiE  )
DToKPiSequence = cms.Sequence(  DToKPi  )

BToDPiMuMuSequence = cms.Sequence( BToDPiMuMu  )
BToDPiEESequence = cms.Sequence( BToDPiEE  )
BToDPiEMuSequence = cms.Sequence( BToDPiEMu  )
BToDPiMuESequence = cms.Sequence( BToDPiMuE  )


#BToKstarEESequence = cms.Sequence(
 
 #  (electronPairsForKstarEE *BToKstarEE ))


#BToKstarLLSequence = cms.Sequence(
 #   ( (muonPairsForKstarMuMu *BToKstarMuMu)
  #   +(electronPairsForKstarEE *BToKstarEE) )   
#)


BToDPiLLTables = cms.Sequence( BToDPiMuMuTable + BToDPiEETable+ BToDPiMuETable + BToDPiEMuTable)

