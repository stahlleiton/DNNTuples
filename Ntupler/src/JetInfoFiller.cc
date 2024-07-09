/*
 * JetInfoAK8.cc
 *
 *  Created on: May 23, 2017
 *      Author: hqu
 */

#include "DeepNTuples/Ntupler/interface/JetInfoFiller.h"

#include <algorithm>


namespace deepntuples {

void JetInfoFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector && cc) {
  minPt_ = iConfig.getUntrackedParameter<double>("jetPtMin", 15);
  maxPt_ = iConfig.getUntrackedParameter<double>("jetPtMax", 1000);
  maxAbsEta_ = iConfig.getUntrackedParameter<double>("jetAbsEtaMax", 2.5);
  btag_discriminators_ = iConfig.getParameter<std::vector<std::string>>("bDiscriminators");

  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  puToken_ = cc.consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("puInfo"));
  rhoToken_ = cc.consumes<double>(iConfig.getParameter<edm::InputTag>("rhoInfo"));
  genParticlesToken_ = cc.consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
}

void JetInfoFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(puToken_, puInfo);
  iEvent.getByToken(rhoToken_, rhoInfo);
  event_ = iEvent.id().event();
  iEvent.getByToken(genParticlesToken_, genParticlesHandle);
  flavorDef.setGenParticles(*genParticlesHandle);

}

bool JetInfoFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {
  // pv selection
  if (vertices->empty()) return false;

  // jet selection
  if (jet.pt() < minPt_) return false;
  if (maxPt_ > 0 && jet.pt() > maxPt_) return false;
  if (std::abs(jet.eta()) > maxAbsEta_) return false;

  // event information
  data.fill<unsigned>("event_no", event_);
  data.fill<unsigned>("jet_no", jetidx);
  data.fill<float>("npv", vertices->size());
  data.fill<float>("rho", *rhoInfo);
  for (const auto &v : *puInfo) {
    int bx = v.getBunchCrossing();
    if (bx == 0) {
      data.fill<float>("ntrueInt", v.getTrueNumInteractions());
    }
  }

  // truth labels
  float gen_pt = jet.genJet() ? jet.genJet()->pt() : 0;
  data.fill<float>("gen_pt", gen_pt);

  auto flavor = flavorDef.jet_flavour(jet);
  data.fill<int>("isB", flavor==JetFlavor::B);
  data.fill<int>("isBB", flavor==JetFlavor::BB);
  data.fill<int>("isLeptonicB", flavor==JetFlavor::LeptonicB);
  data.fill<int>("isLeptonicB_C", flavor==JetFlavor::LeptonicB_C);
  data.fill<int>("isC", flavor==JetFlavor::C);
  data.fill<int>("isUD", flavor==JetFlavor::UD);
  data.fill<int>("isS", flavor==JetFlavor::S);
  data.fill<int>("isG", flavor==JetFlavor::G);
  data.fill<int>("isUndefined", flavor==JetFlavor::UNDEFINED);

  // jet variables
  data.fill<float>("jet_pt", jet.correctedJet("Uncorrected").pt());
  data.fill<float>("jet_corr_pt", jet.pt());
  data.fill<float>("jet_eta", jet.eta());
  data.fill<float>("jet_phi", jet.phi());

  // jet id
  data.fill<float>("jet_tightId", jetIdTight(jet));
  data.fill<float>("jet_tightIdLepVeto", jetIdTightLepVeto(jet));

  // qgl
  data.fill<float>("jet_qgl", jet.hasUserFloat("qgl") ? jet.userFloat("qgl") : 0);

  for(const auto& disc : btag_discriminators_) {
    std::string name(disc);
    std::replace(name.begin(), name.end(), ':', '_');
    data.fill<float>(name, catchInfs(jet.bDiscriminator(disc), -99));
  }

  return true;
}

void JetInfoFiller::book() {
  // event information
  data.add<float>("npv", 0);
  data.add<float>("rho", 0);
  data.add<float>("ntrueInt", 0);
  data.add<unsigned>("event_no", 0);
  data.add<unsigned>("jet_no", 0);

  // truth labels
  data.add<float>("gen_pt", 0);

  data.add<int>("isB", 0);
  data.add<int>("isBB", 0);
  data.add<int>("isLeptonicB", 0);
  data.add<int>("isLeptonicB_C", 0);
  data.add<int>("isC", 0);
  data.add<int>("isUD", 0);
  data.add<int>("isS", 0);
  data.add<int>("isG", 0);
  data.add<int>("isUndefined", 0);

  // jet variables
  data.add<float>("jet_pt", 0);
  data.add<float>("jet_corr_pt", 0);
  data.add<float>("jet_eta", 0);
  data.add<float>("jet_phi", 0);

  // jet id
  data.add<float>("jet_tightId", 0);
  data.add<float>("jet_tightIdLepVeto", 0);

  data.add<float>("jet_qgl", 0);

  for(auto name : btag_discriminators_) {
    std::replace(name.begin(), name.end(), ':', '_');
    data.add<float>(name, 0);
  }
}


} /* namespace deepntuples */

