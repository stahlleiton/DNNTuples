import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.outputFile = 'output.root'
options.inputFiles = 'root://eoscms.cern.ch//eos/cms//store/group/cmst3/group/hintt/Run3/MC/PbPb2023/Embedded/2024_04_19/POWHEG_5p36TeV_2023Run3/TT_hvq_POWHEG_Hydjet_5p36TeV_TuneCP5_2023Run3_MINIAOD_2024_04_19/240419_231333/0000/POWHEG_TT_hvq_MINIAOD_1.root'
options.maxEvents = 100#-1

options.register('skipEvents', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "skip N events")
options.register('inputDataset',
                 '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input dataset")
options.register('isTrainSample', True, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "if the sample is used for training")

options.parseArguments()

globalTagMap = {
    'auto': '132X_mcRun3_2023_realistic_HI_v10',
}

era = 'auto'
for k in globalTagMap:
    if k in options.inputDataset:
        era = k
# ---------------------------------------------------------
from Configuration.Eras.Era_Run3_pp_on_PbPb_2023_cff import Run3_pp_on_PbPb_2023
process = cms.Process("DNNFiller", Run3_pp_on_PbPb_2023)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
    allowUnscheduled=cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(True)
)

print('Using output file ' + options.outputFile)

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(options.outputFile))

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

process.source = cms.Source('PoolSource',
                            fileNames=cms.untracked.vstring(options.inputFiles),
                            skipEvents=cms.untracked.uint32(options.skipEvents)
                            )
# ---------------------------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTagMap[era], '')
print('Using global tag', process.GlobalTag.globaltag)
process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
                                                         ComponentName=cms.string('TransientTrackBuilder')
                                                         )

# Heavy-ion settings
from PhysicsTools.PatAlgos.producersHeavyIons.heavyIonJets_cff import PackedPFTowers, hiPuRho, hiSignalGenParticles, allPartons
process.PackedPFTowers = PackedPFTowers.clone()
process.hiPuRho = hiPuRho.clone(
    src = 'PackedPFTowers'
)
process.hiSignalGenParticles = hiSignalGenParticles.clone(
    src = "prunedGenParticles"
)
process.allPartons = allPartons.clone()
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJetsNoNu
process.ak4GenJetsNoNu = ak4GenJetsNoNu.clone(
    src = 'packedGenParticlesSignal',
    rParam = 0.4
)

# Create unsubtracted reco jets
from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff import ak4PFJets, patJetCorrFactors, patJetPartonMatch, patJetGenJetMatch, patJetPartons, patJetFlavourAssociation, patJets
process.ak4PFMatchingForakCs0PFJets = ak4PFJets.clone(
    src = 'packedPFCandidates',
    rParam = 0.4
)
process.ak4PFMatchingForakCs0PFpatJetCorrFactors = patJetCorrFactors.clone(
    src = 'ak4PFMatchingForakCs0PFJets',
    payload = 'AK4PF'
)
process.ak4PFMatchingForakCs0PFpatJetPartonMatch = patJetPartonMatch.clone(
    src = 'ak4PFMatchingForakCs0PFJets',
    matched = 'hiSignalGenParticles',
    maxDeltaR = 0.4
)
process.ak4PFMatchingForakCs0PFpatJetGenJetMatch = patJetGenJetMatch.clone(
    src = 'ak4PFMatchingForakCs0PFJets',
    matched = 'ak4GenJetsNoNu',
    maxDeltaR = 0.4
)
process.ak4PFMatchingForakCs0PFpatJetPartons = patJetPartons.clone(
    partonMode = 'Pythia8'
)
process.ak4PFMatchingForakCs0PFpatJetFlavourAssociation =  patJetFlavourAssociation.clone(
    jets = 'ak4PFMatchingForakCs0PFJets',
    rParam = 0.4,
    bHadrons = cms.InputTag("ak4PFMatchingForakCs0PFpatJetPartons","bHadrons"),
    cHadrons = cms.InputTag("ak4PFMatchingForakCs0PFpatJetPartons","cHadrons"),
    partons = cms.InputTag("ak4PFMatchingForakCs0PFpatJetPartons","physicsPartons"),
    leptons = cms.InputTag("ak4PFMatchingForakCs0PFpatJetPartons","leptons")
)
process.ak4PFMatchingForakCs0PFpatJets = patJets.clone(
    JetFlavourInfoSource = 'ak4PFMatchingForakCs0PFpatJetFlavourAssociation',
    JetPartonMapSource = 'ak4PFMatchingForakCs0PFpatJetFlavourAssociation',
    genJetMatch = 'ak4PFMatchingForakCs0PFpatJetGenJetMatch',
    genPartonMatch = 'ak4PFMatchingForakCs0PFpatJetPartonMatch',
    jetCorrFactorsSource = ['ak4PFMatchingForakCs0PFpatJetCorrFactors'],
    jetSource = 'ak4PFMatchingForakCs0PFJets',
    addBTagInfo = False,
    addDiscriminators = False,
    addAssociatedTracks = False,
    useLegacyJetMCFlavour = False
)

# Create HIN subtracted reco jets
from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff import akCs4PFJets, patJetPartonAssociationLegacy, patJetFlavourAssociationLegacy
process.akCs0PFJets = akCs4PFJets.clone(
    src = 'packedPFCandidates',
    useModulatedRho = False,
    rParam = 0.4
)
process.akCs0PFpatJetCorrFactors = patJetCorrFactors.clone(
    src = 'akCs0PFJets',
    payload = 'AK4PF'
)
process.akCs0PFpatJetPartonMatch = patJetPartonMatch.clone(
    src = 'akCs0PFJets',
    matched = 'hiSignalGenParticles',
    maxDeltaR = 0.4
)
process.akCs0PFpatJetGenJetMatch = patJetGenJetMatch.clone(
    src = 'akCs0PFJets',
    matched = 'ak4GenJetsNoNu',
    maxDeltaR = 0.4
)
process.akCs0PFpatJetPartonAssociationLegacy = patJetPartonAssociationLegacy.clone(
    jets = 'akCs0PFJets'
)
process.akCs0PFpatJetFlavourAssociationLegacy = patJetFlavourAssociationLegacy.clone(
    srcByReference = 'akCs0PFpatJetPartonAssociationLegacy'
)
process.akCs0PFpatJetPartons = patJetPartons.clone(
    partonMode = 'Pythia8'
)
from RecoBTag.ImpactParameter.pfImpactParameterTagInfos_cfi import pfImpactParameterTagInfos
process.akCs0PFpfImpactParameterTagInfos = pfImpactParameterTagInfos.clone(
    jets = 'akCs0PFJets',
    candidates = 'packedPFCandidates',
    primaryVertex = 'offlineSlimmedPrimaryVertices'
)
from RecoBTag.SecondaryVertex.pfSecondaryVertexTagInfos_cfi import pfSecondaryVertexTagInfos
process.akCs0PFpfSecondaryVertexTagInfos = pfSecondaryVertexTagInfos.clone(
    trackIPTagInfos = 'akCs0PFpfImpactParameterTagInfos'
)
from RecoBTag.Combined.pfDeepCSVTagInfos_cfi import pfDeepCSVTagInfos
process.akCs0PFpfDeepCSVTagInfos = pfDeepCSVTagInfos.clone(
    svTagInfos = 'akCs0PFpfSecondaryVertexTagInfos'
)
from RecoBTag.Combined.pfDeepCSVJetTags_cfi import pfDeepCSVJetTags
process.akCs0PFpfDeepCSVJetTags = pfDeepCSVJetTags.clone(
    src = 'akCs0PFpfDeepCSVTagInfos'
)
from RecoBTag.ImpactParameter.pfJetProbabilityBJetTags_cfi import pfJetProbabilityBJetTags
process.akCs0PFpfJetProbabilityBJetTags = pfJetProbabilityBJetTags.clone(
    tagInfos = ['akCs0PFpfImpactParameterTagInfos']
)
process.akCs0PFpatJets = patJets.clone(
    JetFlavourInfoSource = 'akCs0PFpatJetFlavourAssociation',
    JetPartonMapSource = 'akCs0PFpatJetFlavourAssociationLegacy',
    genJetMatch = 'akCs0PFpatJetGenJetMatch',
    genPartonMatch = 'akCs0PFpatJetPartonMatch',
    jetCorrFactorsSource = ['akCs0PFpatJetCorrFactors'],
    jetSource = 'akCs0PFJets',
    discriminatorSources = [cms.InputTag('akCs0PFpfDeepCSVJetTags','probb'), cms.InputTag('akCs0PFpfDeepCSVJetTags','probc'), cms.InputTag('akCs0PFpfDeepCSVJetTags','probudsg'), cms.InputTag('akCs0PFpfDeepCSVJetTags','probbb'), cms.InputTag('akCs0PFpfJetProbabilityBJetTags')],
    addAssociatedTracks = False,
)
# -------------
process.ak4PFMatchingForakCs0PFJets.jetPtMin = process.akCs0PFJets.jetPtMin
process.akCs0PFpatJetCorrFactors.levels = ['L2Relative', 'L3Absolute']
process.ak4PFMatchedForakCs0PFpatJets = cms.EDProducer("JetMatcherDR", source = cms.InputTag("akCs0PFpatJets"), matched = cms.InputTag("ak4PFMatchingForakCs0PFpatJets"))
process.akCs0PFpatJets.embedPFCandidates = True

# Start of b-tagging sequence ----------------
from RecoBTag.ImpactParameter.pfImpactParameterTagInfos_cfi import pfImpactParameterTagInfos
process.pfImpactParameterTagInfos = pfImpactParameterTagInfos.clone(
    jets = "akCs0PFpatJets",
    candidates = "packedPFCandidates",
    primaryVertex = "offlineSlimmedPrimaryVertices"
)
from RecoBTag.SecondaryVertex.pfSecondaryVertexTagInfos_cfi import pfSecondaryVertexTagInfos
process.pfSecondaryVertexTagInfos = pfSecondaryVertexTagInfos.clone()
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import inclusiveCandidateVertexFinder, candidateVertexMerger, candidateVertexArbitrator, inclusiveCandidateSecondaryVertices
process.inclusiveCandidateVertexFinder = inclusiveCandidateVertexFinder.clone(
    tracks = "packedPFCandidates",
    primaryVertices = "offlineSlimmedPrimaryVertices"
)
process.candidateVertexMerger = candidateVertexMerger.clone()
process.candidateVertexArbitrator = candidateVertexArbitrator.clone(
    tracks = "packedPFCandidates",
    primaryVertices = "offlineSlimmedPrimaryVertices"
)
process.inclusiveCandidateSecondaryVertices = inclusiveCandidateSecondaryVertices.clone()
from RecoBTag.SecondaryVertex.pfInclusiveSecondaryVertexFinderTagInfos_cfi import pfInclusiveSecondaryVertexFinderTagInfos
process.pfInclusiveSecondaryVertexFinderTagInfos = pfInclusiveSecondaryVertexFinderTagInfos.clone()

from RecoBTag.Combined.pfDeepCSVTagInfos_cfi import pfDeepCSVTagInfos
process.pfDeepCSVTagInfos = pfDeepCSVTagInfos.clone(
    svTagInfos = "pfSecondaryVertexTagInfos"
)
from RecoBTag.FeatureTools.pfDeepFlavourTagInfos_cfi import pfDeepFlavourTagInfos
process.pfDeepFlavourTagInfosSlimmedDeepFlavour = pfDeepFlavourTagInfos.clone(
    fallback_puppi_weight = True,
    fallback_vertex_association = True,
    jets = cms.InputTag("akCs0PFpatJets"),
    unsubjet_map = cms.InputTag("ak4PFMatchedForakCs0PFpatJets"),
    puppi_value_map = cms.InputTag(""),
    secondary_vertices = cms.InputTag("inclusiveCandidateSecondaryVertices"),
    shallow_tag_infos = cms.InputTag("pfDeepCSVTagInfos"),
    vertex_associator = cms.InputTag(""),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)
from RecoBTag.FeatureTools.pfParticleTransformerAK4TagInfos_cfi import pfParticleTransformerAK4TagInfos
process.pfParticleTransformerAK4TagInfosSlimmedDeepFlavour = pfParticleTransformerAK4TagInfos.clone(
    fallback_puppi_weight = True,
    fallback_vertex_association = True,
    jets = cms.InputTag("akCs0PFpatJets"),
    unsubjet_map = cms.InputTag("ak4PFMatchedForakCs0PFpatJets"),
    puppi_value_map = cms.InputTag(""),
    secondary_vertices = cms.InputTag("inclusiveCandidateSecondaryVertices"),
    vertex_associator = cms.InputTag(""),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)

from RecoBTag.ONNXRuntime.pfDeepFlavourJetTags_cfi import pfDeepFlavourJetTags
process.pfDeepFlavourJetTagsSlimmedDeepFlavour = pfDeepFlavourJetTags.clone(
    src = cms.InputTag("pfDeepFlavourTagInfosSlimmedDeepFlavour")
)
from RecoBTag.ONNXRuntime.pfParticleTransformerAK4JetTags_cfi import pfParticleTransformerAK4JetTags
process.pfParticleTransformerAK4JetTagsSlimmedDeepFlavour = pfParticleTransformerAK4JetTags.clone(
    src = cms.InputTag("pfParticleTransformerAK4TagInfosSlimmedDeepFlavour")
)

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
    process,
    jetSource = cms.InputTag('akCs0PFpatJets'),
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    btagDiscriminators = ['pfCombinedSecondaryVertexV2BJetTags', 'pfDeepCSVDiscriminatorsJetTags:BvsAll', 'pfDeepCSVDiscriminatorsJetTags:CvsB', 'pfDeepCSVDiscriminatorsJetTags:CvsL'], ## to add discriminators,
    btagPrefix = 'TEST',
)
process.updatedPatJets.addJetCorrFactors = False
process.updatedPatJets.discriminatorSources = cms.VInputTag(
    cms.InputTag("pfDeepFlavourJetTagsSlimmedDeepFlavour","probb"), cms.InputTag("pfDeepFlavourJetTagsSlimmedDeepFlavour","probbb"),   cms.InputTag("pfDeepFlavourJetTagsSlimmedDeepFlavour","probc"),
    cms.InputTag("pfDeepFlavourJetTagsSlimmedDeepFlavour","probg"), cms.InputTag("pfDeepFlavourJetTagsSlimmedDeepFlavour","problepb"), cms.InputTag("pfDeepFlavourJetTagsSlimmedDeepFlavour","probuds"),
    cms.InputTag("pfParticleTransformerAK4JetTagsSlimmedDeepFlavour","probb"), cms.InputTag("pfParticleTransformerAK4JetTagsSlimmedDeepFlavour","probbb"),   cms.InputTag("pfParticleTransformerAK4JetTagsSlimmedDeepFlavour","probc"),
    cms.InputTag("pfParticleTransformerAK4JetTagsSlimmedDeepFlavour","probg"), cms.InputTag("pfParticleTransformerAK4JetTagsSlimmedDeepFlavour","problepb"), cms.InputTag("pfParticleTransformerAK4JetTagsSlimmedDeepFlavour","probuds"),
)
# End of b-tagging sequence ----------------

bTagDiscriminators = [
    'pfDeepFlavourJetTagsSlimmedDeepFlavour:probb',
    'pfDeepFlavourJetTagsSlimmedDeepFlavour:probbb',
    'pfDeepFlavourJetTagsSlimmedDeepFlavour:probc',
    'pfDeepFlavourJetTagsSlimmedDeepFlavour:probg',
    'pfDeepFlavourJetTagsSlimmedDeepFlavour:problepb',
    'pfDeepFlavourJetTagsSlimmedDeepFlavour:probuds',
    'pfParticleTransformerAK4JetTagsSlimmedDeepFlavour:probb',
    'pfParticleTransformerAK4JetTagsSlimmedDeepFlavour:probbb',
    'pfParticleTransformerAK4JetTagsSlimmedDeepFlavour:probc',
    'pfParticleTransformerAK4JetTagsSlimmedDeepFlavour:probg',
    'pfParticleTransformerAK4JetTagsSlimmedDeepFlavour:problepb',
    'pfParticleTransformerAK4JetTagsSlimmedDeepFlavour:probuds',
]
# ------------------------------------------------------------#

srcJets = cms.InputTag('updatedPatJets')
process.unsubJetMap = cms.EDProducer("JetMatcherDR", source = srcJets, matched = cms.InputTag("ak4PFMatchingForakCs0PFpatJets"))
process.jetTask = cms.Task(
    process.PackedPFTowers,
    process.hiPuRho,
    process.hiSignalGenParticles,
    process.allPartons,
    process.ak4GenJetsNoNu,
    process.ak4PFMatchingForakCs0PFJets,
    process.ak4PFMatchingForakCs0PFpatJetCorrFactors,
    process.ak4PFMatchingForakCs0PFpatJetPartonMatch,
    process.ak4PFMatchingForakCs0PFpatJetGenJetMatch,
    process.ak4PFMatchingForakCs0PFpatJetPartons,
    process.ak4PFMatchingForakCs0PFpatJetFlavourAssociation,
    process.ak4PFMatchingForakCs0PFpatJets,
    process.akCs0PFJets,
    process.akCs0PFpatJetCorrFactors,
    process.akCs0PFpatJetPartonMatch,
    process.akCs0PFpatJetGenJetMatch,
    process.akCs0PFpatJetPartonAssociationLegacy,
    process.akCs0PFpatJetFlavourAssociationLegacy,
    process.akCs0PFpatJetPartons,
    process.akCs0PFpfImpactParameterTagInfos,
    process.akCs0PFpfSecondaryVertexTagInfos,
    process.akCs0PFpfDeepCSVTagInfos,
    process.akCs0PFpfDeepCSVJetTags,
    process.akCs0PFpfJetProbabilityBJetTags,
    process.akCs0PFpatJets,
    process.ak4PFMatchedForakCs0PFpatJets,
    process.ak4PFMatchedForakCs0PFpatJets,
    process.pfImpactParameterTagInfos,
    process.pfSecondaryVertexTagInfos,
    process.inclusiveCandidateVertexFinder,
    process.candidateVertexMerger,
    process.candidateVertexArbitrator,
    process.inclusiveCandidateSecondaryVertices,
    process.pfInclusiveSecondaryVertexFinderTagInfos,
    process.pfDeepCSVTagInfos,
    process.pfDeepFlavourTagInfosSlimmedDeepFlavour,
    process.pfParticleTransformerAK4TagInfosSlimmedDeepFlavour,
    process.pfDeepFlavourJetTagsSlimmedDeepFlavour,
    process.pfParticleTransformerAK4JetTagsSlimmedDeepFlavour,
    process.updatedPatJets,
    process.unsubJetMap
)

# ---------------------------------------------------------
# DeepNtuplizer
process.load("DeepNTuples.Ntupler.DeepNtuplizer_cfi")
process.deepntuplizer.jets = srcJets
process.deepntuplizer.isPuppiJets = False
process.deepntuplizer.bDiscriminators = bTagDiscriminators
process.deepntuplizer.unsubjet_map = cms.InputTag("unsubJetMap")

process.deepntuplizer.isQCDSample = '/QCD_' in options.inputDataset
process.deepntuplizer.isPythia = 'pythia' in options.inputDataset.lower()
process.deepntuplizer.isHerwig = 'herwig' in options.inputDataset.lower()
# note: MG can be interfaced w/ either pythia or herwig
process.deepntuplizer.isMadGraph = 'madgraph' in options.inputDataset.lower()

process.deepntuplizer.isTrainSample = options.isTrainSample
if not options.inputDataset:
    # interactive running
    process.deepntuplizer.isTrainSample = False
#==============================================================================================================================#
process.p = cms.Path(process.deepntuplizer)
process.p.associate(process.jetTask)
