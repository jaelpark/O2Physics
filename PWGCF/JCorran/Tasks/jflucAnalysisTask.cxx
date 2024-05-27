// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
/// \author Jasper Parkkila (jparkkil@cern.ch)
/// \author Dong Jo Kim (djkim@jyu.fi)
/// \since Sep 2022

#include <deque>

#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "ReconstructionDataFormats/V0.h"

// #include "CCDB/BasicCCDBManager.h"

#include "PWGCF/JCorran/DataModel/JCatalyst.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "JFFlucAnalysis.h"
#include "JFFlucAnalysisO2Hist.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct jflucAnalysisTask {
  ~jflucAnalysisTask()
  {
    if (pcf)
      delete pcf;
    if (pcf2Prong)
      delete pcf2Prong;
  }

  O2_DEFINE_CONFIGURABLE(etamin, float, 0.4, "Minimal eta for tracks");
  O2_DEFINE_CONFIGURABLE(etamax, float, 0.8, "Maximal eta for tracks");
  O2_DEFINE_CONFIGURABLE(ptmin, float, 0.2, "Minimal pt for tracks");
  O2_DEFINE_CONFIGURABLE(ptmax, float, 0.5, "Maximal pt for tracks");

  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity / centrality axis for histograms"};

  Filter jtrackFilter = (aod::jtrack::pt > ptmin) && (aod::jtrack::pt < ptmax);    // eta cuts done by jfluc
  Filter cftrackFilter = (aod::cftrack::pt > ptmin) && (aod::cftrack::pt < ptmax); // eta cuts done by jfluc
  Filter cf2pFilter = (aod::cf2prongtrack::pt > ptmin) && (aod::jtrack::pt < ptmax);

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    auto a = AxisSpec(axisMultiplicity);
    if (doprocessJDerived || doprocessJDerivedCorrected || doprocessCFDerived || doprocessCFDerivedCorrected) {
      pcf = new JFFlucAnalysisO2Hist(registry, a, "jfluc");
      pcf->AddFlags(JFFlucAnalysis::kFlucEbEWeighting);
      pcf->UserCreateOutputObjects();
    } else
      pcf = 0;
    if (doprocessCF2ProngDerived) {
      pcf2Prong = new JFFlucAnalysisO2Hist(registry, a, "jfluc2prong"); //<--- TODO: constructor, accept phietaz histogram from the other
      pcf2Prong->AddFlags(JFFlucAnalysis::kFlucEbEWeighting);
      pcf2Prong->UserCreateOutputObjects();
    } else
      pcf2Prong = 0;
  }

  template <class CollisionT, class TrackT>
  void analyze(CollisionT const& collision, TrackT const& tracks)
  {
    pcf->Init();
    pcf->SetEventCentrality(collision.multiplicity());
    pcf->SetEventVertex(collision.posZ());
    pcf->FillQA(tracks);
    qvecs.Calculate(tracks, etamin, etamax);
    pcf->SetJQVectors(&qvecs);
    pcf->UserExec("");
  }

  template <class CollisionT, class POITrackT, class REFTrackT>
  void analyze(CollisionT const& collision, POITrackT const& poiTracks, REFTrackT const& refTracks)
  {
    pcf2Prong->Init();
    pcf2Prong->SetEventCentrality(collision.multiplicity());
    pcf2Prong->SetEventVertex(collision.posZ());
    pcf2Prong->FillQA(poiTracks);
    qvecs.Calculate(poiTracks, etamin, etamax);
    qvecsRef.Calculate(refTracks, etamin, etamax);
    pcf2Prong->SetJQVectors(&qvecs, &qvecsRef);
    pcf2Prong->UserExec("");
  }

  void processJDerived(aod::JCollision const& collision, soa::Filtered<aod::JTracks> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processJDerived, "Process derived data", false);

  void processJDerivedCorrected(aod::JCollision const& collision, soa::Filtered<soa::Join<aod::JTracks, aod::JWeights>> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processJDerivedCorrected, "Process derived data with corrections", false);

  void processCFDerived(aod::CFCollision const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processCFDerived, "Process CF derived data", true);

  void processCFDerivedCorrected(aod::CFCollision const& collision, soa::Filtered<soa::Join<aod::CFTracks, aod::JWeights>> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processCFDerivedCorrected, "Process CF derived data with corrections", false);

  void processCF2ProngDerived(aod::CFCollision const& collision, soa::Filtered<aod::CFTracks> const& tracks, soa::Filtered<aod::CF2ProngTracks> const& p2tracks)
  {
    analyze(collision, p2tracks, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processCF2ProngDerived, "Process CF derived data with 2-prongs as POI and charged particles as REF.", false);

  JFFlucAnalysis::JQVectorsT qvecs;
  JFFlucAnalysis::JQVectorsT qvecsRef;
  JFFlucAnalysisO2Hist* pcf;
  JFFlucAnalysisO2Hist* pcf2Prong;
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<jflucAnalysisTask>(cfgc)};
}
