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
/// \author Dong Jo Kim (djkim@jyu.fi)
/// \author Jasper Parkkila (jparkkil@cern.ch)
/// \since Sep 2022

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF3.h>
#include <TMath.h>
#include <TComplex.h>
#include <algorithm>

#include "JFFlucAnalysis.h"

JFFlucAnalysis::JFFlucAnalysis() : fVertex(0),
                                   fCent(0),
                                   fImpactParameter(-1),
                                   subeventMask(kSubEvent_A | kSubEvent_B),
                                   flags(0),
                                   fEta_min(0),
                                   fEta_max(0)
{
  //
}

//________________________________________________________________________
JFFlucAnalysis::JFFlucAnalysis(const char* /*name*/) : fVertex(0),
                                                       fCent(0),
                                                       fImpactParameter(-1),
                                                       subeventMask(kSubEvent_A | kSubEvent_B),
                                                       flags(0),
                                                       fEta_min(0),
                                                       fEta_max(0)
{
  //
}

//________________________________________________________________________
JFFlucAnalysis::JFFlucAnalysis(const JFFlucAnalysis& a) : fVertex(a.fVertex),
                                                          fCent(a.fCent),
                                                          fImpactParameter(a.fImpactParameter),
                                                          subeventMask(a.subeventMask),
                                                          flags(a.flags),
                                                          fEta_min(a.fEta_min),
                                                          fEta_max(a.fEta_max),
                                                          pqvecs(a.pqvecs)
{
  // copy constructor
}
//________________________________________________________________________
JFFlucAnalysis& JFFlucAnalysis::operator=(const JFFlucAnalysis& ap)
{
  // assignment operator
  this->~JFFlucAnalysis();
  new (this) JFFlucAnalysis(ap);
  return *this;
}
//________________________________________________________________________
void JFFlucAnalysis::Init()
{
  //
}
//________________________________________________________________________
void JFFlucAnalysis::UserCreateOutputObjects()
{
  // const double phiBins[] = {-2
  // new(&fh_phietaz) THnF("phietaz","phietaz",4,numBins,pbins,
#if 0
  fHMG = new JHistManager("JFFlucHistManager", "jfluc");
  // set JBin here //
  fBin_Subset.Set("Sub", "Sub", "Sub:%d", JBin::kSingle).SetBin(2);
  fBin_h.Set("NH", "NH", "NH:%d", JBin::kSingle).SetBin(kNH);
  fBin_k.Set("K", "K", "K:%d", JBin::kSingle).SetBin(nKL);

  fBin_hh.Set("NHH", "NHH", "NHH:%d", JBin::kSingle).SetBin(kcNH);
  fBin_kk.Set("KK", "KK", "KK:%d", JBin::kSingle).SetBin(nKL);

  // TODO: index with binning the array of pointers
  fHistCentBin.Set("CentBin", "CentBin", "Cent:%d", JBin::kSingle).SetBin(numBins);

  fVertexBin.Set("Vtx", "Vtx", "Vtx:%d", JBin::kSingle).SetBin(3);
  fCorrBin.Set("C", "C", "C:%d", JBin::kSingle).SetBin(28);

  fBin_Nptbins.Set("PtBin", "PtBin", "Pt:%d", JBin::kSingle).SetBin(N_ptbins);

  // set JTH1D here //
  fh_cent
    << TH1D("h_cent", "h_cent", 200, 0, 100)
    << "END";

  fh_ImpactParameter
    << TH1D("h_IP", "h_IP", 400, -2, 20)
    << "END";

  fh_TrkQA_TPCvsCent
    << TH2D("h_trk_Cent_vs_TPC", "h_trk_Cent_vs_TPC", 100, 0, 100, 100, 0, 3000)
    << "END";

  fh_TrkQA_TPCvsGlob
    << TH2D("h_trk_Glob_vs_TPC", "h_trk_Glob_vs_TPC", 100, 0, 2000, 100, 0, 3000)
    << "END";

  fh_TrkQA_FB32_vs_FB32TOF
    << TH2D("h_trk_FB32_vs_FB32TOF", "h_trk_FB32_vs_FB32TOF", 200, 0, 4000, 100, 0, 2000)
    << "END";

  fh_vertex
    << TH1D("h_vertex", "h_vertex", 400, -20, 20)
    << fVertexBin
    << "END";
  fh_pt
    << TH1D("hChargedPtJacek", "", JFFlucAnalysis::NpttJacek, JFFlucAnalysis::pttJacek)
    << fHistCentBin
    << "END";

  fh_eta
    << TH1D("h_eta", "h_eta", 40, -2.0, 2.0)
    << fHistCentBin
    << "END";
  fh_phi
    << TH1D("h_phi", "h_phi", 50, -TMath::Pi(), TMath::Pi())
    << fHistCentBin << fBin_Subset
    << "END";
  if (!(flags & kFlucPhiCorrection)) {
    fh_phieta
      << TH2D("h_phieta", "h_phieta", 50, -TMath::Pi(), TMath::Pi(), 40, -2.0, 2.0)
      << fHistCentBin
      << "END";
    fh_phietaz
      << TH3D("h_phietaz", "h_phietaz", 50, -TMath::Pi(), TMath::Pi(), 40, -2.0, 2.0, 20, -10.0, 10.0)
      << fHistCentBin
      << "END";
  }

  fh_psi_n
    << TH1D("h_psi", "h_psi", 50, -0.5 * TMath::Pi(), 0.5 * TMath::Pi())
    << fBin_h
    << fHistCentBin
    << "END";

  fh_cos_n_phi
    << TH1D("h_cos_n_phi", "h_cos_n_phi", 50, -1, 1)
    << fBin_h
    << fHistCentBin
    << "END";

  fh_sin_n_phi
    << TH1D("h_sin_n_phi", "h_sin_n_phi", 50, -1, 1)
    << fBin_h
    << fHistCentBin
    << "END";

  fh_cos_n_psi_n
    << TH1D("h_cos_n_psi", "h_cos_n_psi", 50, -1, 1)
    << fBin_h
    << fHistCentBin
    << "END";

  fh_sin_n_psi_n
    << TH1D("h_sin_n_psi", "h_sin_n_psi", 50, -1, 1)
    << fBin_h
    << fHistCentBin
    << "END";

  fh_ntracks
    << TH1D("h_tracks", "h_tracks", 100, 0, 5000)
    << fHistCentBin
    << "END";

  fh_vn
    << TH1D("hvn", "hvn", 1024, -1.0, 1.0)
    << fBin_h << fBin_k
    << fHistCentBin
    << "END"; // histogram of vn_h^k values for [ih][ik][iCent]
  fh_vna
    << TH1D("hvna", "hvna", 1024, -1.0, 1.0)
    << fBin_h << fBin_k
    << fHistCentBin
    << "END"; // histogram of vn_h^k values for [ih][ik][iCent]
  fh_vn_vn
    << TH1D("hvn_vn", "hvn_vn", 1024, -1.0, 1.0)
    << fBin_h << fBin_k
    << fBin_hh << fBin_kk
    << fHistCentBin
    << "END"; // histo of < vn * vn > for [ih][ik][ihh][ikk][iCent]
  fh_correlator
    << TH1D("h_corr", "h_corr", 1024, -3.0, 3.0)
    << fCorrBin
    << fHistCentBin
    << "END";

  fh_SC_ptdep_4corr
    << TH1D("hvnvm_SC", "hvnvm_SC", 1024, -1.5, 1.5)
    << fBin_h << fBin_k
    << fBin_hh << fBin_kk
    << fHistCentBin << fBin_Nptbins
    << "END";
  fh_SC_ptdep_2corr
    << TH1D("hvn_SC", "hvn_SC", 1024, -1.5, 1.5)
    << fBin_h << fBin_k
    << fHistCentBin << fBin_Nptbins
    << "END";

  fh_SC_with_QC_4corr
    << TH1D("hQC_SC4p", "hQC_SC4p", 1024, -1.5, 1.5)
    << fBin_h << fBin_hh
    << fHistCentBin
    << "END";
  fh_SC_with_QC_2corr
    << TH1D("hQC_SC2p", "hQC_SC2p", 1024, -1.5, 1.5)
    << fBin_h
    << fHistCentBin
    << "END";
  fh_SC_with_QC_2corr_gap
    << TH1D("hQC_SC2p_eta10", "hQC_SC2p_eta10", 1024, -1.5, 1.5)
    << fBin_h
    << fHistCentBin
    << "END";
  /*fh_evt_SP_QC_ratio_4p
    << TH1D("hSPQCratio4p", "hSPQCratio4p", 1024, -100, 100)
    << fBin_h
    << fHistCentBin
    << "END"; // fBin_h > not stand for harmonics, just case(32, 42, 52, 53, 43)

  fh_evt_SP_QC_ratio_2p
    << TH1D("hSPQCratio", "hSPQCratio", 1024, -100, 100)
    << fBin_h
    << fHistCentBin
    << "END"; // fBin_h > not stand for harmonics, only for v2, v3, v4, v5*/
  // JTH1D set done.

  fHMG->Print();
  // fHMG->WriteConfig();
#endif
}

//________________________________________________________________________
JFFlucAnalysis::~JFFlucAnalysis()
{
  //
}

#define A i
#define B (1 - i)
#define C(u) TComplex::Conjugate(u)
// TODO: conjugate macro
inline TComplex TwoGap(const TComplex (*pQq)[JFFlucAnalysis::kNH][JFFlucAnalysis::nKL], uint i, uint a, uint b)
{
  return pQq[A][a][1] * C(pQq[B][b][1]);
}

inline TComplex ThreeGap(const TComplex (*pQq)[JFFlucAnalysis::kNH][JFFlucAnalysis::nKL], uint i, uint a, uint b, uint c)
{
  return pQq[A][a][1] * C(pQq[B][b][1] * pQq[B][c][1] - pQq[B][b + c][2]);
}

inline TComplex FourGap22(const TComplex (*pQq)[JFFlucAnalysis::kNH][JFFlucAnalysis::nKL], uint i, uint a, uint b, uint c, uint d)
{
  return pQq[A][a][1] * pQq[A][b][1] * C(pQq[B][c][1] * pQq[B][d][1]) - pQq[A][a + b][2] * C(pQq[B][c][1] * pQq[B][d][1]) - pQq[A][a][1] * pQq[A][b][1] * C(pQq[B][c + d][2]) + pQq[A][a + b][2] * C(pQq[B][c + d][2]);
}

inline TComplex FourGap13(const TComplex (*pQq)[JFFlucAnalysis::kNH][JFFlucAnalysis::nKL], uint i, uint a, uint b, uint c, uint d)
{
  return pQq[A][a][1] * C(pQq[B][b][1] * pQq[B][c][1] * pQq[B][d][1] - pQq[B][b + c][2] * pQq[B][d][1] - pQq[B][b + d][2] * pQq[B][c][1] - pQq[B][c + d][2] * pQq[B][b][1] + 2.0 * pQq[B][b + c + d][3]);
}

inline TComplex SixGap33(const TComplex (*pQq)[JFFlucAnalysis::kNH][JFFlucAnalysis::nKL], uint i, uint n1, uint n2, uint n3, uint n4, uint n5, uint n6)
{
  return pQq[A][n1][1] * pQq[A][n2][1] * pQq[A][n3][1] * C(pQq[B][n4][1] * pQq[B][n5][1] * pQq[B][n6][1]) - pQq[A][n1][1] * pQq[A][n2][1] * pQq[A][n3][1] * C(pQq[B][n4 + n5][2] * pQq[B][n6][1]) - pQq[A][n1][1] * pQq[A][n2][1] * pQq[A][n3][1] * C(pQq[B][n4 + n6][2] * pQq[B][n5][1]) - pQq[A][n1][1] * pQq[A][n2][1] * pQq[A][n3][1] * C(pQq[B][n5 + n6][2] * pQq[B][n4][1]) + 2.0 * pQq[A][n1][1] * pQq[A][n2][1] * pQq[A][n3][1] * C(pQq[B][n4 + n5 + n6][3]) - pQq[A][n1 + n2][2] * pQq[A][n3][1] * C(pQq[B][n4][1] * pQq[B][n5][1] * pQq[B][n6][1]) + pQq[A][n1 + n2][2] * pQq[A][n3][1] * C(pQq[B][n4 + n5][2] * pQq[B][n6][1]) + pQq[A][n1 + n2][2] * pQq[A][n3][1] * C(pQq[B][n4 + n6][2] * pQq[B][n5][1]) + pQq[A][n1 + n2][2] * pQq[A][n3][1] * C(pQq[B][n5 + n6][2] * pQq[B][n4][1]) - 2.0 * pQq[A][n1 + n2][2] * pQq[A][n3][1] * C(pQq[B][n4 + n5 + n6][3]) - pQq[A][n1 + n3][2] * pQq[A][n2][1] * C(pQq[B][n4][1] * pQq[B][n5][1] * pQq[B][n6][1]) + pQq[A][n1 + n3][2] * pQq[A][n2][1] * C(pQq[B][n4 + n5][2] * pQq[B][n6][1]) + pQq[A][n1 + n3][2] * pQq[A][n2][1] * C(pQq[B][n4 + n6][2] * pQq[B][n5][1]) + pQq[A][n1 + n3][2] * pQq[A][n2][1] * C(pQq[B][n5 + n6][2] * pQq[B][n4][1]) - 2.0 * pQq[A][n1 + n3][2] * pQq[A][n2][1] * C(pQq[B][n4 + n5 + n6][3]) - pQq[A][n2 + n3][2] * pQq[A][n1][1] * C(pQq[B][n4][1] * pQq[B][n5][1] * pQq[B][n6][1]) + pQq[A][n2 + n3][2] * pQq[A][n1][1] * C(pQq[B][n4 + n5][2] * pQq[B][n6][1]) + pQq[A][n2 + n3][2] * pQq[A][n1][1] * C(pQq[B][n4 + n6][2] * pQq[B][n5][1]) + pQq[A][n2 + n3][2] * pQq[A][n1][1] * C(pQq[B][n5 + n6][2] * pQq[B][n4][1]) - 2.0 * pQq[A][n2 + n3][2] * pQq[A][n1][1] * C(pQq[B][n4 + n5 + n6][3]) + 2.0 * pQq[A][n1 + n2 + n3][3] * C(pQq[B][n4][1] * pQq[B][n5][1] * pQq[B][n6][1]) - 2.0 * pQq[A][n1 + n2 + n3][3] * C(pQq[B][n4 + n5][2] * pQq[B][n6][1]) - 2.0 * pQq[A][n1 + n2 + n3][3] * C(pQq[B][n4 + n6][2] * pQq[B][n5][1]) - 2.0 * pQq[A][n1 + n2 + n3][3] * C(pQq[B][n5 + n6][2] * pQq[B][n4][1]) + 4.0 * pQq[A][n1 + n2 + n3][3] * C(pQq[B][n4 + n5 + n6][3]);
}

TComplex JFFlucAnalysis::Q(int n, int p)
{
  // Return QvectorQC
  // Q{-n, p} = Q{n, p}*
  return n >= 0 ? pqvecs->QvectorQC[n][p] : C(pqvecs->QvectorQC[-n][p]);
}

TComplex JFFlucAnalysis::Two(int n1, int n2)
{
  // two-particle correlation <exp[i(n1*phi1 + n2*phi2)]>
  return Q(n1, 1) * Q(n2, 1) - Q(n1 + n2, 2);
}

TComplex JFFlucAnalysis::Four(int n1, int n2, int n3, int n4)
{

  return Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) - Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) - Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) - Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) + 2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) - Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) + Q(n2 + n3, 2) * Q(n1 + n4, 2) - Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) + Q(n1 + n3, 2) * Q(n2 + n4, 2) + 2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) - Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) + Q(n1 + n2, 2) * Q(n3 + n4, 2) + 2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) + 2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) - 6. * Q(n1 + n2 + n3 + n4, 4);
}
#undef C

//________________________________________________________________________
void JFFlucAnalysis::UserExec(Option_t* /*popt*/)
{
#if 0
  for (UInt_t ih = 2; ih < kNH; ih++) {
    fh_cos_n_phi[ih][fCBin]->Fill(pqvecs->QvectorQC[ih][1].Re() / pqvecs->QvectorQC[0][1].Re());
    fh_sin_n_phi[ih][fCBin]->Fill(pqvecs->QvectorQC[ih][1].Im() / pqvecs->QvectorQC[0][1].Re());
    //
    //
    Double_t psi = pqvecs->QvectorQC[ih][1].Theta();
    fh_psi_n[ih][fCBin]->Fill(psi);
    fh_cos_n_psi_n[ih][fCBin]->Fill(TMath::Cos((Double_t)ih * psi));
    fh_sin_n_psi_n[ih][fCBin]->Fill(TMath::Sin((Double_t)ih * psi));
  }
#endif
  Double_t vn2_vn2[kNH][nKL][kNH][nKL];

  TComplex corr[kNH][nKL];
  TComplex ncorr[kNH][nKL];
  TComplex ncorr2[kNH][nKL][kcNH][nKL];

  const TComplex(*pQq)[kNH][nKL] = pqvecs->QvectorQCgap;

  for (UInt_t i = 0; i < 2; ++i) {
    if ((subeventMask & (1 << i)) == 0)
      continue;
    Double_t ref_2p = TwoGap(pQq, i, 0, 0).Re();
    Double_t ref_3p = ThreeGap(pQq, i, 0, 0, 0).Re();
    Double_t ref_4p = FourGap22(pQq, i, 0, 0, 0, 0).Re();
    Double_t ref_4pB = FourGap13(pQq, i, 0, 0, 0, 0).Re();
    Double_t ref_6p = SixGap33(pQq, i, 0, 0, 0, 0, 0, 0).Re();

    Double_t ebe_2p_weight = 1.0;
    Double_t ebe_3p_weight = 1.0;
    Double_t ebe_4p_weight = 1.0;
    Double_t ebe_4p_weightB = 1.0;
    Double_t ebe_6p_weight = 1.0;
    if (flags & kFlucEbEWeighting) {
      ebe_2p_weight = ref_2p;
      ebe_3p_weight = ref_3p;
      ebe_4p_weight = ref_4p;
      ebe_4p_weightB = ref_4pB;
      ebe_6p_weight = ref_6p;
    }
    Double_t ref_2Np[2 * nKL] = {
      ref_2p,
      ref_4p,
      ref_6p};
    Double_t ebe_2Np_weight[2 * nKL] = {
      ebe_2p_weight,
      ebe_4p_weight,
      ebe_6p_weight};
    if (flags & kFlucEbEWeighting) {
      for (UInt_t ik = 3; ik < 2 * nKL; ik++) {
        double dk = static_cast<double>(ik);
        ref_2Np[ik] = ref_2Np[ik - 1] * std::max(pQq[A][0][1].Re() - dk, 1.0) * std::max(pQq[B][0][1].Re() - dk, 1.0);
        ebe_2Np_weight[ik] = ebe_2Np_weight[ik - 1] * std::max(pQq[A][0][1].Re() - dk, 1.0) * std::max(pQq[B][0][1].Re() - dk, 1.0);
      }
    } else {
      for (UInt_t ik = 3; ik < 2 * nKL; ik++) {
        double dk = static_cast<double>(ik);
        ref_2Np[ik] = ref_2Np[ik - 1] * std::max(pQq[A][0][1].Re() - dk, 1.0) * std::max(pQq[B][0][1].Re() - dk, 1.0);
        ebe_2Np_weight[ik] = 1.0;
      }
    }

    for (UInt_t ih = 2; ih < kNH; ih++) {
      corr[ih][1] = TwoGap(pQq, i, ih, ih);
      for (UInt_t ik = 2; ik < nKL; ik++)
        corr[ih][ik] = corr[ih][ik - 1] * corr[ih][1]; // TComplex::Power(corr[ih][1],ik);
      ncorr[ih][1] = corr[ih][1];
      ncorr[ih][2] = FourGap22(pQq, i, ih, ih, ih, ih);
      ncorr[ih][3] = SixGap33(pQq, i, ih, ih, ih, ih, ih, ih);
      for (UInt_t ik = 4; ik < nKL; ik++)
        ncorr[ih][ik] = corr[ih][ik]; // for 8,...-particle correlations, ignore the autocorrelation / weight dependency for now

      for (UInt_t ihh = 2; ihh < kcNH; ihh++) {
        ncorr2[ih][1][ihh][1] = FourGap22(pQq, i, ih, ihh, ih, ihh);
        ncorr2[ih][1][ihh][2] = SixGap33(pQq, i, ih, ihh, ihh, ih, ihh, ihh);
        ncorr2[ih][2][ihh][1] = SixGap33(pQq, i, ih, ih, ihh, ih, ih, ihh);
        for (UInt_t ik = 2; ik < nKL; ik++)
          for (UInt_t ikk = 2; ikk < nKL; ikk++)
            ncorr2[ih][ik][ihh][ikk] = ncorr[ih][ik] * ncorr[ihh][ikk];
      }
    }

    for (UInt_t ih = 2; ih < kNH; ih++) {
      for (UInt_t ik = 1; ik < nKL; ik++) { // 2k(0) =1, 2k(1) =2, 2k(2)=4....
                                            // vn2[ih][ik] = corr[ih][ik].Re() / ref_2Np[ik - 1];
        // fh_vn[ih][ik][fCBin]->Fill(vn2[ih][ik], ebe_2Np_weight[ik - 1]);
        // fh_vna[ih][ik][fCBin]->Fill(ncorr[ih][ik].Re() / ref_2Np[ik - 1], ebe_2Np_weight[ik - 1]);
        pht[HIST_THN_VN]->Fill(fCent, ih, ik, ncorr[ih][ik].Re() / ref_2Np[ik - 1], ebe_2Np_weight[ik - 1]);
        for (UInt_t ihh = 2; ihh < kcNH; ihh++) {
          for (UInt_t ikk = 1; ikk < nKL; ikk++) {
            vn2_vn2[ih][ik][ihh][ikk] = ncorr2[ih][ik][ihh][ikk] / ref_2Np[ik + ikk - 1];
            // fh_vn_vn[ih][ik][ihh][ikk][fCBin]->Fill(vn2_vn2[ih][ik][ihh][ikk], ebe_2Np_weight[ik + ikk - 1]);
            pht[HIST_THN_VN_VN]->Fill(fCent, ih, ik, ihh, ikk, vn2_vn2[ih][ik][ihh][ikk], ebe_2Np_weight[ik + ikk - 1]);
          }
        }
      }
    }

    //************************************************************************
    TComplex V4V2star_2 = pQq[A][4][1] * pQq[B][2][1] * pQq[B][2][1];
    TComplex V4V2starv2_2 = V4V2star_2 * corr[2][1] / ref_2Np[0];                                       // vn[2][1]
    TComplex V4V2starv2_4 = V4V2star_2 * corr[2][2] / ref_2Np[1];                                       // vn2[2][2]
    TComplex V5V2starV3starv2_2 = pQq[A][5][1] * pQq[B][2][1] * pQq[B][3][1] * corr[2][1] / ref_2Np[0]; // vn2[2][1]
    TComplex V5V2starV3star = pQq[A][5][1] * pQq[B][2][1] * pQq[B][3][1];
    TComplex V5V2starV3startv3_2 = V5V2starV3star * corr[3][1] / ref_2Np[0]; // vn2[3][1]
    TComplex V6V2star_3 = pQq[A][6][1] * pQq[B][2][1] * pQq[B][2][1] * pQq[B][2][1];
    TComplex V6V3star_2 = pQq[A][6][1] * pQq[B][3][1] * pQq[B][3][1];
    TComplex V6V2starV4star = pQq[A][6][1] * pQq[B][2][1] * pQq[B][4][1];
    TComplex V7V2star_2V3star = pQq[A][7][1] * pQq[B][2][1] * pQq[B][2][1] * pQq[B][3][1];
    TComplex V7V2starV5star = pQq[A][7][1] * pQq[B][2][1] * pQq[B][5][1];
    TComplex V7V3starV4star = pQq[A][7][1] * pQq[B][3][1] * pQq[B][4][1];
    TComplex V8V2starV3star_2 = pQq[A][8][1] * pQq[B][2][1] * pQq[B][3][1] * pQq[B][3][1];
    TComplex V8V2star_4 = pQq[A][8][1] * TComplex::Power(pQq[B][2][1], 4);

    // New correlators (Modified by You's correction term for self-correlations)
    TComplex nV4V2star_2 = ThreeGap(pQq, i, 4, 2, 2) / ref_3p;
    TComplex nV5V2starV3star = ThreeGap(pQq, i, 5, 2, 3) / ref_3p;
    TComplex nV6V2star_3 = FourGap13(pQq, i, 6, 2, 2, 2) / ref_4pB;
    TComplex nV6V3star_2 = ThreeGap(pQq, i, 6, 3, 3) / ref_3p;
    TComplex nV6V2starV4star = ThreeGap(pQq, i, 6, 2, 4) / ref_3p;
    TComplex nV7V2star_2V3star = FourGap13(pQq, i, 7, 2, 2, 3) / ref_4pB;
    TComplex nV7V2starV5star = ThreeGap(pQq, i, 7, 2, 5) / ref_3p;
    TComplex nV7V3starV4star = ThreeGap(pQq, i, 7, 3, 4) / ref_3p;
    TComplex nV8V2starV3star_2 = FourGap13(pQq, i, 8, 2, 3, 3) / ref_4pB;

    TComplex nV4V4V2V2 = FourGap22(pQq, i, 4, 2, 4, 2) / ref_4p;
    TComplex nV3V3V2V2 = FourGap22(pQq, i, 3, 2, 3, 2) / ref_4p;
    TComplex nV5V5V2V2 = FourGap22(pQq, i, 5, 2, 5, 2) / ref_4p;
    TComplex nV5V5V3V3 = FourGap22(pQq, i, 5, 3, 5, 3) / ref_4p;
    TComplex nV4V4V3V3 = FourGap22(pQq, i, 4, 3, 4, 3) / ref_4p;

    pht[HIST_THN_V4V2starv2_2]->Fill(fCent, V4V2starv2_2.Re());
    pht[HIST_THN_V4V2starv2_4]->Fill(fCent, V4V2starv2_4.Re());
    pht[HIST_THN_V4V2star_2]->Fill(fCent, V4V2star_2.Re(), ebe_3p_weight); // added 2015.3.18
    pht[HIST_THN_V5V2starV3starv2_2]->Fill(fCent, V5V2starV3starv2_2.Re());
    pht[HIST_THN_V5V2starV3star]->Fill(fCent, V5V2starV3star.Re(), ebe_3p_weight);
    pht[HIST_THN_V5V2starV3startv3_2]->Fill(fCent, V5V2starV3startv3_2.Re());
    pht[HIST_THN_V6V2star_3]->Fill(fCent, V6V2star_3.Re(), ebe_4p_weightB);
    pht[HIST_THN_V6V3star_2]->Fill(fCent, V6V3star_2.Re(), ebe_3p_weight);
    pht[HIST_THN_V7V2star_2V3star]->Fill(fCent, V7V2star_2V3star.Re(), ebe_4p_weightB);

    pht[HIST_THN_V4V2star_2]->Fill(fCent, nV4V2star_2.Re(), ebe_3p_weight); // added 2015.6.10
    pht[HIST_THN_V5V2starV3star]->Fill(fCent, nV5V2starV3star.Re(), ebe_3p_weight);
    pht[HIST_THN_V6V3star_2]->Fill(fCent, nV6V3star_2.Re(), ebe_3p_weight);

    // use this to avoid self-correlation 4p correlation (2 particles from A, 2 particles from B) -> MA(MA-1)MB(MB-1) : evt weight..
    pht[HIST_THN_nV4V4V2V2]->Fill(fCent, nV4V4V2V2.Re(), ebe_2Np_weight[1]);
    pht[HIST_THN_nV3V3V2V2]->Fill(fCent, nV3V3V2V2.Re(), ebe_2Np_weight[1]);

    pht[HIST_THN_nV5V5V2V2]->Fill(fCent, nV5V5V2V2.Re(), ebe_2Np_weight[1]);
    pht[HIST_THN_nV5V5V3V3]->Fill(fCent, nV5V5V3V3.Re(), ebe_2Np_weight[1]);
    pht[HIST_THN_nV4V4V3V3]->Fill(fCent, nV4V4V3V3.Re(), ebe_2Np_weight[1]);

    // higher order correlators, added 2017.8.10
    pht[HIST_THN_V8V2starV3star_2]->Fill(fCent, V8V2starV3star_2.Re(), ebe_4p_weightB);
    pht[HIST_THN_V8V2star_4]->Fill(fCent, V8V2star_4.Re()); // 5p weight
    pht[HIST_THN_V6V2star_3]->Fill(fCent, nV6V2star_3.Re(), ebe_4p_weightB);
    pht[HIST_THN_V7V2star_2V3star]->Fill(fCent, nV7V2star_2V3star.Re(), ebe_4p_weightB);
    pht[HIST_THN_V8V2starV3star_2]->Fill(fCent, nV8V2starV3star_2.Re(), ebe_4p_weightB);

    pht[HIST_THN_V6V2starV4star]->Fill(fCent, V6V2starV4star.Re(), ebe_3p_weight);
    pht[HIST_THN_V7V2starV5star]->Fill(fCent, V7V2starV5star.Re(), ebe_3p_weight);
    pht[HIST_THN_V7V3starV4star]->Fill(fCent, V7V3starV4star.Re(), ebe_3p_weight);
    pht[HIST_THN_V6V2starV4star]->Fill(fCent, nV6V2starV4star.Re(), ebe_3p_weight);
    pht[HIST_THN_V7V2starV5star]->Fill(fCent, nV7V2starV5star.Re(), ebe_3p_weight);
    pht[HIST_THN_V7V3starV4star]->Fill(fCent, nV7V3starV4star.Re(), ebe_3p_weight);

#if 0
  enum { kSubA,
         kSubB,
         kNSub };

  Double_t event_weight_four = 1.0;
  Double_t event_weight_two = 1.0;
  Double_t event_weight_two_gap = 1.0;
  if (flags & kFlucEbEWeighting) {
    event_weight_four = Four(0, 0, 0, 0).Re();
    event_weight_two = Two(0, 0).Re();
    event_weight_two_gap = (pqvecs->QvectorQCgap[kSubA][0][1] * pqvecs->QvectorQCgap[kSubB][0][1]).Re();
  }

  for (UInt_t ih = 2; ih < kNH; ih++) {
    for (UInt_t ihh = 2, mm = (ih < kcNH ? ih : kcNH); ihh < mm; ihh++) {
      TComplex scfour = Four(ih, ihh, -ih, -ihh) / Four(0, 0, 0, 0).Re();

      fh_SC_with_QC_4corr[ih][ihh][fCBin]->Fill(scfour.Re(), event_weight_four);
    }

    TComplex sctwo = Two(ih, -ih) / Two(0, 0).Re();
    fh_SC_with_QC_2corr[ih][fCBin]->Fill(sctwo.Re(), event_weight_two);

    TComplex sctwoGap = (pqvecs->QvectorQCgap[kSubA][ih][1] * TComplex::Conjugate(pqvecs->QvectorQCgap[kSubB][ih][1])) / (pqvecs->QvectorQCgap[kSubA][0][1] * pqvecs->QvectorQCgap[kSubB][0][1]).Re();
    fh_SC_with_QC_2corr_gap[ih][fCBin]->Fill(sctwoGap.Re(), event_weight_two_gap);
#endif
  }
}

//________________________________________________________________________
void JFFlucAnalysis::Terminate(Option_t* /*popt*/)
{
  //
}
