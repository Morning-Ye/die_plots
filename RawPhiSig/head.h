#include "TH1.h"

TFile *fin;
TH2D* H2d_Ne;
TH1D* HePlusPt;
TH1D* HeMinusPt;

TH2D* H2d_US;
TH2D* H2d_LSPos;
TH2D* H2d_LSNeg;
TH2D* H2d_LS;
TH2D* H2d_Sig;
TH1D* Hmass_US;
TH1D* Hmass_LS;
TH1D* Hmass_Sig;
TH1D* Hpt_Sig;
TH2D* H2d_PhoeUS;
TH2D* H2d_PhoeLSPos;
TH2D* H2d_PhoeLSNeg;
TH2D* H2d_PhoeLS;
TH2D* H2d_PhoeSig;
TH1D* Hmass_PhoeUS;
TH1D* Hmass_PhoeLS;
TH1D* Hmass_PhoeSig;

TH2D* H2d_PhiVm_US;
TH2D* H2d_PhiVm_LSPos;
TH2D* H2d_PhiVm_LSNeg;
TH2D* H2d_PhiVm_LS;
TH2D* H2d_PhiVm_Sig;

TH2D* H2d_wPhiV_US;
TH2D* H2d_wPhiV_LSPos;
TH2D* H2d_wPhiV_LSNeg;
TH2D* H2d_wPhiV_LS;
TH2D* H2d_wPhiV_Sig;
TH1D* Hmass_wPhiV_US;
TH1D* Hmass_wPhiV_LS;
TH1D* Hmass_wPhiV_Sig;
TH1D* Hpt_wPhiV_Sig;
TH2D* H2d_wPhiV_PhoeUS;
TH2D* H2d_wPhiV_PhoeLSPos;
TH2D* H2d_wPhiV_PhoeLSNeg;
TH2D* H2d_wPhiV_PhoeLS;
TH2D* H2d_wPhiV_PhoeSig;
TH1D* Hmass_wPhiV_PhoeUS;
TH1D* Hmass_wPhiV_PhoeLS;
TH1D* Hmass_wPhiV_PhoeSig;

//rebin
TH1D* Hmass_wPhiV_Rebin_US;
TH1D* Hmass_wPhiV_Rebin_LS;
TH1D* Hmass_wPhiV_Rebin_Sig;
//---------------------------------
//---------------------------------
TH1D* H_cent9;
TH1D* H_cent9_aftwt;
TH1D* H_cent16;
TH1D* H_cent16_aftwt;
//---------------------------------
//event plane information
const int nEvtPsiModes = 4;
TH1D* HEvtPsi[nEvtPsiModes];
const TString EvtPsiMode[nEvtPsiModes] = {"Raw", "Recent", "RejE", "fullCorr"};
const TString EvtPsiType[nEvtPsiModes] = {"Raw", "aft RecentCorr", "aft RejE", "aft ShiftCorr"};

//---------------------------------
//---------------------------------
//---------------------------------
//const TString centTitle[nCentBins]   = {"0-80%", "0-20%", "20-40%", "40-60%","60-80%"};
//const int centBinIndex[nCentBins][2] = { {0,8},   {6,8},   {4,5},    {2,3},    {0,1} };
//const TString centTitle[nCentBins]   = {"0-80%", "0-10%", "10-30%", "30-50%", "50-80%"};
//const int centBinIndex[nCentBins][2] = { {0,8},   {7,8},   {5,6},    {3,4},    {0,2} };
const int nCentBins = 4+4;
const TString centTitle[nCentBins]   = {"0-80%", "0-10%", "10-40%", "40-80%", "0-20%", "20-40%", "40-60%", "60-80%"};
const int centBinIndex[nCentBins][2] = { {0,8},   {7,8},   {4,6},    {0,3},    {6,8},   {4,5},    {2,3},    {0,1}};
//const int nPtBins = 7;
//const double ptBinBounds[nPtBins][2] = { {0.0,5.0}, {0.0,0.5}, {0.5,1.0}, {1.0,1.5}, {1.5,2.0}, {2.0,3.0}, {3.0,5.0} };
//const Int_t nptBins = 6;
//Double_t ptBinBDs[nptBins+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0};


const int nPtBins = 7;
const double ptBinBounds[nPtBins][2] = { {0.0,5.0}, {0.0,0.5}, {0.5,1.0}, {1.0,1.5}, {1.5,2.0}, {2.0,3.0}, {3.0,5.0} };

TH1D* H1d_M_US[nCentBins];
TH1D* H1d_M_LSPos[nCentBins];
TH1D* H1d_M_LSNeg[nCentBins];
TH1D* H1d_M_MixUS[nCentBins];
TH1D* H1d_M_MixLSPos[nCentBins];
TH1D* H1d_M_MixLSNeg[nCentBins];
TH1D* H1d_M_Reb_MixUS[nCentBins];
TH1D* H1d_M_Reb_MixLSPos[nCentBins];
TH1D* H1d_M_Reb_MixLSNeg[nCentBins];
TH1D* H1d_M_Reb_MixGmLS[nCentBins];

TH1D* H1d_M_SignAccFactor[nCentBins];

TH1D* H1d_M_perPt_Reb_MixUS[nCentBins][nPtBins];
TH1D* H1d_M_perPt_Reb_MixLSPos[nCentBins][nPtBins];
TH1D* H1d_M_perPt_Reb_MixLSNeg[nCentBins][nPtBins];
TH1D* H1d_M_perPt_Reb_MixGmLS[nCentBins][nPtBins];
TH1D* H1d_M_perPt_SignAccFactor[nCentBins][nPtBins];

const int nCent9 = 9;
TFile *infile[nCent9];

TH2D* H2d_MvsPt_US0[nCent9];
TH2D* H2d_MvsPt_LSPos0[nCent9];
TH2D* H2d_MvsPt_LSNeg0[nCent9];
TH2D* H2d_MvsMt_US0[nCent9];
TH2D* H2d_MvsMt_LSPos0[nCent9];
TH2D* H2d_MvsMt_LSNeg0[nCent9];
//after merge to interested centralities in analysis
TH2D* H2d_MvsPt_US[nCentBins];
TH2D* H2d_MvsPt_LSPos[nCentBins];
TH2D* H2d_MvsPt_LSNeg[nCentBins];
TH2D* H2d_MvsMt_US[nCentBins];
TH2D* H2d_MvsMt_LSPos[nCentBins];
TH2D* H2d_MvsMt_LSNeg[nCentBins];
//------------------------------------
//------------------------------------
TH2D* H2d_MvsPt_MixUS0[nCent9];
TH2D* H2d_MvsPt_MixLSPos0[nCent9];
TH2D* H2d_MvsPt_MixLSNeg0[nCent9];
TH2D* H2d_MvsMt_MixUS0[nCent9];
TH2D* H2d_MvsMt_MixLSPos0[nCent9];
TH2D* H2d_MvsMt_MixLSNeg0[nCent9];
//after merge to interested centralities in analysis
TH2D* H2d_MvsPt_MixUS[nCentBins];
TH2D* H2d_MvsPt_MixLSPos[nCentBins];
TH2D* H2d_MvsPt_MixLSNeg[nCentBins];
TH2D* H2d_MvsMt_MixUS[nCentBins];
TH2D* H2d_MvsMt_MixLSPos[nCentBins];
TH2D* H2d_MvsMt_MixLSNeg[nCentBins];
//------------------------------------
//------------------------------------
TH2D* H2d_MvsPt_Reb_US[nCentBins];
TH2D* H2d_MvsPt_Reb_LSPos[nCentBins];
TH2D* H2d_MvsPt_Reb_LSNeg[nCentBins];
TH2D* H2d_MvsPt_Reb_MixUS[nCentBins];
TH2D* H2d_MvsPt_Reb_MixLSPos[nCentBins];
TH2D* H2d_MvsPt_Reb_MixLSNeg[nCentBins];

TH2D* H2d_MvsMt_Reb_US[nCentBins];
TH2D* H2d_MvsMt_Reb_LSPos[nCentBins];
TH2D* H2d_MvsMt_Reb_LSNeg[nCentBins];
TH2D* H2d_MvsMt_Reb_MixUS[nCentBins];
TH2D* H2d_MvsMt_Reb_MixLSPos[nCentBins];
TH2D* H2d_MvsMt_Reb_MixLSNeg[nCentBins];

TH2D* H2d_MvsPt_Reb_GmLS[nCentBins];
TH2D* H2d_MvsPt_Reb_GeomMixLS[nCentBins];
TH2D* H2d_MvsPt_Reb_PSAC[nCentBins];
TH2D* H2d_MvsPt_Reb_GmLS_aftPSAC[nCentBins];

TH2D* H2d_MvsMt_Reb_GmLS[nCentBins];
TH2D* H2d_MvsMt_Reb_GeomMixLS[nCentBins];
TH2D* H2d_MvsMt_Reb_PSAC[nCentBins];
TH2D* H2d_MvsMt_Reb_GmLS_aftPSAC[nCentBins];

TH2D* H2d_MvsPt_LSPos_aftPSAC[nCentBins];   //LSPos after PSAC
TH2D* H2d_MvsPt_LSNeg_aftPSAC[nCentBins];   //LSNeg after PSAC

TH2D* H2d_MvsPt_Reb_Sig_Raw[nCentBins];
TH2D* H2d_MvsPt_Reb_Sig_aftPSAC[nCentBins];
TH2D* H2d_MvsPt_Reb_Sig_aftEff[nCentBins];
TH2D* H2d_MvsPt_Reb_Sig_aftAcc[nCentBins];

TH2D* H2d_MvsMt_Reb_Sig_Raw[nCentBins];
TH2D* H2d_MvsMt_Reb_Sig_aftPSAC[nCentBins];
TH2D* H2d_MvsMt_Reb_Sig_aftEff[nCentBins];
TH2D* H2d_MvsMt_Reb_Sig_aftAcc[nCentBins];

TH2D* H2d_MvsPt_Reb_GmLS4Ratio[nCentBins];
TH2D* H2d_MvsPt_Reb_GmLS_PSAC4Ratio[nCentBins];

//------------------------------------
//for Ratio, use cocktail bins
TH2D* H2d_MvsPt_Reb_US4Ratio[nCentBins];
TH2D* H2d_MvsPt_Reb_LSPos4Ratio[nCentBins];
TH2D* H2d_MvsPt_Reb_LSNeg4Ratio[nCentBins];
//------------------------------------
TH1D* H1d_M_Reb_US[nCentBins];
TH1D* H1d_M_Reb_LSPos[nCentBins];
TH1D* H1d_M_Reb_LSNeg[nCentBins];
TH1D* H1d_M_Reb_LS[nCentBins];
TH1D* H1d_M_Reb_GmLS[nCentBins];
TH1D* H1d_M_Reb_Sig_Raw[nCentBins];

TH1D* H1d_M_Reb_LSPos_aftPSAC[nCentBins];
TH1D* H1d_M_Reb_LSNeg_aftPSAC[nCentBins];
TH1D* H1d_M_Reb_GmLS_aftPSAC[nCentBins];
TH1D* H1d_M_Reb_Sig_aftPSAC[nCentBins];
//------------------------------------
//Signal/Background Ratio
TH1D* H1d_SoBvsM_Reb[nCentBins];

TH2D* H2d_MvsPt_LSPos_PSAC[nCentBins];  //PSAC for LSPos
TH2D* H2d_MvsPt_LSNeg_PSAC[nCentBins];  //PSAC for LSNeg


TH1D* H1d_M_Reb_US4Ratio[nCentBins];
TH1D* H1d_M_Reb_GmLS4Ratio[nCentBins];

//after the pair sign acc. correction
TH1D* H1d_M_Reb_GmLS_PSAC4Ratio[nCentBins];
TH1D* H1d_M_Reb_Sig_PSAC4Ratio[nCentBins];
//------------------------------------
//------------------------------------

//read in pure Omega and Phi Signal for Removal
TFile *infile_OmegaPhi;
TH1D* H1d_M_Omega_Pure;
TH1D* H1d_M_Phi_Pure;
TH1D* H1d_M_Omega_Reb_Pure;
TH1D* H1d_M_Phi_Reb_Pure;

TH1D* H1d_M_Reb_Sig_NoOmegaAndPhi;
TH1D* H1d_M_Reb_Sig_NoOmegaAndPhi_aftEff;
TH1D* HRatio_Data2CKT;
TH1D* HRatio_Data2CKT_SysErr;

//for Pair Eff and Acc, default as after Reb to same MvsPt, or MvsMt binning
TH2D* H2d_MvsPt_PairEff[nCentBins];
TH2D* H2d_MvsPt_PairAcc[nCentBins];
TH2D* H2d_MvsMt_PairEff[nCentBins];
TH2D* H2d_MvsMt_PairAcc[nCentBins];

TH1D* H1d_M_Reb_Sig_aftEff[nCentBins];
TH1D* H1d_M_Reb_Sig_SysErr;
TH1D* H1d_Mt_Reb_AftEff_IMR[nCentBins];
TH1D* H1d_Mt_totCKT_Reb_IMR[nCentBins];

const Int_t nCKTs = 6;
const Int_t nPtBins4Phys            = 6; //only for 0-80% centrality
Double_t PtBDs4Phys[nPtBins4Phys+1] = {0., 0.15, 0.50, 1.0, 1.5, 2.5, 5.0};

TH1D* H1d_Mt_CKT_Reb_IMR[nCentBins][nCKTs+2];

TH2D* H2d_MvsPt_Reb4PtPhys_Sig_aftEff; //only 0-80%
TH1D* H1d_DataMinusCKT[nCentBins];
TH1D* H1d_M_Reb4PtPhys_Sig_aftEff[nPtBins4Phys];
//for plot
TH1D* H_M_CKT4PtPhys[nCKTs+2][nPtBins4Phys];//0-80%, different pt bin
TH1D* H_M_totCKT4PtPhys[nPtBins4Phys];
//for ratio
TH2D* H2D_MvsPt_totCKT_Reb4PtPhys;
TH1D* H_M_totCKT_Reb4PtPhys[nPtBins4Phys];
TH1D* H1D_Ratio2CKT_4PtPhys[nPtBins4Phys];

//get Mt-M spectral and fit for the 
TH2D* HExY_MvsMt[nCentBins];         // (Data after PairEff - CKT)
TH2D* HExY_MvsMt_AftAcc[nCentBins];  // (Data after PairEff - CKT) after Pair Acc. Corr.


TH1D* HExY_Mt_IMR[nCentBins];            //Direct Eff. Corred Data-CKT 
TH1D* HExY_Mt_AftAcc_IMR[nCentBins];     //(Data-CKT)/Acc, final mT-M
TH1D* HExY_MtSpec_AftAcc_IMR[nCentBins]; //final mT-M, (1/mT)*(dN/dmT) good to Fitting for Temperature 
TH1D* HExY_MtSpec_AftAcc_Reb_IMR[nCentBins]; //final mT-M, (1/mT)*(dN/dmT) good to Fitting for Temperature 

TH1D* HExY_M_AftAcc_AllMass[nCentBins];



TH1D* HRatioVsPt_IMR;
TH1D* HExYVsPt_IMR;

TFile *infile_pairEff;
TH1D* H1D_PairEff_Reb;
TFile * inf_eff[nCentBins];
const int nCuts = 5;
TH2D* HRc_MvsPt[nCentBins][nCuts];
TH2D* HRc_MvsPt_Reb[nCentBins][nCuts];

TH2D* HRc_MvsMt[nCentBins][nCuts];
TH2D* HRc_MvsMt_Reb[nCentBins][nCuts];

TFile *infile_pairEff4Ratio;
TH1D* H1D_PairEff4Ratio_Reb;
TH1D* H1d_M_Reb_Sig_aftEff4Ratio[nCentBins];

//read Zaochen's Cocktail
const TString cktDecayName[nCKTs+2] = 
{
	"#pi^{0} #rightarrow #gammaee", 
	"#eta #rightarrow #gammaee", 
	"#omega #rightarrowee & #omega #rightarrow #pi^{0}ee",
	"#eta' #rightarrow #gammaee",
	"#phi #rightarrow ee & #phi #rightarrow #etaee",
	"J/#psi #rightarrow ee",
	"c#bar{c} #rightarrow ee",
	"DY #rightarrow ee"
};
const int cktColors[nCKTs+2] = {1,6,8,12,15,20,95, 98};

TH2D* H2D_MvsPt_CKT[nCentBins][nCKTs+2];
TH2D* H2D_MvsPt_totCKT[nCentBins];
TH2D* H2D_MvsPt_CKT_Reb[nCentBins][nCKTs+2];
TH2D* H2D_MvsPt_totCKT_Reb[nCentBins];

TH2D* H2D_MvsMt_CKT[nCentBins][nCKTs+2];
TH2D* H2D_MvsMt_totCKT[nCentBins];
TH2D* H2D_MvsMt_CKT_Reb[nCentBins][nCKTs+2];
TH2D* H2D_MvsMt_totCKT_Reb[nCentBins];

TH1D* H1D_M_CKT[nCentBins][nCKTs+2];
TH1D* H1D_M_totCKT[nCentBins];
TH1D* H1D_M_CKT_Reb[nCentBins][nCKTs+2];
TH1D* H1D_M_totCKT_Reb[nCentBins];

TH1D* H1D_Mt_CKT[nCentBins][nCKTs+2];
TH1D* H1D_Mt_totCKT[nCentBins];
TH1D* H1D_Mt_CKT_Reb[nCentBins][nCKTs+2];
TH1D* H1D_Mt_totCKT_Reb[nCentBins];

//Ratio of Data/totCKT 
TH1D* H1D_Ratio2CKT[nCentBins];

//read in Run11 data
TFile *infile_27Run11;
TH1D* hSpecM_27Run11 ;
TH1D* hSpecM_SysErr  ;
TH1D* gCKT           ;
TH1D* gCKT4Ratio     ;
TH1D* gCKTsys        ;
TH1D* gRCKSys        ;
TH1D* gpion          ;
TH1D* geta           ;
TH1D* getap          ;
TH1D* gomega         ;
TH1D* gphi           ;
TH1D* gccbar1        ;
TH1D* gjpsi          ;
TH1D* gCKT_noOmegaPhi;
TH1D* hRatioRun11;
TH1D* hRatioRun11Sys;

TFile* infile_27Rapp;
TGraph* gTheory_Rapp_CGA;


//0-80, npart = 123.35
const int    nNPartBins = 3;
const double NpartValues[nNPartBins] = {40.2, 166.1, 315.7}; //40-80%, 10-40%, 0-10%

TGraphErrors* ge_ExYVsNpart_Rho;    // Data-CKT
TGraphErrors* ge_ExYVsNpart_IMR;    // Data-CKT
TGraphErrors* ge_YVsNpart_Omega;    // Data
TGraphErrors* ge_YVsNpart_Phi;      // Data
const double scale_ExY_IMR = 1.e2;

