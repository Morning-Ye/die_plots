#include "/Users/zaochenrice/myFunction.h"
#include "head.h"
#include "massBins.h"

const double goodPtLow = 0.00;
const double goodPtHig = 5.00;

const double scaleM_Low[nCentBins] = {1.5, 1.5, 1.5, 1.5};
const double scaleM_Hig[nCentBins] = {2.5, 2.5, 2.5, 2.5};
const double RepM_Low[nCentBins]   = {1.5, 1.5, 1.5, 1.5};
const double RepM_Hig[nCentBins]   = {2.5, 2.5, 2.5, 2.5};

double fun_Rcut(double* x, double* par);
TF1* fRcut;
double fun_MtSpec(  double* x, double* par);
double fun_AccedExY(double* x, double* par);
//---------------------------------------------------------------------------------------------------------------------------
double bg_formula(double* x, double* par);
double BreitWigner_formula(double* x, double* par);
double BreitWigner_BG_formula(double* x, double* par);
TH1D* CombGmLSMixUS( TString outname, TH1D*h1, TH1D*h2, double lowMass, double higMass ); //replace h1 with h2 for mass region [lowMass, higMass]
TH1D* hist1_4fit;
TH1D* hist2_4fit;
TH1D* hist3_4fit;
TGraph* gh_4fit;
double HistGraphAsFitf(double *x, double *par);
//---------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------
//------------------------------------------------------------------------
void prepareRawSig();
void preparePSCA(         int viewFlag);
void preparePairEffAndAcc(int viewFlag);
void corrAndGetSignals();

void readTheory();

void drawEvtPsi();
void drawSigInCent();
void drawSignAccFactor();
void getInitCKTs(); //these are from the traditional scales
void prepareCKT4fit2Data();
void drawPhysicsInPt();
void drawPhysicsInCent();
void drawExYvsTheory();
void fitAndRmVectMeson4ExYield();
void getAllExcessYields();//must run after fitAndRmVectMeson4ExYield()
void drawExYieldsVsNpart();
void drawVectorMesons();
double getDnchDy(int icent);
void saveFiles();

Int_t GetnEvtsInCent(Int_t centBinIdx);
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
void plotDieSignals()
{
	prepareRawSig();
	preparePairEffAndAcc(0); //0:noDraw, 1:draw PairEff, 2: draw PairAcc
	preparePSCA(0);          //0:noDraw, 1:draw
	getInitCKTs();           //get the traditionally scaled CKTs
	prepareCKT4fit2Data();
	corrAndGetSignals();
	readTheory();

	//drawVectorMesons();
	drawExYvsTheory();
	fitAndRmVectMeson4ExYield();

	getAllExcessYields();
	drawExYieldsVsNpart();

	//	//	//////drawEvtPsi();
	//	drawSigInCent();
	//	drawPhysicsInCent();
	//	drawPhysicsInPt();
	//	//	drawPhysics();
	//	//checkBump();
	//	saveFiles();
}
//------------------------------------------------------------------------
//------------------------------------------------------------------------

//------------------------------------------------------------------------
//------------------------------------------------------------------------
void prepareRawSig()
{
	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
	for(int icent9=0; icent9<nCent9; icent9++)
	{
		//------------------------------------------------------------------------
		//------------------------------------------------------------------------
		//infile[icent9] = new TFile(Form("./../inputfiles/splitCent/previousPhiV/dataWMt/wMixMt/outAnaHist_Run18_27GeV_cent%d.root",icent9));
		infile[icent9] = new TFile(Form("./inputfiles/data/outAnaHist_Run18_27GeV_Vz35cm_cent%d.root",icent9));
		cout<<"readin: "<<infile[icent9]->GetName()<<endl;
		//------------------------------------------------------------------------
		//------------------------------------------------------------------------
		if(icent9==0)
		{
			H_cent9        = (TH1D*)infile[icent9]->Get("hcent9");
			H_cent16       = (TH1D*)infile[icent9]->Get("hcent16");
			H_cent9_aftwt  = (TH1D*)infile[icent9]->Get("hcent9_aftwt");
			H_cent16_aftwt = (TH1D*)infile[icent9]->Get("hcent16_aftwt");

			cout<<"------------------------------------------- "<<endl;
			//read in event plane information
			for(int imd=0; imd<nEvtPsiModes; imd++)
			{
				TString temEvtPsiHist = "hEvtPsi_"+EvtPsiMode[imd];
				cout<<"readin: "<<temEvtPsiHist<<endl;
				HEvtPsi[imd] = (TH1D*)infile[icent9]->Get("hEvtPsi_"+EvtPsiMode[imd]);
			}
			cout<<"------------------------------------------- "<<endl;
		}
		else 
		{
			TH1D* htem_cent9        = (TH1D*)infile[icent9]->Get("hcent9");
			TH1D* htem_cent16       = (TH1D*)infile[icent9]->Get("hcent16");
			TH1D* htem_cent9_aftwt  = (TH1D*)infile[icent9]->Get("hcent9_aftwt");
			TH1D* htem_cent16_aftwt = (TH1D*)infile[icent9]->Get("hcent16_aftwt");

			H_cent9        -> Add(htem_cent9);
			H_cent16       -> Add(htem_cent16);
			H_cent9_aftwt  -> Add(htem_cent9_aftwt);
			H_cent16_aftwt -> Add(htem_cent16_aftwt);

			cout<<"------------------------------------------- "<<endl;
			//read in event plane information
			for(int imd=0; imd<nEvtPsiModes; imd++)
			{
				TString temEvtPsiHist = "hEvtPsi_"+EvtPsiMode[imd];
				cout<<"readin: "<<temEvtPsiHist<<endl;

				TH1D* htem_HEvtPsi = (TH1D*)infile[icent9]->Get("hEvtPsi_"+EvtPsiMode[imd]);

				HEvtPsi[imd] ->Add(htem_HEvtPsi);

				delete htem_HEvtPsi;
			}
			cout<<"------------------------------------------- "<<endl;

			delete htem_cent9;
			delete htem_cent16;
			delete htem_cent9_aftwt;
			delete htem_cent16_aftwt;
		}
		//------------------------------------------------------------------------
		//read in the real signal histograms from each of the 9 centrality bins
		//------------------------------------------------------------------------
		H2d_MvsPt_US0[icent9]    = (TH2D*)infile[icent9]->Get("hMvsPt_wPhiV_US");
		H2d_MvsPt_LSPos0[icent9] = (TH2D*)infile[icent9]->Get("hMvsPt_wPhiV_LSPos");
		H2d_MvsPt_LSNeg0[icent9] = (TH2D*)infile[icent9]->Get("hMvsPt_wPhiV_LSNeg");
		//------------------------------------------------------------------------
		H2d_MvsMt_US0[icent9]    = (TH2D*)infile[icent9]->Get("hMvsMt_wPhiV_US");
		H2d_MvsMt_LSPos0[icent9] = (TH2D*)infile[icent9]->Get("hMvsMt_wPhiV_LSPos");
		H2d_MvsMt_LSNeg0[icent9] = (TH2D*)infile[icent9]->Get("hMvsMt_wPhiV_LSNeg");
	}//for Raw Signals from 9 centranlity bins
	//------------------------------------------------------------------------
	//------------------------------------------------------------------------

	//const int centBinIndex[nCentBins][2] = { {0,8},   {7,8},   {5,6},    {3,4},    {0,2} };
	//merge the histogram for the interested centrality bins
	for( int icent=0; icent<nCentBins; icent++ )
	{
		//for( int icentbin = centBinIndex[icent][0]; icentbin<centBinIndex[icent][1]; icentbin++ ) //for QM, WRONG!!!
		for( int icentbin = centBinIndex[icent][0]; icentbin<=centBinIndex[icent][1]; icentbin++ )
		{
			if( icentbin==centBinIndex[icent][0] )
			{
				H2d_MvsPt_US[icent]    = (TH2D*) H2d_MvsPt_US0[icentbin]   ->Clone();
				H2d_MvsPt_LSPos[icent] = (TH2D*) H2d_MvsPt_LSPos0[icentbin]->Clone();
				H2d_MvsPt_LSNeg[icent] = (TH2D*) H2d_MvsPt_LSNeg0[icentbin]->Clone();
				//------------------------------------------------------------------------
				H2d_MvsMt_US[icent]    = (TH2D*) H2d_MvsMt_US0[icentbin]   ->Clone();
				H2d_MvsMt_LSPos[icent] = (TH2D*) H2d_MvsMt_LSPos0[icentbin]->Clone();
				H2d_MvsMt_LSNeg[icent] = (TH2D*) H2d_MvsMt_LSNeg0[icentbin]->Clone();
			}
			else
			{
				H2d_MvsPt_US[icent]    -> Add(H2d_MvsPt_US0[icentbin]);
				H2d_MvsPt_LSPos[icent] -> Add(H2d_MvsPt_LSPos0[icentbin]);
				H2d_MvsPt_LSNeg[icent] -> Add(H2d_MvsPt_LSNeg0[icentbin]);
				//------------------------------------------------------------------------
				H2d_MvsMt_US[icent]    -> Add(H2d_MvsMt_US0[icentbin]);
				H2d_MvsMt_LSPos[icent] -> Add(H2d_MvsMt_LSPos0[icentbin]);
				H2d_MvsMt_LSNeg[icent] -> Add(H2d_MvsMt_LSNeg0[icentbin]);
			}
		}

		H2d_MvsPt_US[icent]   ->SetName( Form("H2d_MvsPt_US_icent%d",   icent) );
		H2d_MvsPt_LSPos[icent]->SetName( Form("H2d_MvsPt_LSPos_icent%d",icent) );
		H2d_MvsPt_LSNeg[icent]->SetName( Form("H2d_MvsPt_LSNeg_icent%d",icent) );
		//------------------------------------------------------------------------
		H2d_MvsMt_US[icent]   ->SetName( Form("H2d_MvsMt_US_icent%d",   icent) );
		H2d_MvsMt_LSPos[icent]->SetName( Form("H2d_MvsMt_LSPos_icent%d",icent) );
		H2d_MvsMt_LSNeg[icent]->SetName( Form("H2d_MvsMt_LSNeg_icent%d",icent) );

		//---------------------------------------------------------------------------
		//scale by nMB events
		double nEvt_icent = 1.*GetnEvtsInCent(icent);
		H2d_MvsPt_US[icent]     ->Scale(1./nEvt_icent);
		H2d_MvsPt_LSPos[icent]  ->Scale(1./nEvt_icent);
		H2d_MvsPt_LSNeg[icent]  ->Scale(1./nEvt_icent);
		H2d_MvsMt_US[icent]     ->Scale(1./nEvt_icent);
		H2d_MvsMt_LSPos[icent]  ->Scale(1./nEvt_icent);
		H2d_MvsMt_LSNeg[icent]  ->Scale(1./nEvt_icent);
	}//get Raw signals for ana cent bins

	//------------------------------------------------------------------------------------------------------------------------------------------
	//Rebin the 2d histograms and calculate the geom LS
	//------------------------------------------------------------------------------------------------------------------------------------------
	for( int icent=0; icent<nCentBins; icent++ )//icent here is the analysis centbins
	{
		//------------------------------------------------------------------------------------------------------------------------------------------
		cout<<"rebin histograms"<<endl;
		cout<<"------------------------------------------- "<<endl;

		if(icent==0)
		{
			H2d_MvsPt_Reb_US[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_US[icent],    Form("H2d_MvsPt_Reb_US_icent%d",    icent), nptBins,ptBds4Corr, nMBins_ct0,mBds4Corr_ct0,"XY","");
			H2d_MvsPt_Reb_LSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_LSPos[icent], Form("H2d_MvsPt_Reb_LSPos_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct0,mBds4Corr_ct0,"XY","");
			H2d_MvsPt_Reb_LSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_LSNeg[icent], Form("H2d_MvsPt_Reb_LSNeg_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct0,mBds4Corr_ct0,"XY","");
			//------------------------------------------------------------------------------------------------------------------------------------------
			H2d_MvsMt_Reb_US[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_US[icent],    Form("H2d_MvsMt_Reb_US_icent%d",    icent), nmtBins,mtBinBDs, nMBins_ct0,mBds4Corr_ct0,"XY","");
			H2d_MvsMt_Reb_LSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_LSPos[icent], Form("H2d_MvsMt_Reb_LSPos_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct0,mBds4Corr_ct0,"XY","");
			H2d_MvsMt_Reb_LSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_LSNeg[icent], Form("H2d_MvsMt_Reb_LSNeg_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct0,mBds4Corr_ct0,"XY","");
		}
		else if(icent==1)
		{
			H2d_MvsPt_Reb_US[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_US[icent],    Form("H2d_MvsPt_Reb_US_icent%d",    icent), nptBins,ptBds4Corr, nMBins_ct1,mBds4Corr_ct1,"XY","");
			H2d_MvsPt_Reb_LSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_LSPos[icent], Form("H2d_MvsPt_Reb_LSPos_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct1,mBds4Corr_ct1,"XY","");
			H2d_MvsPt_Reb_LSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_LSNeg[icent], Form("H2d_MvsPt_Reb_LSNeg_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct1,mBds4Corr_ct1,"XY","");
			//------------------------------------------------------------------------------------------------------------------------------------------
			H2d_MvsMt_Reb_US[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_US[icent],    Form("H2d_MvsMt_Reb_US_icent%d",    icent), nmtBins,mtBinBDs, nMBins_ct1,mBds4Corr_ct1,"XY","");
			H2d_MvsMt_Reb_LSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_LSPos[icent], Form("H2d_MvsMt_Reb_LSPos_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct1,mBds4Corr_ct1,"XY","");
			H2d_MvsMt_Reb_LSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_LSNeg[icent], Form("H2d_MvsMt_Reb_LSNeg_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct1,mBds4Corr_ct1,"XY","");
		}
		else if(icent==2)
		{
			H2d_MvsPt_Reb_US[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_US[icent],    Form("H2d_MvsPt_Reb_US_icent%d",    icent), nptBins,ptBds4Corr, nMBins_ct2,mBds4Corr_ct2,"XY","");
			H2d_MvsPt_Reb_LSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_LSPos[icent], Form("H2d_MvsPt_Reb_LSPos_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct2,mBds4Corr_ct2,"XY","");
			H2d_MvsPt_Reb_LSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_LSNeg[icent], Form("H2d_MvsPt_Reb_LSNeg_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct2,mBds4Corr_ct2,"XY","");
			//------------------------------------------------------------------------------------------------------------------------------------------
			H2d_MvsMt_Reb_US[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_US[icent],    Form("H2d_MvsMt_Reb_US_icent%d",    icent), nmtBins,mtBinBDs, nMBins_ct2,mBds4Corr_ct2,"XY","");
			H2d_MvsMt_Reb_LSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_LSPos[icent], Form("H2d_MvsMt_Reb_LSPos_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct2,mBds4Corr_ct2,"XY","");
			H2d_MvsMt_Reb_LSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_LSNeg[icent], Form("H2d_MvsMt_Reb_LSNeg_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct2,mBds4Corr_ct2,"XY","");
		}
		else if(icent==3)
		{
			H2d_MvsPt_Reb_US[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_US[icent],    Form("H2d_MvsPt_Reb_US_icent%d",    icent), nptBins,ptBds4Corr, nMBins_ct3,mBds4Corr_ct3,"XY","");
			H2d_MvsPt_Reb_LSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_LSPos[icent], Form("H2d_MvsPt_Reb_LSPos_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct3,mBds4Corr_ct3,"XY","");
			H2d_MvsPt_Reb_LSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_LSNeg[icent], Form("H2d_MvsPt_Reb_LSNeg_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct3,mBds4Corr_ct3,"XY","");
			//------------------------------------------------------------------------------------------------------------------------------------------
			H2d_MvsMt_Reb_US[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_US[icent],    Form("H2d_MvsMt_Reb_US_icent%d",    icent), nmtBins,mtBinBDs, nMBins_ct3,mBds4Corr_ct3,"XY","");
			H2d_MvsMt_Reb_LSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_LSPos[icent], Form("H2d_MvsMt_Reb_LSPos_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct3,mBds4Corr_ct3,"XY","");
			H2d_MvsMt_Reb_LSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_LSNeg[icent], Form("H2d_MvsMt_Reb_LSNeg_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct3,mBds4Corr_ct3,"XY","");
		}

		//---------------------------------------------------------------------------------------------------------------------------------------------------------
	}//icent
	//---------------------------------------------------------------------------------------------------------------------------------------
}
//---------------------------------------------------------------------------
void preparePSCA(int viewFlag)
{
	//read in the mixed histograms for the pair sign acc. diff. corr
	TFile* infile_PSAC = new TFile("./out4PSAC_2DCorr/outhist_4PSAC_2DCorr_27GeV.root", "read");
	cout<<"readin: "<<infile_PSAC->GetName()<<endl;
	for(int icent=0; icent<nCentBins; icent++)
	{
		H2d_MvsPt_MixUS[icent]    = (TH2D*) infile_PSAC->Get( Form("H2d_MvsPt_MixUS_icent%d",    icent) );
		H2d_MvsPt_MixLSPos[icent] = (TH2D*) infile_PSAC->Get( Form("H2d_MvsPt_MixLSPos_icent%d", icent) );
		H2d_MvsPt_MixLSNeg[icent] = (TH2D*) infile_PSAC->Get( Form("H2d_MvsPt_MixLSNeg_icent%d", icent) );

		H2d_MvsMt_MixUS[icent]    = (TH2D*) infile_PSAC->Get( Form("H2d_MvsMt_MixUS_icent%d",    icent) );
		H2d_MvsMt_MixLSPos[icent] = (TH2D*) infile_PSAC->Get( Form("H2d_MvsMt_MixLSPos_icent%d", icent) );
		H2d_MvsMt_MixLSNeg[icent] = (TH2D*) infile_PSAC->Get( Form("H2d_MvsMt_MixLSNeg_icent%d", icent) );

		if(icent==0)
		{
			H2d_MvsPt_Reb_MixUS[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_MixUS[icent],    Form("H2d_MvsPt_Reb_MixUS_icent%d",    icent), nptBins,ptBds4Corr, nMBins_ct0,mBds4Corr_ct0,"XY","");
			H2d_MvsPt_Reb_MixLSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_MixLSPos[icent], Form("H2d_MvsPt_Reb_MixLSPos_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct0,mBds4Corr_ct0,"XY","");
			H2d_MvsPt_Reb_MixLSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_MixLSNeg[icent], Form("H2d_MvsPt_Reb_MixLSNeg_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct0,mBds4Corr_ct0,"XY","");
			//---------------------------------------------------------------------------------------------------------------------------------------
			H2d_MvsMt_Reb_MixUS[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_MixUS[icent],    Form("H2d_MvsMt_Reb_MixUS_icent%d",    icent), nmtBins,mtBinBDs, nMBins_ct0,mBds4Corr_ct0,"XY","");
			H2d_MvsMt_Reb_MixLSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_MixLSPos[icent], Form("H2d_MvsMt_Reb_MixLSPos_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct0,mBds4Corr_ct0,"XY","");
			H2d_MvsMt_Reb_MixLSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_MixLSNeg[icent], Form("H2d_MvsMt_Reb_MixLSNeg_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct0,mBds4Corr_ct0,"XY","");
		}
		else if(icent==1)
		{
			H2d_MvsPt_Reb_MixUS[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_MixUS[icent],    Form("H2d_MvsPt_Reb_MixUS_icent%d",    icent), nptBins,ptBds4Corr, nMBins_ct1,mBds4Corr_ct1,"XY","");
			H2d_MvsPt_Reb_MixLSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_MixLSPos[icent], Form("H2d_MvsPt_Reb_MixLSPos_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct1,mBds4Corr_ct1,"XY","");
			H2d_MvsPt_Reb_MixLSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_MixLSNeg[icent], Form("H2d_MvsPt_Reb_MixLSNeg_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct1,mBds4Corr_ct1,"XY","");
			//---------------------------------------------------------------------------------------------------------------------------------------
			H2d_MvsMt_Reb_MixUS[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_MixUS[icent],    Form("H2d_MvsMt_Reb_MixUS_icent%d",    icent), nmtBins,mtBinBDs, nMBins_ct1,mBds4Corr_ct1,"XY","");
			H2d_MvsMt_Reb_MixLSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_MixLSPos[icent], Form("H2d_MvsMt_Reb_MixLSPos_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct1,mBds4Corr_ct1,"XY","");
			H2d_MvsMt_Reb_MixLSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_MixLSNeg[icent], Form("H2d_MvsMt_Reb_MixLSNeg_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct1,mBds4Corr_ct1,"XY","");
		}
		else if(icent==2)
		{
			H2d_MvsPt_Reb_MixUS[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_MixUS[icent],    Form("H2d_MvsPt_Reb_MixUS_icent%d",    icent), nptBins,ptBds4Corr, nMBins_ct2,mBds4Corr_ct2,"XY","");
			H2d_MvsPt_Reb_MixLSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_MixLSPos[icent], Form("H2d_MvsPt_Reb_MixLSPos_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct2,mBds4Corr_ct2,"XY","");
			H2d_MvsPt_Reb_MixLSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_MixLSNeg[icent], Form("H2d_MvsPt_Reb_MixLSNeg_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct2,mBds4Corr_ct2,"XY","");
			//---------------------------------------------------------------------------------------------------------------------------------------
			H2d_MvsMt_Reb_MixUS[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_MixUS[icent],    Form("H2d_MvsMt_Reb_MixUS_icent%d",    icent), nmtBins,mtBinBDs, nMBins_ct2,mBds4Corr_ct2,"XY","");
			H2d_MvsMt_Reb_MixLSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_MixLSPos[icent], Form("H2d_MvsMt_Reb_MixLSPos_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct2,mBds4Corr_ct2,"XY","");
			H2d_MvsMt_Reb_MixLSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_MixLSNeg[icent], Form("H2d_MvsMt_Reb_MixLSNeg_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct2,mBds4Corr_ct2,"XY","");
		}
		else if(icent==3)
		{
			H2d_MvsPt_Reb_MixUS[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_MixUS[icent],    Form("H2d_MvsPt_Reb_MixUS_icent%d",    icent), nptBins,ptBds4Corr, nMBins_ct3,mBds4Corr_ct3,"XY","");
			H2d_MvsPt_Reb_MixLSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_MixLSPos[icent], Form("H2d_MvsPt_Reb_MixLSPos_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct3,mBds4Corr_ct3,"XY","");
			H2d_MvsPt_Reb_MixLSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsPt_MixLSNeg[icent], Form("H2d_MvsPt_Reb_MixLSNeg_icent%d", icent), nptBins,ptBds4Corr, nMBins_ct3,mBds4Corr_ct3,"XY","");
			//---------------------------------------------------------------------------------------------------------------------------------------
			H2d_MvsMt_Reb_MixUS[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_MixUS[icent],    Form("H2d_MvsMt_Reb_MixUS_icent%d",    icent), nmtBins,mtBinBDs, nMBins_ct3,mBds4Corr_ct3,"XY","");
			H2d_MvsMt_Reb_MixLSPos[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_MixLSPos[icent], Form("H2d_MvsMt_Reb_MixLSPos_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct3,mBds4Corr_ct3,"XY","");
			H2d_MvsMt_Reb_MixLSNeg[icent]
				= (TH2D*) rebHist2d(H2d_MvsMt_MixLSNeg[icent], Form("H2d_MvsMt_Reb_MixLSNeg_icent%d", icent), nmtBins,mtBinBDs, nMBins_ct3,mBds4Corr_ct3,"XY","");
		}

		//---------------------------------------------------------------------------------------------------------------------------------------
		//get PSAC_LSPos and PSAC_LSNeg separately
		H2d_MvsPt_LSPos_PSAC[icent] = (TH2D*) calRatioH1toH2(H2d_MvsPt_Reb_MixUS[icent], H2d_MvsPt_Reb_MixLSPos[icent], Form("H2d_MvsPt_LSPos_PSAC_icent%d",icent), "NOCORR");
		H2d_MvsPt_LSNeg_PSAC[icent] = (TH2D*) calRatioH1toH2(H2d_MvsPt_Reb_MixUS[icent], H2d_MvsPt_Reb_MixLSNeg[icent], Form("H2d_MvsPt_LSNeg_PSAC_icent%d",icent), "NOCORR");
		H2d_MvsPt_LSPos_PSAC[icent] ->Scale(1./2.);
		H2d_MvsPt_LSNeg_PSAC[icent] ->Scale(1./2.);
		//---------------------------------------------------------------------------------------------------------------------------------------

		//---------------------------------------------------------------------------------------------------------------------------------------
		//get mixed geom LS, if use the same PSAC for LSPos and LSNeg
		H2d_MvsPt_Reb_GeomMixLS[icent] = (TH2D*) getGeomMeanHist2d( H2d_MvsPt_Reb_MixLSPos[icent], H2d_MvsPt_Reb_MixLSNeg[icent], Form("H2d_MvsPt_Reb_GeomMixLS_icent%d", icent) );
		H2d_MvsMt_Reb_GeomMixLS[icent] = (TH2D*) getGeomMeanHist2d( H2d_MvsMt_Reb_MixLSPos[icent], H2d_MvsMt_Reb_MixLSNeg[icent], Form("H2d_MvsMt_Reb_GeomMixLS_icent%d", icent) );

		//get 2D PSAC
		H2d_MvsPt_Reb_PSAC[icent]      = (TH2D*) calRatioH1toH2(H2d_MvsPt_Reb_MixUS[icent], H2d_MvsPt_Reb_GeomMixLS[icent], Form("H2d_MvsPt_Reb_PSAC_icent%d",icent), "NOCORR");
		H2d_MvsMt_Reb_PSAC[icent]      = (TH2D*) calRatioH1toH2(H2d_MvsMt_Reb_MixUS[icent], H2d_MvsMt_Reb_GeomMixLS[icent], Form("H2d_MvsMt_Reb_PSAC_icent%d",icent), "NOCORR");
		//---------------------------------------------------------------------------------------------------------------------------------------
	}//icent

	//commented before
	if(viewFlag==1)
	{
		//---------------------------------------------------------------------------
		//make plots of the pair sign acc. diff. corr. factors in 2D
		//---------------------------------------------------------------------------
		TCanvas *c1 = new TCanvas("c1");
		gStyle->SetOptStat(0);
		system("mkdir -p outAnaPlots/plots_PSAC");
		
		TPDF* mypdf = new TPDF("./outAnaPlots/plots_PSAC/Run18_27GeV_AuAu_PSAC_factors_Summary.pdf",111);
		mypdf->Off();
		updatePDF(c1,mypdf);

		Int_t nColumns = 3;
		Int_t nRows    = 3;
		Int_t nPads    = nColumns*nRows;
				
		TH1D* Htempt1 = (TH1D*)H2d_MvsPt_Reb_PSAC[0]->ProjectionX("Htempt1");
		TH1D* H1D_PSAC[nCentBins][nptBins];
		
		for(int icent=0; icent<nCentBins; icent++)
		{
			c1->cd();
			H2d_MvsPt_Reb_PSAC[icent] ->SetTitle("Pair Sign Acc. Diff. Corr. Factor in Cent: "+centTitle[icent]);
			H2d_MvsPt_Reb_PSAC[icent] ->Draw("LEGO2Z");
			//c1->SaveAs(Form("./outAnaPlots/plots_PSAC/H2d_MvsPt_Reb_PSAC_icent%d.png", icent));
			//c1->SaveAs(Form("./outAnaPlots/plots_PSAC/H2d_MvsPt_Reb_PSAC_icent%d.pdf", icent));

			updatePDF(c1,mypdf);
			c1->Clear();

			c1->Divide(nColumns,nRows);
			for(Int_t ipad=0; ipad<nColumns*nRows; ipad++)
			{
				c1->cd(ipad+1);
				setPad(0.1,0.02,0.08,0.1);
			}

			for(int ipt=0; ipt<Htempt1->GetNbinsX(); ipt++)
			{
				Int_t idx = ipt;
				c1->cd(idx%nPads+1);

				cout<<Form("%.2f<pt<%.2f", ptBds4Corr[ipt], ptBds4Corr[ipt+1])<<endl;
				TString temptName = Form("%.2f<pt<%.2f", ptBds4Corr[ipt], ptBds4Corr[ipt+1]);

				H1D_PSAC[icent][ipt] = (TH1D*) H2d_MvsPt_Reb_PSAC[icent]   ->ProjectionY(Form("H1D_PSAC_icent%d_ipt%d", icent, ipt), ipt+1, ipt+1);

				H1D_PSAC[icent][ipt]->SetTitle("Pair S.A.C Factor in Cent: "+centTitle[icent]+" "+temptName);
				H1D_PSAC[icent][ipt]->SetYTitle("PSAC Factor");
				H1D_PSAC[icent][ipt]->SetMarkerStyle(24); 
				H1D_PSAC[icent][ipt]->SetMarkerSize(0.3); 
				H1D_PSAC[icent][ipt]->SetMarkerColor(1); 
				H1D_PSAC[icent][ipt]->SetLineColor(1); 
				H1D_PSAC[icent][ipt]->SetAxisRange(0.90, 1.10, "y");
				H1D_PSAC[icent][ipt]->Draw("pe");

				//c1->SaveAs(Form("./outAnaPlots/plots_PSAC/H1D_PSAC_icent%d_ipt%d.png", icent, ipt));
				if( idx%nPads == nPads-1 ) updatePDF(c1,mypdf);
			}
			
			if((Htempt1->GetNbinsX())%nPads!=0) updatePDF(c1,mypdf);

			c1->Clear();
		}//icent

		//---------------------------------------------------------------------------------------------------------------------------------------
		//draw for Mt
		//---------------------------------------------------------------------------------------------------------------------------------------
		TH1D* Htemmt1 = (TH1D*)H2d_MvsMt_Reb_PSAC[0]->ProjectionX("Htemmt1");
		TH1D* H1D_PSAC_4Mt[nCentBins][nmtBins];
		for(int icent=0; icent<nCentBins; icent++)
		{
			c1->cd();
			H2d_MvsMt_Reb_PSAC[icent] ->SetTitle("Pair Sign Acc. Diff. Corr. Factor in Cent: "+centTitle[icent]);
			H2d_MvsMt_Reb_PSAC[icent] ->Draw("LEGO2Z");
			//c1->SaveAs(Form("./outAnaPlots/plots_PSAC/H2d_MvsMt_Reb_PSAC_icent%d.png", icent));
			//c1->SaveAs(Form("./outAnaPlots/plots_PSAC/H2d_MvsMt_Reb_PSAC_icent%d.pdf", icent));

			updatePDF(c1,mypdf);
			c1->Clear();

			c1->Divide(nColumns,nRows);
			for(Int_t ipad=0; ipad<nColumns*nRows; ipad++)
			{
				c1->cd(ipad+1);
				setPad(0.1,0.02,0.08,0.1);
			}

			for(int imt=0; imt<Htemmt1->GetNbinsX(); imt++)
			{
				Int_t idx = imt;
				c1->cd(idx%nPads+1);

				cout<<Form("%.2f<mt<%.2f", mtBinBDs[imt], mtBinBDs[imt+1])<<endl;
				TString temmtName = Form("%.2f<mt<%.2f", mtBinBDs[imt], mtBinBDs[imt+1]);

				H1D_PSAC_4Mt[icent][imt] = (TH1D*) H2d_MvsMt_Reb_PSAC[icent]   ->ProjectionY(Form("H1D_PSAC_4Mt_icent%d_imt%d", icent, imt), imt+1, imt+1);
				H1D_PSAC_4Mt[icent][imt]->SetTitle("Pair S.A.C Factor in Cent: "+centTitle[icent]+" "+temmtName);
				H1D_PSAC_4Mt[icent][imt]->SetYTitle("PSAC Factor");
				H1D_PSAC_4Mt[icent][imt]->SetMarkerStyle(24); 
				H1D_PSAC_4Mt[icent][imt]->SetMarkerSize(0.3); 
				H1D_PSAC_4Mt[icent][imt]->SetMarkerColor(1); 
				H1D_PSAC_4Mt[icent][imt]->SetLineColor(1); 
				H1D_PSAC_4Mt[icent][imt]->SetAxisRange(0.90, 1.10, "y");
				H1D_PSAC_4Mt[icent][imt]->Draw("pe");

				//c1->SaveAs(Form("./outAnaPlots/plots_PSAC/H1D_PSAC_icent%d_imt%d.png", icent, imt));
				if( idx%nPads == nPads-1 ) updatePDF(c1,mypdf);
			}
			
			if((Htemmt1->GetNbinsX())%nPads!=0) updatePDF(c1,mypdf);

			c1->Clear();

		}//icent

		mypdf->On();
		mypdf->Close();
		delete c1;
	}//if viewflag == 1
}
//------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------
void preparePairEffAndAcc(int viewFlag)
{
	for(int icent=0; icent<nCentBins; icent++)
	{
		//inf_eff[icent] = new TFile(Form("./inputfiles/pairEff/out4PairEff_icent%d_Vz35cm.root",icent),"read"); //Vz<35cm
		inf_eff[icent] = new TFile(Form("./inputfiles/pairEff/out4PairEff_icent%d_noCKTSumWT_Vz35cm.root",icent),"read"); //Vz<35cm
		cout<<"read in: "<<inf_eff[icent]->GetName()<<endl;

		for(int icut=0; icut<nCuts; icut++)//need to be very very careful of the icut index when fill the histograms
		{
			HRc_MvsPt[icent][icut] = (TH2D*)inf_eff[icent]->Get( Form("hRc_MvsPt_icut%d", icut) );
			HRc_MvsPt[icent][icut] -> SetName(Form("HRc_MvsPt_icent%d_icut%d", icent, icut));
			HRc_MvsPt[icent][icut] -> SetYTitle("M_{ee} (GeV/c^{2})");

			HRc_MvsMt[icent][icut] = (TH2D*)inf_eff[icent]->Get( Form("hRc_MvsMt_icut%d", icut) );
			HRc_MvsMt[icent][icut] -> SetName(Form("HRc_MvsMt_icent%d_icut%d", icent, icut));
			HRc_MvsMt[icent][icut] -> SetYTitle("M_{ee} (GeV/c^{2})");

			TString ptHistName = Form("HRc_MvsPt_Reb_icent%d_icut%d", icent, icut);
			TString mtHistName = Form("HRc_MvsMt_Reb_icent%d_icut%d", icent, icut);

			if(     icent==0)
			{
				HRc_MvsPt_Reb[icent][icut] 
					= (TH2D*) rebHist2d(HRc_MvsPt[icent][icut], ptHistName, nptBins,ptBds4Corr, nMBins_ct0,mBds4Corr_ct0,"XY","");
				HRc_MvsMt_Reb[icent][icut] 
					= (TH2D*) rebHist2d(HRc_MvsMt[icent][icut], mtHistName, nmtBins,mtBinBDs, nMBins_ct0,mBds4Corr_ct0,"XY","");
			}
			else if(icent==1)
			{
				HRc_MvsPt_Reb[icent][icut] 
					= (TH2D*) rebHist2d(HRc_MvsPt[icent][icut], ptHistName, nptBins,ptBds4Corr, nMBins_ct1,mBds4Corr_ct1,"XY","");
				HRc_MvsMt_Reb[icent][icut] 
					= (TH2D*) rebHist2d(HRc_MvsMt[icent][icut], mtHistName, nmtBins,mtBinBDs, nMBins_ct1,mBds4Corr_ct1,"XY","");
			}
			else if(icent==2)
			{
				HRc_MvsPt_Reb[icent][icut] 
					= (TH2D*) rebHist2d(HRc_MvsPt[icent][icut], ptHistName, nptBins,ptBds4Corr, nMBins_ct2,mBds4Corr_ct2,"XY","");
				HRc_MvsMt_Reb[icent][icut] 
					= (TH2D*) rebHist2d(HRc_MvsMt[icent][icut], mtHistName, nmtBins,mtBinBDs, nMBins_ct2,mBds4Corr_ct2,"XY","");
			}
			else if(icent==3)
			{
				HRc_MvsPt_Reb[icent][icut] 
					= (TH2D*) rebHist2d(HRc_MvsPt[icent][icut], ptHistName, nptBins,ptBds4Corr, nMBins_ct3,mBds4Corr_ct3,"XY","");
				HRc_MvsMt_Reb[icent][icut] 
					= (TH2D*) rebHist2d(HRc_MvsMt[icent][icut], mtHistName, nmtBins,mtBinBDs, nMBins_ct3,mBds4Corr_ct3,"XY","");
			}
		}//icut

		//calculate pair effiency
		H2d_MvsPt_PairEff[icent] = (TH2D*) calRatioH1toH2(HRc_MvsPt_Reb[icent][4], HRc_MvsPt_Reb[icent][2], Form("H2d_MvsPt_PairEff_icent%d",icent), "CORR");
		H2d_MvsMt_PairEff[icent] = (TH2D*) calRatioH1toH2(HRc_MvsMt_Reb[icent][4], HRc_MvsMt_Reb[icent][2], Form("H2d_MvsMt_PairEff_icent%d",icent), "CORR");

		//calculate pair acceptance
		H2d_MvsPt_PairAcc[icent] = (TH2D*) calRatioH1toH2(HRc_MvsPt_Reb[icent][2], HRc_MvsPt_Reb[icent][1], Form("H2d_MvsPt_PairAcc_icent%d",icent), "CORR");
		H2d_MvsMt_PairAcc[icent] = (TH2D*) calRatioH1toH2(HRc_MvsMt_Reb[icent][2], HRc_MvsMt_Reb[icent][1], Form("H2d_MvsMt_PairAcc_icent%d",icent), "CORR");

		//calcualte pairEffAcc
		H2d_MvsPt_PairEffAcc[icent] = (TH2D*) calRatioH1toH2(HRc_MvsPt_Reb[icent][4], HRc_MvsPt_Reb[icent][1], Form("H2d_MvsPt_PairEffAcc_icent%d",icent), "CORR");
		H2d_MvsMt_PairEffAcc[icent] = (TH2D*) calRatioH1toH2(HRc_MvsMt_Reb[icent][4], HRc_MvsMt_Reb[icent][1], Form("H2d_MvsMt_PairEffAcc_icent%d",icent), "CORR");

	}//icent

	if(viewFlag==1)
	{
		//---------------------------------------------------------------------------
		//make plots of the pair reconstruction efficiency in 2D
		//---------------------------------------------------------------------------
		TCanvas *c1 = new TCanvas("c1");
		gStyle->SetOptStat(0);
		system("mkdir -p outAnaPlots/plots_PairEff");

		TPDF* mypdf = new TPDF("./outAnaPlots/plots_PairEff/Run18_27GeV_AuAu_PairEff_factors_Summary.pdf",111);
		mypdf->Off();
		updatePDF(c1,mypdf);

		Int_t nColumns = 3;
		Int_t nRows    = 3;
		Int_t nPads    = nColumns*nRows;

		TH1D* Htempt1 = (TH1D*)H2d_MvsPt_PairEff[0]->ProjectionX("Htempt1");
		TH1D* H1D_PairEff[nCentBins][nptBins];
		for(int icent=0; icent<nCentBins; icent++)
		{
			c1->cd();
			H2d_MvsPt_PairEff[icent] ->SetTitle("Pair Reconstruction Efficiency in Cent: "+centTitle[icent]);
			H2d_MvsPt_PairEff[icent] ->GetYaxis()->CenterTitle();
			H2d_MvsPt_PairEff[icent] ->GetYaxis()->SetTitleSize(0.05);
			H2d_MvsPt_PairEff[icent] ->GetYaxis()->SetTitleOffset(1.25);
			H2d_MvsPt_PairEff[icent] ->GetXaxis()->SetTitleSize(0.05);
			H2d_MvsPt_PairEff[icent] ->GetXaxis()->SetTitleOffset(1.30);
			H2d_MvsPt_PairEff[icent] ->GetXaxis()->CenterTitle();
			H2d_MvsPt_PairEff[icent] ->Draw("LEGO2Z");
			c1->SaveAs(Form("./outAnaPlots/plots_PairEff/H2d_MvsPt_PairEff_icent%d.png", icent));
			c1->SaveAs(Form("./outAnaPlots/plots_PairEff/H2d_MvsPt_PairEff_icent%d.pdf", icent));

			updatePDF(c1,mypdf);
			c1->Clear();

			c1->Divide(nColumns,nRows);
			for(Int_t ipad=0; ipad<nColumns*nRows; ipad++)
			{
				c1->cd(ipad+1);
				setPad(0.1,0.02,0.08,0.1);
			}

			for(int ipt=0; ipt<Htempt1->GetNbinsX(); ipt++)
			{
				Int_t idx = ipt;
				c1->cd(idx%nPads+1);

				//const Int_t nptBins = 6;
				//Double_t ptBds4Corr[nptBins+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0};
				cout<<Form("%.2f<pt<%.2f", ptBds4Corr[ipt], ptBds4Corr[ipt+1])<<endl;
				TString temptName = Form("%.2f<pt<%.2f", ptBds4Corr[ipt], ptBds4Corr[ipt+1]);

				H1D_PairEff[icent][ipt] = (TH1D*) H2d_MvsPt_PairEff[icent]   ->ProjectionY(Form("H1D_PairEff_icent%d_ipt%d", icent, ipt), ipt+1, ipt+1);

				H1D_PairEff[icent][ipt]->SetTitle("Pair RcEff. in Cent: "+centTitle[icent]+" "+temptName);
				H1D_PairEff[icent][ipt]->SetYTitle("Pair Rec. Eff");
				H1D_PairEff[icent][ipt]->SetMarkerStyle(24);
				H1D_PairEff[icent][ipt]->SetMarkerSize(0.3);
				H1D_PairEff[icent][ipt]->SetMarkerColor(1);
				H1D_PairEff[icent][ipt]->SetLineColor(1);
				H1D_PairEff[icent][ipt]->SetAxisRange(0., 0.5, "y");
				H1D_PairEff[icent][ipt]->Draw("pe");

				//c1->SaveAs(Form("./outAnaPlots/plots_PairEff/H1D_PairEff_icent%d_ipt%d.png", icent, ipt));
				if( idx%nPads == nPads-1 ) updatePDF(c1,mypdf);
			}//ipt

			if((Htempt1->GetNbinsX())%nPads!=0) updatePDF(c1,mypdf);

			c1->Clear();

		}//icent
		//---------------------------------------------------------------------------
		//draw for Mt
		//---------------------------------------------------------------------------

		TH1D* Htemmt1 = (TH1D*)H2d_MvsMt_PairEff[0]->ProjectionX("Htemmt1");
		TH1D* H1D_PairEff_4Mt[nCentBins][nmtBins];
		for(int icent=0; icent<nCentBins; icent++)
		{
			c1->cd();
			H2d_MvsMt_PairEff[icent] ->SetTitle("Pair Reconstruction Efficiency in Cent: "+centTitle[icent]);
			H2d_MvsMt_PairEff[icent] ->Draw("LEGO2Z");
			//c1->SaveAs(Form("./outAnaPlots/plots_PairEff/H2d_MvsMt_PairEff_icent%d.png", icent));
			//c1->SaveAs(Form("./outAnaPlots/plots_PairEff/H2d_MvsMt_PairEff_icent%d.pdf", icent));

			updatePDF(c1,mypdf);
			c1->Clear();

			c1->Divide(nColumns,nRows);
			for(Int_t ipad=0; ipad<nColumns*nRows; ipad++)
			{
				c1->cd(ipad+1);
				setPad(0.1,0.02,0.08,0.1);
			}

			for(int imt=0; imt<Htemmt1->GetNbinsX(); imt++)
			{
				Int_t idx = imt;
				c1->cd(idx%nPads+1);

				//const Int_t nmtBins = 6;
				//Double_t mtBinBDs[nmtBins+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0};
				cout<<Form("%.2f<mt<%.2f", mtBinBDs[imt], mtBinBDs[imt+1])<<endl;
				TString temmtName = Form("%.2f<mt<%.2f", mtBinBDs[imt], mtBinBDs[imt+1]);

				H1D_PairEff[icent][imt] = (TH1D*) H2d_MvsMt_PairEff[icent]   ->ProjectionY(Form("H1D_PairEff_icent%d_imt%d", icent, imt), imt+1, imt+1);

				H1D_PairEff[icent][imt]->SetTitle("Pair RcEff. in Cent: "+centTitle[icent]+" "+temmtName);
				H1D_PairEff[icent][imt]->SetYTitle("Pair Rec. Eff");
				H1D_PairEff[icent][imt]->SetMarkerStyle(24);
				H1D_PairEff[icent][imt]->SetMarkerSize(0.3);
				H1D_PairEff[icent][imt]->SetMarkerColor(1);
				H1D_PairEff[icent][imt]->SetLineColor(1);
				H1D_PairEff[icent][imt]->SetAxisRange(0., 0.5, "y");
				H1D_PairEff[icent][imt]->Draw("pe");

				//c1->SaveAs(Form("./outAnaPlots/plots_PairEff/H1D_PairEff_icent%d_imt%d.png", icent, imt));
				if( idx%nPads == nPads-1 ) updatePDF(c1,mypdf);
			}//imt

			if((Htemmt1->GetNbinsX())%nPads!=0) updatePDF(c1,mypdf);

			c1->Clear();

		}//icent

		mypdf->On();
		mypdf->Close();
		delete c1;
	}

	//plot the pair Acc
	if(viewFlag==2)
	{
		//---------------------------------------------------------------------------
		//make plots of the pair reconstruction efficiency in 2D
		//---------------------------------------------------------------------------
		TCanvas *c1 = new TCanvas("c1");
		gStyle->SetOptStat(0);
		system("mkdir -p outAnaPlots/plots_PairAcc");

		TPDF* mypdf = new TPDF("./outAnaPlots/plots_PairAcc/Run18_27GeV_AuAu_PairAcc_factors_Summary.pdf",111);
		mypdf->Off();
		updatePDF(c1,mypdf);

		Int_t nColumns = 3;
		Int_t nRows    = 3;
		Int_t nPads    = nColumns*nRows;

		TH1D* Htempt1 = (TH1D*)H2d_MvsPt_PairAcc[0]->ProjectionX("Htempt1");
		TH1D* H1D_PairAcc[nCentBins][nptBins];
		for(int icent=0; icent<nCentBins; icent++)
		{
			c1->cd();
			H2d_MvsPt_PairAcc[icent] ->SetTitle("Pair Reconstruction Acc. in Cent: "+centTitle[icent]);
			H2d_MvsPt_PairAcc[icent] ->GetYaxis()->CenterTitle();
			H2d_MvsPt_PairAcc[icent] ->GetYaxis()->SetTitleSize(0.05);
			H2d_MvsPt_PairAcc[icent] ->GetYaxis()->SetTitleOffset(1.25);
			H2d_MvsPt_PairAcc[icent] ->GetXaxis()->SetTitleSize(0.05);
			H2d_MvsPt_PairAcc[icent] ->GetXaxis()->SetTitleOffset(1.30);
			H2d_MvsPt_PairAcc[icent] ->GetXaxis()->CenterTitle();
			H2d_MvsPt_PairAcc[icent] ->Draw("LEGO2Z");
			c1->SaveAs(Form("./outAnaPlots/plots_PairAcc/H2d_MvsPt_PairAcc_icent%d.png", icent));
			c1->SaveAs(Form("./outAnaPlots/plots_PairAcc/H2d_MvsPt_PairAcc_icent%d.pdf", icent));
;
			updatePDF(c1,mypdf);
			c1->Clear();

			c1->Divide(nColumns,nRows);
			for(Int_t ipad=0; ipad<nColumns*nRows; ipad++)
			{
				c1->cd(ipad+1);
				setPad(0.1,0.02,0.08,0.1);
			}

			for(int ipt=0; ipt<Htempt1->GetNbinsX(); ipt++)
			{
				Int_t idx = ipt;
				c1->cd(idx%nPads+1);

				//const Int_t nptBins = 6;
				//Double_t ptBds4Corr[nptBins+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0};
				cout<<Form("%.2f<pt<%.2f", ptBds4Corr[ipt], ptBds4Corr[ipt+1])<<endl;
				TString temptName = Form("%.2f<pt<%.2f", ptBds4Corr[ipt], ptBds4Corr[ipt+1]);

				H1D_PairAcc[icent][ipt] = (TH1D*) H2d_MvsPt_PairAcc[icent]   ->ProjectionY(Form("H1D_PairAcc_icent%d_ipt%d", icent, ipt), ipt+1, ipt+1);

				H1D_PairAcc[icent][ipt]->SetTitle("Pair RcAcc. in Cent: "+centTitle[icent]+" "+temptName);
				H1D_PairAcc[icent][ipt]->SetYTitle("Pair Rec. Acc");
				H1D_PairAcc[icent][ipt]->SetMarkerStyle(24);
				H1D_PairAcc[icent][ipt]->SetMarkerSize(0.3);
				H1D_PairAcc[icent][ipt]->SetMarkerColor(1);
				H1D_PairAcc[icent][ipt]->SetLineColor(1);
				H1D_PairAcc[icent][ipt]->SetAxisRange(0., 1.0, "y");
				H1D_PairAcc[icent][ipt]->Draw("pe");

				//c1->SaveAs(Form("./outAnaPlots/plots_PairAcc/H1D_PairAcc_icent%d_ipt%d.png", icent, ipt));
				if( idx%nPads == nPads-1 ) updatePDF(c1,mypdf);
			}//ipt

			if((Htempt1->GetNbinsX())%nPads!=0) updatePDF(c1,mypdf);

			c1->Clear();

		}//icent

		//---------------------------------------------------------------------------
		//draw for Mt
		//---------------------------------------------------------------------------
		TH1D* Htemmt1 = (TH1D*)H2d_MvsMt_PairAcc[0]->ProjectionX("Htemmt1");
		TH1D* H1D_PairAcc_4Mt[nCentBins][nmtBins];
		for(int icent=0; icent<nCentBins; icent++)
		{
			c1->cd();
			H2d_MvsMt_PairAcc[icent] ->SetTitle("Pair Reconstruction Acc. in Cent: "+centTitle[icent]);
			H2d_MvsMt_PairAcc[icent] ->Draw("LEGO2Z");
			//c1->SaveAs(Form("./outAnaPlots/plots_PairAcc/H2d_MvsMt_PairAcc_icent%d.png", icent));
			//c1->SaveAs(Form("./outAnaPlots/plots_PairAcc/H2d_MvsMt_PairAcc_icent%d.pdf", icent));

			updatePDF(c1,mypdf);
			c1->Clear();

			c1->Divide(nColumns,nRows);
			for(Int_t ipad=0; ipad<nColumns*nRows; ipad++)
			{
				c1->cd(ipad+1);
				setPad(0.1,0.02,0.08,0.1);
			}

			for(int imt=0; imt<Htemmt1->GetNbinsX(); imt++)
			{
				Int_t idx = imt;
				c1->cd(idx%nPads+1);

				//const Int_t nmtBins = 6;
				//Double_t mtBinBDs[nmtBins+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0};
				cout<<Form("%.2f<mt<%.2f", mtBinBDs[imt], mtBinBDs[imt+1])<<endl;
				TString temmtName = Form("%.2f<mt<%.2f", mtBinBDs[imt], mtBinBDs[imt+1]);

				H1D_PairAcc[icent][imt] = (TH1D*) H2d_MvsMt_PairAcc[icent]   ->ProjectionY(Form("H1D_PairAcc_icent%d_imt%d", icent, imt), imt+1, imt+1);

				H1D_PairAcc[icent][imt]->SetTitle("Pair RcAcc. in Cent: "+centTitle[icent]+" "+temmtName);
				H1D_PairAcc[icent][imt]->SetYTitle("Pair Rec. Acc");
				H1D_PairAcc[icent][imt]->SetMarkerStyle(24);
				H1D_PairAcc[icent][imt]->SetMarkerSize(0.3);
				H1D_PairAcc[icent][imt]->SetMarkerColor(1);
				H1D_PairAcc[icent][imt]->SetLineColor(1);
				H1D_PairAcc[icent][imt]->SetAxisRange(0., 1.0, "y");
				H1D_PairAcc[icent][imt]->Draw("pe");

				//c1->SaveAs(Form("./outAnaPlots/plots_PairAcc/H1D_PairAcc_icent%d_imt%d.png", icent, imt));
				if( idx%nPads == nPads-1 ) updatePDF(c1,mypdf);
			}//imt

			if((Htemmt1->GetNbinsX())%nPads!=0) updatePDF(c1,mypdf);

			c1->Clear();

		}//icent

		mypdf->On();
		mypdf->Close();
		delete c1;
	}
}
//------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------
void getInitCKTs()
{
	TFile* inf_ckt = new TFile("../../Cocktails/CombineAllCKTs/outCKTplots/out4Vz35cm/outAll_CKTS_AuAu27GeV_scaleFlag0_Vz35cm.root", "read");
	cout<<"readin : "<<inf_ckt->GetName()<<endl;

	for(int icent=0; icent<nCentBins; icent++)
	{
		for(int ickt=0; ickt<nCKTs+2; ickt++)//for different ckt components
		{
			//---------------------------------------------------------------------------
			//get 2d CKT histograms
			H2D_MvsPt_CKT[icent][ickt] = (TH2D*) inf_ckt->Get(Form("H2D_MvsPt_icent%d_ickt%d", icent, ickt));
			H2D_MvsMt_CKT[icent][ickt] = (TH2D*) inf_ckt->Get(Form("H2D_MvsMt_icent%d_ickt%d", icent, ickt));
			//----------------------------------------------------------------------------------------------------------------
			H2D_MvsPt_CKTbfAcc[icent][ickt]  = (TH2D*) inf_ckt->Get(Form("H2D_MvsPt_CKTbfAcc_icent%d_ickt%d",  icent, ickt));
			//----------------------------------------------------------------------------------------------------------------
			//----------------------------------------------------------------------------------------------------------------

			TString PtHistName = Form("H2D_MvsPt_Reb_icent%d_ickt%d",icent, ickt);
			TString MtHistName = Form("H2D_MvsMt_Reb_icent%d_ickt%d",icent, ickt);

			//---------------------------------------------------------------------------
			//rebin the histograms to be same binning as real data
			if(     icent==0)
			{
				H2D_MvsPt_CKT_Reb[icent][ickt] 
					= (TH2D*) rebHist2d(H2D_MvsPt_CKT[icent][ickt], PtHistName, nptBins,ptBds4Corr, nMBins_ct0,mBds4Corr_ct0,"XY","");
				H2D_MvsMt_CKT_Reb[icent][ickt] 
					= (TH2D*) rebHist2d(H2D_MvsMt_CKT[icent][ickt], MtHistName, nmtBins,mtBinBDs,   nMBins_ct0,mBds4Corr_ct0,"XY","");
			}
			else if(icent==1)
			{
				H2D_MvsPt_CKT_Reb[icent][ickt] 
					= (TH2D*) rebHist2d(H2D_MvsPt_CKT[icent][ickt], PtHistName, nptBins,ptBds4Corr, nMBins_ct1,mBds4Corr_ct1,"XY","");
				H2D_MvsMt_CKT_Reb[icent][ickt] 
					= (TH2D*) rebHist2d(H2D_MvsMt_CKT[icent][ickt], MtHistName, nmtBins,mtBinBDs,   nMBins_ct1,mBds4Corr_ct1,"XY","");
			}
			else if(icent==2)
			{
				H2D_MvsPt_CKT_Reb[icent][ickt] 
					= (TH2D*) rebHist2d(H2D_MvsPt_CKT[icent][ickt], PtHistName, nptBins,ptBds4Corr, nMBins_ct2,mBds4Corr_ct2,"XY","");
				H2D_MvsMt_CKT_Reb[icent][ickt] 
					= (TH2D*) rebHist2d(H2D_MvsMt_CKT[icent][ickt], MtHistName, nmtBins,mtBinBDs,   nMBins_ct2,mBds4Corr_ct2,"XY","");
			}
			else if(icent==3)
			{
				H2D_MvsPt_CKT_Reb[icent][ickt] 
					= (TH2D*) rebHist2d(H2D_MvsPt_CKT[icent][ickt], PtHistName, nptBins,ptBds4Corr, nMBins_ct3,mBds4Corr_ct3,"XY","");
				H2D_MvsMt_CKT_Reb[icent][ickt] 
					= (TH2D*) rebHist2d(H2D_MvsMt_CKT[icent][ickt], MtHistName, nmtBins,mtBinBDs,   nMBins_ct3,mBds4Corr_ct3,"XY","");
			}
			//---------------------------------------------------------------------------
		}//ickt, 6+2 cocktails

		//the CKTSum
		H2D_MvsPt_totCKT[icent]        = (TH2D*) inf_ckt->Get(Form("H2D_MvsPt_totCKT_icent%d",       icent));
		H2D_MvsMt_totCKT[icent]        = (TH2D*) inf_ckt->Get(Form("H2D_MvsMt_totCKT_icent%d",       icent));
		
		//---------------------------------------------------------------------------
		//just for check the pairAccEff on the CKT
		H2D_MvsPt_totCKTbfAcc[icent]   = (TH2D*) inf_ckt->Get(Form("H2D_MvsPt_totCKTbfAcc_icent%d",  icent));
		H2D_MvsPt_totCKTwAccEff[icent] = (TH2D*) inf_ckt->Get(Form("H2D_MvsPt_totCKTwAccEff_icent%d",icent));
		//---------------------------------------------------------------------------

		TString totPtHistName  = Form("H2D_MvsPt_totCKT_Reb_icent%d",       icent);
		TString totMtHistName  = Form("H2D_MvsMt_totCKT_Reb_icent%d",       icent);
		TString totPtHistName2 = Form("H2D_MvsPt_totCKTbfAcc_Reb_icent%d",  icent);
		TString totPtHistName3 = Form("H2D_MvsPt_totCKTwAccEff_Reb_icent%d",icent);

		if(     icent==0)
		{
			H2D_MvsPt_totCKT_Reb[icent]
				= (TH2D*) rebHist2d(H2D_MvsPt_totCKT[icent], totPtHistName, nptBins,ptBds4Corr, nMBins_ct0,mBds4Corr_ct0,"XY","");
			H2D_MvsMt_totCKT_Reb[icent]
				= (TH2D*) rebHist2d(H2D_MvsMt_totCKT[icent], totMtHistName, nmtBins,mtBinBDs,   nMBins_ct0,mBds4Corr_ct0,"XY","");

			H2D_MvsPt_totCKTbfAcc_Reb[icent]
				= (TH2D*) rebHist2d(H2D_MvsPt_totCKTbfAcc[icent],   totPtHistName2, nptBins,ptBds4Corr, nMBins_ct0,mBds4Corr_ct0,"XY","");
			H2D_MvsPt_totCKTwAccEff_Reb[icent]
				= (TH2D*) rebHist2d(H2D_MvsPt_totCKTwAccEff[icent], totPtHistName3, nptBins,ptBds4Corr, nMBins_ct0,mBds4Corr_ct0,"XY","");
		}
		else if(icent==1)
		{
			H2D_MvsPt_totCKT_Reb[icent] 
				= (TH2D*) rebHist2d(H2D_MvsPt_totCKT[icent], totPtHistName, nptBins,ptBds4Corr, nMBins_ct1,mBds4Corr_ct1,"XY","");
			H2D_MvsMt_totCKT_Reb[icent] 
				= (TH2D*) rebHist2d(H2D_MvsMt_totCKT[icent], totMtHistName, nmtBins,mtBinBDs,   nMBins_ct1,mBds4Corr_ct1,"XY","");

			H2D_MvsPt_totCKTbfAcc_Reb[icent] 
				= (TH2D*) rebHist2d(H2D_MvsPt_totCKTbfAcc[icent],   totPtHistName2, nptBins,ptBds4Corr, nMBins_ct1,mBds4Corr_ct1,"XY","");
			H2D_MvsPt_totCKTwAccEff_Reb[icent] 
				= (TH2D*) rebHist2d(H2D_MvsPt_totCKTwAccEff[icent], totPtHistName3, nptBins,ptBds4Corr, nMBins_ct1,mBds4Corr_ct1,"XY","");
		}
		else if(icent==2)
		{
			H2D_MvsPt_totCKT_Reb[icent] 
				= (TH2D*) rebHist2d(H2D_MvsPt_totCKT[icent], totPtHistName, nptBins,ptBds4Corr, nMBins_ct2,mBds4Corr_ct2,"XY","");
			H2D_MvsMt_totCKT_Reb[icent] 
				= (TH2D*) rebHist2d(H2D_MvsMt_totCKT[icent], totMtHistName, nmtBins,mtBinBDs,   nMBins_ct2,mBds4Corr_ct2,"XY","");

			H2D_MvsPt_totCKTbfAcc_Reb[icent] 
				= (TH2D*) rebHist2d(H2D_MvsPt_totCKTbfAcc[icent],   totPtHistName2, nptBins,ptBds4Corr, nMBins_ct2,mBds4Corr_ct2,"XY","");
			H2D_MvsPt_totCKTwAccEff_Reb[icent] 
				= (TH2D*) rebHist2d(H2D_MvsPt_totCKTwAccEff[icent], totPtHistName3, nptBins,ptBds4Corr, nMBins_ct2,mBds4Corr_ct2,"XY","");
		}
		else if(icent==3)
		{
			H2D_MvsPt_totCKT_Reb[icent] 
				= (TH2D*) rebHist2d(H2D_MvsPt_totCKT[icent], totPtHistName, nptBins,ptBds4Corr, nMBins_ct3,mBds4Corr_ct3,"XY","");
			H2D_MvsMt_totCKT_Reb[icent] 
				= (TH2D*) rebHist2d(H2D_MvsMt_totCKT[icent], totMtHistName, nmtBins,mtBinBDs,   nMBins_ct3,mBds4Corr_ct3,"XY","");

			H2D_MvsPt_totCKTbfAcc_Reb[icent] 
				= (TH2D*) rebHist2d(H2D_MvsPt_totCKTbfAcc[icent],   totPtHistName2, nptBins,ptBds4Corr, nMBins_ct3,mBds4Corr_ct3,"XY","");
			H2D_MvsPt_totCKTwAccEff_Reb[icent] 
				= (TH2D*) rebHist2d(H2D_MvsPt_totCKTwAccEff[icent], totPtHistName3, nptBins,ptBds4Corr, nMBins_ct3,mBds4Corr_ct3,"XY","");
		}
	}//icent
	//------------------------------------------------------------------------------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------
void prepareCKT4fit2Data()
{
	cout<<"----------------------------------------------------------"<<endl;
	cout<<"prepare CKT for CKT_noOmegaPhi+a*Omega+b*Phi+c*Theory"<<endl;
	cout<<"----------------------------------------------------------"<<endl;

	for(int icent=0; icent<nCentBins; icent++)
	{
		//Add other components except the omega and phi contributions
		for(int ickt=0; ickt<nCKTs+2; ickt++)//for different ckt components
		{
			if(ickt==0)      //Pi0
			{
				H2D_MvsPt_CKTnoVect[icent] = (TH2D*) H2D_MvsPt_CKT_Reb[icent][ickt]->Clone();
			}
			else if(ickt==2) //Omega
			{
				H2D_MvsPt_CKTw[icent] = (TH2D*) H2D_MvsPt_CKT_Reb[icent][ickt]->Clone();
			}
			else if(ickt==4) //Phi
			{
				H2D_MvsPt_CKTphi[icent] = (TH2D*) H2D_MvsPt_CKT_Reb[icent][ickt]->Clone();
			}
			else             //if ickt!=2&ickt!=4, then add together to get CktSum without omega and phi components
			{
				H2D_MvsPt_CKTnoVect[icent] ->Add(H2D_MvsPt_CKT_Reb[icent][ickt]);
			}

			if(ickt==6) H2D_MvsPt_CKTcc[icent] = (TH2D*) H2D_MvsPt_CKT_Reb[icent][ickt]->Clone();
			if(ickt==7) H2D_MvsPt_CKTdy[icent] = (TH2D*) H2D_MvsPt_CKT_Reb[icent][ickt]->Clone();
		}//ickt
	}//icent

}
//------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------
void corrAndGetSignals()
{
	fRcut = new TF1("fRcut", fun_Rcut, 0., 0.415, 2);
	fRcut -> SetParameters(-1.02, 0.17);
	//Apply the Arched cut: Rcut on the 2D hists to reject the edge effects
	//apply the Rcut reject the edge effects
	for(int icent=0; icent<nCentBins; icent++)
	{
		H2d_MvsPt_PairEffAcc[icent]        = (TH2D*)applyRcut2D(H2d_MvsPt_PairEffAcc[icent],        Form("H2d_MvsPt_PairEffAccRcut_icent%d",       icent), fRcut);
		H2d_MvsPt_PairEff[icent]           = (TH2D*)applyRcut2D(H2d_MvsPt_PairEff[icent],           Form("H2d_MvsPt_PairEffRcut_icent%d",          icent), fRcut);
		H2d_MvsPt_PairAcc[icent]           = (TH2D*)applyRcut2D(H2d_MvsPt_PairAcc[icent],           Form("H2d_MvsPt_PairAccRcut_icent%d",          icent), fRcut);

		H2d_MvsPt_Reb_US[icent]            = (TH2D*)applyRcut2D(H2d_MvsPt_Reb_US[icent],            Form("H2d_MvsPt_Reb_USRcut_icent%d",           icent), fRcut);
		H2d_MvsPt_Reb_LSPos[icent]         = (TH2D*)applyRcut2D(H2d_MvsPt_Reb_LSPos[icent],         Form("H2d_MvsPt_Reb_LSPosRcut_icent%d",        icent), fRcut);
		H2d_MvsPt_Reb_LSNeg[icent]         = (TH2D*)applyRcut2D(H2d_MvsPt_Reb_LSNeg[icent],         Form("H2d_MvsPt_Reb_LSNegRcut_icent%d",        icent), fRcut);
		H2d_MvsPt_Reb_MixUS[icent]         = (TH2D*)applyRcut2D(H2d_MvsPt_Reb_MixUS[icent],         Form("H2d_MvsPt_Reb_MixUSRcut_icent%d",        icent), fRcut);

		H2d_MvsPt_LSPos_PSAC[icent]        = (TH2D*)applyRcut2D(H2d_MvsPt_LSPos_PSAC[icent],        Form("H2d_MvsPt_LSPos_PSACRcut_icent%d",       icent), fRcut);
		H2d_MvsPt_LSNeg_PSAC[icent]        = (TH2D*)applyRcut2D(H2d_MvsPt_LSNeg_PSAC[icent],        Form("H2d_MvsPt_LSNeg_PSACRcut_icent%d",       icent), fRcut);

		H2D_MvsPt_totCKT_Reb[icent]        = (TH2D*)applyRcut2D(H2D_MvsPt_totCKT_Reb[icent],        Form("H2D_MvsPt_totCKT_RebRcut_icent%d",       icent), fRcut);
		H2D_MvsPt_totCKTbfAcc_Reb[icent]   = (TH2D*)applyRcut2D(H2D_MvsPt_totCKTbfAcc_Reb[icent],   Form("H2D_MvsPt_totCKTbfAcc_RebRcut_icent%d",  icent), fRcut);
		H2D_MvsPt_totCKTwAccEff_Reb[icent] = (TH2D*)applyRcut2D(H2D_MvsPt_totCKTwAccEff_Reb[icent], Form("H2D_MvsPt_totCKTwAccEff_RebRcut_icent%d",icent), fRcut);
		
		//apply on sub cocktails
		H2D_MvsPt_CKTnoVect[icent]         = (TH2D*)applyRcut2D(H2D_MvsPt_CKTnoVect[icent],         Form("H2D_MvsPt_CKTnoVect_Rcut_icent%d",icent), fRcut);
		H2D_MvsPt_CKTw[icent]              = (TH2D*)applyRcut2D(H2D_MvsPt_CKTw[icent],              Form("H2D_MvsPt_CKTw_Rcut_icent%d",     icent), fRcut);
		H2D_MvsPt_CKTphi[icent]            = (TH2D*)applyRcut2D(H2D_MvsPt_CKTphi[icent],            Form("H2D_MvsPt_CKTphi_Rcut_icent%d",   icent), fRcut);
		H2D_MvsPt_CKTcc[icent]             = (TH2D*)applyRcut2D(H2D_MvsPt_CKTcc[icent],             Form("H2D_MvsPt_CKTcc_Rcut_icent%d",    icent), fRcut);
		H2D_MvsPt_CKTdy[icent]             = (TH2D*)applyRcut2D(H2D_MvsPt_CKTdy[icent],             Form("H2D_MvsPt_CKTdy_Rcut_icent%d",    icent), fRcut);
	}//icent

	//do all corrections
	for(int icent=0; icent<nCentBins; icent++)
	{
		//---------------------------------------------------------------------------
		//apply pairEff and PairAcc on US
		//---------------------------------------------------------------------------
		//H2d_MvsPt_US_aftEff[icent]    = (TH2D*) calRatioH1toH2( H2d_MvsPt_Reb_US[icent],    H2d_MvsPt_PairEff[icent], Form("H2d_MvsPt_US_aftEff_icent%d",    icent), "NOCORR" );
		//H2d_MvsPt_US_aftEffAcc[icent] = (TH2D*) calRatioH1toH2( H2d_MvsPt_US_aftEff[icent], H2d_MvsPt_PairAcc[icent], Form("H2d_MvsPt_US_aftEffAcc_icent%d", icent), "NOCORR" );
		H2d_MvsPt_US_aftEffAcc[icent]
			= (TH2D*) calRatioH1toH2( H2d_MvsPt_Reb_US[icent],    H2d_MvsPt_PairEffAcc[icent], Form("H2d_MvsPt_US_aftEffAcc_icent%d",    icent), "NOCORR" );
		H2d_MvsPt_MixUS_aftEffAcc[icent]
			= (TH2D*) calRatioH1toH2( H2d_MvsPt_Reb_MixUS[icent], H2d_MvsPt_PairEffAcc[icent], Form("H2d_MvsPt_MixUS_aftEffAcc_icent%d", icent), "NOCORR" );

		//for tem check, don't forget to change it back to normal, just comment this out!!!!!!!
		//H2d_MvsPt_US_aftEffAcc[icent] = (TH2D*) calRatioH1toH2( H2d_MvsPt_Reb_US[icent], H2d_MvsPt_PairEff[icent], Form("H2d_MvsPt_US_aftEffAcc_icent%d", icent), "NOCORR" );
		//---------------------------------------------------------------------------

		//---------------------------------------------------------------------------
		//Apply PSAC_LSPos, PSAC_LSNeg 
		//---------------------------------------------------------------------------
		//H2d_MvsPt_Reb_GmLS_aftPSAC[icent] = calH1xH2(H2d_MvsPt_Reb_GmLS[icent],H2d_MvsPt_Reb_PSAC[icent], Form("H2d_MvsPt_Reb_GmLS_aftPSAC_icent%d", icent));
		H2d_MvsPt_LSPos_aftPSAC[icent]
			= (TH2D*) calH1xH2( H2d_MvsPt_Reb_LSPos[icent],  H2d_MvsPt_LSPos_PSAC[icent],  Form("H2d_MvsPt_LSPos_aftPSAC_icent%d", icent) );
		H2d_MvsPt_LSNeg_aftPSAC[icent]
			= (TH2D*) calH1xH2( H2d_MvsPt_Reb_LSNeg[icent],  H2d_MvsPt_LSNeg_PSAC[icent],  Form("H2d_MvsPt_LSNeg_aftPSAC_icent%d", icent) );

		////---------------------------------------------------------------------------
		////Apply same PSAC for LSPos and LSNeg
		////---------------------------------------------------------------------------
		//H2d_MvsPt_LSPos_aftPSAC[icent] 
		//	= (TH2D*) calH1xH2( H2d_MvsPt_Reb_LSPos[icent],        H2d_MvsPt_Reb_PSAC[icent],  Form("H2d_MvsPt_LSPos_aftPSAC_icent%d",    icent) );
		//H2d_MvsPt_LSNeg_aftPSAC[icent] 
		//	= (TH2D*) calH1xH2( H2d_MvsPt_Reb_LSNeg[icent],        H2d_MvsPt_Reb_PSAC[icent],  Form("H2d_MvsPt_LSNeg_aftPSAC_icent%d",    icent) );
		//---------------------------------------------------------------------------

		////TEST, if Don't use PSAC 
		//H2d_MvsPt_LSPos_aftPSAC[icent] = (TH2D*) H2d_MvsPt_Reb_LSPos[icent]->Clone(Form("H2d_MvsPt_LSPos_aftPSAC_icent%d",    icent));
		//H2d_MvsPt_LSNeg_aftPSAC[icent] = (TH2D*) H2d_MvsPt_Reb_LSNeg[icent]->Clone(Form("H2d_MvsPt_LSNeg_aftPSAC_icent%d",    icent));

		//---------------------------------------------------------------------------
		//Apply pairEff, pairAcc on LSPos, LSNeg
		//---------------------------------------------------------------------------
		//H2d_MvsPt_LSPos_aftPSACEff[icent] //aplly pairEff
		//	= (TH2D*) calRatioH1toH2( H2d_MvsPt_LSPos_aftPSAC[icent],    H2d_MvsPt_PairEff[icent],    Form("H2d_MvsPt_LSPos_aftPSACEff_icent%d",icent), "NOCORR" );
		//H2d_MvsPt_LSPos_aftAllCorr[icent] //apply pairAcc
		//	= (TH2D*) calRatioH1toH2( H2d_MvsPt_LSPos_aftPSACEff[icent], H2d_MvsPt_PairAcc[icent],    Form("H2d_MvsPt_LSPos_aftAllCorr_icent%d", icent), "NOCORR" );
		H2d_MvsPt_LSPos_aftAllCorr[icent]  //if directly apply PairEffAcc
			= (TH2D*) calRatioH1toH2( H2d_MvsPt_LSPos_aftPSAC[icent], H2d_MvsPt_PairEffAcc[icent],    Form("H2d_MvsPt_LSPos_aftAllCorr_icent%d", icent), "NOCORR" );

		//H2d_MvsPt_LSNeg_aftPSACEff[icent] //aplly pairEff
		//	= (TH2D*) calRatioH1toH2( H2d_MvsPt_LSNeg_aftPSAC[icent],    H2d_MvsPt_PairEff[icent],    Form("H2d_MvsPt_LSNeg_aftPSACEff_icent%d",icent), "NOCORR" );
		//H2d_MvsPt_LSNeg_aftAllCorr[icent] //apply pairAcc
		//	= (TH2D*) calRatioH1toH2( H2d_MvsPt_LSNeg_aftPSACEff[icent], H2d_MvsPt_PairAcc[icent],    Form("H2d_MvsPt_LSNeg_aftAllCorr_icent%d", icent), "NOCORR" );
		H2d_MvsPt_LSNeg_aftAllCorr[icent] //if directly apply PairEffAcc
			= (TH2D*) calRatioH1toH2( H2d_MvsPt_LSNeg_aftPSAC[icent], H2d_MvsPt_PairEffAcc[icent],    Form("H2d_MvsPt_LSNeg_aftAllCorr_icent%d", icent), "NOCORR" );
		//---------------------------------------------------------------------------

		//---------------------------------------------------------------------------
		//Apply PairAcc on CKTSum
		//---------------------------------------------------------------------------
		H2D_MvsPt_totCKT_Reb_aftAcc[icent]
			= (TH2D*) calRatioH1toH2( H2D_MvsPt_totCKT_Reb[icent], H2d_MvsPt_PairAcc[icent],   Form("H2D_MvsPt_totCKT_Reb_aftAcc_icent%d", icent), "NOCORR" );

		H2D_MvsPt_CKTnoVect_aftAcc[icent]
			= (TH2D*) calRatioH1toH2( H2D_MvsPt_CKTnoVect[icent],  H2d_MvsPt_PairAcc[icent],   Form("H2D_MvsPt_CKTnoVect_aftAcc_icent%d",  icent), "NOCORR" );
		H2D_MvsPt_CKTw_aftAcc[icent]
			= (TH2D*) calRatioH1toH2( H2D_MvsPt_CKTw[icent],       H2d_MvsPt_PairAcc[icent],   Form("H2D_MvsPt_CKTw_aftAcc_icent%d",       icent), "NOCORR" );
		H2D_MvsPt_CKTphi_aftAcc[icent]
			= (TH2D*) calRatioH1toH2( H2D_MvsPt_CKTphi[icent],     H2d_MvsPt_PairAcc[icent],   Form("H2D_MvsPt_CKTphi_aftAcc_icent%d",     icent), "NOCORR" );
		H2D_MvsPt_CKTcc_aftAcc[icent]
			= (TH2D*) calRatioH1toH2( H2D_MvsPt_CKTcc[icent],      H2d_MvsPt_PairAcc[icent],   Form("H2D_MvsPt_CKTcc_aftAcc_icent%d",       icent), "NOCORR" );
		H2D_MvsPt_CKTdy_aftAcc[icent]
			= (TH2D*) calRatioH1toH2( H2D_MvsPt_CKTdy[icent],      H2d_MvsPt_PairAcc[icent],   Form("H2D_MvsPt_CKTdy_aftAcc_icent%d",     icent), "NOCORR" );
		//---------------------------------------------------------------------------
		//Apply PairEff on the CKTSum_wAccEff, which then could be compared to the CKTSum within STAR Acc ?
		//check whether pairEff could correct the CKT sum to be what they should be
		H2D_MvsPt_totCKT4PairEffTest_Reb[icent]
			= (TH2D*) calRatioH1toH2( H2D_MvsPt_totCKTwAccEff_Reb[icent], H2d_MvsPt_PairEff[icent], Form("H2D_MvsPt_totCKT4PairEffTest_icent%d", icent), "NOCORR" );
	}//icent

	//---------------------------------------------------------------------------
	//project to 1D, the Mass dimension
	//---------------------------------------------------------------------------
	for(int icent=0; icent<nCentBins; icent++)
	{
		const int LowPtBin = H2d_MvsPt_US_aftEffAcc[icent]->GetXaxis()->FindBin(goodPtLow);
		const int HigPtBin = H2d_MvsPt_US_aftEffAcc[icent]->GetXaxis()->FindBin(goodPtHig);

		//project 1D US
		H1d_M_US_aftAllCorr[icent]     = (TH1D*)H2d_MvsPt_US_aftEffAcc[icent]    ->ProjectionY(Form("H1d_M_US_aftAllCorr_icent%d",   icent), LowPtBin, HigPtBin);
		H1d_M_US_aftAllCorr[icent]     ->SetName(Form("H1d_M_US_aftAllCorr_icent%d",    icent));
		H1d_M_MixUS_aftAllCorr[icent]  = (TH1D*)H2d_MvsPt_MixUS_aftEffAcc[icent] ->ProjectionY(Form("H1d_M_MixUS_aftAllCorr_icent%d",icent), LowPtBin, HigPtBin);
		H1d_M_MixUS_aftAllCorr[icent]  ->SetName(Form("H1d_M_MixUS_aftAllCorr_icent%d", icent));

		//project 1D LSPos, LSNeg
		H1d_M_LSPos_aftAllCorr[icent]  = (TH1D*)H2d_MvsPt_LSPos_aftAllCorr[icent]  ->ProjectionY(Form("H1d_M_LSPos_aftAllCorr_icent%d", icent), LowPtBin, HigPtBin);
		H1d_M_LSNeg_aftAllCorr[icent]  = (TH1D*)H2d_MvsPt_LSNeg_aftAllCorr[icent]  ->ProjectionY(Form("H1d_M_LSNeg_aftAllCorr_icent%d", icent), LowPtBin, HigPtBin);
		H1d_M_LSPos_aftAllCorr[icent]  ->SetName(Form("H1d_M_LSPos_aftAllCorr_icent%d",  icent));
		H1d_M_LSNeg_aftAllCorr[icent]  ->SetName(Form("H1d_M_LSNeg_aftAllCorr_icent%d",  icent));

		//project 1D CKTSum
		H1d_M_CKTSum_aftAllCorr[icent] = (TH1D*)H2D_MvsPt_totCKT_Reb_aftAcc[icent] ->ProjectionY(Form("H1d_M_CKTSum_aftAllCorr_icent%d", icent), LowPtBin, HigPtBin);
		H1d_M_CKTSum_aftAllCorr[icent] ->SetName(Form("H1d_M_CKTSum_aftAllCorr_icent%d", icent));
		//in principle: does H1d_M_CKTSum_Expect == H1d_M_CKTSum_aftAllCorr ???
		H1d_M_CKTSum_Expect[icent] = (TH1D*)H2D_MvsPt_totCKTbfAcc_Reb[icent] ->ProjectionY(Form("H1d_M_CKTSum_Expect_icent%d", icent), LowPtBin, HigPtBin);
		H1d_M_CKTSum_Expect[icent] ->SetName(Form("H1d_M_CKTSum_Expect_icent%d",         icent));

		H1d_M_CKTnoVect_aftAllCorr[icent] 
			= (TH1D*)H2D_MvsPt_CKTnoVect_aftAcc[icent]->ProjectionY(Form("H1d_M_CKTnoVect_aftAllCorr_icent%d", icent), LowPtBin, HigPtBin);
		H1d_M_CKTnoVect_aftAllCorr[icent] ->SetName(Form("H1d_M_CKTnoVect_aftAllCorr_icent%d", icent));
		H1d_M_CKTw_aftAllCorr[icent] 
			= (TH1D*)H2D_MvsPt_CKTw_aftAcc[icent]    ->ProjectionY(Form("H1d_M_CKTw_aftAllCorr_icent%d",       icent), LowPtBin, HigPtBin);
		H1d_M_CKTw_aftAllCorr[icent]      ->SetName(Form("H1d_M_CKTw_aftAllCorr_icent%d", icent));
		H1d_M_CKTphi_aftAllCorr[icent] 
			= (TH1D*)H2D_MvsPt_CKTphi_aftAcc[icent]  ->ProjectionY(Form("H1d_M_CKTphi_aftAllCorr_icent%d",     icent), LowPtBin, HigPtBin);
		H1d_M_CKTphi_aftAllCorr[icent]    ->SetName(Form("H1d_M_CKTphi_aftAllCorr_icent%d", icent));
		H1d_M_CKTcc_aftAllCorr[icent] 
			= (TH1D*)H2D_MvsPt_CKTcc_aftAcc[icent]    ->ProjectionY(Form("H1d_M_CKTcc_aftAllCorr_icent%d",     icent), LowPtBin, HigPtBin);
		H1d_M_CKTcc_aftAllCorr[icent]      ->SetName(Form("H1d_M_CKTcc_aftAllCorr_icent%d", icent));
		H1d_M_CKTdy_aftAllCorr[icent] 
			= (TH1D*)H2D_MvsPt_CKTdy_aftAcc[icent]  ->ProjectionY(Form("H1d_M_CKTdy_aftAllCorr_icent%d",     icent), LowPtBin, HigPtBin);
		H1d_M_CKTdy_aftAllCorr[icent]    ->SetName(Form("H1d_M_CKTdy_aftAllCorr_icent%d", icent));
		

		//project 1D US, LSPos, LSNeg for given Pt bin
		for(int ipt=0; ipt<nPtBins4Vm; ipt++)
		{
			double iptBinLow = H2d_MvsPt_US_aftEffAcc[icent]->GetXaxis()->FindBin(ptBds4Vm[ipt]);
			double iptBinHig = H2d_MvsPt_US_aftEffAcc[icent]->GetXaxis()->FindBin(ptBds4Vm[ipt+1]);

			H1d_MinPt_US_aftAllCorr[icent][ipt]    = (TH1D*) H2d_MvsPt_US_aftEffAcc[icent]    ->ProjectionY(Form("H1d_MinPt_US_aftAllCorr_icent%d_ipt%d",    icent, ipt), iptBinLow, iptBinHig);
			H1d_MinPt_LSPos_aftAllCorr[icent][ipt] = (TH1D*) H2d_MvsPt_LSPos_aftAllCorr[icent]->ProjectionY(Form("H1d_MinPt_LSPos_aftAllCorr_icent%d_ipt%d", icent, ipt), iptBinLow, iptBinHig);
			H1d_MinPt_LSNeg_aftAllCorr[icent][ipt] = (TH1D*) H2d_MvsPt_LSNeg_aftAllCorr[icent]->ProjectionY(Form("H1d_MinPt_LSNeg_aftAllCorr_icent%d_ipt%d", icent, ipt), iptBinLow, iptBinHig);
			
			H1d_MinPt_US_aftAllCorr[icent][ipt]    ->SetName(Form("H1d_MinPt_US_aftAllCorr_icent%d_ipt%d",       icent, ipt));
			H1d_MinPt_LSPos_aftAllCorr[icent][ipt] ->SetName(Form("H1d_MinPt_LSPos_aftAllCorr_icent%d_ipt%d",    icent, ipt));
			H1d_MinPt_LSNeg_aftAllCorr[icent][ipt] ->SetName(Form("H1d_MinPt_LSNeg_aftAllCorr_icent%d_ipt%d",    icent, ipt));
		
			//project different ckt compoents into seperate bins
		}

	}//icent
	//---------------------------------------------------------------------------
	//---------------------------------------------------------------------------
	//Get 1D fully corrected "US-GeomLS-CKTSum" for final Excess Yield
	//replace GeomLS with MixUS for given mass region
	for(int icent=0; icent<nCentBins; icent++)
	{
		int item_nMassB = nMBins4Phys_ct0;
		
		if     (icent==0) item_nMassB = nMBins4Phys_ct0;
		else if(icent==1) item_nMassB = nMBins4Phys_ct1;
		else if(icent==2) item_nMassB = nMBins4Phys_ct2;
		else if(icent==3) item_nMassB = nMBins4Phys_ct3;
		
		double item_massBds[item_nMassB];

		if     (icent==0) for(int im=0; im<item_nMassB+1; im++) item_massBds[im] = mBds4Phys_ct0[im]; 
		else if(icent==1) for(int im=0; im<item_nMassB+1; im++) item_massBds[im] = mBds4Phys_ct1[im]; 
		else if(icent==2) for(int im=0; im<item_nMassB+1; im++) item_massBds[im] = mBds4Phys_ct2[im]; 
		else if(icent==3) for(int im=0; im<item_nMassB+1; im++) item_massBds[im] = mBds4Phys_ct3[im]; 
		
		//Rebin small bins to wide bins for Physics figures
		H1d_M_Reb_US_aftAllCorr[icent] 
			= (TH1D*) rebHist1d(H1d_M_US_aftAllCorr[icent],        Form("H1d_M_Reb_US_aftAllCorr_icent%d",   icent),     item_nMassB, item_massBds, "NO");
		H1d_M_Reb_MixUS_aftAllCorr[icent] 
			= (TH1D*) rebHist1d(H1d_M_MixUS_aftAllCorr[icent],     Form("H1d_M_Reb_MixUS_aftAllCorr_icent%d",icent),     item_nMassB, item_massBds, "NO");
		H1d_M_Reb_LSPos_aftAllCorr[icent] 
			= (TH1D*) rebHist1d(H1d_M_LSPos_aftAllCorr[icent],     Form("H1d_M_Reb_LSPos_aftAllCorr_icent%d",icent),     item_nMassB, item_massBds, "NO");
		H1d_M_Reb_LSNeg_aftAllCorr[icent] 
			= (TH1D*) rebHist1d(H1d_M_LSNeg_aftAllCorr[icent],     Form("H1d_M_Reb_LSNeg_aftAllCorr_icent%d",icent),     item_nMassB, item_massBds, "NO");
		//Rebin CKTSum
		H1d_M_Reb_CKTSum_aftAllCorr[icent]
			= (TH1D*) rebHist1d(H1d_M_CKTSum_aftAllCorr[icent],    Form("H1d_M_Reb_CKTSum_aftAllCorr_icent%d",   icent), item_nMassB, item_massBds, "NO"); 
		H1d_M_Reb_CKTnoVect_aftAllCorr[icent]
			= (TH1D*) rebHist1d(H1d_M_CKTnoVect_aftAllCorr[icent], Form("H1d_M_Reb_CKTnoVect_aftAllCorr_icent%d",icent), item_nMassB, item_massBds, "NO"); 
		H1d_M_Reb_CKTw_aftAllCorr[icent]
			= (TH1D*) rebHist1d(H1d_M_CKTw_aftAllCorr[icent],      Form("H1d_M_Reb_CKTw_aftAllCorr_icent%d",     icent), item_nMassB, item_massBds, "NO"); 
		H1d_M_Reb_CKTphi_aftAllCorr[icent]
			= (TH1D*) rebHist1d(H1d_M_CKTphi_aftAllCorr[icent],    Form("H1d_M_Reb_CKTphi_aftAllCorr_icent%d",   icent), item_nMassB, item_massBds, "NO"); 
		H1d_M_Reb_CKTcc_aftAllCorr[icent]
			= (TH1D*) rebHist1d(H1d_M_CKTcc_aftAllCorr[icent],     Form("H1d_M_Reb_CKTcc_aftAllCorr_icent%d",    icent), item_nMassB, item_massBds, "NO"); 
		H1d_M_Reb_CKTdy_aftAllCorr[icent]
			= (TH1D*) rebHist1d(H1d_M_CKTdy_aftAllCorr[icent],     Form("H1d_M_Reb_CKTdy_aftAllCorr_icent%d",    icent), item_nMassB, item_massBds, "NO"); 
		//---------------------------------------------------------------------------

		//	//---------------------------------------------------------------------------
		//	//Rebin small bins to wide bins for Physics figures
		//	H1d_M_Reb_US_aftAllCorr[icent] 
		//		= (TH1D*) rebHist1d(H1d_M_US_aftAllCorr[icent],        Form("H1d_M_Reb_US_aftAllCorr_icent%d",   icent),     nMBins4Phys_ct0, mBds4Phys_ct0, "NO");
		//	H1d_M_Reb_MixUS_aftAllCorr[icent] 
		//		= (TH1D*) rebHist1d(H1d_M_MixUS_aftAllCorr[icent],     Form("H1d_M_Reb_MixUS_aftAllCorr_icent%d",icent),     nMBins4Phys_ct0, mBds4Phys_ct0, "NO");
		//	H1d_M_Reb_LSPos_aftAllCorr[icent] 
		//		= (TH1D*) rebHist1d(H1d_M_LSPos_aftAllCorr[icent],     Form("H1d_M_Reb_LSPos_aftAllCorr_icent%d",icent),     nMBins4Phys_ct0, mBds4Phys_ct0, "NO");
		//	H1d_M_Reb_LSNeg_aftAllCorr[icent] 
		//		= (TH1D*) rebHist1d(H1d_M_LSNeg_aftAllCorr[icent],     Form("H1d_M_Reb_LSNeg_aftAllCorr_icent%d",icent),     nMBins4Phys_ct0, mBds4Phys_ct0, "NO");
		//	//Rebin CKTSum
		//	H1d_M_Reb_CKTSum_aftAllCorr[icent]
		//		= (TH1D*) rebHist1d(H1d_M_CKTSum_aftAllCorr[icent],    Form("H1d_M_Reb_CKTSum_aftAllCorr_icent%d",   icent), nMBins4Phys_ct0, mBds4Phys_ct0, "NO"); 
		//	H1d_M_Reb_CKTnoVect_aftAllCorr[icent]
		//		= (TH1D*) rebHist1d(H1d_M_CKTnoVect_aftAllCorr[icent], Form("H1d_M_Reb_CKTnoVect_aftAllCorr_icent%d",icent), nMBins4Phys_ct0, mBds4Phys_ct0, "NO"); 
		//	H1d_M_Reb_CKTw_aftAllCorr[icent]
		//		= (TH1D*) rebHist1d(H1d_M_CKTw_aftAllCorr[icent],      Form("H1d_M_Reb_CKTw_aftAllCorr_icent%d",     icent), nMBins4Phys_ct0, mBds4Phys_ct0, "NO"); 
		//	H1d_M_Reb_CKTphi_aftAllCorr[icent]
		//		= (TH1D*) rebHist1d(H1d_M_CKTphi_aftAllCorr[icent],    Form("H1d_M_Reb_CKTphi_aftAllCorr_icent%d",   icent), nMBins4Phys_ct0, mBds4Phys_ct0, "NO"); 
		//	H1d_M_Reb_CKTcc_aftAllCorr[icent]
		//		= (TH1D*) rebHist1d(H1d_M_CKTcc_aftAllCorr[icent],     Form("H1d_M_Reb_CKTcc_aftAllCorr_icent%d",    icent), nMBins4Phys_ct0, mBds4Phys_ct0, "NO"); 
		//	H1d_M_Reb_CKTdy_aftAllCorr[icent]
		//		= (TH1D*) rebHist1d(H1d_M_CKTdy_aftAllCorr[icent],     Form("H1d_M_Reb_CKTdy_aftAllCorr_icent%d",    icent), nMBins4Phys_ct0, mBds4Phys_ct0, "NO"); 
		//	//---------------------------------------------------------------------------
		
		//get GeomLS
		H1d_M_Reb_GmLS_aftAllCorr[icent] 
			= (TH1D*) getGeomMeanHist1d( H1d_M_Reb_LSPos_aftAllCorr[icent], H1d_M_Reb_LSNeg_aftAllCorr[icent], Form("H1d_M_Reb_GmLS_aftAllCorr_icent%d", icent) );
		//---------------------------------------------------------------------------
		//Scale MixUS to GeomLS for given mass region
		double scf = ScaleFactorH1toH2(scaleM_Low[icent], scaleM_Hig[icent], H1d_M_Reb_MixUS_aftAllCorr[icent], H1d_M_Reb_GmLS_aftAllCorr[icent]);
		H1d_M_Reb_MixUS_aftAllCorr[icent] ->Scale(scf);
		
		H1d_M_Reb_CombBG_aftAllCorr[icent] = (TH1D*)CombGmLSMixUS( Form("H1d_M_Reb_ComBG_aftAllCorr_icent%d", icent), H1d_M_Reb_GmLS_aftAllCorr[icent], H1d_M_Reb_MixUS_aftAllCorr[icent], RepM_Low[icent], RepM_Hig[icent] );

		//---------------------------------------------------------------------------
		//US-GeomLS
		H1d_M_Reb_Sig_aftAllCorr[icent]  = (TH1D*)H1d_M_Reb_US_aftAllCorr[icent]->Clone(Form("H1d_M_Reb_Sig_aftAllCorr_icent%d", icent));
		H1d_M_Reb_Sig_aftAllCorr[icent]  ->Add(H1d_M_Reb_GmLS_aftAllCorr[icent],   -1);
		//H1d_M_Reb_Sig_aftAllCorr[icent]->Add(H1d_M_Reb_CombBG_aftAllCorr[icent], -1); //use GeomLS+MixUS as background, replace 1.5-2.5 GeV/c2

		//---------------------------------------------------------------------------
		H1d_M_Reb_US_aftAllCorr[icent]    ->SetName(Form("H1d_M_Reb_US_aftAllCorr_icent%d",     icent));
		H1d_M_Reb_GmLS_aftAllCorr[icent]  ->SetName(Form("H1d_M_Reb_GmLS_aftAllCorr_icent%d",   icent));
		H1d_M_Reb_Sig_aftAllCorr[icent]   ->SetName(Form("H1d_M_Reb_Sig_aftAllCorr_icent%d",    icent));
		H1d_M_Reb_CKTSum_aftAllCorr[icent]->SetName(Form("H1d_M_Reb_CKTSum_aftAllCorr_icent%d", icent));

		//---------------------------------------------------------------------------
		//(US-GeomLS)/nMB - CKTSum use Bin for Physics
		H1d_M_Reb_ExY_aftAllCorr[icent]   = (TH1D*) H1d_M_Reb_Sig_aftAllCorr[icent]->Clone(Form("H1d_M_Reb_ExY_aftAllCorr_icent%d", icent));
		H1d_M_Reb_ExY_aftAllCorr[icent]   ->Add(H1d_M_Reb_CKTSum_aftAllCorr[icent], -1);
		H1d_M_Reb_ExY_aftAllCorr[icent]   ->SetName(Form("H1d_M_Reb_ExY_aftAllCorr_icent%d",    icent));
	}//icent
	//---------------------------------------------------------------------------

	//---------------------------------------------------------------------------
	//Get fully correct 1D "US-GeomLS" for Vector Meson: Omega, Phi, J/psi
	for(int icent=0; icent<nCentBins; icent++)
	{
		//---------------------------------------------------------------------------
		//Rebin small bins to wide bins for Physics figures
		H1d_M_4Vm_US_aftAllCorr[icent] 
			= (TH1D*) rebHist1d(H1d_M_US_aftAllCorr[icent],     Form("H1d_M_4Vm_US_aftAllCorr_icent%d",   icent), nMBins4Vm_ct0, mBds4Vm_ct0, "NO");
		H1d_M_4Vm_LSPos_aftAllCorr[icent] 
			= (TH1D*) rebHist1d(H1d_M_LSPos_aftAllCorr[icent],  Form("H1d_M_4Vm_LSPos_aftAllCorr_icent%d",icent), nMBins4Vm_ct0, mBds4Vm_ct0, "NO");
		H1d_M_4Vm_LSNeg_aftAllCorr[icent] 
			= (TH1D*) rebHist1d(H1d_M_LSNeg_aftAllCorr[icent],  Form("H1d_M_4Vm_LSNeg_aftAllCorr_icent%d",icent), nMBins4Vm_ct0, mBds4Vm_ct0, "NO");
		//---------------------------------------------------------------------------
		//get GeomLS
		H1d_M_4Vm_GmLS_aftAllCorr[icent] 
			= (TH1D*) getGeomMeanHist1d( H1d_M_4Vm_LSPos_aftAllCorr[icent], H1d_M_4Vm_LSNeg_aftAllCorr[icent], Form("H1d_M_4Vm_GmLS_aftAllCorr_icent%d", icent) );
		//---------------------------------------------------------------------------
		//US-GeomLS
		H1d_M_4Vm_Sig_aftAllCorr[icent]  = (TH1D*)H1d_M_4Vm_US_aftAllCorr[icent]->Clone(Form("H1d_M_4Vm_Sig_aftAllCorr_icent%d", icent));
		H1d_M_4Vm_Sig_aftAllCorr[icent]  ->Add(H1d_M_4Vm_GmLS_aftAllCorr[icent],   -1);

		//---------------------------------------------------------------------------
		H1d_M_4Vm_US_aftAllCorr[icent]    ->SetName(Form("H1d_M_4Vm_US_aftAllCorr_icent%d",   icent));
		H1d_M_4Vm_GmLS_aftAllCorr[icent]  ->SetName(Form("H1d_M_4Vm_GmLS_aftAllCorr_icent%d", icent));
		H1d_M_4Vm_Sig_aftAllCorr[icent]   ->SetName(Form("H1d_M_4Vm_Sig_aftAllCorr_icent%d",  icent));

		//-----------------------------------------------------------------------------------------------------------
		//do above thre procedures for each pt bins
		//-----------------------------------------------------------------------------------------------------------
		for(int ipt=0; ipt<nPtBins4Vm; ipt++)
		{
			//Rebin small bins to wide bins for Physics figures
			H1d_M4VminPt_US_aftAllCorr[icent][ipt] 
				= (TH1D*) rebHist1d(H1d_MinPt_US_aftAllCorr[icent][ipt],    Form("H1d_M4VminPt_US_aftAllCorr_icent%d_ipt%d",    icent, ipt), nMBins4Vm_ct0, mBds4Vm_ct0, "NO");
			H1d_M4VminPt_LSPos_aftAllCorr[icent][ipt] 
				= (TH1D*) rebHist1d(H1d_MinPt_LSPos_aftAllCorr[icent][ipt], Form("H1d_M4VminPt_LSPos_aftAllCorr_icent%d_ipt%d", icent, ipt), nMBins4Vm_ct0, mBds4Vm_ct0, "NO");
			H1d_M4VminPt_LSNeg_aftAllCorr[icent][ipt] 
				= (TH1D*) rebHist1d(H1d_MinPt_LSNeg_aftAllCorr[icent][ipt], Form("H1d_M4VminPt_LSNeg_aftAllCorr_icent%d_ipt%d", icent, ipt), nMBins4Vm_ct0, mBds4Vm_ct0, "NO");

			//---------------------------------------------------------------------------
			//get GeomLS
			H1d_M4VminPt_GmLS_aftAllCorr[icent][ipt] 
				= (TH1D*) getGeomMeanHist1d( H1d_M4VminPt_LSPos_aftAllCorr[icent][ipt], H1d_M4VminPt_LSNeg_aftAllCorr[icent][ipt], Form("H1d_M4VminPt_GmLS_aftAllCorr_icent%d_ipt%d", icent, ipt) );
			//---------------------------------------------------------------------------
			//US-GeomLS
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]  = (TH1D*)H1d_M4VminPt_US_aftAllCorr[icent][ipt]->Clone(Form("H1d_M4VminPt_Sig_aftAllCorr_icent%d_ipt%d", icent, ipt));
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]  ->Add(H1d_M4VminPt_GmLS_aftAllCorr[icent][ipt],   -1);

			//---------------------------------------------------------------------------
			H1d_M4VminPt_US_aftAllCorr[icent][ipt]    ->SetName(Form("H1d_M4VminPt_US_aftAllCorr_icent%d_ipt%d",   icent,ipt));
			H1d_M4VminPt_GmLS_aftAllCorr[icent][ipt]  ->SetName(Form("H1d_M4VminPt_GmLS_aftAllCorr_icent%d_ipt%d", icent,ipt));
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]   ->SetName(Form("H1d_M4VminPt_Sig_aftAllCorr_icent%d_ipt%d",  icent,ipt));
		}//ipt
		//-----------------------------------------------------------------------------------------------------------
		//-----------------------------------------------------------------------------------------------------------
	}//icent

	TFile* outtem = new TFile("outtem.root", "recreate");
	outtem->cd();
	for(int icent=0; icent<nCentBins; icent++)
	{
		H2d_MvsPt_PairEff[icent]     ->Write();
		H2d_MvsPt_PairAcc[icent]     ->Write();
		H2d_MvsPt_PairEffAcc[icent]     ->Write();
		H2d_MvsPt_Reb_US[icent]         ->Write();
		H2d_MvsPt_Reb_LSPos[icent]      ->Write();
		H2d_MvsPt_Reb_LSNeg[icent]      ->Write();
		H2d_MvsPt_LSPos_aftPSAC[icent]  ->Write();
		H2d_MvsPt_LSNeg_aftPSAC[icent]  ->Write();

		H2d_MvsPt_US_aftEffAcc[icent]     ->Write();
		H2d_MvsPt_LSPos_aftAllCorr[icent] ->Write();
		H2d_MvsPt_LSNeg_aftAllCorr[icent] ->Write();

		H1d_M_US_aftAllCorr[icent]      ->Write();
		H1d_M_CKTSum_aftAllCorr[icent]  ->Write();
		H1d_M_CKTSum_Expect[icent]      ->Write();
		H1d_M_Reb_US_aftAllCorr[icent]  ->Write();
		H1d_M_Reb_GmLS_aftAllCorr[icent]->Write();
		H1d_M_Reb_Sig_aftAllCorr[icent] ->Write();
		H1d_M_Reb_CKTSum_aftAllCorr[icent] ->Write();
		H1d_M_Reb_ExY_aftAllCorr[icent] ->Write();
	}
	outtem->Close();
}
//------------------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------------------
void readTheory()
{
	TFile *infileTheory = new TFile("../RemovePhoE/Signals/prepareTheoryCurve/FullMassTheoryCurve/output/out_TheoryCurve_Ex_diffCent.root", "read");
	cout<<"infileTheory: "<<infileTheory->GetName()<<endl;

	for(int it=0; it<nTs; it++)
	{
		for(int icent=0; icent<nCentBins; icent++)
		{//g_Ex_MediumRho_it0_icent
			g_Ex_MediumRho_inCent[it][icent]   = (TGraph*)infileTheory->Get(Form("g_Ex_MediumRho_it%d_icent%d",   it, icent));
			g_Ex_QGPEmission_inCent[it][icent] = (TGraph*)infileTheory->Get(Form("g_Ex_QGPEmission_it%d_icent%d", it, icent));
			g_Ex_TheorySum_inCent[it][icent]   = (TGraph*)infileTheory->Get(Form("g_Ex_Sum_it%d_icent%d",         it, icent));
		}//icent
	}//it

}
//------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------
void drawVectorMesons()
{
	double Ratio_Omega2Phi[nNPartBins];
	double RatioErr_Omega2Phi[nNPartBins];

	TCanvas *c1 = new TCanvas("c1");
	gStyle->SetOptStat(0);
	gStyle->SetLineWidth(2);
	gPad->SetTopMargin(0.08);
	gPad->SetBottomMargin(0.11);
	gPad->SetLeftMargin(0.11);
	gPad->SetRightMargin(0.05);
	//---------------------------------------------------------------------------
	//---------------------------------------------------------------------------
	TLatex tl;
	tl.SetTextSize(0.06);
	tl.SetNDC();
	c1->SetLogy(0);

	//---------------------------------------------------------------------------
	//plot Phi mass region
	//---------------------------------------------------------------------------
	for(int icent=0; icent<nCentBins; icent++)
	{
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->Scale(1, "width");
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetMarkerStyle(20);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetMarkerColor(1);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetLineColor(1);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetLineColor(1);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetYTitle("dN/dM");
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetTitle("Fully Corrected");
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->GetYaxis()->SetTitleSize(0.06);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->GetYaxis()->SetTitleOffset(0.80);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->GetXaxis()->SetTitleSize(0.05);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->GetXaxis()->SetTitleOffset(0.85);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetAxisRange(0.85, 1.15, "x");
		//H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetAxisRange(0.65, 0.85, "x");
		//if(icent==1)H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetMinimum(1.e-5);
		int itemBin_phi    = H1d_M_4Vm_Sig_aftAllCorr[icent]->FindBin(1.015);
		double itemContent = H1d_M_4Vm_Sig_aftAllCorr[icent]->GetBinContent(itemBin_phi);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetMaximum(itemContent*1.2);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->Draw("pe");
	
		const double LM_Phi = 0.84; double HM_Phi = 1.200;

		TF1 *fitf_Phi = new TF1("fitf_Phi", BreitWigner_BG_formula, LM_Phi, HM_Phi, 5);
		fitf_Phi->SetParNames("NS",   "#Gamma", "mass",  "a",    "b" );
		fitf_Phi->SetParameters(1.e-5,  0.01,    1.020,   1.e5,    5. );
		fitf_Phi->SetParLimits( 0, 0.0,   1.e2 );
		fitf_Phi->SetParLimits( 1, 0.0,   0.06 );
		fitf_Phi->SetParLimits( 2, 1.00,  1.05 );
		fitf_Phi->SetParLimits( 3, -1.e6, 1.e6 );
		fitf_Phi->SetParLimits( 4, -50, 50 );
		fitf_Phi->SetLineColor(1);
		fitf_Phi->SetLineWidth(2);
		fitf_Phi->SetNpx(20000);

		H1d_M_4Vm_Sig_aftAllCorr[icent]->Fit(fitf_Phi, "E", "", 1.04, 1.2 );
		if(icent==1) H1d_M_4Vm_Sig_aftAllCorr[icent]->Fit(fitf_Phi, "E", "", 0.85, 0.92 );
		
		H1d_M_4Vm_Sig_aftAllCorr[icent]->Fit(fitf_Phi, "ER", "", LM_Phi, HM_Phi+0.005 );
		
		double chi2 = fitf_Phi ->GetChisquare();
		double ndf  = fitf_Phi ->GetNDF();
		double Chi2overNDF = chi2/ndf;

		double parm[5]; double parmerr[5];

		for(int iparm=0; iparm<5; iparm++)
		{
			parm[iparm]    =  fitf_Phi->GetParameter( iparm );
			parmerr[iparm] =  fitf_Phi->GetParError(  iparm );
		}

		TF1 *fsg_Phi = new TF1("fsg_Phi", BreitWigner_formula,  LM_Phi, HM_Phi+0.005, 3);
		TF1* fbg_Phi = new TF1("fbg",     bg_formula,           LM_Phi, HM_Phi+0.005, 2);
		fsg_Phi->SetNpx(20000);

		for(int isg=0; isg<3; isg++) 
		{
			cout<<parm[isg]<<endl;
			fsg_Phi->SetParameter(isg,  parm[isg]   );
			fsg_Phi->SetParError( isg,  parmerr[isg]);
		}
		for(int ibg=0; ibg<2; ibg++) 
		{
			fbg_Phi->SetParameter(ibg, parm[ibg+3]   );
			fbg_Phi->SetParError( ibg, parmerr[ibg+3]);
		}

		fsg_Phi->SetLineColor(4);
		fsg_Phi->SetLineWidth(2);
		fsg_Phi->Draw("lsame");
		fbg_Phi->SetLineStyle(5);
		fbg_Phi->SetLineColor(2);
		fbg_Phi->SetLineWidth(2);
		fbg_Phi->Draw("lsame");

		TLegend* leg_Phi = new TLegend(0.15, 0.60, 0.45, 0.80);
		leg_Phi->SetBorderSize(0);
		leg_Phi->SetFillColor(0);
		leg_Phi->SetTextSize(0.035);
		leg_Phi->AddEntry( H1d_M_4Vm_Sig_aftAllCorr[icent], "Data",                 "lp" );
		leg_Phi->AddEntry( fitf_Phi,         "Breit-Wigner+Expo.",    "l"  );
		leg_Phi->AddEntry( fsg_Phi,          "Breit-Wigner",         "l"  );
		leg_Phi->AddEntry( fbg_Phi,          "Expo.",                 "l"  );
		leg_Phi->Draw();

		tl.SetTextSize(0.050);
		tl.SetTextColor(1);
		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
		tl.SetTextSize(0.040);
		tl.SetTextColor(1);
		
		tl.DrawLatex(0.63, 0.82, Form("#chi^{2}/ndf = %.2f",      Chi2overNDF) );
		tl.DrawLatex(0.63, 0.76, Form("M = %.1f #pm %.1f MeV/c^{2}",      parm[2]*1.e3, parmerr[2]*1.e3) );
		tl.DrawLatex(0.63, 0.70, Form("#Gamma = %.1f #pm %.1f MeV/c^{2}", parm[1]*1.e3, parmerr[1]*1.e3) );

		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Phi_aftAllCorr_icent%d.png", icent) );
		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Phi_aftAllCorr_icent%d.pdf", icent) );
		//---------------------------------------------------------------------------
		//---------------------------------------------------------------------------
		//Draw Vector mass in different pt bins

		for(int ipt=0; ipt<nPtBins4Vm; ipt++)
		{
			//---------------------------------------------------------------------------
			//plot Phi mass region
			//---------------------------------------------------------------------------
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->Scale(1, "width");
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetMarkerStyle(20);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetMarkerColor(1);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetLineColor(1);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetLineColor(1);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetYTitle("dN/dM");
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetTitle("Fully Corrected");
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->GetYaxis()->SetTitleSize(0.06);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->GetYaxis()->SetTitleOffset(0.80);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->GetXaxis()->SetTitleSize(0.05);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->GetXaxis()->SetTitleOffset(0.85);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetAxisRange(0.85, 1.15, "x");
			
			int    itemBin_phi = H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->FindBin(1.015);
			double itemContent = H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->GetBinContent(itemBin_phi);
			//H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetMaximum(itemContent*1.2);
			
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->Draw("pe");
			//H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->Draw("histsame");
			//---------------------------------------------------------------------------

			//---------------------------------------------------------------------------
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Phi, "E", "", 1.04, 1.2 );
			if(icent==1)H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Phi, "E", "", 0.85, 0.92 );
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Phi, "E", "", 1.0, 1.05 );
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Phi, "E", "", 1.0, 1.20 );
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Phi, "ER", "", LM_Phi, HM_Phi+0.005 );
			//---------------------------------------------------------------------------

			double chi2 = fitf_Phi ->GetChisquare();
			double ndf  = fitf_Phi ->GetNDF();
			double Chi2overNDF = chi2/ndf;

			double parm[5]; double parmerr[5];

			for(int iparm=0; iparm<5; iparm++)
			{
				parm[iparm]    =  fitf_Phi->GetParameter( iparm );
				parmerr[iparm] =  fitf_Phi->GetParError(  iparm );
			}

			TF1 *fsg_Phi = new TF1("fsg_Phi", BreitWigner_formula,  LM_Phi, HM_Phi+0.005, 3);
			TF1* fbg_Phi = new TF1("fbg",     bg_formula,           LM_Phi, HM_Phi+0.005, 2);
			fsg_Phi->SetNpx(20000);

			for(int isg=0; isg<3; isg++) 
			{
				cout<<parm[isg]<<endl;
				fsg_Phi->SetParameter(isg,  parm[isg]   );
				fsg_Phi->SetParError( isg,  parmerr[isg]);
			}
			for(int ibg=0; ibg<2; ibg++) 
			{
				fbg_Phi->SetParameter(ibg, parm[ibg+3]   );
				fbg_Phi->SetParError( ibg, parmerr[ibg+3]);
			}

			fsg_Phi->SetLineColor(4);
			fsg_Phi->SetLineWidth(2);
			fsg_Phi->Draw("lsame");
			fbg_Phi->SetLineStyle(5);
			fbg_Phi->SetLineColor(2);
			fbg_Phi->SetLineWidth(2);
			fbg_Phi->Draw("lsame");

			TLegend* leg_Phi = new TLegend(0.15, 0.60, 0.45, 0.80);
			leg_Phi->SetBorderSize(0);
			leg_Phi->SetFillColor(0);
			leg_Phi->SetTextSize(0.035);
			leg_Phi->AddEntry( H1d_M4VminPt_Sig_aftAllCorr[icent][ipt], Form("Data: %.1f<p_{T}<%.1f GeV/c", ptBds4Vm[ipt], ptBds4Vm[ipt+1]),  "lp" );
			leg_Phi->AddEntry( fitf_Phi,         "Breit-Wigner+Expo.",   "l"  );
			leg_Phi->AddEntry( fsg_Phi,          "Breit-Wigner",         "l"  );
			leg_Phi->AddEntry( fbg_Phi,          "Expo.",                "l"  );
			leg_Phi->Draw();

			tl.SetTextSize(0.050);
			tl.SetTextColor(1);
			tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
			tl.SetTextSize(0.040);
			tl.SetTextColor(1);

			tl.DrawLatex(0.63, 0.82, Form("#chi^{2}/ndf = %.2f",      Chi2overNDF) );
			tl.DrawLatex(0.63, 0.76, Form("M = %.1f #pm %.1f MeV/c^{2}",      parm[2]*1.e3, parmerr[2]*1.e3) );
			tl.DrawLatex(0.63, 0.70, Form("#Gamma = %.1f #pm %.1f MeV/c^{2}", parm[1]*1.e3, parmerr[1]*1.e3) );

			c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Phi_aftAllCorr_icent%d_ipt%d.png", icent, ipt) );
		}//ipt
	}//icent

	//---------------------------------------------------------------------------
	//plot Omega mass region
	//---------------------------------------------------------------------------
	const double LM_Omega = 0.62; double HM_Omega = 0.86;
	for(int icent=0; icent<nCentBins; icent++)
	{
		//H1d_M_4Vm_Sig_aftAllCorr[icent]    ->Scale(1, "width");
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetMarkerStyle(20);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetMarkerColor(1);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetLineColor(1);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetLineWidth(2);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetYTitle("dN/dM");
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetTitle("Fully Corrected");
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->GetYaxis()->SetTitleSize(0.06);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->GetYaxis()->SetTitleOffset(0.80);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->GetXaxis()->SetTitleSize(0.05);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->GetXaxis()->SetTitleOffset(0.85);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetAxisRange(LM_Omega, HM_Omega, "x");
		//if(icent==1)H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetMinimum(1.e-5);

		int    itemBin_Omega  = H1d_M_4Vm_Sig_aftAllCorr[icent]->FindBin(0.78);
		double itemContent    = H1d_M_4Vm_Sig_aftAllCorr[icent]->GetBinContent(itemBin_Omega);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetMaximum(itemContent*1.3);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->Draw("pe");
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->Draw("histsame");
	
		TF1 *fitf_Omega = new TF1("fitf_Omega", BreitWigner_BG_formula, LM_Omega, HM_Omega, 5);
		fitf_Omega->SetParNames("NS",   "#Gamma", "mass",  "a",    "b" );
		fitf_Omega->SetParameters(1.e-5,  0.02,    0.780,   1.e5,    5. );
		fitf_Omega->SetParLimits( 0, 0.0,   1.e2 );
		fitf_Omega->SetParLimits( 1, 0.0,   0.06 );
		fitf_Omega->SetParLimits( 2, 0.75,  0.80 );
		fitf_Omega->SetParLimits( 3, -1.e6, 1.e6 );
		fitf_Omega->SetParLimits( 4, -50,   50   );
		fitf_Omega->SetLineColor(1);
		fitf_Omega->SetLineWidth(2);
		fitf_Omega->SetNpx(20000);

		H1d_M_4Vm_Sig_aftAllCorr[icent]->Fit(fitf_Omega, "E",  "", 0.80, 0.88 );
		H1d_M_4Vm_Sig_aftAllCorr[icent]->Fit(fitf_Omega, "E",  "", 0.77, 0.88 );
		H1d_M_4Vm_Sig_aftAllCorr[icent]->Fit(fitf_Omega, "ER", "", LM_Omega, HM_Omega+0.03 );

		double chi2 = fitf_Omega ->GetChisquare();
		double ndf  = fitf_Omega ->GetNDF();
		double Chi2overNDF = chi2/ndf;

		double parm[5]; double parmerr[5];

		for(int iparm=0; iparm<5; iparm++)
		{
			parm[iparm]    =  fitf_Omega->GetParameter( iparm );
			parmerr[iparm] =  fitf_Omega->GetParError(  iparm );
		}

		TF1 *fsg_Omega = new TF1("fsg_Omega", BreitWigner_formula,  LM_Omega, HM_Omega+0.005, 3);
		TF1* fbg_Omega = new TF1("fbg",       bg_formula,           LM_Omega, HM_Omega+0.005, 2);
		fsg_Omega->SetNpx(20000);

		for(int isg=0; isg<3; isg++) 
		{
			cout<<parm[isg]<<endl;
			fsg_Omega->SetParameter(isg,  parm[isg]   );
			fsg_Omega->SetParError( isg,  parmerr[isg]);
		}
		for(int ibg=0; ibg<2; ibg++) 
		{
			fbg_Omega->SetParameter(ibg, parm[ibg+3]   );
			fbg_Omega->SetParError( ibg, parmerr[ibg+3]);
		}

		fsg_Omega->SetLineColor(4);
		fsg_Omega->SetLineWidth(2);
		fsg_Omega->Draw("lsame");
		fbg_Omega->SetLineStyle(5);
		fbg_Omega->SetLineColor(2);
		fbg_Omega->SetLineWidth(2);
		fbg_Omega->Draw("lsame");

		TLegend* leg_Omega = new TLegend(0.15, 0.60, 0.45, 0.80);
		leg_Omega->SetBorderSize(0);
		leg_Omega->SetFillColor(0);
		leg_Omega->SetTextSize(0.035);
		leg_Omega->AddEntry( H1d_M_4Vm_Sig_aftAllCorr[icent], "Data",                 "lp" );
		leg_Omega->AddEntry( fitf_Omega,                      "Breit-Wigner+Expo.",   "l"  );
		leg_Omega->AddEntry( fsg_Omega,                       "Breit-Wigner",         "l"  );
		leg_Omega->AddEntry( fbg_Omega,                       "Expo.",                "l"  );
		leg_Omega->Draw();

		tl.SetTextSize(0.050);
		tl.SetTextColor(1);
		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
		tl.SetTextSize(0.040);
		tl.SetTextColor(1);
		
		tl.DrawLatex(0.65, 0.82, Form("#chi^{2}/ndf = %.2f",              Chi2overNDF) );
		tl.DrawLatex(0.65, 0.76, Form("M = %.1f #pm %.1f MeV/c^{2}",      parm[2]*1.e3, parmerr[2]*1.e3) );
		tl.DrawLatex(0.65, 0.70, Form("#Gamma = %.1f #pm %.1f MeV/c^{2}", parm[1]*1.e3, parmerr[1]*1.e3) );

		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Omega_aftAllCorr_icent%d.png", icent) );
		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Omega_aftAllCorr_icent%d.pdf", icent) );
		//---------------------------------------------------------------------------
		

		//---------------------------------------------------------------------------
		//Draw Omega and phi Signals on same page
		//---------------------------------------------------------------------------
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetAxisRange(0.62, 1.18, "x");
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->GetFunction("fitf_Omega")->SetBit(TF1::kNotDraw);
		
		int    itemBin_OmegaPhi    = H1d_M_4Vm_Sig_aftAllCorr[icent]->FindBin(1.015);
		double itemCont= H1d_M_4Vm_Sig_aftAllCorr[icent]->GetBinContent(itemBin_OmegaPhi);
		
		if(icent==3) H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetMaximum(itemCont*2.2);
		else         H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetMaximum(itemCont*1.5);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->Draw("pe");
		//H1d_M_4Vm_Sig_aftAllCorr[icent]    ->Draw("histsame");
		tl.SetTextSize(0.050);
		tl.SetTextColor(1);
		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_OmegaPhi_aftAllCorr_icent%d.png", icent) );
	
		//---------------------------------------------------------------------------
		//Draw Vector mass in different pt bins
		//---------------------------------------------------------------------------
		for(int ipt=0; ipt<nPtBins4Vm; ipt++)
		{
			//---------------------------------------------------------------------------
			//plot Omega mass region
			//---------------------------------------------------------------------------
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetMarkerStyle(20);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetMarkerColor(1);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetLineColor(1);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetLineWidth(2);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetYTitle("dN/dM");
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetTitle("Fully Corrected");
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->GetYaxis()->SetTitleSize(0.06);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->GetYaxis()->SetTitleOffset(0.80);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->GetXaxis()->SetTitleSize(0.05);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->GetXaxis()->SetTitleOffset(0.85);
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetAxisRange(LM_Omega, HM_Omega, "x");

			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->Draw("pe");
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->Draw("histsame");
			//---------------------------------------------------------------------------

			//---------------------------------------------------------------------------
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Omega, "E", "", 0.80, 0.88 );
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Omega, "E", "", 0.78, 0.88 );
			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Omega, "ER", "", LM_Omega, HM_Omega+0.03 );
			//			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Omega, "E", "", 1.04, 1.2 );
			//			if(icent==1)H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Omega, "E", "", 0.85, 0.92 );
			//			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Omega, "E", "", 1.0, 1.05 );
			//			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Omega, "E", "", 1.0, 1.20 );
			//			H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Omega, "ER", "", LM_Omega, HM_Omega+0.005 );
			//			//---------------------------------------------------------------------------
			//
			double chi2 = fitf_Omega ->GetChisquare();
			double ndf  = fitf_Omega ->GetNDF();
			double Chi2overNDF = chi2/ndf;

			double parm[5]; double parmerr[5];

			for(int iparm=0; iparm<5; iparm++)
			{
				parm[iparm]    =  fitf_Omega->GetParameter( iparm );
				parmerr[iparm] =  fitf_Omega->GetParError(  iparm );
			}

			TF1 *fsg_Omega = new TF1("fsg_Omega", BreitWigner_formula,  LM_Omega, HM_Omega+0.005, 3);
			TF1* fbg_Omega = new TF1("fbg",       bg_formula,           LM_Omega, HM_Omega+0.005, 2);
			fsg_Omega->SetNpx(20000);

			for(int isg=0; isg<3; isg++) 
			{
				cout<<parm[isg]<<endl;
				fsg_Omega->SetParameter(isg,  parm[isg]   );
				fsg_Omega->SetParError( isg,  parmerr[isg]);
			}
			for(int ibg=0; ibg<2; ibg++) 
			{
				fbg_Omega->SetParameter(ibg, parm[ibg+3]   );
				fbg_Omega->SetParError( ibg, parmerr[ibg+3]);
			}

			fsg_Omega->SetLineColor(4);
			fsg_Omega->SetLineWidth(2);
			fsg_Omega->Draw("lsame");
			fbg_Omega->SetLineStyle(5);
			fbg_Omega->SetLineColor(2);
			fbg_Omega->SetLineWidth(2);
			fbg_Omega->Draw("lsame");

			TLegend* leg_Omega = new TLegend(0.15, 0.60, 0.45, 0.80);
			leg_Omega->SetBorderSize(0);
			leg_Omega->SetFillColor(0);
			leg_Omega->SetTextSize(0.035);
			leg_Omega->AddEntry( H1d_M4VminPt_Sig_aftAllCorr[icent][ipt], Form("Data: %.1f<p_{T}<%.1f GeV/c", ptBds4Vm[ipt], ptBds4Vm[ipt+1]), "lp" );
			leg_Omega->AddEntry( fitf_Omega,     "Breit-Wigner+Expo.",    "l" );
			leg_Omega->AddEntry( fsg_Omega,      "Breit-Wigner",          "l" );
			leg_Omega->AddEntry( fbg_Omega,      "Expo.",                 "l" );
			leg_Omega->Draw();

			tl.SetTextSize(0.050);
			tl.SetTextColor(1);
			tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
			tl.SetTextSize(0.040);
			tl.SetTextColor(1);

			tl.DrawLatex(0.63, 0.82, Form("#chi^{2}/ndf = %.2f",              Chi2overNDF) );
			tl.DrawLatex(0.63, 0.76, Form("M = %.1f #pm %.1f MeV/c^{2}",      parm[2]*1.e3, parmerr[2]*1.e3) );
			tl.DrawLatex(0.63, 0.70, Form("#Gamma = %.1f #pm %.1f MeV/c^{2}", parm[1]*1.e3, parmerr[1]*1.e3) );

			c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Omega_aftAllCorr_icent%d_ipt%d.png", icent, ipt) );
		}//ipt
	}//icent
	
	//---------------------------------------------------------------------------------------------------------------
	//Fit the Phi Signals with the phi-ckt + expo.
	//---------------------------------------------------------------------------------------------------------------
	//---------------------------------------------------------------------------
	//plot Phi mass region
	//---------------------------------------------------------------------------
	cout<<"---------------------------------------------------------------------------"<<endl;
	cout<<"fit phi signals with phi-ckt"<<endl;
	cout<<"---------------------------------------------------------------------------"<<endl;
	for(int icent=0; icent<nCentBins; icent++)
	{
		//H1d_M_4Vm_Sig_aftAllCorr[icent]    ->Scale(1, "width");
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetMarkerStyle(20);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetMarkerColor(1);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetLineColor(1);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetLineColor(1);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetYTitle("dN/dM");
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetTitle("Fully Corrected");
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->GetYaxis()->SetTitleSize(0.06);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->GetYaxis()->SetTitleOffset(0.80);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->GetXaxis()->SetTitleSize(0.05);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->GetXaxis()->SetTitleOffset(0.85);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetAxisRange(0.85, 1.15, "x");
		//H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetAxisRange(0.65, 0.85, "x");
		//if(icent==1)H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetMinimum(1.e-5);
		int itemBin_phi    = H1d_M_4Vm_Sig_aftAllCorr[icent]->FindBin(1.015);
		double itemContent = H1d_M_4Vm_Sig_aftAllCorr[icent]->GetBinContent(itemBin_phi);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->SetMaximum(itemContent*1.2);
		H1d_M_4Vm_Sig_aftAllCorr[icent]    ->Draw("pe");
	
		const double LM_Phi = 0.84; double HM_Phi = 1.200;
		
		hist4fit = (TH1D*) H2D_MvsPt_CKT[icent][4]->ProjectionY();
		hist4fit ->Scale(1., "width");

		TF1 *fitf_Phi = new TF1("fitf_Phi", histAsFitf_SigBg, LM_Phi, HM_Phi, 3);
		fitf_Phi->SetParNames("NS",  "a",    "b" );
		fitf_Phi->SetParameters(1.e-5,  1.e5,    5. );
		fitf_Phi->SetParLimits( 0, 0.0,   1.e2 );
		fitf_Phi->SetParLimits( 1, -1.e6, 1.e6 );
		fitf_Phi->SetParLimits( 2, -50, 50 );
		fitf_Phi->SetLineColor(1);
		fitf_Phi->SetLineWidth(2);
		//fitf_Phi->SetNpx(200);

		H1d_M_4Vm_Sig_aftAllCorr[icent]->Fit(fitf_Phi, "E", "", 1.04, 1.2 );
		if(icent==1) H1d_M_4Vm_Sig_aftAllCorr[icent]->Fit(fitf_Phi, "E", "", 0.85, 0.92 );
		
		TFitResultPtr FitReS = H1d_M_4Vm_Sig_aftAllCorr[icent]->Fit(fitf_Phi, "ERS", "", LM_Phi, HM_Phi+0.005 );
		
		TMatrixDSym covTotal = FitReS -> GetCovarianceMatrix();
		covTotal.Print();
		TMatrixDSym covSig;
		covTotal.GetSub(0,0,0,0, covSig);
		covSig.Print();

		

		//exclude the 0.94-1.0 GeV/c2 mass region, then fit the signals

		double chi2 = fitf_Phi ->GetChisquare();
		double ndf  = fitf_Phi ->GetNDF();
		double Chi2overNDF = chi2/ndf;

		double parm[3]; double parmerr[3];

		for(int iparm=0; iparm<3; iparm++)
		{
			parm[iparm]    =  fitf_Phi->GetParameter( iparm );
			parmerr[iparm] =  fitf_Phi->GetParError(  iparm );
		}

		TF1 *fsg_Phi = new TF1("fsg_Phi", histAsFitf_Sig,  LM_Phi, HM_Phi+0.005, 1);
		TF1* fbg_Phi = new TF1("fbg",     bg_formula,      LM_Phi, HM_Phi+0.005, 2);
		fsg_Phi->SetNpx(500);

		for(int isg=0; isg<1; isg++) 
		{
			cout<<parm[isg]<<endl;
			fsg_Phi->SetParameter(isg,  parm[isg]   );
			fsg_Phi->SetParError( isg,  parmerr[isg]);
		}
		for(int ibg=0; ibg<2; ibg++) 
		{
			fbg_Phi->SetParameter(ibg, parm[ibg+1]   );
			fbg_Phi->SetParError( ibg, parmerr[ibg+1]);
		}

		
		double nPhi_fromFit    = fsg_Phi->Integral(0.95, 1.05);
		double nPhiErr_fromFit = fsg_Phi->IntegralError(0.95, 1.05, fsg_Phi->GetParameters(), covSig.GetMatrixArray());



		fsg_Phi->SetLineColor(4);
		fsg_Phi->SetLineWidth(2);
		fsg_Phi->Draw("csame");
		fbg_Phi->SetLineStyle(5);
		fbg_Phi->SetLineColor(2);
		fbg_Phi->SetLineWidth(2);
		fbg_Phi->Draw("lsame");

		TLegend* leg_Phi = new TLegend(0.15, 0.60, 0.45, 0.80);
		leg_Phi->SetBorderSize(0);
		leg_Phi->SetFillColor(0);
		leg_Phi->SetTextSize(0.035);
		leg_Phi->AddEntry( H1d_M_4Vm_Sig_aftAllCorr[icent], "Data",         "lp" );
		leg_Phi->AddEntry( fitf_Phi,                        "#phi+Expo.",   "l"  );
		leg_Phi->AddEntry( fsg_Phi,                         "#phi",         "l"  );
		leg_Phi->AddEntry( fbg_Phi,                         "Expo.",        "l"  );
		leg_Phi->Draw();

		tl.SetTextSize(0.050);
		tl.SetTextColor(1);
		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
		tl.SetTextSize(0.040);
		tl.SetTextColor(1);
		
		tl.DrawLatex(0.63, 0.82, Form("#chi^{2}/ndf = %.2f",      Chi2overNDF) );

		//---------------------------------------------------------------------------
		//calculate the phi yield from integral of ckt fit function
		//double yield_phi = fsg_Phi->Integral();
		
		int    Phi_LowMBin =  H1d_M_4Vm_Sig_aftAllCorr[icent] ->FindBin(0.95);
		int    Phi_HigMBin =  H1d_M_4Vm_Sig_aftAllCorr[icent] ->FindBin(1.05);
		
		double nYerr_Phi   = 0.0;
		double nYvalue_Phi = H1d_M_4Vm_Sig_aftAllCorr[icent]->IntegralAndError(Phi_LowMBin, Phi_HigMBin, nYerr_Phi, "width");
		
		nYvalue_Phi  *= 1.e4;
		nYerr_Phi    *= 1.e4;
		nPhi_fromFit *= 1.e4;
		nPhiErr_fromFit *= 1.e4;

		tl.DrawLatex(0.65, 0.70, Form("Yield of #phi = %.2f#pm%.2f (10^{-4})", nYvalue_Phi, nYerr_Phi) );
		tl.DrawLatex(0.65, 0.60, Form("YdFit of #phi = %.2f#pm%.2f (10^{-4})", nPhi_fromFit,nPhiErr_fromFit) );
		//---------------------------------------------------------------------------

		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Phi_cktFit_aftAllCorr_icent%d.png", icent) );
		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Phi_cktFit_aftAllCorr_icent%d.pdf", icent) );
		//---------------------------------------------------------------------------
		//---------------------------------------------------------------------------
		//Draw Vector mass in different pt bins
		
		//for(int ipt=0; ipt<nPtBins4Vm; ipt++)
		//{
		//	//---------------------------------------------------------------------------
		//	//plot Phi mass region
		//	//---------------------------------------------------------------------------
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->Scale(1, "width");
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetMarkerStyle(20);
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetMarkerColor(1);
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetLineColor(1);
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetLineColor(1);
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetYTitle("dN/dM");
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetTitle("Fully Corrected");
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->GetYaxis()->SetTitleSize(0.06);
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->GetYaxis()->SetTitleOffset(0.80);
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->GetXaxis()->SetTitleSize(0.05);
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->GetXaxis()->SetTitleOffset(0.85);
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetAxisRange(0.85, 1.15, "x");
		//	
		//	int    itemBin_phi = H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->FindBin(1.015);
		//	double itemContent = H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->GetBinContent(itemBin_phi);
		//	//H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->SetMaximum(itemContent*1.2);
		//	
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->Draw("pe");
		//	//H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]    ->Draw("histsame");
		//	//---------------------------------------------------------------------------

		//	//---------------------------------------------------------------------------
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Phi, "E", "", 1.04, 1.2 );
		//	if(icent==1)H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Phi, "E", "", 0.85, 0.92 );
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Phi, "E", "", 1.0, 1.05 );
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Phi, "E", "", 1.0, 1.20 );
		//	H1d_M4VminPt_Sig_aftAllCorr[icent][ipt]->Fit(fitf_Phi, "ER", "", LM_Phi, HM_Phi+0.005 );
		//	//---------------------------------------------------------------------------

		//	double chi2 = fitf_Phi ->GetChisquare();
		//	double ndf  = fitf_Phi ->GetNDF();
		//	double Chi2overNDF = chi2/ndf;

		//	double parm[5]; double parmerr[5];

		//	for(int iparm=0; iparm<5; iparm++)
		//	{
		//		parm[iparm]    =  fitf_Phi->GetParameter( iparm );
		//		parmerr[iparm] =  fitf_Phi->GetParError(  iparm );
		//	}

		//	TF1 *fsg_Phi = new TF1("fsg_Phi", BreitWigner_formula,  LM_Phi, HM_Phi+0.005, 3);
		//	TF1* fbg_Phi = new TF1("fbg",     bg_formula,           LM_Phi, HM_Phi+0.005, 2);
		//	fsg_Phi->SetNpx(20000);

		//	for(int isg=0; isg<3; isg++) 
		//	{
		//		cout<<parm[isg]<<endl;
		//		fsg_Phi->SetParameter(isg,  parm[isg]   );
		//		fsg_Phi->SetParError( isg,  parmerr[isg]);
		//	}
		//	for(int ibg=0; ibg<2; ibg++) 
		//	{
		//		fbg_Phi->SetParameter(ibg, parm[ibg+3]   );
		//		fbg_Phi->SetParError( ibg, parmerr[ibg+3]);
		//	}

		//	fsg_Phi->SetLineColor(4);
		//	fsg_Phi->SetLineWidth(2);
		//	fsg_Phi->Draw("lsame");
		//	fbg_Phi->SetLineStyle(5);
		//	fbg_Phi->SetLineColor(2);
		//	fbg_Phi->SetLineWidth(2);
		//	fbg_Phi->Draw("lsame");

		//	TLegend* leg_Phi = new TLegend(0.15, 0.60, 0.45, 0.80);
		//	leg_Phi->SetBorderSize(0);
		//	leg_Phi->SetFillColor(0);
		//	leg_Phi->SetTextSize(0.035);
		//	leg_Phi->AddEntry( H1d_M4VminPt_Sig_aftAllCorr[icent][ipt], Form("Data: %.1f<p_{T}<%.1f GeV/c", ptBds4Vm[ipt], ptBds4Vm[ipt+1]),                 "lp" );
		//	leg_Phi->AddEntry( fitf_Phi,         "Breit-Wigner+Expo.",   "l"  );
		//	leg_Phi->AddEntry( fsg_Phi,          "Breit-Wigner",         "l"  );
		//	leg_Phi->AddEntry( fbg_Phi,          "Expo.",                "l"  );
		//	leg_Phi->Draw();

		//	tl.SetTextSize(0.050);
		//	tl.SetTextColor(1);
		//	tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
		//	tl.SetTextSize(0.040);
		//	tl.SetTextColor(1);

		//	tl.DrawLatex(0.63, 0.82, Form("#chi^{2}/ndf = %.2f",      Chi2overNDF) );
		//	tl.DrawLatex(0.63, 0.76, Form("M = %.1f #pm %.1f MeV/c^{2}",      parm[2]*1.e3, parmerr[2]*1.e3) );
		//	tl.DrawLatex(0.63, 0.70, Form("#Gamma = %.1f #pm %.1f MeV/c^{2}", parm[1]*1.e3, parmerr[1]*1.e3) );

		//	c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Phi_aftAllCorr_icent%d_ipt%d.png", icent, ipt) );
		//}//ipt
	}
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
void drawExYvsTheory()
{
	double Ratio_Omega2Phi[nNPartBins];
	double RatioErr_Omega2Phi[nNPartBins];

	TCanvas *c1 = new TCanvas("c1");
	gStyle->SetOptStat(0);
	gStyle->SetLineWidth(2);
	gPad->SetTopMargin(0.08);
	gPad->SetBottomMargin(0.11);
	gPad->SetLeftMargin(0.11);
	gPad->SetRightMargin(0.05);
	//---------------------------------------------------------------------------
	//---------------------------------------------------------------------------

	TLatex tl;
	tl.SetTextSize(0.06);
	tl.SetNDC();
	c1->SetLogy(1);

	for(int icent=0; icent<nCentBins; icent++)
	{
		//---------------------------------------------------------------------------
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->Scale(1, "width");
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetMarkerStyle(20);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetMarkerColor(4);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetLineColor(4);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetYTitle("dN/dM");
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetTitle("After All Corrections");
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->GetYaxis()->SetTitleSize(0.06);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->GetYaxis()->SetTitleOffset(0.80);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->GetXaxis()->SetTitleSize(0.05);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->GetXaxis()->SetTitleOffset(0.85);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetAxisRange(0.,3.4, "x");

		if(icent==1)H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetMinimum(1.e-5);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->Draw("pe");

		//rebined CKTSum
		H1d_M_Reb_CKTSum_aftAllCorr[icent] ->Scale(1, "width");
		H1d_M_Reb_CKTSum_aftAllCorr[icent] ->SetLineColor(1);
		H1d_M_Reb_CKTSum_aftAllCorr[icent] ->SetLineWidth(2);
		//H1d_M_Reb_CKTSum_aftAllCorr[icent] ->Draw("histsame");
		H1d_M_Reb_CKTSum_aftAllCorr[icent] ->Draw("chistsame");
		//	
		//	//finner bin CKTSum
		//	H1d_M_CKTSum_aftAllCorr[icent] ->SetLineColor(1);
		//	H1d_M_CKTSum_aftAllCorr[icent] ->SetLineWidth(2);
		//	H1d_M_CKTSum_aftAllCorr[icent] ->Scale(1, "width");
		//	//H1d_M_CKTSum_aftAllCorr[icent] ->Draw("histsame");
		//	H1d_M_CKTSum_aftAllCorr[icent] ->Draw("histsame");


		//	//Expected CKTSum
		//	H1d_M_CKTSum_Expect[icent] ->SetLineColor(2);
		//	H1d_M_CKTSum_Expect[icent] ->SetLineWidth(2);
		//  H1d_M_CKTSum_Expect[icent] ->Scale(1, "width");
		//	//H1d_M_CKTSum_Expect[icent] ->Draw("histsame");
		//	H1d_M_CKTSum_Expect[icent] ->Draw("histsame");

		tl.SetTextSize(0.050);
		tl.SetTextColor(1);
		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
		//tl.DrawLatex(0.15, 0.78, Form("%.2f<p_{T}^{ee}<%.1f GeV/c", goodPtLow, goodPtHig));
		tl.SetTextSize(0.040);
		tl.SetTextSize(0.050);
		tl.SetTextColor(1);

		TLegend* leg0 = new TLegend(0.62, 0.69, 0.92, 0.90);
		leg0->SetBorderSize(0);
		leg0->SetFillColor(0);
		leg0->AddEntry(H1d_M_Reb_Sig_aftAllCorr[icent], "Data",        "lp" );
		leg0->AddEntry(H1d_M_CKTSum_aftAllCorr[icent],  "CKTSum",      "lp" );
		//leg0->AddEntry(H1d_M_CKTSum_aftAllCorr[icent],  "CKTSum aftCorr",      "lp" );
		//leg0->AddEntry(H1d_M_CKTSum_Expect[icent],      "CKTSum Direct",       "lp" );
		leg0->SetTextSize(0.045);
		leg0->Draw();

		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Reb_SigVsCKTSum_aftAllCorr_icent%d.png", icent) );
		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Reb_SigVsCKTSum_aftAllCorr_icent%d.pdf", icent) );
		
		//----------------------------------------------------------------------------------------------------------
		//----------------------------------------------------------------------------------------------------------
		//draw in linear scale to see the Omega and Phi signals: Data vs. CKT
		c1->SetLogy(0);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetAxisRange(0.70, 1.10, "x");
		//if(icent==1)H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetMinimum(1.e-5);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->Draw("pe");

		//rebined CKTSum
		H1d_M_Reb_CKTSum_aftAllCorr[icent] ->Draw("histsame");
		//H1d_M_Reb_CKTSum_aftAllCorr[icent] ->Draw("chistsame");
		//	
		//	//finner bin CKTSum
		//	H1d_M_CKTSum_aftAllCorr[icent] ->SetLineColor(1);
		//	H1d_M_CKTSum_aftAllCorr[icent] ->SetLineWidth(2);
		//	H1d_M_CKTSum_aftAllCorr[icent] ->Scale(1, "width");
		//	//H1d_M_CKTSum_aftAllCorr[icent] ->Draw("histsame");
		//	H1d_M_CKTSum_aftAllCorr[icent] ->Draw("histsame");

		//	//Expected CKTSum
		//	H1d_M_CKTSum_Expect[icent] ->SetLineColor(2);
		//	H1d_M_CKTSum_Expect[icent] ->SetLineWidth(2);
		//	H1d_M_CKTSum_Expect[icent] ->Scale(1, "width");
		//	//H1d_M_CKTSum_Expect[icent] ->Draw("histsame");
		//  H1d_M_CKTSum_Expect[icent] ->Draw("histsame");

		tl.SetTextSize(0.050);
		tl.SetTextColor(1);
		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
		tl.SetTextSize(0.040);
		tl.SetTextSize(0.050);
		tl.SetTextColor(1);

		TLegend* leg4omegaphi = new TLegend(0.75, 0.75, 0.99, 0.99);
		leg4omegaphi->SetBorderSize(0);
		leg4omegaphi->SetFillColor(0);
		leg4omegaphi->AddEntry(H1d_M_Reb_Sig_aftAllCorr[icent], "Data",                "lp" );
		leg4omegaphi->AddEntry(H1d_M_CKTSum_aftAllCorr[icent],  "CKTSum",      "lp" );
		//leg4omegaphi->AddEntry(H1d_M_CKTSum_Expect[icent],      "CKTSum Direct",       "lp" );
		leg4omegaphi->SetTextSize(0.045);
		leg4omegaphi->Draw();

		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_OmegaPhi_SigVsCKT_icent%d.png", icent) );
		
		//calculate the omega/phi ratio
		int    Omega_LowMBin   =  H1d_M_Reb_Sig_aftAllCorr[icent] ->FindBin(OmegaMassWindow[0]);
		int    Omega_HigMBin   =  H1d_M_Reb_Sig_aftAllCorr[icent] ->FindBin(OmegaMassWindow[1]);
		double iOmega_Yerr     =  0.;
		double iOmega_Y        =  H1d_M_Reb_Sig_aftAllCorr[icent] ->IntegralAndError(Omega_LowMBin, Omega_HigMBin, iOmega_Yerr, "width");

		int    Phi_LowMBin     =  H1d_M_Reb_Sig_aftAllCorr[icent] ->FindBin(PhiMassWindow[0]);
		int    Phi_HigMBin     =  H1d_M_Reb_Sig_aftAllCorr[icent] ->FindBin(PhiMassWindow[1]);
		double iPhi_Yerr       =  0.;
		double iPhi_Y          =  H1d_M_Reb_Sig_aftAllCorr[icent] ->IntegralAndError(Phi_LowMBin,    Phi_HigMBin,   iPhi_Yerr,  "width");
		
		//ratioErr = sqrt( (pow(nDen*nNumErr,2)+pow(nNum*nDenErr,2))/pow(nDen,4) );
		double iRatio_Omega2Phi    = iOmega_Y/iPhi_Y;
		double iRatioErr_Omega2Phi = sqrt( (pow(iPhi_Y*iOmega_Yerr,2)+pow(iOmega_Y*iPhi_Yerr,2))/pow(iPhi_Y,4) );
		
		int iNpartIndex = -1;
		if     (icent==1) iNpartIndex = 2;
		else if(icent==2) iNpartIndex = 1;
		else if(icent==3) iNpartIndex = 0;

		Ratio_Omega2Phi[iNpartIndex]    = iRatio_Omega2Phi;
		RatioErr_Omega2Phi[iNpartIndex] = iRatioErr_Omega2Phi;

		//Compare CKTSum_wAcc/PairAcc to CKTSum
		c1->SetLogy(0);
		TH1D* HRatio_CKTSum = (TH1D*) H1d_M_CKTSum_Expect[icent] ->Clone(Form("HRatio_CKTSum_icent%d", icent));
		HRatio_CKTSum -> Divide(H1d_M_CKTSum_aftAllCorr[icent]);

		HRatio_CKTSum ->SetTitle(centTitle[icent]);
		HRatio_CKTSum ->GetYaxis()->SetTitleSize(0.06);
		HRatio_CKTSum ->GetYaxis()->SetTitleOffset(0.80);
		HRatio_CKTSum ->GetXaxis()->SetTitleSize(0.05);
		HRatio_CKTSum ->GetXaxis()->SetTitleOffset(0.85);
		HRatio_CKTSum ->SetYTitle("CKTSum_Direct/CKTSum_aftCorr");
		HRatio_CKTSum ->Draw("pe");
		HRatio_CKTSum ->Draw("histsame");

		c1->SaveAs( Form("./outAnaPlots/Physics/HRatio_CKTSum_icent%d.png", icent) );
		c1->SetLogy(1);
		//---------------------------------------------------------------------------
		
		//compare the CKT within STAR Acc. to the one coming from the CKT_wAccEff/PairEff, check whether they are same
		//or have the large difference?
		
		TH1D* H1D_M_totCKT  = (TH1D*) H2D_MvsPt_totCKT_Reb[icent]            ->ProjectionY(Form("H1D_M_totCKT_icent%d",            icent));
		TH1D* H1D_M_totCKT2 = (TH1D*) H2D_MvsPt_totCKT4PairEffTest_Reb[icent]->ProjectionY(Form("H1D_M_totCKT4PairEffTest_icent%d",icent));

		H1D_M_totCKT ->Scale(1., "width");
		H1D_M_totCKT2->Scale(1., "width");

		H1D_M_totCKT    ->SetYTitle("dN/dM");
		H1D_M_totCKT    ->SetTitle("check PairEff Performance on CKTs in " + centTitle[icent]);
		H1D_M_totCKT    ->GetYaxis()->SetTitleSize(0.06);
		H1D_M_totCKT    ->GetYaxis()->SetTitleOffset(0.80);
		H1D_M_totCKT    ->GetXaxis()->SetTitleSize(0.05);
		H1D_M_totCKT    ->GetXaxis()->SetTitleOffset(0.85);
		
		H1D_M_totCKT    ->SetLineColor(2);
		H1D_M_totCKT    ->SetLineWidth(2);
		H1D_M_totCKT    ->Draw("hist");
		H1D_M_totCKT2   ->SetLineColor(4);
		H1D_M_totCKT2   ->SetLineWidth(2);
		H1D_M_totCKT2   ->Draw("histsame");
	
		TLegend* leg4checkPairEff = new TLegend(0.58, 0.69, 0.89, 0.89);
		leg4checkPairEff->SetBorderSize(0);
		leg4checkPairEff->SetFillColor(0);
		leg4checkPairEff->AddEntry(H1D_M_totCKT,   "Direct CKT",            "lp" );
		leg4checkPairEff->AddEntry(H1D_M_totCKT2,  "CKT_wAccEff/PairEff",   "lp" );
		leg4checkPairEff->SetTextSize(0.045);
		leg4checkPairEff->Draw();

		c1->SaveAs( Form("./outAnaPlots/Physics/H_CKTSumComp_4PairEff_icent%d.png", icent) );
	
		//take ratio: H1D_M_totCKT2/H1D_M_totCKT
		TH1D* HRatio_btwTwoCKT = (TH1D*) H1D_M_totCKT2->Clone();
		HRatio_btwTwoCKT -> Divide(H1D_M_totCKT);
		c1->SetLogy(0);
		HRatio_btwTwoCKT ->SetTitle(centTitle[icent]);
		HRatio_btwTwoCKT ->GetYaxis()->SetTitleSize(0.06);
		HRatio_btwTwoCKT ->GetYaxis()->SetTitleOffset(0.80);
		HRatio_btwTwoCKT ->GetXaxis()->SetTitleSize(0.05);
		HRatio_btwTwoCKT ->GetXaxis()->SetTitleOffset(0.85);
		HRatio_btwTwoCKT ->SetYTitle("(CKT_wAccEff/PairEff)/(Direct CKT)");
		HRatio_btwTwoCKT ->Draw("pe");
		HRatio_btwTwoCKT ->Draw("histsame");

		c1->SaveAs( Form("./outAnaPlots/Physics/HRatio_btwTwoCKT_icent%d.png", icent) );
		c1->SetLogy(1);

		delete H1D_M_totCKT;
		delete H1D_M_totCKT2;
		delete HRatio_btwTwoCKT;
		//---------------------------------------------------------------------------
	}//icent

	//---------------------------------------------------------------------------
	c1->SetLogy(0);
	//---------------------------------------------------------------------------
	double NpartXerr[nNPartBins] = {0. ,0., 0.};
	ge_Ratio_Omega2Phi = new TGraphErrors(nNPartBins, NpartValues, Ratio_Omega2Phi, NpartXerr, RatioErr_Omega2Phi);
	//---------------------------------------------------------------------------
	//---------------------------------------------------------------------------
	TH2D* Htem2d_4Omega2PhiRatio = new TH2D("Htem2d_4Omega2PhiRatio", "", 10, 0, 400, 10, 0.5, 1.70);
	Htem2d_4Omega2PhiRatio->SetTitle("");
	Htem2d_4Omega2PhiRatio->SetYTitle("##omega/##phi");
	Htem2d_4Omega2PhiRatio->SetXTitle("<N_{part}>");
	Htem2d_4Omega2PhiRatio->GetYaxis()->SetTitleSize(0.06);
	Htem2d_4Omega2PhiRatio->GetYaxis()->SetTitleOffset(0.75);
	Htem2d_4Omega2PhiRatio->GetXaxis()->SetTitleSize(0.05);
	Htem2d_4Omega2PhiRatio->GetXaxis()->SetTitleOffset(0.85);
	Htem2d_4Omega2PhiRatio->Draw();

	ge_Ratio_Omega2Phi->SetMarkerStyle(21);
	ge_Ratio_Omega2Phi->SetMarkerSize(1.5);
	ge_Ratio_Omega2Phi->Draw("pesame");

	tl.DrawLatex(0.50, 0.82, Form("#omega-like: %.2f<m<%.2f GeV/c^{2}", OmegaMassWindow[0], OmegaMassWindow[1]));
	tl.DrawLatex(0.50, 0.75, Form("#phi-like:   %.2f<m<%.2f GeV/c^{2}", PhiMassWindow[0],   PhiMassWindow[1])  );

	c1->SaveAs( "./outAnaPlots/Physics/HRatio_Omega2Phi_vsNpart.png" );


	//---------------------------------------------------------------------------
	//---------------------------------------------------------------------------
	const double lowM4AccedExY = 1.05;
	const double higM4AccedExY = 2.5;
	TF1* fit_AccedExY = new TF1("fit_AccedExY", fun_AccedExY, lowM4AccedExY, higM4AccedExY, 2);
	fit_AccedExY->SetLineColor(4);
	fit_AccedExY->SetLineWidth(2);
	fit_AccedExY->SetLineStyle(5);
	fit_AccedExY->SetParameters(1.e-6, 295.);

	const double dY = 2.0;
	for(int icent=0; icent<nCentBins; icent++)
	{
		c1->SetLogy(1);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent] = (TH1D*) H1d_M_Reb_ExY_aftAllCorr[icent]->Clone(Form("H1d_M_Reb_ExY_aftAllCorr_scl1_icent%d",icent));
		//---------------------------------------------------------------------------
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->Scale(1, "width");
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->Scale(1./dY);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->Scale(20/1000.);//scale to unit 20MeV
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->Scale(1./getDnchDy(icent));
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetYTitle("(dN_{Excess}^{2}/dM/dy)/(dN_{ch}/dy) (20 MeV/c^{2})^{-1}");
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetTitle("Excess Yields");
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetMarkerStyle(20);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetMarkerSize(1.3);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetMarkerColor(2);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetLineColor(2);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetLineWidth(2);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetMaximum(1.e-4);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetMinimum(1.e-11);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetAxisRange(0.01, 3.4, "x");
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->GetYaxis()->SetTitleSize(0.05);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->GetYaxis()->SetTitleOffset(0.99);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->GetXaxis()->SetTitleSize(0.05);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->GetXaxis()->SetTitleOffset(0.85);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->Draw("pe");

		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->Fit("fit_AccedExY", "ER", "", lowM4AccedExY, higM4AccedExY);
		fit_AccedExY->Draw("lsame");

		double T_Value = fit_AccedExY->GetParameter(1)*1000.; //here 1000 is the the unit GeV to MeV
		double T_Error = fit_AccedExY->GetParError(1) *1000.;

		tl.SetTextSize(0.050);
		tl.SetTextColor(1);
		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
		//tl.DrawLatex(0.15, 0.78, Form("%.2f<p_{T}^{ee}<%.1f GeV/c", goodPtLow, goodPtHig));
		tl.SetTextColor(1);
		tl.SetTextSize(0.040);
		tl.DrawLatex(0.60, 0.60, Form("%.2f<M_{ee}<%.2f GeV/c^{2}", lowM4AccedExY, higM4AccedExY));
		tl.DrawLatex(0.60, 0.55, "fit by: A*M_{ee}^{3/2}*exp(-M_{ee}/T_{med})");
		tl.SetTextColor(4);
		tl.DrawLatex(0.60, 0.50, Form("T_{med} = %.1f #pm %.1f MeV", T_Value, T_Error));

		//g_Ex_TheorySum_inCent[0][icent]  -> SetLineColor(4);
		//g_Ex_TheorySum_inCent[0][icent]  -> Draw("lsame");
		g_Ex_TheorySum_inCent[1][icent]  -> SetLineColor(1);
		g_Ex_TheorySum_inCent[1][icent]  -> Draw("lsame");
		g_Ex_MediumRho_inCent[1][icent]  -> SetLineColor(kYellow+2);
		g_Ex_MediumRho_inCent[1][icent]  -> Draw("lsame");
		g_Ex_QGPEmission_inCent[1][icent]-> SetLineColor(2);
		g_Ex_QGPEmission_inCent[1][icent]-> Draw("lsame"); 

		TLegend* leg1 = new TLegend(0.65, 0.65, 0.92, 0.90);
		leg1->SetBorderSize(0);
		leg1->SetFillColor(0);
		leg1->AddEntry(H1d_M_Reb_ExY_aftAllCorr_scl1[icent],    "Data(aft ArchCut)",   "lp" );
		leg1->AddEntry(g_Ex_TheorySum_inCent[1][icent],         "Broaden #rho + QGP",  "l"  );
		leg1->AddEntry(g_Ex_MediumRho_inCent[1][icent],         "Broaden #rho",        "lp" );
		leg1->AddEntry(g_Ex_QGPEmission_inCent[1][icent],       "QGP Emission",        "lp" );
		leg1->SetTextSize(0.045);
		leg1->Draw();

		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Reb_ExY_aftAllCorr_icent%d.png", icent) );
		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Reb_ExY_aftAllCorr_icent%d.pdf", icent) );
		//---------------------------------------------------------------------------

		//---------------------------------------------------------------------------
		//Draw same plots in the RhoMass region, in linear sclae
		c1->SetLogy(0);
		//---------------------------------------------------------------------------
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetYTitle("(dN_{Excess}^{2}/dM/dy)/(dN_{ch}/dy) (20 MeV/c^{2})^{-1}");
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetTitle("Excess Yields");
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetMarkerStyle(20);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetMarkerSize(1.3);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetMarkerColor(2);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetLineColor(2);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetLineWidth(2);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetMaximum(2.e-6);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetMinimum(2.e-9);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->SetAxisRange(0.1, 1.1, "x");
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->GetYaxis()->SetTitleSize(0.06);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->GetYaxis()->SetTitleOffset(0.80);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->GetXaxis()->SetTitleSize(0.05);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->GetXaxis()->SetTitleOffset(0.85);
		H1d_M_Reb_ExY_aftAllCorr_scl1[icent]->Draw("pe");

		tl.SetTextSize(0.050);
		tl.SetTextColor(1);
		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
		//tl.DrawLatex(0.15, 0.78, Form("%.2f<p_{T}^{ee}<%.1f GeV/c", goodPtLow, goodPtHig));
		tl.SetTextColor(1);
		tl.SetTextSize(0.040);

		//g_Ex_TheorySum_inCent[0][icent]  -> SetLineColor(4);
		//g_Ex_TheorySum_inCent[0][icent]  -> Draw("lsame");
		g_Ex_TheorySum_inCent[1][icent]  -> SetLineColor(1);
		g_Ex_TheorySum_inCent[1][icent]  -> Draw("lsame");
		g_Ex_MediumRho_inCent[1][icent]  -> SetLineColor(kYellow+2);
		g_Ex_MediumRho_inCent[1][icent]  -> Draw("lsame");
		g_Ex_QGPEmission_inCent[1][icent]-> SetLineColor(2);
		g_Ex_QGPEmission_inCent[1][icent]-> Draw("lsame"); 

		TLegend* leg_RhoM = new TLegend(0.65, 0.65, 0.92, 0.90);
		leg_RhoM->SetBorderSize(0);
		leg_RhoM->SetFillColor(0);
		leg_RhoM->AddEntry(H1d_M_Reb_ExY_aftAllCorr_scl1[icent],    "Data(aft ArchCut)",   "lp" );
		leg_RhoM->AddEntry(g_Ex_TheorySum_inCent[1][icent],         "Broaden #rho + QGP",  "l"  );
		leg_RhoM->AddEntry(g_Ex_MediumRho_inCent[1][icent],         "Broaden #rho",        "lp" );
		leg_RhoM->AddEntry(g_Ex_QGPEmission_inCent[1][icent],       "QGP Emission",        "lp" );
		leg_RhoM->SetTextSize(0.045);
		leg_RhoM->Draw();

		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Reb_RhoMassExY_aftAllCorr_icent%d.png", icent) );
		//c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Reb_RhoMassExY_aftAllCorr_icent%d.pdf", icent) );
		//---------------------------------------------------------------------------

	}//icent

	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//	TH2D* Htem2d_4AccedExY = new TH2D("Htem2d_4AccedExY", "", 10, 0.0, 3.0, 10, 1.e-7, 5.e-1);
	//	Htem2d_4AccedExY->SetTitle("");
	//	Htem2d_4AccedExY->SetYTitle("dN/dM (Data-CKT) after Acc. Corr");
	//	Htem2d_4AccedExY->SetXTitle("M_{ee} (GeV/c^{2})");
	//	Htem2d_4AccedExY->GetYaxis()->SetTitleSize(0.06);
	//	Htem2d_4AccedExY->GetYaxis()->SetTitleOffset(0.75);
	//	Htem2d_4AccedExY->GetXaxis()->SetTitleSize(0.05);
	//	Htem2d_4AccedExY->GetXaxis()->SetTitleOffset(0.85);
	//
	//	//const double lowM4AccedExY = 1.26;
	//	//const double higM4AccedExY = 2.60;
	//	const double lowM4AccedExY = 1.20;
	//	const double higM4AccedExY = 3.0;
	//	TF1* fit_AccedExY = new TF1("fit_AccedExY", fun_AccedExY, lowM4AccedExY, higM4AccedExY, 2);
	//	fit_AccedExY->SetLineColor(4);
	//	fit_AccedExY->SetLineWidth(2);
	//	fit_AccedExY->SetParameters(1.e-2, 250.);
	//
	//	//read in Rapp Predictions
	//	TFile* inf_theory    = new TFile("./inputfiles/Histograms4Yi.root", "read");
	//	TGraph* ge_RappSum   = (TGraph*) inf_theory->Get("gRappFull");
	//	TGraph* ge_RappHGMed = (TGraph*) inf_theory->Get("gRappHGMedFull");
	//	TGraph* ge_RappQGP   = (TGraph*) inf_theory->Get("gRappQGPFull");
	//
	//	ge_RappSum  ->SetLineWidth(2);
	//	ge_RappHGMed->SetLineWidth(2);
	//	ge_RappQGP  ->SetLineWidth(2);
	//	ge_RappSum  ->SetLineColor(1);
	//	ge_RappHGMed->SetLineColor(8);
	//	ge_RappQGP  ->SetLineColor(2);
	//
	//	c1->SetLogy(1);
	//	for(int icent=0; icent<nCentBins; icent++)
	//	{
	//		Htem2d_4AccedExY->Draw();
	//
	//		HExY_M_AftAcc_AllMass[icent]->SetMarkerStyle(20);
	//		HExY_M_AftAcc_AllMass[icent]->SetMarkerSize(2);
	//		HExY_M_AftAcc_AllMass[icent]->SetMarkerColor(2);
	//		HExY_M_AftAcc_AllMass[icent]->SetLineColor(2);
	//		HExY_M_AftAcc_AllMass[icent]->SetLineWidth(2);
	//		HExY_M_AftAcc_AllMass[icent]->Draw("pesame");
	//
	//		ge_RappSum  ->Draw("lsame");
	//		ge_RappHGMed->Draw("lsame");
	//		ge_RappQGP  ->Draw("lsame");
	//		
	//		HExY_M_AftAcc_AllMass[icent]->Fit("fit_AccedExY", "ER", "", lowM4AccedExY, higM4AccedExY);
	//
	//		fit_AccedExY->Draw("lsame");
	//
	//		double T_Value = fit_AccedExY->GetParameter(1)*1000.; //here 1000 is the the unit GeV to MeV
	//		double T_Error = fit_AccedExY->GetParError(1)*1000.;
	//		
	//		tl.SetTextSize(0.050);
	//		tl.SetTextColor(1);
	//		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
	//		tl.SetTextSize(0.040);
	//		tl.SetTextSize(0.050);
	//		tl.SetTextColor(1);
	//		tl.DrawLatex(0.60, 0.80, Form("%.2f<M_{ee}<%.2f GeV/c^{2}", lowM4AccedExY, higM4AccedExY));
	//		tl.DrawLatex(0.60, 0.72, "fit by: M_{ee}^{3/2}*exp(-M_{ee}/T)");
	//		tl.SetTextColor(2);
	//		tl.DrawLatex(0.60, 0.65, Form("T = %.1f #pm %.1f MeV", T_Value, T_Error));
	//
	//		c1->SaveAs( Form("./outAnaPlots/Physics/HExY_M_AftAcc_AllMass_icent%d.png", icent) );
	//		c1->SaveAs( Form("./outAnaPlots/Physics/HExY_M_AftAcc_AllMass_icent%d.pdf", icent) );
	//	}
	//	delete Htem2d;
	//	delete Htem2d_ratio;
	//	delete Htem2d_ratio2;
	//	delete Htem2d_4MtSpec;
	//	delete Htem2d4Mt;
	//	delete Htem2d_4AccedExY;
	//	delete c1;
}
void fitAndRmVectMeson4ExYield()
{
	TCanvas *c1 = new TCanvas("c1");
	gStyle->SetOptStat(0);
	gStyle->SetLineWidth(2);
	gPad->SetTopMargin(0.08);
	gPad->SetBottomMargin(0.11);
	gPad->SetLeftMargin(0.11);
	gPad->SetRightMargin(0.05);
	//---------------------------------------------------------------------------
	//---------------------------------------------------------------------------
	TLatex tl;
	tl.SetTextSize(0.06);
	tl.SetNDC();
	c1->SetLogy(1);
	
	//---------------------------------------------------------------------------
	TGraph* gh_RappSum;
	const double lowM4AccedExY = 1.10;//1.05;
	const double higM4AccedExY = 2.9;
	TF1* fit_AccedExY = new TF1("fit_AccedExY", fun_AccedExY, lowM4AccedExY, higM4AccedExY, 2);
	fit_AccedExY->SetLineColor(4);
	fit_AccedExY->SetLineWidth(2);
	fit_AccedExY->SetLineStyle(5);
	fit_AccedExY->SetParameters(1.e-6, 295.);
	const double dY = 2.0;

	for(int icent=0; icent<nCentBins; icent++)
	{
		//---------------------------------------------------------------------------
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetMarkerStyle(20);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetMarkerColor(4);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetLineColor(4);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetYTitle("dN/dM");
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetTitle("After All Corrections");
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->GetYaxis()->SetTitleSize(0.06);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->GetYaxis()->SetTitleOffset(0.80);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->GetXaxis()->SetTitleSize(0.05);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->GetXaxis()->SetTitleOffset(0.85);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetAxisRange(0.,3.4, "x");

		if(icent==1)H1d_M_Reb_Sig_aftAllCorr[icent]    ->SetMinimum(1.e-5);
		H1d_M_Reb_Sig_aftAllCorr[icent]    ->Draw("pe");

		H1d_M_Reb_CKTnoVect_aftAllCorr[icent] ->Scale(1.,"width");
		H1d_M_Reb_CKTnoVect_aftAllCorr[icent] ->SetLineColor(1);
		H1d_M_Reb_CKTnoVect_aftAllCorr[icent] ->SetLineWidth(2);
		H1d_M_Reb_CKTnoVect_aftAllCorr[icent] ->SetLineStyle(2);
		H1d_M_Reb_CKTnoVect_aftAllCorr[icent] ->Draw("chistsame");

		H1d_M_Reb_CKTw_aftAllCorr[icent] ->Scale(1.,"width");
		H1d_M_Reb_CKTw_aftAllCorr[icent] ->SetLineColor(1);
		H1d_M_Reb_CKTw_aftAllCorr[icent] ->SetLineWidth(2);
		H1d_M_Reb_CKTw_aftAllCorr[icent] ->SetLineStyle(2);
		//H1d_M_Reb_CKTw_aftAllCorr[icent] ->Draw("chistsame");
		
		H1d_M_Reb_CKTphi_aftAllCorr[icent] ->Scale(1.,"width");
		H1d_M_Reb_CKTphi_aftAllCorr[icent] ->SetLineColor(1);
		H1d_M_Reb_CKTphi_aftAllCorr[icent] ->SetLineWidth(2);
		H1d_M_Reb_CKTphi_aftAllCorr[icent] ->SetLineStyle(2);
		//H1d_M_Reb_CKTphi_aftAllCorr[icent] ->Draw("chistsame");

		gh_RappSum = (TGraph*)g_Ex_TheorySum_inCent[1][icent]->Clone();
		for(int ip=0; ip<gh_RappSum->GetN(); ip++) gh_RappSum->GetY()[ip] *= (getDnchDy(icent)*1000./20.);
		
		//gh_RappSum ->SetLineColor(2);
		//gh_RappSum ->Draw("lsame");

		//-------------------------------------------------------------------------------------
		//construct the fitting function: fitf = ckt_noVect + a*omega + b*phi + c*RappSum
		hist1_4fit = (TH1D*)H1d_M_Reb_CKTnoVect_aftAllCorr[icent]->Clone();
		hist2_4fit = (TH1D*)H1d_M_Reb_CKTw_aftAllCorr[icent]     ->Clone();
		hist3_4fit = (TH1D*)H1d_M_Reb_CKTphi_aftAllCorr[icent]   ->Clone();
		gh_4fit    = (TGraph*) gh_RappSum->Clone();
		
		TF1* fitf  = new TF1("fitf", HistGraphAsFitf, 0.40, 1.2, 3);
		fitf ->SetParNames(  "a",  "b",    "c" );
		fitf ->SetParameters(1.0,   1.,    1.  );
		fitf ->SetParLimits( 0, 0.6, 1.5 );
		fitf ->SetParLimits( 1, 0.6, 1.5 );
		fitf ->SetParLimits( 2, 0.6, 2.5 );
		
		const double fitLowM = 0.40;
		const double fitHigM = 1.20;

		TFitResultPtr FitReS = H1d_M_Reb_Sig_aftAllCorr[icent]->Fit(fitf, "ERS0", "", fitLowM, fitHigM );
		double chi2        = fitf->GetChisquare();
		double ndf         = fitf->GetNDF();
		double Chi2overNDF = chi2/ndf;
		const double a = fitf->GetParameter(0);
		const double b = fitf->GetParameter(1);
		const double c = fitf->GetParameter(2);
		
		tl.SetTextSize(0.035);
		tl.DrawLatex(0.25, 0.79, Form("fit range: %.2f<M_{ee}<%.2f GeV/c^{2}", fitLowM, fitHigM) );
		tl.DrawLatex(0.25, 0.75, Form("a = %.2f, b = %.2f, c = %.2f", a, b, c) );
		tl.DrawLatex(0.33, 0.70, Form("#chi^{2}/ndf = %.2f",      Chi2overNDF) );

		H1d_M_Reb_CKTw_aftAllCorr[icent]  ->Scale(a);
		H1d_M_Reb_CKTphi_aftAllCorr[icent]->Scale(b);
		H1d_M_Reb_CKTw_aftAllCorr[icent]  ->Draw("chistsame");
		H1d_M_Reb_CKTphi_aftAllCorr[icent]->Draw("chistsame");
		for(int ip=0; ip<gh_RappSum->GetN(); ip++) gh_RappSum->GetY()[ip] *= c;
		gh_RappSum ->SetLineColor(2);
		gh_RappSum ->SetLineStyle(2);
		gh_RappSum ->Draw("lsame");

		//-------------------------------------------------------------------------------------
		//Add CKT_noVect + a*Omega + b*Phi to get the new total cocktail
		TH1D* H1d_M_Reb_CKTAllAdd_aftAllCorr = (TH1D*)H1d_M_Reb_CKTnoVect_aftAllCorr[icent]->Clone();
		H1d_M_Reb_CKTAllAdd_aftAllCorr ->Add(H1d_M_Reb_CKTw_aftAllCorr[icent]  );
		H1d_M_Reb_CKTAllAdd_aftAllCorr ->Add(H1d_M_Reb_CKTphi_aftAllCorr[icent]);
		H1d_M_Reb_CKTAllAdd_aftAllCorr ->SetLineColor(1);
		H1d_M_Reb_CKTAllAdd_aftAllCorr ->SetLineWidth(2);
		H1d_M_Reb_CKTAllAdd_aftAllCorr ->SetLineStyle(1);
		H1d_M_Reb_CKTAllAdd_aftAllCorr ->Draw("chistsame");
		//-------------------------------------------------------------------------------------

		tl.SetTextSize(0.050);
		tl.SetTextColor(1);
		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
		//tl.DrawLatex(0.15, 0.78, Form("%.2f<p_{T}^{ee}<%.1f GeV/c", goodPtLow, goodPtHig));
		tl.SetTextSize(0.040);
		tl.SetTextSize(0.050);
		tl.SetTextColor(1);

		TLegend* leg0 = new TLegend(0.60, 0.60, 0.92, 0.90);
		leg0->SetBorderSize(0);
		leg0->SetFillColor(0);
		leg0->AddEntry(H1d_M_Reb_Sig_aftAllCorr[icent],        "Data",             "lp" );
		leg0->AddEntry(H1d_M_Reb_CKTnoVect_aftAllCorr[icent],  "CKTnoOmegaPhi",    "lp" );
		leg0->AddEntry(H1d_M_Reb_CKTw_aftAllCorr[icent],       "a*#omega",           "lp" );
		leg0->AddEntry(H1d_M_Reb_CKTphi_aftAllCorr[icent],     "b*#phi",             "lp" );
		leg0->AddEntry(gh_RappSum,                             "c*RappCurve",        "lp" );
		leg0->AddEntry(H1d_M_Reb_CKTAllAdd_aftAllCorr,         "CKTSUM",             "lp" );
		leg0->SetTextSize(0.045);
		leg0->Draw();

		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Reb_SigVsCKTwFit_aftAllCorr_icent%d.png", icent) );
		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Reb_SigVsCKTwFit_aftAllCorr_icent%d.pdf", icent) );
		
		//------------------------------------------------------------------------------------------------------------------------------------------------
		//calculate the Excess Yield = Data - CKTSUM, after scaled to unit then compare to the theory curve
		//------------------------------------------------------------------------------------------------------------------------------------------------
		TH1D* H1d_M_ExY_new_aftAllCorr = (TH1D*) H1d_M_Reb_Sig_aftAllCorr[icent]->Clone();
		H1d_M_ExY_new_aftAllCorr       -> Add(H1d_M_Reb_CKTAllAdd_aftAllCorr, -1);
		
		TH1D* H1d_M_ExY_new_4ccdyRatio = (TH1D*) H1d_M_ExY_new_aftAllCorr->Clone();
		
		//make a copy for the later Excess Yields calculations
		H1d_M_ExY_new_aftAllCorr_copy[icent] = (TH1D*) H1d_M_ExY_new_aftAllCorr->Clone(Form("H1d_M_ExY_new_aftAllCorr_icent%d",icent));

		//H1d_M_ExY_new_aftAllCorr = (TH1D*) H1d_M_Reb_ExY_aftAllCorr[icent]->Clone(Form("H1d_M_Reb_ExY_aftAllCorr_scl1_icent%d",icent));
		//---------------------------------------------------------------------------
		//H1d_M_ExY_new_aftAllCorr->Scale(1, "width");
		H1d_M_ExY_new_aftAllCorr->Scale(1./dY);
		H1d_M_ExY_new_aftAllCorr->Scale(20/1000.);//scale to unit 20MeV
		H1d_M_ExY_new_aftAllCorr->Scale(1./getDnchDy(icent));
		H1d_M_ExY_new_aftAllCorr->SetYTitle("(dN_{Excess}^{2}/dM/dy)/(dN_{ch}/dy) (20 MeV/c^{2})^{-1}");
		H1d_M_ExY_new_aftAllCorr->SetTitle("Excess Yields");
		H1d_M_ExY_new_aftAllCorr->SetMarkerStyle(20);
		H1d_M_ExY_new_aftAllCorr->SetMarkerSize(1.3);
		H1d_M_ExY_new_aftAllCorr->SetMarkerColor(2);
		H1d_M_ExY_new_aftAllCorr->SetLineColor(2);
		H1d_M_ExY_new_aftAllCorr->SetLineWidth(2);
		H1d_M_ExY_new_aftAllCorr->SetMaximum(1.e-4);
		H1d_M_ExY_new_aftAllCorr->SetMinimum(1.e-11);
		H1d_M_ExY_new_aftAllCorr->SetAxisRange(0.01, 2.90, "x");
		H1d_M_ExY_new_aftAllCorr->GetYaxis()->SetTitleSize(0.05);
		H1d_M_ExY_new_aftAllCorr->GetYaxis()->SetTitleOffset(0.99);
		H1d_M_ExY_new_aftAllCorr->GetXaxis()->SetTitleSize(0.05);
		H1d_M_ExY_new_aftAllCorr->GetXaxis()->SetTitleOffset(0.85);
		H1d_M_ExY_new_aftAllCorr->Draw("pe");

		H1d_M_ExY_new_aftAllCorr->Fit("fit_AccedExY", "ER", "", lowM4AccedExY, higM4AccedExY);
		fit_AccedExY->Draw("lsame");

		double T_Value = fit_AccedExY->GetParameter(1)*1000.; //here 1000 is the the unit GeV to MeV
		double T_Error = fit_AccedExY->GetParError(1) *1000.;

		tl.SetTextSize(0.050);
		tl.SetTextColor(1);
		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
		tl.DrawLatex(0.18, 0.75, "Subtract floating #omega, #phi");
		//tl.DrawLatex(0.15, 0.78, Form("%.2f<p_{T}^{ee}<%.1f GeV/c", goodPtLow, goodPtHig));
		tl.SetTextColor(1);
		tl.SetTextSize(0.040);
		tl.DrawLatex(0.60, 0.60, Form("%.2f<M_{ee}<%.2f GeV/c^{2}", lowM4AccedExY, higM4AccedExY));
		tl.DrawLatex(0.60, 0.55, "fit by: A*M_{ee}^{3/2}*exp(-M_{ee}/T_{med})");
		tl.SetTextColor(4);
		tl.DrawLatex(0.60, 0.50, Form("T_{med} = %.1f #pm %.1f MeV", T_Value, T_Error));
		tl.SetTextColor(1);

		//g_Ex_TheorySum_inCent[0][icent]  -> SetLineColor(4);
		//g_Ex_TheorySum_inCent[0][icent]  -> Draw("lsame");
		g_Ex_TheorySum_inCent[1][icent]  -> SetLineColor(1);
		g_Ex_TheorySum_inCent[1][icent]  -> Draw("lsame");
		g_Ex_MediumRho_inCent[1][icent]  -> SetLineColor(kYellow+2);
		g_Ex_MediumRho_inCent[1][icent]  -> Draw("lsame");
		g_Ex_QGPEmission_inCent[1][icent]-> SetLineColor(2);
		g_Ex_QGPEmission_inCent[1][icent]-> Draw("lsame"); 

		TLegend* leg1 = new TLegend(0.65, 0.65, 0.92, 0.90);
		leg1->SetBorderSize(0);
		leg1->SetFillColor(0);
		leg1->AddEntry(H1d_M_ExY_new_aftAllCorr,                "Data",                "lp" );
		leg1->AddEntry(g_Ex_TheorySum_inCent[1][icent],         "Broaden #rho + QGP",  "l"  );
		leg1->AddEntry(g_Ex_MediumRho_inCent[1][icent],         "Broaden #rho",        "lp" );
		leg1->AddEntry(g_Ex_QGPEmission_inCent[1][icent],       "QGP Emission",        "lp" );
		leg1->SetTextSize(0.045);
		leg1->Draw();

		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Reb_ExYnew_icent%d.png", icent) );
		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Reb_ExYnew_icent%d.pdf", icent) );
		//---------------------------------------------------------------------------

		//---------------------------------------------------------------------------
		//---------------------------------------------------------------------------
		//calculate the ExY/cc and ExY/dy ratios to estimate the sensitivity of Temperature on cc, dy scale factors
		//---------------------------------------------------------------------------
		//---------------------------------------------------------------------------
		//c1->SetLogy(0);
		tl.SetTextSize(0.050);
		gPad->SetGridy();
		TH1D* HRatio_ExY2cc = (TH1D*)H1d_M_ExY_new_4ccdyRatio->Clone("HRatio_ExY2cc");
		H1d_M_Reb_CKTcc_aftAllCorr[icent] ->Scale(1., "width");
		HRatio_ExY2cc              ->Divide(H1d_M_Reb_CKTcc_aftAllCorr[icent]);
		
		HRatio_ExY2cc->SetYTitle("Ratio (ExcessYield/c#bar{c})");
		HRatio_ExY2cc->SetTitle("");
		HRatio_ExY2cc->SetMarkerStyle(20);
		HRatio_ExY2cc->SetMarkerSize(1.3);
		HRatio_ExY2cc->SetMarkerColor(2);
		HRatio_ExY2cc->SetLineColor(2);
		HRatio_ExY2cc->SetLineWidth(2);
		//HRatio_ExY2cc->SetMaximum(1.e-4);
		HRatio_ExY2cc->SetMinimum(1.e-1);
		HRatio_ExY2cc->SetAxisRange(0.05, 3.4, "x");
		HRatio_ExY2cc->GetYaxis()->SetTitleSize(0.06);
		HRatio_ExY2cc->GetYaxis()->SetTitleOffset(0.85);
		HRatio_ExY2cc->GetXaxis()->SetTitleSize(0.05);
		HRatio_ExY2cc->GetXaxis()->SetTitleOffset(0.85);

		HRatio_ExY2cc->SetMarkerColor(1);
		HRatio_ExY2cc->SetLineColor(1);
		HRatio_ExY2cc->Draw("pe");
		
		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
		
		c1->SaveAs( Form("./outAnaPlots/Physics/HRatio_ExY2cc_icent%d.png", icent) );
		c1->SetLogy(1);

		TH1D* HRatio_ExY2dy = (TH1D*)H1d_M_ExY_new_4ccdyRatio->Clone("HRatio_ExY2dy");
		H1d_M_Reb_CKTdy_aftAllCorr[icent] ->Scale(1., "width");
		HRatio_ExY2dy              ->Divide(H1d_M_Reb_CKTdy_aftAllCorr[icent]);
		
		HRatio_ExY2dy->SetYTitle("Ratio (ExcessYield/Drell-Yan)");
		HRatio_ExY2dy->SetTitle("");
		HRatio_ExY2dy->SetMarkerStyle(20);
		HRatio_ExY2dy->SetMarkerSize(1.3);
		HRatio_ExY2dy->SetMarkerColor(2);
		HRatio_ExY2dy->SetLineColor(2);
		HRatio_ExY2dy->SetLineWidth(2);
		//HRatio_ExY2dy->SetMaximum(1.e-4);
		HRatio_ExY2dy->SetMinimum(1.e-1);
		HRatio_ExY2dy->SetAxisRange(0.05, 3.4, "x");
		HRatio_ExY2dy->GetYaxis()->SetTitleSize(0.06);
		HRatio_ExY2dy->GetYaxis()->SetTitleOffset(0.85);
		HRatio_ExY2dy->GetXaxis()->SetTitleSize(0.05);
		HRatio_ExY2dy->GetXaxis()->SetTitleOffset(0.85);

		HRatio_ExY2dy->SetMarkerColor(1);
		HRatio_ExY2dy->SetLineColor(1);
		HRatio_ExY2dy->Draw("pe");
		
		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
		
		c1->SaveAs( Form("./outAnaPlots/Physics/HRatio_ExY2dy_icent%d.png", icent) );
		
		TH1D* H1d_M_Reb_CKTccdy_aftAllCorr = (TH1D*)H1d_M_Reb_CKTcc_aftAllCorr[icent] -> Clone();
		H1d_M_Reb_CKTccdy_aftAllCorr ->Add(H1d_M_Reb_CKTdy_aftAllCorr[icent]);

		TH1D* HRatio_ExY2ccdy = (TH1D*)H1d_M_ExY_new_4ccdyRatio->Clone("HRatio_ExY2ccdy");
		HRatio_ExY2ccdy              ->Divide(H1d_M_Reb_CKTccdy_aftAllCorr);
		
		HRatio_ExY2ccdy->SetYTitle("Ratio (ExcessYield/(c#bar{c}+Drell-Yan))");
		HRatio_ExY2ccdy->SetTitle("");
		HRatio_ExY2ccdy->SetMarkerStyle(20);
		HRatio_ExY2ccdy->SetMarkerSize(1.3);
		HRatio_ExY2ccdy->SetMarkerColor(2);
		HRatio_ExY2ccdy->SetLineColor(2);
		HRatio_ExY2ccdy->SetLineWidth(2);
		//HRatio_ExY2ccdy->SetMaximum(1.e-4);
		HRatio_ExY2ccdy->SetMinimum(1.e-1);
		HRatio_ExY2ccdy->SetAxisRange(0.05, 3.4, "x");
		HRatio_ExY2ccdy->GetYaxis()->SetTitleSize(0.06);
		HRatio_ExY2ccdy->GetYaxis()->SetTitleOffset(0.85);
		HRatio_ExY2ccdy->GetXaxis()->SetTitleSize(0.05);
		HRatio_ExY2ccdy->GetXaxis()->SetTitleOffset(0.85);

		HRatio_ExY2ccdy->SetMarkerColor(1);
		HRatio_ExY2ccdy->SetLineColor(1);
		HRatio_ExY2ccdy->Draw("pe");
		
		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
		
		c1->SaveAs( Form("./outAnaPlots/Physics/HRatio_ExY2ccdy_icent%d.png", icent) );

	}//icent
}	
//------------------------------------------------------------------------------------------------------------------------------------------------------
void getAllExcessYields( )
{
	double Rho_ExY[nNPartBins];
	double Rho_ExYerr[nNPartBins];
	double IMR_ExY[nNPartBins];
	double IMR_ExYerr[nNPartBins];
	double Omega_Y[nNPartBins];
	double Omega_Yerr[nNPartBins];
	double Phi_Y[nNPartBins];
	double Phi_Yerr[nNPartBins];
	//------------------------------------------------------------------------------------------------------------------------------------------------------
	for(int icent=0; icent<nCentBins; icent++)
	{
		H1d_DataMinusCKT[icent] = (TH1D*) H1d_M_ExY_new_aftAllCorr_copy[icent]->Clone( Form("H1d_DataMinusCKT_icent%d", icent) );
		H1d_DataMinusCKT[icent] -> SetName(Form("H1d_DataMinusCKT_icent%d", icent));
	
		//get the excess yield for rho mass region for 0-10, 10-40, 40-80% centrality, but need to reverse the order
		if(icent>0) //let 1, 2, 3 pass
		{
			int iNpartIndex = -1;
			if     (icent==1) iNpartIndex = 2;
			else if(icent==2) iNpartIndex = 1;
			else if(icent==3) iNpartIndex = 0;
			
			//--------------------------------------------------------------------------------------------------------------------
			//get Rho ExY
			int    Rho_LowMBin     =  H1d_DataMinusCKT[icent]->FindBin(RhoMassWindow[0]);
			int    Rho_HigMBin     =  H1d_DataMinusCKT[icent]->FindBin(RhoMassWindow[1]);
			double iRho_ExYerr     = 0.;
			double iRho_ExY        = H1d_DataMinusCKT[icent] ->IntegralAndError(Rho_LowMBin, Rho_HigMBin, iRho_ExYerr, "width");
			
			cout<<"icent: "<<icent<<" # Rho Yield = "<<iRho_ExY<<" +/- "<<iRho_ExYerr<<endl;
			
			//get IMR ExY
			int    IMR_LowMBin     =  H1d_DataMinusCKT[icent]->FindBin(IMRMassWindow[0]);
			int    IMR_HigMBin     =  H1d_DataMinusCKT[icent]->FindBin(IMRMassWindow[1]);
			double iIMR_ExYerr     = 0.;
			double iIMR_ExY        = H1d_DataMinusCKT[icent] ->IntegralAndError(IMR_LowMBin, IMR_HigMBin, iIMR_ExYerr, "width");
			//--------------------------------------------------------------------------------------------------------------------
			//get Omega Y
			int    Omega_LowMBin   =  H1d_M_Reb_Sig_aftAllCorr[icent]->FindBin(OmegaMassWindow[0]);
			int    Omega_HigMBin   =  H1d_M_Reb_Sig_aftAllCorr[icent]->FindBin(OmegaMassWindow[1]);
			double iOmega_Yerr     = 0.;
			double iOmega_Y        = H1d_M_Reb_Sig_aftAllCorr[icent] ->IntegralAndError(Omega_LowMBin, Omega_HigMBin, iOmega_Yerr, "width");
			
			//get Phi Y
			int    Phi_LowMBin     =  H1d_M_Reb_Sig_aftAllCorr[icent]->FindBin(PhiMassWindow[0]);
			int    Phi_HigMBin     =  H1d_M_Reb_Sig_aftAllCorr[icent]->FindBin(PhiMassWindow[1]);
			double iPhi_Yerr       = 0.;
			double iPhi_Y          = H1d_M_Reb_Sig_aftAllCorr[icent] ->IntegralAndError( Phi_LowMBin,   Phi_HigMBin,   iPhi_Yerr,  "width");
			//--------------------------------------------------------------------------------------------------------------------
			
			double iNpartValue      = NpartValues[iNpartIndex];
			
			Rho_ExY[iNpartIndex]    = iRho_ExY/iNpartValue;
			Rho_ExYerr[iNpartIndex] = iRho_ExYerr/iNpartValue;
			
			IMR_ExY[iNpartIndex]    = (iIMR_ExY/iNpartValue);
			IMR_ExYerr[iNpartIndex] = (iIMR_ExYerr/iNpartValue);
			
			Omega_Y[iNpartIndex]    = iOmega_Y/iNpartValue;
			Omega_Yerr[iNpartIndex] = iOmega_Yerr/iNpartValue;
			Phi_Y[iNpartIndex]      = iPhi_Y/iNpartValue;
			Phi_Yerr[iNpartIndex]   = iPhi_Yerr/iNpartValue;
		}//only for icent=1,2,3, no 0
		else
		{
			//get Rho ExY
			int    Rho_LowMBin     =  H1d_DataMinusCKT[icent]->FindBin(RhoMassWindow[0]);
			int    Rho_HigMBin     =  H1d_DataMinusCKT[icent]->FindBin(RhoMassWindow[1]);
			double iRho_ExYerr     = 0.;
			double iRho_ExY        = H1d_DataMinusCKT[icent] ->IntegralAndError(Rho_LowMBin, Rho_HigMBin, iRho_ExYerr, "width");

			cout<<"0-80 centrality, # Rho Yield = "<<iRho_ExY<<" +/- "<<iRho_ExYerr<<endl;
		}
	}//icent

	double NpartXerr[nNPartBins] = {0. ,0., 0.};
	ge_ExYVsNpart_Rho = new TGraphErrors(nNPartBins, NpartValues, Rho_ExY, NpartXerr, Rho_ExYerr);
	ge_ExYVsNpart_IMR = new TGraphErrors(nNPartBins, NpartValues, IMR_ExY, NpartXerr, IMR_ExYerr);

	ge_YVsNpart_Omega = new TGraphErrors(nNPartBins, NpartValues, Omega_Y, NpartXerr, Omega_Yerr);
	ge_YVsNpart_Phi   = new TGraphErrors(nNPartBins, NpartValues, Phi_Y,   NpartXerr, Phi_Yerr  );

}//tem
//------------------------------------------------------------------------------------------------------------------------------------------------------
	//HRatioVsPt_IMR   = new TH1D("HRatioVsPt_IMR",   "Ratio of Data/CKT in 0-80% vs. pt",       nPtBins4Phys, PtBDs4Phys);
	//HExYVsPt_IMR     = new TH1D("HExYVsPt_IMR",     "Excess Yields: Data-CKT in 0-80% vs. pt", nPtBins4Phys, PtBDs4Phys);
	//
	//	for(int ipt=0; ipt<nPtBins4Phys; ipt++)
	//	{
	//		H1D_Ratio2CKT_4PtPhys[ipt] = (TH1D*) calRatioH1toH2( H1d_M_Reb4PtPhys_Sig_aftEff[ipt], H_M_totCKT_Reb4PtPhys[ipt], Form("H1D_Ratio2CKT_4PtPhys_ipt%d", ipt) );
	//	
	//		//calculate the Ratio within mass window(1.2-2.8 GeV/c2) as function of pt
	//		double data_value = 0.;
	//		double data_err   = 0.;
	//		double ckt_value  = 0.;
	//		double ckt_err    = 0.;
	//		
	//		int LowMBin =  H1d_M_Reb4PtPhys_Sig_aftEff[ipt]->FindBin(IMRMassWindow[0]);
	//		int HigMBin =  H1d_M_Reb4PtPhys_Sig_aftEff[ipt]->FindBin(IMRMassWindow[1]);
	//
	//		data_value  = H1d_M_Reb4PtPhys_Sig_aftEff[ipt]->IntegralAndError(LowMBin, HigMBin, data_err, "width");
	//		ckt_value   = H_M_totCKT_Reb4PtPhys[ipt]      ->IntegralAndError(LowMBin, HigMBin, ckt_err,  "width");
	//
	//		double Ratio_Value = data_value/ckt_value;
	//		double Ratio_Error = Ratio_Value*(data_err/data_value);
	//		
	//		HRatioVsPt_IMR->SetBinContent(ipt+1, Ratio_Value);
	//		HRatioVsPt_IMR->SetBinError(  ipt+1, Ratio_Error);
	//			
	//		double ExY_Value = data_value - ckt_value;
	//		double ExY_Error = data_err;
	//
	//		HExYVsPt_IMR->SetBinContent(ipt+1, ExY_Value);
	//		HExYVsPt_IMR->SetBinError(  ipt+1, ExY_Error);
	//	}//ipt
	//HExYVsPt_IMR->Scale(1., "width");
//	//Get Excess Yield vs (mT-M) for the IMR mass region
//	for(int icent=0; icent<nCentBins; icent++)
//	{
//		//----------------------------------------------------------------------
//		//ExY = (Eff. Corr. Data) - CKT
//		//----------------------------------------------------------------------
//		HExY_MvsMt[icent] =  (TH2D*) H2d_MvsMt_Reb_Sig_aftEff[icent]->Clone(Form("HExY_MvsMt_icent%d", icent));
//		HExY_MvsMt[icent] -> Add(H2D_MvsMt_totCKT_Reb[icent], -1);
//		
//		HExY_MvsPt[icent] =  (TH2D*) H2d_MvsPt_Reb_Sig_aftEff[icent]->Clone(Form("HExY_MvsPt_icent%d", icent));
//		HExY_MvsPt[icent] -> Add(H2D_MvsPt_totCKT_Reb[icent], -1);
//
//		//----------------------------------------------------------------------
//		//apply the Acc. Corr on (Data - CKT)
//		//----------------------------------------------------------------------
//		HExY_MvsMt_AftAcc[icent]      // this is the acceptance corrected access yield (Data-CKT)
//			= (TH2D*)calRatioH1toH2(HExY_MvsMt[icent], H2d_MvsMt_PairAcc[icent], Form("HExY_MvsMt_AftAcc_icent%d", icent) );
//		HExY_MvsPt_AftAcc[icent]      // this is the acceptance corrected access yield (Data-CKT)
//			= (TH2D*)calRatioH1toH2(HExY_MvsPt[icent], H2d_MvsPt_PairAcc[icent], Form("HExY_MvsPt_AftAcc_icent%d", icent) );
//
//		//Projection to Yaxis for the dN/dM vs M of the whole mass region, then fit the IMR to extract the Temperature
//		//HExY_M_AftAcc_AllMass[icent] = (TH1D*) HExY_MvsMt_AftAcc[icent] -> ProjectionY( Form("HExY_M_AftAcc_AllMass_icent%d", icent) );
//		//HExY_M_AftAcc_AllMass[icent] ->Scale(1., "width"); //dN/dM vs M of acceptance corrected access yield (Data-CKT)
//		HExY_M_AftAcc_AllMass[icent] = (TH1D*) HExY_MvsPt_AftAcc[icent] -> ProjectionY( Form("HExY_M_AftAcc_AllMass_icent%d", icent) );
//		HExY_M_AftAcc_AllMass[icent] ->Scale(1., "width"); //dN/dM vs M of acceptance corrected access yield (Data-CKT)
//
//		//----------------------------------------------------------------------
//		//Projection to Mt Axis for the IMR region
//		//----------------------------------------------------------------------
//		int    IMR_LowMBin     =  HExY_MvsMt[icent]->GetYaxis()->FindBin(IMRMassWindow[0]);
//		int    IMR_HigMBin     =  HExY_MvsMt[icent]->GetYaxis()->FindBin(IMRMassWindow[1]);
//		
//		//----------------------------------------------------------------------
//		//ExY
//		HExY_Mt_IMR[icent]            = (TH1D*) HExY_MvsMt[icent]        ->ProjectionX(Form("HExY_Mt_IMR_icent%d",icent),        IMR_LowMBin, IMR_HigMBin);
//		
//		//----------------------------------------------------------------------
//		//ExY aft Acc. Corr.
//		HExY_Mt_AftAcc_IMR[icent]     = (TH1D*) HExY_MvsMt_AftAcc[icent] ->ProjectionX(Form("HExY_Mt_AftAcc_IMR_icent%d",icent), IMR_LowMBin, IMR_HigMBin);
//		
//		//----------------------------------------------------------------------
//		//ExY Mt Spectral, dN/dMt/Mt/2Pi
//		HExY_MtSpec_AftAcc_IMR[icent] = (TH1D*) HExY_Mt_AftAcc_IMR[icent]->Clone(Form("HExY_MtSpec_AftAcc_IMR_icent%d",icent));
//		
//		//Rebin for Physics results
//		HExY_MtSpec_AftAcc_Reb_IMR[icent]  = (TH1D*) rebHisto( HExY_MtSpec_AftAcc_IMR[icent], Form("HExY_MtSpec_AftAcc_Reb_IMR_icent%d",icent), nmtBins4Phys, mtBinBDs4Phys, "NO");
//
//		HExY_MtSpec_AftAcc_Reb_IMR[icent] -> Scale(1., "width"); //ExY/dmt
//		
//		for(int ib=1; ib<= HExY_MtSpec_AftAcc_Reb_IMR[icent]->GetNbinsX(); ib++)
//		{
//			double ibContent = HExY_MtSpec_AftAcc_Reb_IMR[icent]->GetBinContent(ib)/HExY_MtSpec_AftAcc_Reb_IMR[icent]->GetBinCenter(ib);
//			double ibError   = HExY_MtSpec_AftAcc_Reb_IMR[icent]->GetBinError(ib)/HExY_MtSpec_AftAcc_Reb_IMR[icent]->GetBinCenter(ib);
//			HExY_MtSpec_AftAcc_Reb_IMR[icent]->SetBinContent(ib, ibContent);
//			HExY_MtSpec_AftAcc_Reb_IMR[icent]->SetBinError(  ib, ibError  );
//		} //ExY/mt
//		//----------------------------------------------------------------------
//	}//icent
//
//	//Get the Data-CKT for whole mass region
//	//Correct the acceptance for (Data-CKT)
//	//Fit the acceptance corrected d(Data-CKT)/dM vs M in the intermedium mass region to extrac the pure temperature
//
//}
////------------------------------------------------------------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------------------------------------------------------------
//void drawEvtPsi()
//{
//	TCanvas *c1 = new TCanvas("c1");
//	gStyle->SetOptStat(0);
//
//	TLatex tl;
//	tl.SetTextSize(0.06);
//	tl.SetNDC();
//
//	//plot events in different centralities
//	//=========================================================================
//	H_cent9->SetTitle("nEvents in 9 Centrality Bins, Before Reweight");
//	H_cent9->SetMaximum(H_cent9->GetMaximum()*1.1);
//	H_cent9->SetXTitle("Centrality Index");
//	H_cent9->SetYTitle("#Events");
//	H_cent9->GetYaxis()->SetTitleSize(0.05);
//	H_cent9->GetYaxis()->SetTitleOffset(0.75);
//	H_cent9->GetXaxis()->SetTitleSize(0.05);
//	H_cent9->GetXaxis()->SetTitleOffset(0.75);
//	H_cent9->SetLineWidth(2);
//	H_cent9->Draw("hist");
//	c1->SaveAs("./outAnaPlots/H_cent9.png");
//	c1->SaveAs("./outAnaPlots/H_cent9.pdf");
//
//	H_cent16->SetTitle("nEvents in 16 Centrality Bins, Before Reweight");
//	H_cent16->SetMaximum(H_cent16->GetMaximum()*1.1);
//	H_cent16->SetXTitle("Centrality Index");
//	H_cent16->SetYTitle("#Events");
//	H_cent16->GetYaxis()->SetTitleSize(0.05);
//	H_cent16->GetYaxis()->SetTitleOffset(0.75);
//	H_cent16->GetXaxis()->SetTitleSize(0.05);
//	H_cent16->GetXaxis()->SetTitleOffset(0.75);
//	H_cent16->SetLineWidth(2);
//	H_cent16->Draw("hist");
//	c1->SaveAs("./outAnaPlots/H_cent16.png");
//	c1->SaveAs("./outAnaPlots/H_cent16.pdf");
//
//	H_cent9_aftwt->SetTitle("nEvents in 9 Centrality Bins, After Reweight");
//	H_cent9_aftwt->SetMaximum(H_cent9_aftwt->GetMaximum()*1.1);
//	H_cent9_aftwt->SetXTitle("Centrality Index");
//	H_cent9_aftwt->SetYTitle("#Events");
//	H_cent9_aftwt->GetYaxis()->SetTitleSize(0.05);
//	H_cent9_aftwt->GetYaxis()->SetTitleOffset(0.75);
//	H_cent9_aftwt->GetXaxis()->SetTitleSize(0.05);
//	H_cent9_aftwt->GetXaxis()->SetTitleOffset(0.75);
//	H_cent9_aftwt->SetLineWidth(2);
//	H_cent9_aftwt->Draw("hist");
//	c1->SaveAs("./outAnaPlots/H_cent9_aftwt.png");
//	c1->SaveAs("./outAnaPlots/H_cent9_aftwt.pdf");
//
//	H_cent16_aftwt->SetTitle("nEvents in 16 Centrality Bins, After Reweight");
//	H_cent16_aftwt->SetMaximum(H_cent16_aftwt->GetMaximum()*1.1);
//	H_cent16_aftwt->SetXTitle("Centrality Index");
//	H_cent16_aftwt->SetYTitle("#Events");
//	H_cent16_aftwt->GetYaxis()->SetTitleSize(0.05);
//	H_cent16_aftwt->GetYaxis()->SetTitleOffset(0.75);
//	H_cent16_aftwt->GetXaxis()->SetTitleSize(0.05);
//	H_cent16_aftwt->GetXaxis()->SetTitleOffset(0.75);
//	H_cent16_aftwt->SetLineWidth(2);
//	H_cent16_aftwt->Draw("hist");
//	c1->SaveAs("./outAnaPlots/H_cent16_aftwt.png");
//	c1->SaveAs("./outAnaPlots/H_cent16_aftwt.pdf");
//
//	//plot event plane information
//	for(int imd=0; imd<nEvtPsiModes; imd++)
//	{
//		HEvtPsi[imd]->SetTitle("");
//		HEvtPsi[imd]->SetAxisRange(2.4e6, 3.8e6, "y");
//		HEvtPsi[imd]->SetAxisRange(0.00, 3.13, "x");
//		//HEvtPsi[imd]->Draw("pe");
//		//c1->SaveAs( Form("outAnaPlots/HEvtPsi_imd%d.png",imd) );
//		//c1->SaveAs( Form("outAnaPlots/HEvtPsi_imd%d.pdf",imd) );
//
//		//overlay them
//		HEvtPsi[imd]->SetLineColor(imd+1);
//		HEvtPsi[imd]->SetLineWidth(2);
//
//		if(imd==0) HEvtPsi[imd]->Draw("hist");
//		else 
//		{
//			HEvtPsi[imd]->Draw("histsame");	
//		}
//	}
//
//	TLegend * leg_evtpsi = new TLegend(0.6,0.6,0.89,0.89);
//	for(int imd=0; imd<nEvtPsiModes; imd++) leg_evtpsi->AddEntry(HEvtPsi[imd], EvtPsiType[imd],  "lp");
//	leg_evtpsi->Draw();
//	c1->SaveAs( "outAnaPlots/HEvtPsi_Comp.png" );
//	c1->SaveAs( "outAnaPlots/HEvtPsi_Comp.pdf" );
//	
//	delete c1;
//	delete leg_evtpsi;
//}
//
////----------------------------------------------------------------------
//void drawSigInCent( )
//{
//	TCanvas *c1 = new TCanvas("c1");
//	gStyle->SetOptStat(0);
//	c1->SetLogy();
//	
//	gPad->SetTopMargin(0.08);
//	gPad->SetBottomMargin(0.11);
//	gPad->SetLeftMargin(0.11);
//	gPad->SetRightMargin(0.05);
//
//	TLatex tl;
//	tl.SetTextSize(0.06);
//	tl.SetNDC();
//
//	for(int icent=0; icent<nCentBins; icent++)
//	{
//		cout<<"H1d_M_Reb_US[icent]->GetBinContent(10): "<<H1d_M_Reb_US[icent]->GetBinContent(10)<<endl;
//		
//		H1d_M_Reb_US[icent]    ->SetTitle(""); 
//		H1d_M_Reb_US[icent]    ->SetYTitle("dN/dM_{ee} (c^{2}/GeV)");
//		H1d_M_Reb_US[icent]    ->SetAxisRange(-0.2,  3.5, "x");
//		H1d_M_Reb_US[icent]    ->SetAxisRange(5.e-8, 80., "y"); 
//		H1d_M_Reb_US[icent]    ->GetYaxis()->SetTitleSize(0.06);
//		H1d_M_Reb_US[icent]    ->GetYaxis()->SetTitleOffset(0.75);
//		H1d_M_Reb_US[icent]    ->GetXaxis()->SetTitleSize(0.05);
//		H1d_M_Reb_US[icent]    ->GetXaxis()->SetTitleOffset(0.85);
//		H1d_M_Reb_US[icent]    ->GetYaxis()->SetNdivisions(4);
//		H1d_M_Reb_US[icent]    ->SetMarkerStyle(24); 
//		H1d_M_Reb_US[icent]    ->SetMarkerColor(1); 
//		H1d_M_Reb_US[icent]    ->SetLineColor(1); 
//		H1d_M_Reb_US[icent]    ->SetLineWidth(2); 
//		H1d_M_Reb_US[icent]    ->Draw("pe"); 
//
//		//H1d_M_Reb_US[icent]    ->SetTitle(""); 
//		//H1d_M_Reb_US[icent]    ->SetYTitle("dN/dM_{ee} (c^{2}/GeV)"); 
//		//H1d_M_Reb_US[icent]    ->SetAxisRange(5.e-8, 80., "y"); 
//		//H1d_M_Reb_US[icent]    ->SetMarkerStyle(24); 
//		//H1d_M_Reb_US[icent]    ->SetMarkerColor(1); 
//		//H1d_M_Reb_US[icent]    ->SetLineColor(1); 
//		//H1d_M_Reb_US[icent]    ->Draw("pe"); 
//		H1d_M_Reb_GmLS[icent]  ->SetLineColor(1);
//		H1d_M_Reb_GmLS[icent]  ->SetLineWidth(2);
//		H1d_M_Reb_GmLS[icent]  ->Draw("histsame");
//		H1d_M_Reb_Sig_Raw[icent] ->SetMarkerStyle(20);
//		H1d_M_Reb_Sig_Raw[icent] ->SetMarkerColor(4);
//		H1d_M_Reb_Sig_Raw[icent] ->SetLineColor(4);
//		H1d_M_Reb_Sig_Raw[icent] ->SetLineWidth(2);
//		H1d_M_Reb_Sig_Raw[icent] ->Draw("pesame");
//
//		TLegend* leg1 = new TLegend(0.70, 0.65, 0.92, 0.90);
//		leg1->SetBorderSize(0);
//		leg1->SetFillColor(0);
//		leg1->AddEntry((TObject *)0,             "Cent: "+centTitle[icent], " " );
//		leg1->AddEntry(H1d_M_Reb_US[icent],      "UnlikeSign",             "lp" );
//		leg1->AddEntry(H1d_M_Reb_GmLS[icent],    "GeomLS",                 "l"  );
//		leg1->AddEntry(H1d_M_Reb_Sig_Raw[icent], "Signal",                 "lp" );
//		leg1->SetTextSize(0.045);
//		leg1->Draw();
//
//
//		tl.SetTextSize(0.050);
//		tl.DrawLatex(0.14, 0.86, "2018 Au+Au #sqrt{s_{NN}} = 27 GeV (MinBias)");
//		tl.SetTextSize(0.040);
//		tl.DrawLatex(0.17, 0.80, "p_{T}^{e}>0.2 GeV/c, |#eta^{e}|<1, |y^{ee}|<1");
//		
//		drawLatex(0.68, 0.55, "Raw",                42, 0.06, 1);
//		
//		c1->SaveAs(Form("outAnaPlots/H1d_M_Reb_Raw_icent%d.png",icent));
//		c1->SaveAs(Form("outAnaPlots/H1d_M_Reb_Raw_icent%d.pdf",icent));
//		
//		//----------------------------------------------------------------------
//		//after the pair sign acc. correction
//		//----------------------------------------------------------------------
//		H1d_M_Reb_US[icent]    ->SetTitle(""); 
//		H1d_M_Reb_US[icent]    ->SetYTitle("dN/dM_{ee} (c^{2}/GeV)");
//		H1d_M_Reb_US[icent]    ->SetAxisRange(-0.2,  3.5, "x");
//		H1d_M_Reb_US[icent]    ->SetAxisRange(5.e-8, 80., "y"); 
//		H1d_M_Reb_US[icent]    ->GetYaxis()->SetTitleSize(0.06);
//		H1d_M_Reb_US[icent]    ->GetYaxis()->SetTitleOffset(0.75);
//		H1d_M_Reb_US[icent]    ->GetXaxis()->SetTitleSize(0.05);
//		H1d_M_Reb_US[icent]    ->GetXaxis()->SetTitleOffset(0.85);
//		H1d_M_Reb_US[icent]    ->GetYaxis()->SetNdivisions(4);
//		H1d_M_Reb_US[icent]    ->SetMarkerStyle(24); 
//		H1d_M_Reb_US[icent]    ->SetMarkerColor(1); 
//		H1d_M_Reb_US[icent]    ->SetLineColor(1); 
//		H1d_M_Reb_US[icent]    ->SetLineWidth(2); 
//		H1d_M_Reb_US[icent]    ->Draw("pe"); 
//		H1d_M_Reb_GmLS_aftPSAC[icent]->SetLineColor(1);
//		H1d_M_Reb_GmLS_aftPSAC[icent]->SetLineWidth(2);
//		H1d_M_Reb_GmLS_aftPSAC[icent]->Draw("histsame");
//
//		H1d_M_Reb_Sig_aftPSAC[icent]   ->SetMarkerStyle(20);
//		H1d_M_Reb_Sig_aftPSAC[icent]   ->SetMarkerColor(4);
//		H1d_M_Reb_Sig_aftPSAC[icent]   ->SetLineColor(4);
//		H1d_M_Reb_Sig_aftPSAC[icent]   ->SetLineWidth(2);
//		H1d_M_Reb_Sig_aftPSAC[icent]   ->Draw("pesame");
//
//		TLegend* leg2 = new TLegend(0.70, 0.65, 0.92, 0.90);
//		leg2->SetBorderSize(0);
//		leg2->SetFillColor(0);
//		leg2->AddEntry((TObject *)0,                  "Cent: "+centTitle[icent], ""  );
//		leg2->AddEntry(H1d_M_Reb_US[icent],           "UnlikeSign",              "lp");
//		leg2->AddEntry(H1d_M_Reb_GmLS_aftPSAC[icent], "GeomLS*PSAC",             "l" );
//		leg2->AddEntry(H1d_M_Reb_Sig_aftPSAC[icent],  "Signal",                  "lp");
//		leg2->SetTextSize(0.041);
//		leg2->Draw();
//
//		tl.SetTextSize(0.050);
//		tl.DrawLatex(0.14, 0.86, "2018 Au+Au #sqrt{s_{NN}} = 27 GeV (MinBias)");
//		tl.SetTextSize(0.040);
//		tl.DrawLatex(0.17, 0.80, "p_{T}^{e}>0.2 GeV/c, |#eta^{e}|<1, |y^{ee}|<1");
//		
//		//drawLatex(0.16, 0.19, "#font[22]{STAR Preliminary}",      42, 0.06, 2);
//		drawLatex(0.68, 0.55, "before eff. corr.",                42, 0.06, 1);
//		
//		c1->SaveAs(Form("outAnaPlots/H1d_M_aftPSAC_icent%d.png",icent));
//		c1->SaveAs(Form("outAnaPlots/H1d_M_aftPSAC_icent%d.pdf",icent));
//		
//		//draw Signal/Background
//		H1d_SoBvsM_Reb[icent]   ->SetTitle(""); 
//		H1d_SoBvsM_Reb[icent]   ->SetYTitle("Signal/geomLS");
//		H1d_SoBvsM_Reb[icent]   ->SetAxisRange(-0.2,  3.5, "x");
//		H1d_SoBvsM_Reb[icent]   ->SetAxisRange(5.e-4, 50., "y"); 
//		H1d_SoBvsM_Reb[icent]   ->GetYaxis()->SetTitleSize(0.06);
//		H1d_SoBvsM_Reb[icent]   ->GetYaxis()->SetTitleOffset(0.75);
//		H1d_SoBvsM_Reb[icent]   ->GetXaxis()->SetTitleSize(0.05);
//		H1d_SoBvsM_Reb[icent]   ->GetXaxis()->SetTitleOffset(0.85);
//		H1d_SoBvsM_Reb[icent]   ->SetMarkerStyle(20);
//		H1d_SoBvsM_Reb[icent]   ->SetMarkerColor(4);
//		H1d_SoBvsM_Reb[icent]   ->SetLineColor(4);
//		H1d_SoBvsM_Reb[icent]   ->SetLineWidth(2);
//		H1d_SoBvsM_Reb[icent]   ->Draw("pe");
//
//		tl.SetTextSize(0.050);
//		tl.DrawLatex(0.14, 0.86, "2018 Au+Au #sqrt{s_{NN}} = 27 GeV Cent: "+centTitle[icent]);
//		tl.SetTextSize(0.040);
//		tl.DrawLatex(0.17, 0.80, "p_{T}^{e}>0.2 GeV/c, |#eta^{e}|<1, |y^{ee}|<1");
//		
//
//		c1->SaveAs(Form("outAnaPlots/H1d_SoBvsM_Reb_icent%d.png",icent));
//		c1->SaveAs(Form("outAnaPlots/H1d_SoBvsM_Reb_icent%d.pdf",icent));
//
//
//		//if(icent==0)
//		//{
//		//	c1->SaveAs("./MakeFigures4QM19/figures4QM/RawSignal_27GeV.pdf");
//		//}
//	}
//	delete c1;
//}
//
//----------------------------------------------------------------------
//----------------------------------------------------------------------
Int_t GetnEvtsInCent(Int_t centBinIdx)
{
	Int_t cent_low = centBinIndex[centBinIdx][0]+2;
	Int_t cent_hig = centBinIndex[centBinIdx][1]+2;

	Int_t nEvts    = H_cent9_aftwt->Integral(cent_low, cent_hig);

	cout<<"centBinIdx"<<endl;
	cout<<"centBinIndex[centBinIdx][0]+2: "<<centBinIndex[centBinIdx][0]+2<<endl;
	cout<<"centBinIndex[centBinIdx][1]+2: "<<centBinIndex[centBinIdx][1]+2<<endl;

	cout<<"icent: "<<centBinIdx<<", number of MinBias Events: "<<nEvts/1.e6<<"M, Go to the final signals "<<endl;

	return nEvts;
}
//----------------------------------------------------------------------
//void saveFiles()
//{
//	TFile *fout = new TFile("outAnaPlots/signals_dielectron_run18_27GeV.root", "recreate");
//	cout<<"outAnaPlots rootfile: "<<fout->GetName()<<endl;
//	fout->cd();
//	//	
//	//	for(int icent=0; icent<nCentBins; icent++)
//	//	{
//	//		H2d_MvsPt_US[icent]       ->Write();
//	//		H2d_MvsPt_LSPos[icent]    ->Write();
//	//		H2d_MvsPt_LSNeg[icent]    ->Write();
//	//		H2D_MvsPt_Sig[icent]      ->Write();
//	//		H1d_M_Reb_Sig_aftPSAC[icent] ->Write();
//	//	}
//	//	
//	//	H1D_PairEff_Reb->Write();
//	//	
//	//	H1d_M_Reb_Sig_aftEff[0]->SetTitle("Run18 Final Mass Spectra");
//	//	H1d_M_Reb_Sig_aftEff[0]->Write();
//	//	H1d_M_Reb_Sig_SysErr->SetName("H1d_M_Reb_Sig_SysErr");
//	//	H1d_M_Reb_Sig_SysErr->SetTitle("H1d_M_Reb_Sig_SysErr");
//	//	H1d_M_Reb_Sig_SysErr->Write();
//	//	
//	//	hSpecM_27Run11->SetTitle("Run11 Final Mass Spectra Data");
//	//	hSpecM_27Run11->SetName("hSpecM_27Run11");
//	//	hSpecM_27Run11->Write();
//	//	
//	//	//H1d_M_Reb_Sig_NoOmegaAndPhi_aftEff->SetName("H1d_M_Reb_Sig_NoOmegaAndPhi_aftEff");
//	//	//H1d_M_Reb_Sig_NoOmegaAndPhi_aftEff->SetTitle("Run18 Final Mass Spectra, no Omega and Phi");
//	//	//H1d_M_Reb_Sig_NoOmegaAndPhi_aftEff->Write();
//	//	
//	//	gCKT->SetName("gCKT");
//	//	gCKT->SetTitle("gCKT");
//	//	gCKT->Write();
//	//
//	//	HRatio_Data2CKT->Write();
//	//	HRatio_Data2CKT_SysErr->Write();
//	//	
//	fout->Close();
//
//	//TFile *outfile = new TFile("outfile_data_Run18_27GeV_aftRmPhoE.root", "recreate");
//	//TFile *outfile = new TFile("outfile_data_Run18_27GeV_4highptpure_aftRmPhoE.root", "recreate");
//	//TFile *outfile = new TFile("outfile_data_Run18_27GeV_4allpure_aftRmPhoE.root", "recreate");
//	//TFile *outfile = new TFile("outfile_data_Run18_27GeV_pt0.3_aftRmPhoE.root", "recreate");
//	//TFile *outfile = new TFile("outfile_data_Run18_27GeV_pt0.4_aftRmPhoE.root", "recreate");
//	TFile *outfile = new TFile("outfile_data_Run18_27GeV_aftRmPhoE.root", "recreate");
//	outfile->cd();
//	for(int icent=0; icent<nCentBins; icent++)
//	{
//		H2d_MvsPt_US[icent]       ->Write();
//		H2d_MvsPt_LSPos[icent]    ->Write();
//		H2d_MvsPt_LSNeg[icent]    ->Write();
//		
//		H1d_SoBvsM_Reb[icent]          ->Write();
//		H1d_M_Reb_Sig_Raw[icent]       ->Write();
//		H1d_M_Reb_Sig_aftPSAC[icent]   ->Write();
//		H1d_M_Reb_Sig_aftEff[icent]    ->Write();
//		H2d_MvsPt_Reb_Sig_Raw[icent]   ->Write();
//		
//		H2d_MvsPt_Reb_Sig_aftEff[icent]->Write();
//		H2D_MvsMt_totCKT_Reb[icent]->Write();
//
//		H1D_M_totCKT[icent]     ->Write();
//		H1D_M_totCKT_Reb[icent] ->Write();
//		
//		for(int ickt=0; ickt<nCKTs+2; ickt++)
//		{
//			H1D_M_CKT[icent][ickt] -> Write();
//		}
//
//		H1D_Ratio2CKT[icent]      ->Write();
//		H2d_MvsPt_Reb_US[icent]   ->Write();
//		H2d_MvsPt_Reb_GmLS[icent] ->Write();
//		H2d_MvsPt_Reb_GmLS_aftPSAC[icent] -> Write();
//		
//		H2d_MvsPt_Reb_PSAC[icent]->Write();
//		H2d_MvsPt_PairEff[icent] ->Write();
//	}
//	for(int icent=0; icent<nCentBins; icent++)
//	{
//		HExY_MvsMt[icent] ->Write();
//		HExY_MvsPt[icent] ->Write();
//		HExY_MvsMt_AftAcc[icent] -> Write();
//		HExY_MvsPt_AftAcc[icent] -> Write();
//	}
//	
//	//write down physics in pt bins
//	//H2d_MvsPt_Reb4PtPhys_Sig_aftEff->Write();
//
//	for(int ipt=0; ipt<nPtBins4Phys; ipt++)
//	{
//		//H1d_M_Reb4PtPhys_Sig_aftEff[ipt]->Write();
//	}
//	for(int ipt=0; ipt<nPtBins4Phys; ipt++)
//	{
//		H_M_totCKT_Reb4PtPhys[ipt]->Write();
//	}
//	for(int ipt=0; ipt<nPtBins4Phys; ipt++)
//	{
//		H_M_totCKT4PtPhys[ipt]->Write();
//	}
//	for(int ipt=0; ipt<nPtBins4Phys; ipt++)
//	{
//		H1D_Ratio2CKT_4PtPhys[ipt] ->Write();
//	}
//	cout<<"tem outfile: "<<outfile->GetName()<<endl;
//	outfile->Close();
//
//	//-----------------------------------------------------------------------------	
//	//-----------------------------------------------------------------------------	
//	TFile *outfile4EXY = new TFile("outfile4EXY_data_Run18_27GeV_aftRmPhoE.root", "recreate");
//	outfile4EXY->cd();
//	
//	H_cent9_aftwt->Write();
//
//	for(int icent=0; icent<nCentBins; icent++)
//	{
//		H2d_MvsPt_Reb_US[icent]      ->Write(); //US
//		H2d_MvsPt_Reb_GmLS[icent]    ->Write(); //GeomLS
//		H2d_MvsPt_Reb_GmLS_aftPSAC[icent] ->Write();
//		H1d_M_Reb_Sig_Raw[icent]     ->Write(); //US-GeomLS
//		H1d_M_Reb_Sig_aftPSAC[icent] ->Write(); //US-GeomLS*PSAC
//		H1d_M_Reb_Sig_aftEff[icent]  ->Write();
//	}
//
//	//
//	//	for(int icent=0; icent<nCentBins; icent++)
//	//	{
//	//		H1D_M_totCKT[icent]            ->Write();
//	//		H1D_M_totCKT_Reb[icent]        ->Write();
//	//
//	//		H2d_MvsPt_Reb_Sig_aftEff[icent]->Write();
//	//		H2D_MvsMt_totCKT_Reb[icent]    ->Write();
//	//
//	//		H1D_Ratio2CKT[icent]           ->Write();
//	//	}
//	//	
//	//	for(int icent=0; icent<nCentBins; icent++)
//	//	{
//	//		HExY_MvsMt[icent] ->Write();
//	//		HExY_MvsPt[icent] ->Write();
//	//		HExY_MvsMt_AftAcc[icent] -> Write();
//	//		HExY_MvsPt_AftAcc[icent] -> Write();
//	//	}
//
//	cout<<"tem outfile4EXY: "<<outfile4EXY->GetName()<<endl;
//	outfile4EXY->Close();
//
//
//
//
//}
////----------------------------------------------------------------------
//void drawPhysicsInCent()
//{
//	TCanvas *c1 = new TCanvas("c1");
//	gStyle->SetOptStat(0);
//	gPad->SetTopMargin(0.08);
//	gPad->SetBottomMargin(0.11);
//	gPad->SetLeftMargin(0.11);
//	gPad->SetRightMargin(0.05);
//	//---------------------------------------------------------------------------
//	//---------------------------------------------------------------------------
//
//	TLatex tl;
//	tl.SetTextSize(0.06);
//	tl.SetNDC();
//	
//	//----------------------------------------------------------------------
//	//Draw Corrected Data vs CKT Figures
//	//----------------------------------------------------------------------
//	TH2D* Htem2d = new TH2D("Htem2d", "", 10, -0.1, 3.6, 10, 1.e-7, 30.);
//	Htem2d->SetTitle(""); 
//	Htem2d->SetYTitle("dN/dM_{ee} (c^{2}/GeV)");
//	Htem2d->SetXTitle("M_{ee} (GeV/c^{2})");
//	//Htem2d->SetYTitle("dN/dM_{ee}/nEvents (c^{2}/GeV)");
//	Htem2d->GetYaxis()->SetTitleSize(0.06);
//	Htem2d->GetYaxis()->SetTitleOffset(0.75);
//	Htem2d->GetXaxis()->SetTitleSize(0.05);
//	Htem2d->GetXaxis()->SetTitleOffset(0.80);
//	Htem2d->SetMarkerStyle(24); 
//	Htem2d->SetMarkerStyle(20);
//	Htem2d->SetMarkerColor(4);
//	Htem2d->SetLineColor(4);
//	Htem2d->SetLineWidth(2);
//	//Htem2d->GetYaxis()->SetNdivisions(4);
//
//	c1->SetLogy();
//	for(int icent=0; icent<nCentBins; icent++)
//	{
//		Htem2d->Draw("");
//	
//		for(int ickt=0; ickt<nCKTs+2; ickt++)
//		{
//			//H1D_M_CKT_Reb[icent][ickt]->SetLineColor(cktColors[ickt]);
//			//H1D_M_CKT_Reb[icent][ickt]->SetLineWidth(2);
//			//H1D_M_CKT_Reb[icent][ickt]->SetLineStyle(7);
//			//H1D_M_CKT_Reb[icent][ickt]->Draw("chistsame");
//			
//			H1D_M_CKT[icent][ickt]->SetLineColor(cktColors[ickt]);
//			H1D_M_CKT[icent][ickt]->SetLineWidth(2);
//			H1D_M_CKT[icent][ickt]->SetLineStyle(7);
//			H1D_M_CKT[icent][ickt]->Draw("chistsame");
//		}//ickt 6+1
//		
//		//	H1D_M_totCKT_Reb[icent]       ->SetLineWidth(2);
//		//	H1D_M_totCKT_Reb[icent]       ->SetLineColor(1);
//		//	H1D_M_totCKT_Reb[icent]       ->Draw("chistsame");
//		H1D_M_totCKT[icent]       ->SetLineWidth(2);
//		H1D_M_totCKT[icent]       ->SetLineColor(1);
//		H1D_M_totCKT[icent]       ->Draw("chistsame");
//		
//		H1d_M_Reb_Sig_aftEff[icent]   ->SetMarkerStyle(20);
//		H1d_M_Reb_Sig_aftEff[icent]   ->SetMarkerSize(0.5);
//		H1d_M_Reb_Sig_aftEff[icent]   ->SetMarkerColor(4);
//		H1d_M_Reb_Sig_aftEff[icent]   ->SetLineColor(4);
//		H1d_M_Reb_Sig_aftEff[icent]   ->Draw("pesame");
//
//		tl.SetTextSize(0.050);
//		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
//		tl.SetTextSize(0.040);
//		tl.DrawLatex(0.15, 0.79, "p_{T}^{e}>0.2 GeV/c, |#eta^{e}|<1, |y^{ee}|<1");
//		//drawLatex(0.60, 0.50, "#font[22]{STAR Preliminary}",      42, 0.06, 2);
//
//		TLegend *leg = new TLegend(0.17,0.65,0.93,0.78);
//		leg->SetNColumns(3);
//		leg->SetFillStyle(0);
//		leg->SetBorderSize(0);
//		leg->SetTextSize(0.035);
//		for(int ickt=0; ickt<nCKTs+2; ickt++)
//		{
//			//leg->AddEntry(H1D_M_CKT_Reb[icent][ickt], cktDecayName[ickt], "l");
//			leg->AddEntry(H1D_M_CKT[icent][ickt], cktDecayName[ickt], "l");
//		}
//		leg->AddEntry(H1D_M_totCKT[icent],  "Cocktail Sum",          "lf");
//		leg->Draw("same");
//
//		c1->SaveAs( Form("./outAnaPlots/Physics/H1D_DataVsCKT_icent%d.png", icent ));
//		c1->SaveAs( Form("./outAnaPlots/Physics/H1D_DataVsCKT_icent%d.pdf", icent ));
//		//H1D_Ratio2CKT[icent]    ->Write();
//		
//		delete leg;
//	}//icent
//
//	//----------------------------------------------------------------------
//	//Draw Ratio Figures
//	//----------------------------------------------------------------------
//	infile_27Rapp     = new TFile("prepareTheoryCurve/forZaochen/27Cocktail_Theory.root", "read");//updated with new CKT
//	gTheory_Rapp_CGA = (TGraph*)infile_27Rapp->Get("gTheroy_Cocktail_Ratio");
//	//infile_27Rapp     = new TFile("../../inputfiles/previous/Theory_Rapp_27.root", "read");
//	//gTheory_Rapp_CGA = (TGraph*)infile_27Rapp->Get("gRCKThe");
//
//	double lowx_Ratio = -0.05;
//	double higx_Ratio = 3.35;
//	TH2D* Htem2d_ratio = new TH2D("Htem2d_ratio", "", 10, lowx_Ratio, higx_Ratio, 10, 0, 15);
//	Htem2d_ratio->SetTitle("");
//	Htem2d_ratio->SetYTitle("Data/Cocktail");
//	Htem2d_ratio->SetXTitle("M_{ee} (GeV/c^{2})");
//	//Htem2d_ratio->SetAxisRange(0.0,   10., "y");
//	Htem2d_ratio->GetYaxis()->SetTitleSize(0.07);
//	Htem2d_ratio->GetYaxis()->SetTitleOffset(0.55);
//	Htem2d_ratio->GetXaxis()->SetTitleSize(0.05);
//	Htem2d_ratio->GetXaxis()->SetTitleOffset(0.75);
//
//	c1->SetLogy(0);
//	for(int icent=0; icent<nCentBins; icent++)
//	{
//		Htem2d_ratio->Draw("");
//		drawLine(lowx_Ratio, 1, higx_Ratio, 1, 2,1,1);
//	
//		gTheory_Rapp_CGA->SetMarkerStyle(0);
//		gTheory_Rapp_CGA->SetLineColor(2);
//		gTheory_Rapp_CGA->SetLineWidth(3);
//		gTheory_Rapp_CGA->Draw("csame");
//
//		H1D_Ratio2CKT[icent]->SetMarkerStyle(20);
//		H1D_Ratio2CKT[icent]->SetMarkerColor(4);
//		H1D_Ratio2CKT[icent]->SetMarkerSize(0.8);
//		H1D_Ratio2CKT[icent]->SetLineColor(4);
//		H1D_Ratio2CKT[icent]->SetLineWidth(2);
//		H1D_Ratio2CKT[icent] ->Draw("pesame");
//		
//		TLegend* leg_Ratio = new TLegend(0.14, 0.70, 0.45, 0.82);
//		leg_Ratio->SetBorderSize(0);
//		leg_Ratio->SetFillColor(0);
//		//leg_Ratio->AddEntry( HRatio_Data2CKT,    "2018",                 "lp" );
//		leg_Ratio->AddEntry( H1D_Ratio2CKT[icent], "Ratio (2018)",       "lp" );
//		leg_Ratio->AddEntry( gTheory_Rapp_CGA,     "Rapp/(Rapp+CKTSum) (0-80%)",           "lp" );
//		leg_Ratio->SetTextSize(0.045);
//		leg_Ratio->Draw();
//
//		tl.SetTextSize(0.050);
//		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
//
//
//		c1->SaveAs( Form("./outAnaPlots/Physics/H1D_Ratio2CKT_icent%d.png", icent ));
//		c1->SaveAs( Form("./outAnaPlots/Physics/H1D_Ratio2CKT_icent%d.pdf", icent ));
//
//		delete leg_Ratio;
//	}//icent
//
//	//----------------------------------------------------------------------
//	//compare Ratio between centralities, 0-10%, 10-40%, 40-80%
//	//----------------------------------------------------------------------
//	TH2D* Htem2d_ratio2 = new TH2D("Htem2d_ratio2", "", 10, lowx_Ratio, higx_Ratio, 10, 0, 20);
//	Htem2d_ratio2->SetTitle("");
//	Htem2d_ratio2->SetYTitle("Data/Cocktail");
//	Htem2d_ratio2->SetXTitle("M_{ee} (GeV/c^{2})");
//	Htem2d_ratio2->GetYaxis()->SetTitleSize(0.07);
//	Htem2d_ratio2->GetYaxis()->SetTitleOffset(0.55);
//	Htem2d_ratio2->GetXaxis()->SetTitleSize(0.05);
//	Htem2d_ratio2->GetXaxis()->SetTitleOffset(0.75);
//
//	TLegend* leg_Ratio2 = new TLegend(0.14, 0.60, 0.45, 0.82);
//	leg_Ratio2->SetBorderSize(0);
//	leg_Ratio2->SetFillColor(0);
//
//	c1->SetLogy(0);
//	Htem2d_ratio2->Draw("");
//	for(int icent=1; icent<nCentBins; icent++)
//	{
//		if(icent==2) continue;
//		drawLine(lowx_Ratio, 1, higx_Ratio, 1, 2,1,1);
//
//		H1D_Ratio2CKT[icent]->SetMarkerStyle(20);
//		H1D_Ratio2CKT[icent]->SetMarkerSize(0.7);
//		H1D_Ratio2CKT[icent]->SetMarkerColor(cktColors[icent]);
//		H1D_Ratio2CKT[icent]->SetLineColor(cktColors[icent]);
//		H1D_Ratio2CKT[icent]->SetLineWidth(2);
//		H1D_Ratio2CKT[icent] ->Draw("pesame");
//		
//		leg_Ratio2->AddEntry( H1D_Ratio2CKT[icent], "Cent: "+centTitle[icent],       "lp" );
//	}//icent:1,3
//	
//	leg_Ratio2->Draw("same");
//		
//	tl.SetTextSize(0.050);
//	tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV");
//
//	c1->SaveAs("./outAnaPlots/Physics/H1D_Ratio2CKT_compBTWCents.png");
//	c1->SaveAs("./outAnaPlots/Physics/H1D_Ratio2CKT_compBTWCents.pdf");
//
//
//	//----------------------------------------------------------------------
//	//compare Ratio between centralities, 0-80%, 0-10%, 10-40%, 40-80% to check the low mass bump structure
//	//----------------------------------------------------------------------
//	TH2D* Htem2d_ratio3 = new TH2D("Htem2d_ratio3", "", 10, lowx_Ratio, 0.4, 10, 0, 3);
//	Htem2d_ratio3->SetTitle("");
//	Htem2d_ratio3->SetYTitle("Data/Cocktail");
//	Htem2d_ratio3->SetXTitle("M_{ee} (GeV/c^{2})");
//	Htem2d_ratio3->GetYaxis()->SetTitleSize(0.07);
//	Htem2d_ratio3->GetYaxis()->SetTitleOffset(0.55);
//	Htem2d_ratio3->GetXaxis()->SetTitleSize(0.05);
//	Htem2d_ratio3->GetXaxis()->SetTitleOffset(0.75);
//
//	TLegend* leg_Ratio3 = new TLegend(0.14, 0.65, 0.45, 0.89);
//	leg_Ratio3->SetBorderSize(0);
//	leg_Ratio3->SetFillColor(0);
//
//	c1->SetLogy(0);
//	Htem2d_ratio3->Draw("");
//	for(int icent=0; icent<nCentBins; icent++)
//	{
//		drawLine(lowx_Ratio, 1, 0.4, 1, 2,1,1);
//
//		H1D_Ratio2CKT[icent]->SetMarkerStyle(20);
//		H1D_Ratio2CKT[icent]->SetMarkerSize(0.7);
//		H1D_Ratio2CKT[icent]->SetMarkerColor(cktColors[icent]);
//		H1D_Ratio2CKT[icent]->SetLineColor(cktColors[icent]);
//		H1D_Ratio2CKT[icent]->SetLineWidth(2);
//		H1D_Ratio2CKT[icent] ->Draw("pesame");
//		
//		leg_Ratio3->AddEntry( H1D_Ratio2CKT[icent], "Cent: "+centTitle[icent],       "lp" );
//	}//icent
//	
//	leg_Ratio3->Draw("same");
//
//	system("mkdir -p outAnaPlots/plots_Bump");
//	c1->SaveAs("./outAnaPlots/plots_Bump/H1D_Ratio2CKT_compBTWCents_CheckBump.png");
//	c1->SaveAs("./outAnaPlots/plots_Bump/H1D_Ratio2CKT_compBTWCents_CheckBump.pdf");
//
//	//----------------------------------------------------------------------
//	//Draw Rho Mass Shape
//	//----------------------------------------------------------------------
//	double lowx_Rho = 0.3;
//	double higx_Rho = 1.1;
//	for(int icent=0; icent<nCentBins; icent++)
//	{
//		TH1D* H1d_M_Rho = (TH1D*) H1d_DataMinusCKT[icent]->Clone( Form("H1d_M_Rho_icent%d", icent) );
//		H1d_M_Rho->SetYTitle("dN/dM_{ee} (c^{2}/GeV) x 10^{3}");
//		H1d_M_Rho->SetXTitle("M_{ee} (GeV/c^{2})");
//		H1d_M_Rho->SetTitle("#rho Mass Shape in "+centTitle[icent]);
//		H1d_M_Rho->GetYaxis()->SetTitleSize(0.06);
//		H1d_M_Rho->GetYaxis()->SetTitleOffset(0.70);
//		H1d_M_Rho->GetXaxis()->SetTitleSize(0.05);
//		H1d_M_Rho->GetXaxis()->SetTitleOffset(0.80);
//		H1d_M_Rho->SetAxisRange(lowx_Rho, higx_Rho, "x");
//		H1d_M_Rho->SetMarkerStyle(20);
//		H1d_M_Rho->SetMarkerColor(4);
//		H1d_M_Rho->SetLineColor(4);
//		H1d_M_Rho->SetLineWidth(2);
//		H1d_M_Rho->Scale(1.e3);
//		H1d_M_Rho->Draw("pe");
//		H1d_M_Rho->Draw("chistsame");
//		drawLine(lowx_Rho, 0, higx_Rho, 0, 2,1,1);
//
//		tl.SetTextSize(0.050);
//		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV");
//		tl.DrawLatex(0.14, 0.78, "Data - Cocktail");
//		
//		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Rho_icent%d.png", icent) );
//		c1->SaveAs( Form("./outAnaPlots/Physics/H1d_M_Rho_icent%d.pdf", icent) );
//	}//icent
//	
//	//----------------------------------------------------------------------
//	//Data vs CKTs in Mt Axis
//	//----------------------------------------------------------------------
//	TH2D* Htem2d4Mt = new TH2D("Htem2d4Mt", "", 10, 0.0, 1.8, 10, 1.e-8, 1.);
//	Htem2d4Mt->SetTitle(""); 
//	Htem2d4Mt->SetYTitle("dN/dM_{T} (c^{2}/GeV)");
//	Htem2d4Mt->SetXTitle("M_{T}-M (GeV/c^{2})");
//	//Htem2d4Mt->SetYTitle("dN/dM_{ee}/nEvents (c^{2}/GeV)");
//	Htem2d4Mt->GetYaxis()->SetTitleSize(0.06);
//	Htem2d4Mt->GetYaxis()->SetTitleOffset(0.75);
//	Htem2d4Mt->GetXaxis()->SetTitleSize(0.05);
//	Htem2d4Mt->GetXaxis()->SetTitleOffset(0.80);
//	//Htem2d4Mt->GetYaxis()->SetNdivisions(4);
//
//	c1->SetLogy(1);
//	for(int icent=0; icent<nCentBins; icent++)
//	{
//		Htem2d4Mt->Draw("");
//	
//		for(int ickt=0; ickt<nCKTs+2; ickt++)
//		{
//			H1d_Mt_CKT_Reb_IMR[icent][ickt]->SetLineColor(cktColors[ickt]);
//			H1d_Mt_CKT_Reb_IMR[icent][ickt]->SetLineWidth(2);
//			H1d_Mt_CKT_Reb_IMR[icent][ickt]->SetLineStyle(7);
//			H1d_Mt_CKT_Reb_IMR[icent][ickt]->Draw("chistsame");
//		}//ickt 6+2
//		
//		H1d_Mt_totCKT_Reb_IMR[icent]  ->SetLineWidth(2);
//		H1d_Mt_totCKT_Reb_IMR[icent]  ->SetLineColor(1);
//		H1d_Mt_totCKT_Reb_IMR[icent]  ->Draw("histsame");
//		//H1d_Mt_totCKT_Reb_IMR[icent]  ->Draw("chistsame");
//
//		H1d_Mt_Reb_AftEff_IMR[icent]   ->SetMarkerStyle(20);
//		H1d_Mt_Reb_AftEff_IMR[icent]   ->SetMarkerSize(0.5);
//		H1d_Mt_Reb_AftEff_IMR[icent]   ->SetMarkerColor(4);
//		H1d_Mt_Reb_AftEff_IMR[icent]   ->SetLineColor(4);
//		H1d_Mt_Reb_AftEff_IMR[icent]   ->Draw("pesame");
//
//		cout<<H1d_Mt_Reb_AftEff_IMR[icent]->GetEntries()<<endl;
//		cout<<"H1d_Mt_Reb_AftEff_IMR[icent]->GetBinContent(3): "<<H1d_Mt_Reb_AftEff_IMR[icent]->GetBinContent(3)<<endl;
//
//		tl.SetTextSize(0.050);
//		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
//		tl.DrawLatex(0.65, 0.86, Form("%.1f<M_{ee}<%.1f GeV/c^{2}", IMRMassWindow[0], IMRMassWindow[1]));
//		tl.SetTextSize(0.040);
//		tl.DrawLatex(0.15, 0.81, "p_{T}^{e}>0.2 GeV/c, |#eta^{e}|<1, |y^{ee}|<1");
//
//		TLegend *leg = new TLegend(0.14,0.65,0.93,0.78);
//		leg->SetNColumns(3);
//		leg->SetFillStyle(0);
//		leg->SetBorderSize(0);
//		leg->SetTextSize(0.035);
//		for(int ickt=0; ickt<nCKTs+2; ickt++)
//		{
//			leg->AddEntry(H1d_Mt_CKT_Reb_IMR[icent][ickt], cktDecayName[ickt],      "l" );
//		}
//		leg->AddEntry(H1d_Mt_totCKT_Reb_IMR[icent],        "Cocktail Sum",          "lf");
//		leg->AddEntry(H1d_Mt_Reb_AftEff_IMR[icent],        "Data",                  "lf");
//		leg->Draw("same");
//
//		c1->SaveAs( Form("./outAnaPlots/Physics/H1D_Mt_DataVsCKT_IMR_icent%d.png", icent ));
//		c1->SaveAs( Form("./outAnaPlots/Physics/H1D_Mt_DataVsCKT_IMR_icent%d.pdf", icent ));
//		
//		delete leg;
//	}//icent
//
//	//Now, Let's Draw the (mt-m) and Extract the Temperature
//	TH2D* Htem2d_4MtSpec = new TH2D("Htem2d_4MtSpec", "", 10, 0.0, 1.8, 10, 1.e-7, 1.e-0);
//	Htem2d_4MtSpec->SetTitle("");
//	Htem2d_4MtSpec->SetYTitle("1/m_{T}dN/dm_{T} ((c^{2}/GeV)^2)");
//	Htem2d_4MtSpec->SetXTitle("m_{T}-M (GeV/c^{2})");
//	Htem2d_4MtSpec->GetYaxis()->SetTitleSize(0.06);
//	Htem2d_4MtSpec->GetYaxis()->SetTitleOffset(0.75);
//	Htem2d_4MtSpec->GetXaxis()->SetTitleSize(0.05);
//	Htem2d_4MtSpec->GetXaxis()->SetTitleOffset(0.85);
//
//	c1->SetLogy(1);
//
//	const double lowM4MtSpec = 0.0;
//	const double higM4MtSpec = 1.8;
//	TF1* fit_MtSpec = new TF1("fit_MtSpec", fun_MtSpec, lowM4MtSpec, higM4MtSpec, 2);
//	fit_MtSpec->SetLineColor(4);
//	fit_MtSpec->SetLineWidth(2);
//	fit_MtSpec->SetParameters(1.e-3, 200.);
//
//	TH1D* HExY_MtSpec_AftAcc_Reb_IMR_rmNegPoint[nCentBins];
//
//	for(int icent=0; icent<nCentBins; icent++)
//	{
//		Htem2d_4MtSpec->Draw("");
//		HExY_MtSpec_AftAcc_Reb_IMR[icent]->SetMarkerStyle(20);
//		HExY_MtSpec_AftAcc_Reb_IMR[icent]->SetMarkerSize(2);
//		HExY_MtSpec_AftAcc_Reb_IMR[icent]->SetMarkerColor(2);
//		HExY_MtSpec_AftAcc_Reb_IMR[icent]->SetLineColor(2);
//		HExY_MtSpec_AftAcc_Reb_IMR[icent]->SetLineWidth(2);
//		HExY_MtSpec_AftAcc_Reb_IMR[icent]->Draw("pesame");
//
//		fit_MtSpec->SetLineStyle(1);
//		HExY_MtSpec_AftAcc_Reb_IMR[icent]->Fit("fit_MtSpec", "EMR0", "", lowM4MtSpec, higM4MtSpec);
//		fit_MtSpec->Draw("lsame");
//
//		double Teff_Value = fit_MtSpec->GetParameter(1)*1000.; //here 1000 is the the unit GeV to MeV
//		double Teff_Error = fit_MtSpec->GetParError(1)*1000.;
//		
//		tl.SetTextSize(0.050);
//		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
//		tl.SetTextSize(0.040);
//		tl.DrawLatex(0.15, 0.80, Form("%.1f<M_{ee}<%.1f GeV/c^{2}", IMRMassWindow[0], IMRMassWindow[1]));
//		tl.DrawLatex(0.60, 0.65, Form("T_{eff} = %.2f #pm %.2f", Teff_Value, Teff_Error));
//
//		TLegend* leg_4MtSpec = new TLegend(0.60, 0.75, 0.89, 0.80);
//		leg_4MtSpec->SetBorderSize(0);
//		leg_4MtSpec->SetFillColor(0);
//		leg_4MtSpec->AddEntry( fit_MtSpec, "A*exp(-m_{T}/T_{eff}) ",       "lp" );
//		leg_4MtSpec->Draw("same");
//		c1->SaveAs( Form("./outAnaPlots/Physics/HExY_MtSpec_AftAcc_Reb_IMR_icent%d.png", icent) );
//		c1->SaveAs( Form("./outAnaPlots/Physics/HExY_MtSpec_AftAcc_Reb_IMR_icent%d.pdf", icent) );
//		//c1->SaveAs( Form("./outAnaPlots/Physics/HExY_MtSpec_AftAcc_Reb_IMR_icent%d.pdf", icent) );
//
//
//		//remove the negative bins
//		HExY_MtSpec_AftAcc_Reb_IMR_rmNegPoint[icent] = (TH1D*)removeNegPoint(HExY_MtSpec_AftAcc_Reb_IMR[icent], Form("HExY_MtSpec_AftAcc_Reb_IMR_rmNegPoint_icent%d", icent));
//		Htem2d_4MtSpec->Draw("");
//		HExY_MtSpec_AftAcc_Reb_IMR_rmNegPoint[icent]->SetMarkerStyle(24);
//		HExY_MtSpec_AftAcc_Reb_IMR_rmNegPoint[icent]->Draw("pesame");
//
//		fit_MtSpec->SetLineStyle(5);
//		HExY_MtSpec_AftAcc_Reb_IMR_rmNegPoint[icent]->Fit("fit_MtSpec", "EMR", "", lowM4MtSpec, higM4MtSpec);
//		fit_MtSpec->Draw("lsame");
//
//		Teff_Value = fit_MtSpec->GetParameter(1)*1000.; //here 1000 is the the unit GeV to MeV
//		Teff_Error = fit_MtSpec->GetParError(1)*1000.;
//		
//		tl.SetTextSize(0.050);
//		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
//		tl.SetTextSize(0.040);
//		tl.DrawLatex(0.15, 0.80, Form("%.1f<M_{ee}<%.1f GeV/c^{2}", IMRMassWindow[0], IMRMassWindow[1]));
//		tl.DrawLatex(0.15, 0.75, "Reject Negative Points");
//		tl.DrawLatex(0.60, 0.65, Form("T_{eff} = %.2f #pm %.2f", Teff_Value, Teff_Error));
//
//		TLegend* leg_4MtSpec_rmNegPoint = new TLegend(0.60, 0.75, 0.89, 0.80);
//		leg_4MtSpec_rmNegPoint->SetBorderSize(0);
//		leg_4MtSpec_rmNegPoint->SetFillColor(0);
//		leg_4MtSpec_rmNegPoint->AddEntry( fit_MtSpec, "A*exp(-m_{T}/T_{eff}) ",       "lp" );
//		leg_4MtSpec_rmNegPoint->Draw("same");
//
//		c1->SaveAs( Form("./outAnaPlots/Physics/HExY_MtSpec_AftAcc_Reb_IMR_rmNegPoint_icent%d.png", icent) );
//		c1->SaveAs( Form("./outAnaPlots/Physics/HExY_MtSpec_AftAcc_Reb_IMR_rmNegPoint_icent%d.pdf", icent) );
//
//
//
//
//		delete leg_4MtSpec;
//		delete leg_4MtSpec_rmNegPoint;
//	}//icent
//
//	//----------------------------------------------------------------------
//	//Draw Acceptance corrected Data-CKT to extract Temperature (dN/dM vs M) 
//	//----------------------------------------------------------------------
//	TH2D* Htem2d_4AccedExY = new TH2D("Htem2d_4AccedExY", "", 10, 0.0, 3.0, 10, 1.e-7, 5.e-1);
//	Htem2d_4AccedExY->SetTitle("");
//	Htem2d_4AccedExY->SetYTitle("dN/dM (Data-CKT) after Acc. Corr");
//	Htem2d_4AccedExY->SetXTitle("M_{ee} (GeV/c^{2})");
//	Htem2d_4AccedExY->GetYaxis()->SetTitleSize(0.06);
//	Htem2d_4AccedExY->GetYaxis()->SetTitleOffset(0.75);
//	Htem2d_4AccedExY->GetXaxis()->SetTitleSize(0.05);
//	Htem2d_4AccedExY->GetXaxis()->SetTitleOffset(0.85);
//
//	//const double lowM4AccedExY = 1.26;
//	//const double higM4AccedExY = 2.60;
//	const double lowM4AccedExY = 1.20;
//	const double higM4AccedExY = 3.0;
//	TF1* fit_AccedExY = new TF1("fit_AccedExY", fun_AccedExY, lowM4AccedExY, higM4AccedExY, 2);
//	fit_AccedExY->SetLineColor(4);
//	fit_AccedExY->SetLineWidth(2);
//	fit_AccedExY->SetParameters(1.e-2, 250.);
//
//	//read in Rapp Predictions
//	TFile* inf_theory    = new TFile("./inputfiles/Histograms4Yi.root", "read");
//	TGraph* ge_RappSum   = (TGraph*) inf_theory->Get("gRappFull");
//	TGraph* ge_RappHGMed = (TGraph*) inf_theory->Get("gRappHGMedFull");
//	TGraph* ge_RappQGP   = (TGraph*) inf_theory->Get("gRappQGPFull");
//
//	ge_RappSum  ->SetLineWidth(2);
//	ge_RappHGMed->SetLineWidth(2);
//	ge_RappQGP  ->SetLineWidth(2);
//	ge_RappSum  ->SetLineColor(1);
//	ge_RappHGMed->SetLineColor(8);
//	ge_RappQGP  ->SetLineColor(2);
//
//	c1->SetLogy(1);
//	for(int icent=0; icent<nCentBins; icent++)
//	{
//		Htem2d_4AccedExY->Draw();
//
//		HExY_M_AftAcc_AllMass[icent]->SetMarkerStyle(20);
//		HExY_M_AftAcc_AllMass[icent]->SetMarkerSize(2);
//		HExY_M_AftAcc_AllMass[icent]->SetMarkerColor(2);
//		HExY_M_AftAcc_AllMass[icent]->SetLineColor(2);
//		HExY_M_AftAcc_AllMass[icent]->SetLineWidth(2);
//		HExY_M_AftAcc_AllMass[icent]->Draw("pesame");
//
//		ge_RappSum  ->Draw("lsame");
//		ge_RappHGMed->Draw("lsame");
//		ge_RappQGP  ->Draw("lsame");
//		
//		HExY_M_AftAcc_AllMass[icent]->Fit("fit_AccedExY", "ER", "", lowM4AccedExY, higM4AccedExY);
//
//		fit_AccedExY->Draw("lsame");
//
//		double T_Value = fit_AccedExY->GetParameter(1)*1000.; //here 1000 is the the unit GeV to MeV
//		double T_Error = fit_AccedExY->GetParError(1)*1000.;
//		
//		tl.SetTextSize(0.050);
//		tl.SetTextColor(1);
//		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[icent]+")");
//		tl.SetTextSize(0.040);
//		tl.SetTextSize(0.050);
//		tl.SetTextColor(1);
//		tl.DrawLatex(0.60, 0.80, Form("%.2f<M_{ee}<%.2f GeV/c^{2}", lowM4AccedExY, higM4AccedExY));
//		tl.DrawLatex(0.60, 0.72, "fit by: M_{ee}^{3/2}*exp(-M_{ee}/T)");
//		tl.SetTextColor(2);
//		tl.DrawLatex(0.60, 0.65, Form("T = %.1f #pm %.1f MeV", T_Value, T_Error));
//
//		c1->SaveAs( Form("./outAnaPlots/Physics/HExY_M_AftAcc_AllMass_icent%d.png", icent) );
//		c1->SaveAs( Form("./outAnaPlots/Physics/HExY_M_AftAcc_AllMass_icent%d.pdf", icent) );
//	}
//	delete Htem2d;
//	delete Htem2d_ratio;
//	delete Htem2d_ratio2;
//	delete Htem2d_4MtSpec;
//	delete Htem2d4Mt;
//	delete Htem2d_4AccedExY;
//	delete c1;
//}
////----------------------------------------------------------------------
////----------------------------------------------------------------------
//void drawPhysicsInPt()
//{
//	TCanvas *c1 = new TCanvas("c1");
//	gStyle->SetOptStat(0);
//	gPad->SetTopMargin(0.08);
//	gPad->SetBottomMargin(0.11);
//	gPad->SetLeftMargin(0.11);
//	gPad->SetRightMargin(0.05);
//	//---------------------------------------------------------------------------
//	//---------------------------------------------------------------------------
//
//	TLatex tl;
//	tl.SetTextSize(0.06);
//	tl.SetNDC();
//	
//	//--------------------------------------------------------------------------------------------------------------
//	//Draw Cocktail in different pt bins
//	//--------------------------------------------------------------------------------------------------------------
//	TH2D* Htem2d2;	//Htem2d2->GetYaxis()->SetNdivisions(4);
//
//	//const int cktColors[nCKTs+2] = {1,6,8,12,15,20,95};
//	
//	c1->SetLogy();
//	for(int ipt=0; ipt<nPtBins4Phys; ipt++)
//	{
//		cout<<"Draw cocktail in "<<Form("%.2f<pt<%.2f", PtBDs4Phys[ipt], PtBDs4Phys[ipt+1])<<endl;
//		//double yMax = H_M_totCKT4PtPhys[ipt]->GetMaximum()*50.;
//		//double yMin = yMax*5.e-9;
//
//		//if(ipt==0)
//		//{
//		//	yMax *= 100.;
//		//	yMin *= 100.;
//		//}
//		//else if(ipt==3)
//		//{
//		//	yMax *= 0.1;
//		//	yMin *= 10;
//		//}
//		////Htem2d2 = new TH2D("Htem2d2", "", 10, -0.1, 3.6, 10, yMin, yMax);
//
//		Htem2d2 = new TH2D("Htem2d2", "", 10, -0.1, 3.6, 10, 1.e-8, 20.);
//		Htem2d2->SetTitle(Form("Dielectron in %.2f<p_{T}<%.2f GeV/c", PtBDs4Phys[ipt], PtBDs4Phys[ipt+1])); 
//		Htem2d2->SetYTitle("dN/dM_{ee} (c^{2}/GeV)");
//		Htem2d2->SetXTitle("M_{ee} (GeV/c^{2})");
//		Htem2d2->GetYaxis()->SetTitleSize(0.06);
//		Htem2d2->GetYaxis()->SetTitleOffset(0.75);
//		Htem2d2->GetXaxis()->SetTitleSize(0.05);
//		Htem2d2->GetXaxis()->SetTitleOffset(0.80);
//		Htem2d2->SetMarkerStyle(24); 
//		Htem2d2->SetMarkerStyle(20);
//		Htem2d2->SetMarkerColor(4);
//		Htem2d2->SetLineColor(4);
//		Htem2d2->SetLineWidth(2);
//		Htem2d2 -> Draw("");
//
//		for(int ickt=0; ickt<nCKTs+2; ickt++)
//		{
//			H_M_CKT4PtPhys[ickt][ipt]->SetLineColor(cktColors[ickt]);
//			H_M_CKT4PtPhys[ickt][ipt]->SetLineWidth(2);
//			H_M_CKT4PtPhys[ickt][ipt]->SetLineStyle(7);
//			H_M_CKT4PtPhys[ickt][ipt]->Draw("chistsame");
//		}//ickt
//
//		H_M_totCKT4PtPhys[ipt]->SetLineWidth(2);;
//		H_M_totCKT4PtPhys[ipt]->SetLineColor(1);;
//		H_M_totCKT4PtPhys[ipt]-> Draw("chistsame");
//		
//		//H1d_M_Reb4PtPhys_Sig_aftEff[ipt]->SetMarkerStyle(20);
//		//H1d_M_Reb4PtPhys_Sig_aftEff[ipt]->SetMarkerSize(0.5);
//		//H1d_M_Reb4PtPhys_Sig_aftEff[ipt]->SetMarkerColor(4);
//		//H1d_M_Reb4PtPhys_Sig_aftEff[ipt]->SetLineColor(4);
//		//H1d_M_Reb4PtPhys_Sig_aftEff[ipt]->Draw("pesame");
//
//		TLegend *leg = new TLegend(0.16,0.62,0.92,0.76);
//		leg->SetNColumns(3);
//		leg->SetFillStyle(0);
//		leg->SetBorderSize(0);
//		leg->SetTextSize(0.035);
//		for(int ickt=0; ickt<nCKTs+2; ickt++)
//		{
//			leg->AddEntry(H_M_CKT4PtPhys[ickt][ipt], cktDecayName[ickt], "l");
//		}
//		leg->AddEntry(H_M_totCKT4PtPhys[ipt],        "Cocktail Sum",     "l");
//		leg->Draw("same");
//
//		tl.SetTextSize(0.050);
//		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV (0-80%)");
//		tl.SetTextSize(0.040);
//		tl.DrawLatex(0.15, 0.79, "p_{T}^{e}>0.2 GeV/c, |#eta^{e}|<1, |y^{ee}|<1");
//
//		c1->SaveAs(Form("./outAnaPlots/Physics/H1D_DataVsCKT_Cent080_ipt%d.png",ipt));
//		c1->SaveAs(Form("./outAnaPlots/Physics/H1D_DataVsCKT_Cent080_ipt%d.pdf",ipt));
//
//		delete Htem2d2;
//		delete leg;
//	}//ipt4phys
//	
//	
//	//Draw Ratio of Data/CKT
//	double lowx_Ratio = 0.05;
//	//double lowx_Ratio = -0.05;
//	double higx_Ratio = 3.35;
//	TH2D* Htem2d_ratio;
//	c1->SetLogy(0);
//	for(int ipt=0; ipt<nPtBins4Phys; ipt++)
//	{
//		if(ipt==5)Htem2d_ratio = new TH2D("Htem2d_ratio", "", 10, lowx_Ratio, higx_Ratio, 10, -10, 30);
//		else      Htem2d_ratio = new TH2D("Htem2d_ratio", "", 10, lowx_Ratio, higx_Ratio, 10, -10, 30);
//		//if(ipt==5)Htem2d_ratio = new TH2D("Htem2d_ratio", "", 10, lowx_Ratio, higx_Ratio, 10, -50, 600);
//		//else      Htem2d_ratio = new TH2D("Htem2d_ratio", "", 10, lowx_Ratio, higx_Ratio, 10, -10, 100);
//		Htem2d_ratio->SetTitle("");
//		Htem2d_ratio->SetYTitle("Data/Cocktail");
//		Htem2d_ratio->SetXTitle("M_{ee} (GeV/c^{2})");
//		Htem2d_ratio->GetYaxis()->SetTitleSize(0.07);
//		Htem2d_ratio->GetYaxis()->SetTitleOffset(0.55);
//		Htem2d_ratio->GetXaxis()->SetTitleSize(0.05);
//		Htem2d_ratio->GetXaxis()->SetTitleOffset(0.75);
//
//		Htem2d_ratio->SetTitle(Form("Ratio in %.2f<p_{T}<%.2f GeV/c", PtBDs4Phys[ipt], PtBDs4Phys[ipt+1])); 
//		Htem2d_ratio->Draw("");
//		drawLine(lowx_Ratio, 1, higx_Ratio, 1, 2,1,1);
//
//		H1D_Ratio2CKT_4PtPhys[ipt]->SetMarkerStyle(20);
//		H1D_Ratio2CKT_4PtPhys[ipt]->SetMarkerColor(4);
//		H1D_Ratio2CKT_4PtPhys[ipt]->SetMarkerSize(0.5);
//		H1D_Ratio2CKT_4PtPhys[ipt]->SetLineColor(4);
//		H1D_Ratio2CKT_4PtPhys[ipt]->SetLineWidth(2);
//		H1D_Ratio2CKT_4PtPhys[ipt] ->Draw("pesame");
//		H1D_Ratio2CKT_4PtPhys[ipt] ->Draw("chistsame");
//		
//		TLegend* leg_Ratio = new TLegend(0.14, 0.70, 0.45, 0.82);
//		leg_Ratio->SetBorderSize(0);
//		leg_Ratio->SetFillColor(0);
//		leg_Ratio->AddEntry( H1D_Ratio2CKT_4PtPhys[ipt], "Ratio (2018)",       "lp" );
//		leg_Ratio->SetTextSize(0.045);
//		//leg_Ratio->Draw();
//
//		tl.SetTextSize(0.050);
//		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[0]+")");
//
//		c1->SaveAs( Form("./outAnaPlots/Physics/H1D_Ratio2CKT_4PtPhys_ipt%d.png", ipt ));
//		c1->SaveAs( Form("./outAnaPlots/Physics/H1D_Ratio2CKT_4PtPhys_ipt%d.pdf", ipt ));
//
//		delete leg_Ratio;
//		delete Htem2d_ratio;
//	}//ipt
//
//	//Draw Ratio as a function of Pt
//	HRatioVsPt_IMR->SetTitle("Ratio of Data/CKT in 0-80% vs. p_{T}");
//	HRatioVsPt_IMR->SetYTitle("Data/Cocktail");
//	HRatioVsPt_IMR->SetXTitle("p_{T} (GeV/c)");
//	HRatioVsPt_IMR->GetYaxis()->SetTitleSize(0.07);
//	HRatioVsPt_IMR->GetYaxis()->SetTitleOffset(0.55);
//	HRatioVsPt_IMR->GetXaxis()->SetTitleSize(0.05);
//	HRatioVsPt_IMR->GetXaxis()->SetTitleOffset(0.85);
//	HRatioVsPt_IMR->SetAxisRange(-1,   10., "y");
//	//HRatioVsPt_IMR->SetAxisRange(-3,   22., "y");
//	HRatioVsPt_IMR->SetMarkerStyle(20);
//	HRatioVsPt_IMR->SetMarkerColor(4);
//	HRatioVsPt_IMR->SetMarkerSize(1.5);
//	HRatioVsPt_IMR->SetLineColor(4);
//	HRatioVsPt_IMR->SetLineWidth(2);
//	HRatioVsPt_IMR->Draw("pe");
//
//	tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[0]+")");
//	tl.DrawLatex(0.14, 0.76, Form("%.1f<M_{ee}<%.1f GeV/c^{2}", IMRMassWindow[0], IMRMassWindow[1]));
//	c1->SaveAs( "./outAnaPlots/Physics/HRatioVsPt_IMR_Cent080.png");
//	c1->SaveAs( "./outAnaPlots/Physics/HRatioVsPt_IMR_Cent080.pdf");
//	
//	//c1->SetLogy();
//	//Draw ExY as a function of Pt
//	HExYVsPt_IMR->SetTitle("Yield Excess: (Data-CKT) in 0-80% vs. p_{T}");
//	HExYVsPt_IMR->SetYTitle("Data-Cocktail (dN/dP_{T})");
//	HExYVsPt_IMR->SetXTitle("p_{T} (GeV/c)");
//	HExYVsPt_IMR->GetYaxis()->SetTitleSize(0.07);
//	HExYVsPt_IMR->GetYaxis()->SetTitleOffset(0.65);
//	HExYVsPt_IMR->GetXaxis()->SetTitleSize(0.05);
//	HExYVsPt_IMR->GetXaxis()->SetTitleOffset(0.85);
//	//HExYVsPt_IMR->SetAxisRange(1.e-7,   1.e-3, "y");
//	HExYVsPt_IMR->SetMarkerStyle(20);
//	HExYVsPt_IMR->SetMarkerColor(4);
//	HExYVsPt_IMR->SetMarkerSize(1.5);
//	HExYVsPt_IMR->SetLineColor(4);
//	HExYVsPt_IMR->SetLineWidth(2);
//	HExYVsPt_IMR->Draw("pe");
//
//	tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV ("+centTitle[0]+")");
//	tl.DrawLatex(0.14, 0.76, Form("%.1f<M_{ee}<%.1f GeV/c^{2}", IMRMassWindow[0], IMRMassWindow[1]));
//	c1->SaveAs( "./outAnaPlots/Physics/HExYVsPt_IMR_Cent080.png");
//	c1->SaveAs( "./outAnaPlots/Physics/HExYVsPt_IMR_Cent080.pdf");
//
//	delete c1;
//}
//
////----------------------------------------------------------------------
//void checkBump()
//{
//	TCanvas *c1 = new TCanvas("c1");
//	gStyle->SetOptStat(0);
//	gPad->SetTopMargin(0.08);
//	gPad->SetBottomMargin(0.11);
//	gPad->SetLeftMargin(0.11);
//	gPad->SetRightMargin(0.05);
//
//	TLatex tl;
//	tl.SetTextSize(0.06);
//	tl.SetNDC();
//
//	TH2D* Htem2d = new TH2D("Htem2d", "", 10, 0.08, 0.7, 10, 1.e-4, 1.);
//	Htem2d->SetTitle(""); 
//	Htem2d->SetYTitle("dN/dM_{ee} (c^{2}/GeV)");
//	Htem2d->SetXTitle("M_{ee} (GeV/c^{2})");
//	Htem2d->GetYaxis()->SetTitleSize(0.06);
//	Htem2d->GetYaxis()->SetTitleOffset(0.75);
//	Htem2d->GetXaxis()->SetTitleSize(0.05);
//	Htem2d->GetXaxis()->SetTitleOffset(0.80);
//	Htem2d->SetMarkerStyle(24); 
//	Htem2d->SetMarkerStyle(20);
//	Htem2d->SetMarkerColor(4);
//	Htem2d->SetLineColor(4);
//	Htem2d->SetLineWidth(2);
//	//Htem2d->GetYaxis()->SetNdivisions(4);
//
//	//	TH2D* Htem2d2 = new TH2D("Htem2d2", "", 10, -0.01, 0.4, 10, 5.e-4, 5.);
//	TH2D* Htem2d2 = new TH2D("Htem2d2", "", 10, 0.08, 0.7, 10, 2.e-4, 1.);
//	Htem2d2->SetTitle(""); 
//	Htem2d2->SetYTitle("dN/dM_{ee} (c^{2}/GeV)");
//	Htem2d2->SetXTitle("M_{ee} (GeV/c^{2})");
//	Htem2d2->GetYaxis()->SetTitleSize(0.06);
//	Htem2d2->GetYaxis()->SetTitleOffset(0.75);
//	Htem2d2->GetXaxis()->SetTitleSize(0.05);
//	Htem2d2->GetXaxis()->SetTitleOffset(0.80);
//	Htem2d2->SetMarkerStyle(24); 
//	Htem2d2->SetMarkerStyle(20);
//	Htem2d2->SetMarkerColor(4);
//	Htem2d2->SetLineColor(4);
//	Htem2d2->SetLineWidth(2);
//
//	c1->SetLogy();
//
//	for(int icent=0; icent<nCentBins; icent++)
//	{
//		//Htem2d->Draw("");
//		H1d_M_Reb_Sig_Raw[icent]->SetTitle(""); 
//		H1d_M_Reb_Sig_Raw[icent]->SetAxisRange(0.08, 0.70, "x");
//		H1d_M_Reb_Sig_Raw[icent]->SetYTitle("dN/dM_{ee} (c^{2}/GeV)");
//		H1d_M_Reb_Sig_Raw[icent]->SetXTitle("M_{ee} (GeV/c^{2})");
//		H1d_M_Reb_Sig_Raw[icent]->GetYaxis()->SetTitleSize(0.06);
//		H1d_M_Reb_Sig_Raw[icent]->GetYaxis()->SetTitleOffset(0.75);
//		H1d_M_Reb_Sig_Raw[icent]->GetXaxis()->SetTitleSize(0.05);
//		H1d_M_Reb_Sig_Raw[icent]->GetXaxis()->SetTitleOffset(0.80);
//		H1d_M_Reb_Sig_Raw[icent]->SetMarkerStyle(20);
//		H1d_M_Reb_Sig_Raw[icent]->SetMarkerSize(0.5);
//		H1d_M_Reb_Sig_Raw[icent]->SetMarkerColor(2);
//		H1d_M_Reb_Sig_Raw[icent]->SetLineColor(2);
//		H1d_M_Reb_Sig_Raw[icent]->SetLineWidth(2);
//		H1d_M_Reb_Sig_Raw[icent]-> Draw("pe");
//		H1d_M_Reb_Sig_aftPSAC[icent]  -> SetMarkerStyle(20);
//		H1d_M_Reb_Sig_aftPSAC[icent]  -> SetMarkerSize(0.5);
//		H1d_M_Reb_Sig_aftPSAC[icent]  -> SetMarkerColor(4);
//		H1d_M_Reb_Sig_aftPSAC[icent]  -> SetLineColor(4);
//		H1d_M_Reb_Sig_aftPSAC[icent]  -> Draw("pesame");
//		
//		tl.SetTextSize(0.050);
//		tl.DrawLatex(0.14, 0.86, "LowMass in Cent: "+centTitle[icent]);
//
//		TLegend* leg_bump = new TLegend(0.65, 0.65, 0.89, 0.85);
//		leg_bump->SetBorderSize(0);
//		leg_bump->SetFillColor(0);
//		leg_bump->AddEntry( H1d_M_Reb_Sig_Raw[icent],     "Raw",       "lp" );
//		leg_bump->AddEntry( H1d_M_Reb_Sig_aftPSAC[icent], "aftPSAC",   "lp" );
//		leg_bump->Draw("same");
//
//		c1->SaveAs(Form("./outAnaPlots/plots_Bump/H1D_M_CheckBump_RawAftPSAC_icent%d.png",icent));
//		delete leg_bump;
//
//		//Htem2d2->Draw("");
//		
//
//		H1d_M_Reb_Sig_aftEff[icent]->SetTitle(""); 
//		H1d_M_Reb_Sig_aftEff[icent]->SetAxisRange(0.08, 0.70, "x");
//		H1d_M_Reb_Sig_aftEff[icent]->SetYTitle("dN/dM_{ee} (c^{2}/GeV)");
//		H1d_M_Reb_Sig_aftEff[icent]->SetXTitle("M_{ee} (GeV/c^{2})");
//		H1d_M_Reb_Sig_aftEff[icent]->GetYaxis()->SetTitleSize(0.06);
//		H1d_M_Reb_Sig_aftEff[icent]->GetYaxis()->SetTitleOffset(0.75);
//		H1d_M_Reb_Sig_aftEff[icent]->GetXaxis()->SetTitleSize(0.05);
//		H1d_M_Reb_Sig_aftEff[icent]->GetXaxis()->SetTitleOffset(0.80);
//		H1d_M_Reb_Sig_aftEff[icent]->SetLineWidth(2);
//		H1d_M_Reb_Sig_aftEff[icent]-> SetMarkerStyle(20);
//		H1d_M_Reb_Sig_aftEff[icent]-> SetMarkerSize(0.5);
//		H1d_M_Reb_Sig_aftEff[icent]-> SetMarkerColor(2);
//		H1d_M_Reb_Sig_aftEff[icent]-> SetLineColor(2);
//		H1d_M_Reb_Sig_aftEff[icent]-> Draw("pe");
//		H1D_M_totCKT_Reb[icent]    -> Draw("histsame");
//
//		tl.SetTextSize(0.050);
//		tl.DrawLatex(0.14, 0.86, "LowMass in Cent: "+centTitle[icent]);
//
//		TLegend* leg_bump2 = new TLegend(0.65, 0.65, 0.89, 0.85);
//		leg_bump2->SetBorderSize(0);
//		leg_bump2->SetFillColor(0);
//		leg_bump2->AddEntry( H1d_M_Reb_Sig_aftEff[icent],     "aft PairEff",  "lp" );
//		leg_bump2->AddEntry( H1D_M_totCKT_Reb[icent],         "totCKT",       "lp" );
//		leg_bump2->Draw("same");
//
//		c1->SaveAs(Form("./outAnaPlots/plots_Bump/H1D_M_CheckBump_AftEff_icent%d.png",icent));
//		delete leg_bump2;
//	}
//
//
//	//compare the corrected Data vs CKT
//	for(int ipt=0; ipt<nPtBins4Phys; ipt++)
//	{
//		cout<<"Draw cocktail in "<<Form("%.2f<pt<%.2f", PtBDs4Phys[ipt], PtBDs4Phys[ipt+1])<<endl;
//		//Htem2d2 = new TH2D("Htem2d2", "", 10, -0.1, 3.6, 10, 1.e-8, 20.);
//		H_M_totCKT4PtPhys[ipt]->SetTitle(Form("Dielectron in %.2f<p_{T}<%.2f GeV/c", PtBDs4Phys[ipt], PtBDs4Phys[ipt+1]));
//		H_M_totCKT4PtPhys[ipt]->SetYTitle("dN/dM_{ee} (c^{2}/GeV)");
//		H_M_totCKT4PtPhys[ipt]->SetXTitle("M_{ee} (GeV/c^{2})");
//		H_M_totCKT4PtPhys[ipt]->GetYaxis()->SetTitleSize(0.06);
//		H_M_totCKT4PtPhys[ipt]->GetYaxis()->SetTitleOffset(0.75);
//		H_M_totCKT4PtPhys[ipt]->GetXaxis()->SetTitleSize(0.05);
//		H_M_totCKT4PtPhys[ipt]->GetXaxis()->SetTitleOffset(0.80);
//		H_M_totCKT4PtPhys[ipt]->SetLineWidth(2);
//		H_M_totCKT4PtPhys[ipt]->SetLineColor(1);
//		H_M_totCKT4PtPhys[ipt]->SetAxisRange(0.08, 0.70, "x");
//		if(ipt==0) H_M_totCKT4PtPhys[ipt]->SetAxisRange(1.e-10, 1.e-3, "y");
//		if(ipt==5) H_M_totCKT4PtPhys[ipt]->SetAxisRange(1.e-7,  1.e-3, "y");
//		H_M_totCKT4PtPhys[ipt]->Draw("chist");
//
//		//H1d_M_Reb4PtPhys_Sig_aftEff[ipt]->SetMarkerStyle(20);
//		//H1d_M_Reb4PtPhys_Sig_aftEff[ipt]->SetMarkerSize(0.5);
//		//H1d_M_Reb4PtPhys_Sig_aftEff[ipt]->SetMarkerColor(4);
//		//H1d_M_Reb4PtPhys_Sig_aftEff[ipt]->SetLineColor(4);
//		//H1d_M_Reb4PtPhys_Sig_aftEff[ipt]->Draw("pesame");
//
//		TLegend *leg = new TLegend(0.60,0.72,0.92,0.76);
//		leg->SetNColumns(3);
//		leg->SetFillStyle(0);
//		leg->SetBorderSize(0);
//		leg->SetTextSize(0.050);
//		leg->AddEntry(H_M_totCKT4PtPhys[ipt],        "Cocktail Sum",     "l");
//		leg->Draw("same");
//
//		tl.SetTextSize(0.050);
//		tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV (0-80%)");
//		//tl.SetTextSize(0.040);
//		//tl.DrawLatex(0.15, 0.79, "p_{T}^{e}>0.2 GeV/c, |#eta^{e}|<1, |y^{ee}|<1");
//
//		c1->SaveAs(Form("./outAnaPlots/plots_Bump/H1D_CheckBump_DataVsCKT_Cent080_ipt%d.png",ipt));
//		c1->SaveAs(Form("./outAnaPlots/plots_Bump/H1D_CheckBump_DataVsCKT_Cent080_ipt%d.pdf",ipt));
//
//		delete leg;
//	}//ipt4phys
//
//	delete c1;
//}
//
////----------------------------------------------------------------------
//----------------------------------------------------------------------
void drawExYieldsVsNpart()
{
	TCanvas *c1 = new TCanvas("c1");
	gStyle->SetOptStat(0);
	gPad->SetTopMargin(0.08);
	gPad->SetBottomMargin(0.11);
	gPad->SetLeftMargin(0.11);
	gPad->SetRightMargin(0.05);

	TLatex tl;
	tl.SetTextSize(0.06);
	tl.SetNDC();
	
	//c1->SetLogy();
	TH2D* Htem2d = new TH2D("Htem2d", "", 10, 0, 400, 10, 1.e-6, 0.60e-4);
	Htem2d->SetTitle(""); 
	Htem2d->SetYTitle("Yield/N_{part}");
	Htem2d->SetXTitle("N_{part}");
	Htem2d->GetYaxis()->SetTitleSize(0.07);
	Htem2d->GetYaxis()->SetTitleOffset(0.60);
	Htem2d->GetXaxis()->SetTitleSize(0.05);
	Htem2d->GetXaxis()->SetTitleOffset(0.80);
	Htem2d->SetMarkerStyle(20);
	Htem2d->SetMarkerColor(4);
	Htem2d->SetLineColor(4);
	Htem2d->SetLineWidth(2);
	//Htem2d->GetYaxis()->SetNdivisions(4);
	Htem2d->Draw("");

	ge_ExYVsNpart_Rho->SetMarkerStyle(20);
	ge_ExYVsNpart_Rho->SetMarkerSize(1.5);
	ge_ExYVsNpart_Rho->SetMarkerColor(1);
	ge_ExYVsNpart_Rho->SetLineColor(1);
	ge_ExYVsNpart_Rho->SetLineWidth(2);
	ge_ExYVsNpart_Rho->Draw("pesame");
	
	ge_YVsNpart_Omega->SetMarkerStyle(24);
	ge_YVsNpart_Omega->SetMarkerColor(4);
	ge_YVsNpart_Omega->SetMarkerSize(1.5);
	ge_YVsNpart_Omega->SetLineColor(4);
	ge_YVsNpart_Omega->SetLineWidth(2);
	ge_YVsNpart_Omega->Draw("pesame");
	
	ge_YVsNpart_Phi->SetMarkerStyle(25);
	ge_YVsNpart_Phi->SetMarkerSize(1.5);
	ge_YVsNpart_Phi->SetMarkerColor(4);
	ge_YVsNpart_Phi->SetLineColor(4);
	ge_YVsNpart_Phi->SetLineWidth(2);
	ge_YVsNpart_Phi->Draw("pesame");

	TF1 *f_ExYvsnpart_Rho = new TF1("f_ExYvsnpart_Rho", "[0]*x[0]^[1]", 0., 400.);
	f_ExYvsnpart_Rho->SetParNames("b", "a");
	f_ExYvsnpart_Rho->SetParameters(1.e-5, 0.24);

	ge_ExYVsNpart_Rho->Fit(f_ExYvsnpart_Rho, "ER", "", 0, 400);
	double aValue_Rho = f_ExYvsnpart_Rho->GetParameter(1);
	double aError_Rho = f_ExYvsnpart_Rho->GetParError(1);

	f_ExYvsnpart_Rho->SetLineColor(4);
	f_ExYvsnpart_Rho->Draw("lsame");

	//TF1* f_const = new TF1("f_const", "[0]", 0, 400);
	//ge_ExYVsNpart_Rho->Fit(f_const,  "ER", "", 0, 400);
	//double constValue_Rho = f_const->GetParameter(0);
	//double constError_Rho = f_const->GetParError(0);
	//f_const->SetLineStyle(5);
	//f_const->Draw("lsame");


	TF1 *f_ExYvsnpart_Omega = new TF1("f_ExYvsnpart_Omega", "[0]*x[0]^[1]", 0., 400.);
	f_ExYvsnpart_Omega->SetParNames("b", "a");
	f_ExYvsnpart_Omega->SetParameters(1.e-6, 0.44);

	ge_YVsNpart_Omega->Fit(f_ExYvsnpart_Omega, "ER", "", 0, 400);
	double aValue_Omega = f_ExYvsnpart_Omega->GetParameter(1);
	double aError_Omega = f_ExYvsnpart_Omega->GetParError(1);

	f_ExYvsnpart_Omega->SetLineColor(6);
	f_ExYvsnpart_Omega->Draw("lsame");


	TF1 *f_ExYvsnpart_Phi = new TF1("f_ExYvsnpart_Phi", "[0]*x[0]^[1]", 0., 400.);
	f_ExYvsnpart_Phi->SetParNames("b", "a");
	f_ExYvsnpart_Phi->SetParameters(1.e-6, 0.44);

	ge_YVsNpart_Phi->Fit(f_ExYvsnpart_Phi, "ER", "", 0, 400);
	double aValue_Phi = f_ExYvsnpart_Phi->GetParameter(1);
	double aError_Phi = f_ExYvsnpart_Phi->GetParError(1);

	f_ExYvsnpart_Phi->SetLineColor(2);
	f_ExYvsnpart_Phi->Draw("lsame");

	tl.SetTextSize(0.04);
	tl.DrawLatex(0.158, 0.75, "fit by: b #times N_{part}^{a}");
	tl.DrawLatex(0.158, 0.70, Form("#rho-like: fitted a = %.2f #pm %.2f",   aValue_Rho,      aError_Rho));
	tl.DrawLatex(0.158, 0.65, Form("#omega-like: fitted a = %.2f #pm %.2f", aValue_Omega,    aError_Omega));
	tl.DrawLatex(0.158, 0.60, Form("#phi-like: fitted a = %.2f #pm %.2f",   aValue_Phi,      aError_Phi));
	//tl.DrawLatex(0.28,  0.30, Form("#rho-like: const Fit = %.2f #pm %.2f (10^{-3})",  constValue_Rho*1.e3,  constError_Rho*1.e3));

	TLegend* leg_yd = new TLegend(0.55, 0.72, 0.92, 0.90);
	leg_yd->SetBorderSize(0);
	leg_yd->SetFillColor(0);
	leg_yd->AddEntry( ge_ExYVsNpart_Rho,     Form("Data-CKT: %.2f<M_{ee}<%.2f GeV/c^{2}", RhoMassWindow[0],   RhoMassWindow[1]),   "lp" );
	leg_yd->AddEntry( ge_YVsNpart_Omega,     Form("Data:     %.2f<M_{ee}<%.2f GeV/c^{2}", OmegaMassWindow[0], OmegaMassWindow[1]), "lp" );
	leg_yd->AddEntry( ge_YVsNpart_Phi,       Form("Data:     %.2f<M_{ee}<%.2f GeV/c^{2}", PhiMassWindow[0],   PhiMassWindow[1]),   "lp" );
	leg_yd->Draw("same");

	tl.SetTextSize(0.06);
	tl.DrawLatex(0.14, 0.85, "Au+Au #sqrt{s_{NN}} = 27 GeV");

	c1->SaveAs( "./outAnaPlots/Physics/H_YieldVsNpart_RhoOmegaPhi.png");
	c1->SaveAs( "./outAnaPlots/Physics/H_YieldVsNpart_RhoOmegaPhi.pdf");

	//----------------------------------------------------------------------------------------------------------------------------------------
	TH2D* Htem2d_IMR = new TH2D("Htem2d_IMR", "", 10, 0, 400, 10, 1.e-8, 1.2e-5);
	Htem2d_IMR->SetTitle(""); 
	Htem2d_IMR->SetYTitle("Yield/N_{part}");
	Htem2d_IMR->SetXTitle("N_{part}");
	Htem2d_IMR->GetYaxis()->SetTitleSize(0.07);
	Htem2d_IMR->GetYaxis()->SetTitleOffset(0.60);
	Htem2d_IMR->GetXaxis()->SetTitleSize(0.05);
	Htem2d_IMR->GetXaxis()->SetTitleOffset(0.80);
	Htem2d_IMR->SetMarkerStyle(20);
	Htem2d_IMR->SetMarkerColor(4);
	Htem2d_IMR->SetLineColor(4);
	Htem2d_IMR->SetLineWidth(2);
	//Htem2d_IMR->GetYaxis()->SetNdivisions(4);
	Htem2d_IMR->Draw("");

	ge_ExYVsNpart_IMR->SetMarkerStyle(20);
	ge_ExYVsNpart_IMR->SetMarkerSize(1.5);
	ge_ExYVsNpart_IMR->SetMarkerColor(2);
	ge_ExYVsNpart_IMR->SetLineColor(2);
	ge_ExYVsNpart_IMR->SetLineWidth(2);
	ge_ExYVsNpart_IMR->Draw("pesame");
	
	TF1 *f_ExYvsnpart_IMR = new TF1("f_ExYvsnpart_IMR", "[0]*x[0]^[1]", 0., 400.);
	f_ExYvsnpart_IMR->SetParNames("b", "a");
	f_ExYvsnpart_IMR->SetParameters(1.e-6, 0.44);

	ge_ExYVsNpart_IMR->Fit(f_ExYvsnpart_IMR, "ER", "", 0, 400);
	double aValue_IMR = f_ExYvsnpart_IMR->GetParameter(1);
	double aError_IMR = f_ExYvsnpart_IMR->GetParError(1);

	f_ExYvsnpart_IMR->SetLineColor(4);
	//f_ExYvsnpart_IMR->SetLineStyle(5);
	f_ExYvsnpart_IMR->Draw("lsame");

	TLegend* leg_yd_IMR = new TLegend(0.14, 0.60, 0.65, 0.80);
	leg_yd_IMR->SetBorderSize(0);
	leg_yd_IMR->SetFillColor(0);
	leg_yd_IMR->AddEntry( ge_ExYVsNpart_IMR, Form("Data-CKT: %.1f<M_{ee}<%.1f GeV/c^{2}",           IMRMassWindow[0], IMRMassWindow[1]),  "lp" );
	leg_yd_IMR->AddEntry( f_ExYvsnpart_IMR,  Form("b #times N_{part}^{a}: fitted a = %.2f #pm %.2f", aValue_IMR,      aError_IMR),         "l" );
	leg_yd_IMR->Draw("same");

	tl.DrawLatex(0.14, 0.85, "Au+Au #sqrt{s_{NN}} = 27 GeV");
	
	c1->SaveAs( "./outAnaPlots/Physics/H_ExYieldVsNpart_IMR.png");
	c1->SaveAs( "./outAnaPlots/Physics/H_ExYieldVsNpart_IMR.pdf");
	delete leg_yd;
	delete c1;
}
////----------------------------------------------------------------------
////----------------------------------------------------------------------
//void drawPhysics()
//{
//	//read in Previous data
//	infile_27Run11  = new TFile("inputfiles/previous/Data_CK_27.root","read");
//	hSpecM_27Run11  = (TH1D*)infile_27Run11->Get("hData");
//	hSpecM_SysErr   = (TH1D*)infile_27Run11->Get("hSys");
//	gCKT            = (TH1D*)infile_27Run11->Get("gCK");
//	gRCKSys         = (TH1D*)infile_27Run11->Get("gRCKSys");
//	gCKTsys         = (TH1D*)infile_27Run11->Get("gCKSys");
//	gpion           = (TH1D*)infile_27Run11->Get("gpion");
//	geta            = (TH1D*)infile_27Run11->Get("geta");
//	getap           = (TH1D*)infile_27Run11->Get("getap");
//	gomega          = (TH1D*)infile_27Run11->Get("gomega");
//	gphi            = (TH1D*)infile_27Run11->Get("gphi");
//	gccbar1         = (TH1D*)infile_27Run11->Get("gccbar1");
//	gccbar1         -> SetLineStyle(5);
//	gccbar1         -> SetLineColor(1);
//	gjpsi           = (TH1D*)infile_27Run11->Get("gjpsi");
//
//	hRatioRun11     = (TH1D*)infile_27Run11->Get("hData2CK");
//	hRatioRun11     ->SetMarkerStyle(24);
//	hRatioRun11     ->SetMarkerColor(28);
//	hRatioRun11     ->SetLineColor(28);
//
//	hRatioRun11Sys  = (TH1D*) infile_27Run11->Get("hRSys");
//	hRatioRun11Sys->SetMarkerStyle(24);
//	hRatioRun11Sys->SetMarkerColor(28);
//	hRatioRun11Sys->SetFillColor(0);
//	hRatioRun11Sys->SetFillStyle(0);
//	hRatioRun11Sys->SetLineColor(28);
//	hRatioRun11Sys->SetLineWidth(2);
//
//	//----------------------------------------------------------------------------------
//	//----------------------------------------------------------------------------------
//	infile_27Rapp     = new TFile("TheoryOverCKT/ForZaochen/27Cocktail_Theory.root", "read");//updated with new CKT
//	gTheory_Rapp_CGA = (TGraph*)infile_27Rapp->Get("gTheroy_Cocktail_Ratio");
//	//infile_27Rapp     = new TFile("inputfiles/previous/Theory_Rapp_27.root", "read");
//	//gTheory_Rapp_CGA = (TGraph*)infile_27Rapp->Get("gRCKThe");
//	//----------------------------------------------------------------------------------
//	//----------------------------------------------------------------------------------
//	TCanvas *c1 = new TCanvas("c1");
//	gStyle->SetOptStat(0);
//	gPad->SetTopMargin(0.08);
//	gPad->SetBottomMargin(0.11);
//	gPad->SetLeftMargin(0.11);
//	gPad->SetRightMargin(0.05);
//	//---------------------------------------------------------------------------
//	//---------------------------------------------------------------------------
//
//	TLatex tl;
//	tl.SetTextSize(0.06);
//	tl.SetNDC();
//
//	c1->SetLogy();
//
//	//Draw Run18 vs Cocktails
//	H1d_M_Reb_Sig_aftEff[0]->SetTitle(""); 
//	H1d_M_Reb_Sig_aftEff[0]->SetYTitle("dN/dM_{ee} (c^{2}/GeV)");
//	//H1d_M_Reb_Sig_aftEff[0]->SetYTitle("dN/dM_{ee}/nEvents (c^{2}/GeV)");
//	H1d_M_Reb_Sig_aftEff[0]->SetAxisRange(-0.2,  3.60, "x");
//	H1d_M_Reb_Sig_aftEff[0]->SetAxisRange(1.e-7, 20., "y");
//	H1d_M_Reb_Sig_aftEff[0]->GetYaxis()->SetTitleSize(0.06);
//	H1d_M_Reb_Sig_aftEff[0]->GetYaxis()->SetTitleOffset(0.75);
//	H1d_M_Reb_Sig_aftEff[0]->GetXaxis()->SetTitleSize(0.05);
//	H1d_M_Reb_Sig_aftEff[0]->GetXaxis()->SetTitleOffset(0.75);
//	H1d_M_Reb_Sig_aftEff[0]->SetMarkerStyle(24); 
//	H1d_M_Reb_Sig_aftEff[0]->SetMarkerStyle(20);
//	H1d_M_Reb_Sig_aftEff[0]->SetMarkerColor(4);
//	H1d_M_Reb_Sig_aftEff[0]->SetLineColor(4);
//	H1d_M_Reb_Sig_aftEff[0]->SetLineWidth(2);
//	H1d_M_Reb_Sig_aftEff[0]->GetYaxis()->SetNdivisions(4);
//	H1d_M_Reb_Sig_aftEff[0]->Draw("pe");
//
//	//---------------------------------------------------------------------------
//	//---------------------------------------------------------------------------
//	//Assign the Sys.Err from Run11 to Run18
//	H1d_M_Reb_Sig_SysErr = (TH1D*) H1d_M_Reb_Sig_aftEff[0]->Clone(Form("H1d_M_Reb_Sig_SysErr"));
//	H1d_M_Reb_Sig_SysErr -> Reset();
//	for(int ib=0; ib<H1d_M_Reb_Sig_aftEff[0]->GetNbinsX(); ib++)
//	{
//		double ibcter = H1d_M_Reb_Sig_aftEff[0]->GetBinCenter(  ib+1 );
//		double ibcont = H1d_M_Reb_Sig_aftEff[0]->GetBinContent( ib+1 );
//
//		double ibContRun11  = hSpecM_SysErr->GetBinContent( hSpecM_SysErr->FindBin(ibcter) );
//		double iSysErrRun11 = hSpecM_SysErr->GetBinError(   hSpecM_SysErr->FindBin(ibcter) );
//
//		double iSysErrRun18 = iSysErrRun11 * (ibcont/ibContRun11);
//
//		H1d_M_Reb_Sig_SysErr->SetBinContent( ib+1, ibcont       );
//		H1d_M_Reb_Sig_SysErr->SetBinError(   ib+1, iSysErrRun18 );
//
//		//for Bins higher than Run11 Data
//		if(ibcter>3.26)
//		{
//
//			ibContRun11  = hSpecM_SysErr->GetBinContent( hSpecM_SysErr->FindBin(3.26) );
//			iSysErrRun11 = hSpecM_SysErr->GetBinError(   hSpecM_SysErr->FindBin(3.26) );
//			iSysErrRun18 = iSysErrRun11 * (ibcont/ibContRun11);
//
//			H1d_M_Reb_Sig_SysErr->SetBinContent( ib+1, ibcont       );
//			H1d_M_Reb_Sig_SysErr->SetBinError(   ib+1, iSysErrRun18 );
//		}
//	}
//
//	hSpecM_SysErr->SetMarkerStyle(24);
//	hSpecM_SysErr->SetMarkerColor(28);
//	hSpecM_SysErr->SetFillColor(0);
//	hSpecM_SysErr->SetFillStyle(0);
//	hSpecM_SysErr->SetLineColor(28);
//	hSpecM_SysErr->SetLineWidth(2);
//
//	H1d_M_Reb_Sig_SysErr->SetFillColor(0);
//	H1d_M_Reb_Sig_SysErr->SetFillStyle(0);
//	H1d_M_Reb_Sig_SysErr->SetMarkerStyle(20);
//	H1d_M_Reb_Sig_SysErr->SetMarkerColor(4);
//	H1d_M_Reb_Sig_SysErr->SetLineColor(4);
//	//H1d_M_Reb_Sig_SysErr->SetLineColor(13);
//	H1d_M_Reb_Sig_SysErr->SetLineWidth(2);
//
//	gpion  ->Draw("csame");
//	geta   ->Draw("csame");
//	getap  ->Draw("csame");
//	gomega ->Draw("csame");
//	gphi   ->Draw("csame");
//	gccbar1->Draw("csame");
//	gjpsi  ->Draw("csame");
//	gCKTsys->Draw("e3same");
//	gCKT   ->Draw("csame");
//
//	H1d_M_Reb_Sig_SysErr->Draw("e2same");
//	H1d_M_Reb_Sig_aftEff[0]->Draw("pesame");
//
//	tl.SetTextSize(0.050);
//	tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV (MinBias)");
//	tl.SetTextSize(0.040);
//	tl.DrawLatex(0.15, 0.79, "p_{T}^{e}>0.2 GeV/c, |#eta^{e}|<1, |y^{ee}|<1");
//	//drawLatex(0.60, 0.50, "#font[22]{STAR Preliminary}",      42, 0.06, 2);
//
//	TLegend *leg = new TLegend(0.18,0.60,0.92,0.78);
//	leg->SetNColumns(3);
//	leg->SetFillStyle(0);
//	leg->SetBorderSize(0);
//	leg->SetTextSize(0.035);
//	//leg->SetTextFont(mTextFont);
//	leg->AddEntry(gpion,    "#pi^{0} #rightarrow #gammaee",                       "l");
//	leg->AddEntry(geta,     "#eta #rightarrow #gammaee",                          "l");
//	leg->AddEntry(getap,    "#eta' #rightarrow #gammaee",                         "l");
//	leg->AddEntry(gomega,   "#omega #rightarrowee & #omega #rightarrow #pi^{0}ee","f");
//	leg->AddEntry(gphi,     "#phi #rightarrow ee & #phi #rightarrow #etaee",      "f");
//	leg->AddEntry(gjpsi,    "J/#psi #rightarrow ee",                              "f");
//	leg->AddEntry(gccbar1,  "c#bar{c} #rightarrow ee (PYTHIA)",                   "l");
//	leg->AddEntry(gCKTsys,  "Cocktail Sum",                                       "lf");
//	leg->Draw();
//
//	c1->SaveAs("./figures/Run18vsCK.png");
//	c1->SaveAs("./figures/Run18vsCK.pdf");
//	//compare Run18 vs Run11
//	//H1d_M_Reb_Sig_aftEff[0]->SetTitle(""); 
//	//H1d_M_Reb_Sig_aftEff[0]->SetYTitle("dN/dM_{ee} (c^{2}/GeV)");
//	////H1d_M_Reb_Sig_aftEff[0]->SetYTitle("dN/dM_{ee}/nEvents (c^{2}/GeV)");
//	//H1d_M_Reb_Sig_aftEff[0]->SetAxisRange(-0.2,  3.5, "x");
//	//H1d_M_Reb_Sig_aftEff[0]->GetYaxis()->SetTitleSize(0.06);
//	//H1d_M_Reb_Sig_aftEff[0]->GetYaxis()->SetTitleOffset(0.75);
//	//H1d_M_Reb_Sig_aftEff[0]->GetXaxis()->SetTitleSize(0.05);
//	//H1d_M_Reb_Sig_aftEff[0]->GetXaxis()->SetTitleOffset(0.75);
//	////H1d_M_Reb_Sig_aftEff[0]->SetMarkerStyle(20);
//	////H1d_M_Reb_Sig_aftEff[0]->SetMarkerColor(1);
//	////H1d_M_Reb_Sig_aftEff[0]->SetLineColor(1);
//	//H1d_M_Reb_Sig_aftEff[0]->SetLineWidth(2);
//	H1d_M_Reb_Sig_aftEff[0]->Draw("pe");
//	hSpecM_27Run11->SetMarkerStyle(24);
//	hSpecM_27Run11->SetMarkerSize(1.00);
//	hSpecM_27Run11->SetMarkerColor(28);
//	hSpecM_27Run11->SetLineColor(28);
//	hSpecM_27Run11->Draw("pesame");
//	gCKTsys->Draw("e3same");
//	gCKT->Draw("csame");
//
//	hSpecM_SysErr->Draw("e2same");
//	H1d_M_Reb_Sig_SysErr->Draw("e2same");
//
//	hSpecM_27Run11->Draw("pesame");
//	H1d_M_Reb_Sig_aftEff[0]->Draw("pesame");
//
//	TLegend* leg_Run11Run18 = new TLegend(0.70, 0.65, 0.92, 0.90);
//	leg_Run11Run18->SetBorderSize(0);
//	leg_Run11Run18->SetFillColor(0);
//	//leg_Run11Run18->AddEntry(H1d_M_Reb_Sig_aftEff[0],    "Run18",                 "lp" );
//	//leg_Run11Run18->AddEntry(hSpecM_27Run11,             "Run11",                 "lp" );
//	leg_Run11Run18->AddEntry(H1d_M_Reb_Sig_SysErr,    "2018",                 "lpf" );
//	leg_Run11Run18->AddEntry(hSpecM_SysErr,           "2011",                 "lpf" );
//	leg_Run11Run18->AddEntry(gCKTsys,                 "Cocktail Sum",          "lf" );
//	leg_Run11Run18->SetTextSize(0.045);
//	leg_Run11Run18->Draw();
//
//	tl.SetTextSize(0.050);
//	tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV (MinBias)");
//	tl.SetTextSize(0.040);
//	tl.DrawLatex(0.15, 0.79, "p_{T}^{e}>0.2 GeV/c, |#eta^{e}|<1, |y^{ee}|<1");
//
//	//drawLatex(0.60, 0.50, "#font[22]{STAR Preliminary}",      42, 0.06, 2);
//
//	c1->SaveAs("./figures/compRun11Run18_Data.png");
//	c1->SaveAs("./figures/compRun11Run18_Data.pdf");
//
//
//	return;
//
//	////-------------------------------------------------------------------------------------------------------------------------------------
//	////-------------------------------------------------------------------------------------------------------------------------------------
//	//HRatio_Data2CKT = (TH1D*)H1d_M_Reb_Sig_aftEff4Ratio[0]->Clone("HRatio_Data2CKT");
//	////rebin for Ratio plots
//	//HRatio_Data2CKT = (TH1D*) rebHisto( HRatio_Data2CKT, "HRatio_Data2CKT_Reb", nMBins_Ratio2, BinBounds_Ratio2, "NO");
//	//gCKT4Ratio      = (TH1D*) rebHisto( gCKT,            "gCKT4Ratio",          nMBins_Ratio2, BinBounds_Ratio2, "NO");
//	//HRatio_Data2CKT -> Divide(gCKT4Ratio);
//
//	//HRatio_Data2CKT_SysErr = (TH1D*)HRatio_Data2CKT->Clone("HRatio_Data2CKT_SysErr");
//	//HRatio_Data2CKT_SysErr -> Reset();
//
//	//for(int ib=0; ib<HRatio_Data2CKT_SysErr->GetNbinsX(); ib++)
//	//{
//	//	double ibcter = HRatio_Data2CKT->GetBinCenter(  ib+1 );
//	//	double ibcont = HRatio_Data2CKT->GetBinContent( ib+1 );
//
//	//	double ibContRun11  = hRatioRun11Sys->GetBinContent( hRatioRun11Sys->FindBin(ibcter) );
//	//	double iSysErrRun11 = hRatioRun11Sys->GetBinError(   hRatioRun11Sys->FindBin(ibcter) );
//	//	
//	//	double iSysErrRun18 = iSysErrRun11 * (ibcont/ibContRun11);
//	//	
//	//	HRatio_Data2CKT_SysErr->SetBinContent( ib+1, ibcont       );
//	//	HRatio_Data2CKT_SysErr->SetBinError(   ib+1, iSysErrRun18 );
//
//	//	//for Bins higher than Run11 Data
//	//	if(ibcter>3.26)
//	//	{
//
//	//		ibContRun11  = hRatioRun11Sys->GetBinContent( hRatioRun11Sys->FindBin(3.26) );
//	//		iSysErrRun11 = hRatioRun11Sys->GetBinError(   hRatioRun11Sys->FindBin(3.26) );
//	//		iSysErrRun18 = iSysErrRun11 * (ibcont/ibContRun11);
//	//		
//	//		HRatio_Data2CKT_SysErr->SetBinContent( ib+1, ibcont       );
//	//		HRatio_Data2CKT_SysErr->SetBinError(   ib+1, iSysErrRun18 );
//	//	}
//	//}
//	//
//	//HRatio_Data2CKT_SysErr->SetFillColor(0);
//	//HRatio_Data2CKT_SysErr->SetFillStyle(0);
//	//HRatio_Data2CKT_SysErr->SetMarkerStyle(20);
//	//HRatio_Data2CKT_SysErr->SetMarkerColor(4);
//	//HRatio_Data2CKT_SysErr->SetLineColor(4);
//	////HRatio_Data2CKT_SysErr->SetLineColor(13);
//	//HRatio_Data2CKT_SysErr->SetLineWidth(2);
//
//	//double lowx_Ratio = -0.05;
//	//double higx_Ratio = 3.35;
//	//c1->SetLogy(0);
//	//TH2D* Htem2d_ratio = new TH2D("Htem2d_ratio", "", 10, lowx_Ratio, higx_Ratio, 10, 0, 8.);
//	//Htem2d_ratio->SetTitle("");
//	//Htem2d_ratio->SetYTitle("Data/Cocktail");
//	//Htem2d_ratio->SetXTitle("M_{ee} (GeV/c^{2})");
//	//Htem2d_ratio->SetAxisRange(0.0,   10., "y");
//	//Htem2d_ratio->GetYaxis()->SetTitleSize(0.07);
//	//Htem2d_ratio->GetYaxis()->SetTitleOffset(0.55);
//	//Htem2d_ratio->GetXaxis()->SetTitleSize(0.05);
//	//Htem2d_ratio->GetXaxis()->SetTitleOffset(0.75);
//	//Htem2d_ratio->Draw("");
//	//
//	//gRCKSys->Draw("e3same");
//	//drawLine(lowx_Ratio, 1, higx_Ratio, 1, 2,1,1);
//
//	//HRatio_Data2CKT->SetMarkerStyle(20);
//	//HRatio_Data2CKT->SetMarkerColor(4);
//	//HRatio_Data2CKT->SetLineColor(4);
//	//HRatio_Data2CKT->SetLineWidth(2);
//	//HRatio_Data2CKT->Draw("pesame");
//	//HRatio_Data2CKT_SysErr->Draw("e2same");
//
//	//hRatioRun11->Draw("pesame");
//	//hRatioRun11Sys->Draw("e2same");
//	////HRatio_Data2CKT->Draw("pesame");
//	//
//	//TLegend* leg_Ratio = new TLegend(0.70, 0.70, 0.94, 0.91);
//	//leg_Ratio->SetBorderSize(0);
//	//leg_Ratio->SetFillColor(0);
//	//leg_Ratio->AddEntry( HRatio_Data2CKT,    "2018",                 "lp" );
//	//leg_Ratio->AddEntry( hRatioRun11,        "2011",                 "lp" );
//	////leg_Ratio->AddEntry( gTheory_Rapp_CGA,  "Rapp Model",           "lp" );
//	//leg_Ratio->SetTextSize(0.045);
//	//leg_Ratio->Draw();
//
//	//tl.SetTextSize(0.050);
//	//tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV (MinBias)");
//	//tl.SetTextSize(0.040);
//	//drawLatex(0.15, 0.75, "#font[22]{STAR Preliminary}",      42, 0.06, 2);
//
//	//c1->SaveAs("./figures/compRun11Run18_Data2CKT.png");
//	//c1->SaveAs("./figures/compRun11Run18_Data2CKT.pdf");
//
//
//	////draw Ratio 2018 vs Theory
//	//Htem2d_ratio->Draw("");
//
//	//gRCKSys->Draw("e3same");
//	//drawLine(lowx_Ratio, 1, higx_Ratio, 1, 2,1,1);
//
//	//HRatio_Data2CKT->Draw("pesame");
//	//HRatio_Data2CKT_SysErr->Draw("e2same");
//	//gTheory_Rapp_CGA->SetLineColor(2);
//	//gTheory_Rapp_CGA->SetLineWidth(3);
//	//gTheory_Rapp_CGA->Draw("csame");
//
//	//TLegend* leg2_Ratio = new TLegend(0.70, 0.70, 0.94, 0.91);
//	//leg2_Ratio->SetBorderSize(0);
//	//leg2_Ratio->SetFillColor(0);
//	//leg2_Ratio->AddEntry( HRatio_Data2CKT,    "Data",                 "lp" );
//	//leg2_Ratio->AddEntry( gTheory_Rapp_CGA,  "Rapp Model",           "lp" );
//	//leg2_Ratio->SetTextSize(0.045);
//	//leg2_Ratio->Draw();
//
//	//tl.SetTextSize(0.050);
//	//tl.DrawLatex(0.14, 0.86, "Au+Au #sqrt{s_{NN}} = 27 GeV (MinBias)");
//	//tl.SetTextSize(0.040);
//	//drawLatex(0.15, 0.75, "#font[22]{STAR Preliminary}",      42, 0.06, 2);
//
//	//c1->SaveAs("./figures/compRun18vsTheory_Data2CKT.png");
//	//c1->SaveAs("./figures/compRun18vsTheory_Data2CKT.pdf");
//}


double fun_MtSpec(double* x, double* par)
{
	double xcur = x[0];
	return par[0]*exp(-xcur/par[1]);
}


double fun_AccedExY(double* x, double* par)
{
	double xcur = x[0];
	return par[0]*pow(xcur, 3./2.)*exp(-xcur/par[1]);
}
//-----------------------------------------------------------------------------------------
double getDnchDy(int icent)
{
	const double dndy_0005 = 172.9+177.1+31.1+22.6+31.7+6.0;
	const double dndy_0510 = 144.3+147.5+25.8+18.7+26.5+5.1;
	const double dndy_1020 = 109.4+111.6+19.4+14.5+19.4+4.0;
	const double dndy_2030 = 74.3 +75.9 +12.9+9.8 +12.9+2.9;
	const double dndy_3040 = 48.8 +49.9 +8.3 +6.2 +8.9 +2.0;
	const double dndy_4050 = 30.7 +31.5 +5.2 +3.9 +5.6 +1.4;
	const double dndy_5060 = 18.6 +18.9 +2.9 +2.2 +3.2 +0.8;
	const double dndy_6070 = 10.4 +10.6 +1.5 +1.1 +1.7 +0.49;
	const double dndy_7080 = 5.1  +5.3  +0.68+0.51+0.8 +0.23;

	const double dndy_0010 = (dndy_0005*0.05+dndy_0510*0.05)/0.10;
	const double dndy_1040 = (dndy_1020*0.10+dndy_2030*0.10+dndy_3040*0.10)/0.30;
	const double dndy_4080 = (dndy_4050*0.1 +dndy_5060*0.1 +dndy_6070*0.1+dndy_7080*0.1)/0.40;
	const double dndy_0080 = (dndy_0005*0.05+dndy_0510*0.05
			+dndy_1020*0.10+dndy_2030*0.10+dndy_3040*0.10
			+dndy_4050*0.1 +dndy_5060*0.1 +dndy_6070*0.1+dndy_7080*0.1)/0.80;

	//	cout<<"dndy_0010: "<<dndy_0010<<endl;
	//	cout<<"dndy_1040: "<<dndy_1040<<endl;
	//	cout<<"dndy_4080: "<<dndy_4080<<endl;
	//	cout<<"dndy_0080: "<<dndy_0080<<endl;

	if(     icent==0) return dndy_0080;
	else if(icent==1) return dndy_0010;
	else if(icent==2) return dndy_1040;
	else if(icent==3) return dndy_4080;
	else
	{
		cout<<"wrong icent for getDnchDy!!!!"<<endl;
		return -999.;
	}
}

//---------------------------------------------------------------------------
double fun_Rcut(double* x, double* par)
{
	double xcur = x[0];
	double a = par[0];
	double b = par[1];
	double y = a*pow(xcur,2) + b;//a*x^2+b

	if(y>0.) y = sqrt(y);
	else y = 0.0;
	return y; 
}

//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------
double bg_formula(double* x, double* par)
{
	//if use pol2
	//double xcur = x[0];
	//double outValue = par[0]+ par[1]*xcur + par[2]*xcur*xcur;
	//return outValue; 
	
	double xcur = x[0];
	return par[0]*exp(-par[1]*xcur);
}

//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------
double BreitWigner_formula(double* x, double* par)
{
	Double_t arg1 = 14.0/22.0; // 2 over pi
	Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma = par[1],  Mass = par[2]
	Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
	Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
	
	return par[0]*arg1*arg2/(arg3 + arg4);
}

//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------
double BreitWigner_BG_formula(double* x, double* par)
{
	Double_t arg1 = 14.0/22.0; // 2 over pi
	Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
	Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
	Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
	
	Double_t value_BW = par[0]*arg1*arg2/(arg3 + arg4);

	//Double_t value_BG = par[3] + par[4]*x[0]+par[5]*x[0]*x[0]; //if use pol2 for bg
	Double_t value_BG = par[3]*exp(-par[4]*x[0]);
	
	Double_t BreitWignerAndBG = value_BW + value_BG;
	
	return BreitWignerAndBG;
}
//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------
TH1D* CombGmLSMixUS( TString outname, TH1D*h1, TH1D*h2, double lowMass, double higMass ) //replace h1 with h2 for mass region [lowMass, higMass]
{
	TH1D* hout = (TH1D*)h1->Clone("hout");
	hout->Reset();
	hout->SetTitle(outname);
	hout->SetName( outname);

	for(int ib=0; ib<hout->GetNbinsX(); ib++)
	{
		double ibCenter   = h1->GetBinCenter(ib+1);
		
		double ibContent1 = h1->GetBinContent(ib+1);
		double ibError1   = h1->GetBinError(  ib+1);
		
		double ibContent2 = h2->GetBinContent(ib+1);
		double ibError2   = h2->GetBinError(  ib+1);

		if(ibCenter<lowMass || ibCenter>higMass)
		{
			hout->SetBinContent(ib+1, ibContent1);
			hout->SetBinError(  ib+1, ibError1  );
		}
		else
		{
			hout->SetBinContent(ib+1, ibContent2);
			hout->SetBinError(  ib+1, ibError2  );
		}
	}

	return hout;
}

double HistGraphAsFitf(double *x, double *par)
{
	TH1D* h1   = (TH1D*)   hist1_4fit->Clone(); //ckt_noVect
	TH1D* h2   = (TH1D*)   hist2_4fit->Clone(); //ckt_omega
	TH1D* h3   = (TH1D*)   hist3_4fit->Clone(); //ckt_phi
	TGraph* gh = (TGraph*) gh_4fit   ->Clone();

	int    ibin1   = h1 ->FindBin(x[0]);
	double ibcont1 = h1 ->GetBinContent(ibin1);
	int    ibin2   = h2 ->FindBin(x[0]);
	double ibcont2 = h2 ->GetBinContent(ibin2);
	int    ibin3   = h3 ->FindBin(x[0]);
	double ibcont3 = h3 ->GetBinContent(ibin3);

	double ibcont4 = gh ->Eval(x[0]);

	double yValue = ibcont1 + par[0]*ibcont2 + par[1]*ibcont3 + par[2]*ibcont4;

	return yValue;
}
