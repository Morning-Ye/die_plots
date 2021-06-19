#include "/Users/zaochenrice/myFunction.h"
#include "head.h"
#include "massBins.h"
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void readFiles();
void saveSignAccFactor();
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
void getPairSignAccCorr()
{
	readFiles();
	saveSignAccFactor();
}
//------------------------------------------------------------------------
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void readFiles()
{
	for(int icent9=0; icent9<nCent9; icent9++)
	{
		infile[icent9] = new TFile(Form("./inputfiles/data/outAnaHist_Run18_27GeV_Vz35cm_cent%d.root",icent9));
		
		cout<<"readin: "<<infile[icent9]->GetName()<<endl;

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

		//---------------------------------------------------------------------------
		//get the mixed event histograms for the pair sign acc. corr. calculation
		H2d_MvsPt_MixUS0[icent9]    = (TH2D*)infile[icent9]->Get("hMvsPt_wPhiV_MixUS");
		H2d_MvsPt_MixLSPos0[icent9] = (TH2D*)infile[icent9]->Get("hMvsPt_wPhiV_MixLSPos");
		H2d_MvsPt_MixLSNeg0[icent9] = (TH2D*)infile[icent9]->Get("hMvsPt_wPhiV_MixLSNeg");
		//---------------------------------------------------------------------------
		H2d_MvsMt_MixUS0[icent9]    = (TH2D*)infile[icent9]->Get("hMvsMt_wPhiV_MixUS");
		H2d_MvsMt_MixLSPos0[icent9] = (TH2D*)infile[icent9]->Get("hMvsMt_wPhiV_MixLSPos");
		H2d_MvsMt_MixLSNeg0[icent9] = (TH2D*)infile[icent9]->Get("hMvsMt_wPhiV_MixLSNeg");
	}//icent9

	//const int centBinIndex[nCentBins][2] = { {0,8},   {7,8},   {5,6},    {3,4},    {0,2} };
	//merge the histogram for the interested centrality bins
	for( int icent=0; icent<nCentBins; icent++ )
	{
		//for( int icentbin = centBinIndex[icent][0]; icentbin<centBinIndex[icent][1]; icentbin++ ) //used for QM
		for( int icentbin = centBinIndex[icent][0]; icentbin<=centBinIndex[icent][1]; icentbin++ )
		{
			if( icentbin==centBinIndex[icent][0] )
			{
				H2d_MvsPt_MixUS[icent]    = (TH2D*) H2d_MvsPt_MixUS0[icentbin]   ->Clone();
				H2d_MvsPt_MixLSPos[icent] = (TH2D*) H2d_MvsPt_MixLSPos0[icentbin]->Clone();
				H2d_MvsPt_MixLSNeg[icent] = (TH2D*) H2d_MvsPt_MixLSNeg0[icentbin]->Clone();
				//---------------------------------------------------------------------------
				H2d_MvsMt_MixUS[icent]    = (TH2D*) H2d_MvsMt_MixUS0[icentbin]   ->Clone();
				H2d_MvsMt_MixLSPos[icent] = (TH2D*) H2d_MvsMt_MixLSPos0[icentbin]->Clone();
				H2d_MvsMt_MixLSNeg[icent] = (TH2D*) H2d_MvsMt_MixLSNeg0[icentbin]->Clone();
			}
			else
			{
				H2d_MvsPt_MixUS[icent]    -> Add(H2d_MvsPt_MixUS0[icentbin]);
				H2d_MvsPt_MixLSPos[icent] -> Add(H2d_MvsPt_MixLSPos0[icentbin]);
				H2d_MvsPt_MixLSNeg[icent] -> Add(H2d_MvsPt_MixLSNeg0[icentbin]);
				//---------------------------------------------------------------------------
				H2d_MvsMt_MixUS[icent]    -> Add(H2d_MvsMt_MixUS0[icentbin]);
				H2d_MvsMt_MixLSPos[icent] -> Add(H2d_MvsMt_MixLSPos0[icentbin]);
				H2d_MvsMt_MixLSNeg[icent] -> Add(H2d_MvsMt_MixLSNeg0[icentbin]);
			}
		}
		
		H2d_MvsPt_MixUS[icent]        ->SetName(Form("H2d_MvsPt_MixUS_icent%d",   icent));
		H2d_MvsPt_MixLSPos[icent]     ->SetName(Form("H2d_MvsPt_MixLSPos_icent%d",icent));
		H2d_MvsPt_MixLSNeg[icent]     ->SetName(Form("H2d_MvsPt_MixLSNeg_icent%d",icent));
		//---------------------------------------------------------------------------
		H2d_MvsMt_MixUS[icent]        ->SetName(Form("H2d_MvsMt_MixUS_icent%d",   icent));
		H2d_MvsMt_MixLSPos[icent]     ->SetName(Form("H2d_MvsMt_MixLSPos_icent%d",icent));
		H2d_MvsMt_MixLSNeg[icent]     ->SetName(Form("H2d_MvsMt_MixLSNeg_icent%d",icent));
	}//icent
}
//---------------------------------------------------------------------------
void saveSignAccFactor()
{
	//save for the new 2D corrections
	TFile *foutfnew = new  TFile("./out4PSAC_2DCorr/outhist_4PSAC_2DCorr_27GeV.root", "recreate");

	for(int icent=0; icent<nCentBins; icent++)
	{
		H2d_MvsPt_MixUS[icent]    -> SetTitle(centTitle[icent]);
		H2d_MvsPt_MixUS[icent]    -> Write();
		H2d_MvsPt_MixLSPos[icent] -> SetTitle(centTitle[icent]);
		H2d_MvsPt_MixLSPos[icent] -> Write();
		H2d_MvsPt_MixLSNeg[icent] -> SetTitle(centTitle[icent]);
		H2d_MvsPt_MixLSNeg[icent] -> Write();

		H2d_MvsMt_MixUS[icent]    -> SetTitle(centTitle[icent]);
		H2d_MvsMt_MixUS[icent]    -> Write();
		H2d_MvsMt_MixLSPos[icent] -> SetTitle(centTitle[icent]);
		H2d_MvsMt_MixLSPos[icent] -> Write();
		H2d_MvsMt_MixLSNeg[icent] -> SetTitle(centTitle[icent]);
		H2d_MvsMt_MixLSNeg[icent] -> Write();
	}

	cout<<"outplots rootfile: "<<foutfnew->GetName()<<endl;

	foutfnew->Close();
}
