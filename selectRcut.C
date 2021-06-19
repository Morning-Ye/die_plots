//this codes is to find the best Rcut to reject the edge effects
#include "/Users/zaochenrice/myFunction.h"
#include "head.h"

double fun_Rcut(double* x, double* par);
TF1* fRcut;

void readFiles();
void draw_bfRcut();
void applyRcut();
void draw_aftRcut();

void selectRcut()
{
	fRcut = new TF1("fRcut", fun_Rcut, 0., 0.415, 2);
	fRcut -> SetParameters(-1.02, 0.17);
	
	readFiles();
	draw_bfRcut();
	applyRcut();
	draw_aftRcut();
}

void readFiles()
{
	//be careful, this outtem.root need to be produced before apply Rcut in the plotDieSig.C
	TFile *infile = new TFile("outtem.root", "read");
	cout<<"readin: "<<infile->GetName()<<endl;

	for(int icent=0; icent<nCentBins; icent++)
	{
		H2d_MvsPt_PairEff[icent]     = (TH2D*) infile->Get(Form("H2d_MvsPt_PairEff_icent%d",icent));
		H2d_MvsPt_PairAcc[icent]     = (TH2D*) infile->Get(Form("H2d_MvsPt_PairAcc_icent%d",icent));
		H2d_MvsPt_PairEffAcc[icent]  = (TH2D*) infile->Get(Form("H2d_MvsPt_PairEffAcc_icent%d",icent));
		H2d_MvsPt_Reb_US[icent]      = (TH2D*) infile->Get(Form("H2d_MvsPt_Reb_US_icent%d",icent));
		H2d_MvsPt_Reb_LSPos[icent]   = (TH2D*) infile->Get(Form("H2d_MvsPt_Reb_LSPos_icent%d",icent));
		H2d_MvsPt_Reb_LSNeg[icent]   = (TH2D*) infile->Get(Form("H2d_MvsPt_Reb_LSNeg_icent%d",icent));
	}//icent
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void applyRcut()
{
	//apply the Rcut reject the edge effects
	for(int icent=0; icent<nCentBins; icent++)
	{
		H2d_MvsPt_PairEff[icent]   = (TH2D*)applyRcut2D(H2d_MvsPt_PairEff[icent],   Form("H2d_MvsPt_PairEffRcut_icent%d",  icent), fRcut);
		H2d_MvsPt_PairAcc[icent]   = (TH2D*)applyRcut2D(H2d_MvsPt_PairAcc[icent],   Form("H2d_MvsPt_PairAccRcut_icent%d",  icent), fRcut);
		H2d_MvsPt_Reb_US[icent]    = (TH2D*)applyRcut2D(H2d_MvsPt_Reb_US[icent],    Form("H2d_MvsPt_Reb_USRcut_icent%d",   icent), fRcut);
		H2d_MvsPt_Reb_LSPos[icent] = (TH2D*)applyRcut2D(H2d_MvsPt_Reb_LSPos[icent], Form("H2d_MvsPt_Reb_LSPosRcut_icent%d",icent), fRcut);
		H2d_MvsPt_Reb_LSNeg[icent] = (TH2D*)applyRcut2D(H2d_MvsPt_Reb_LSNeg[icent], Form("H2d_MvsPt_Reb_LSNegRcut_icent%d",icent), fRcut);
	}//icent
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void draw_bfRcut()
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
	
	for(int icent=0; icent<nCentBins; icent++)
	{
		//---------------------------------------------------------------------------
		H2d_MvsPt_PairEff[icent]    ->GetYaxis()->SetTitleSize(0.06);
		H2d_MvsPt_PairEff[icent]    ->GetYaxis()->SetTitleOffset(0.80);
		H2d_MvsPt_PairEff[icent]    ->GetXaxis()->SetTitleSize(0.05);
		H2d_MvsPt_PairEff[icent]    ->GetXaxis()->SetTitleOffset(0.85);
		H2d_MvsPt_PairEff[icent]    ->SetAxisRange(0., 0.5, "x");
		H2d_MvsPt_PairEff[icent]    ->SetAxisRange(0., 0.5, "y");
		H2d_MvsPt_PairEff[icent]    ->Draw("TEXT10");

		fRcut -> Draw("lsame");

		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_PairEff_icent%d.png", icent) );
		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_PairEff_icent%d.pdf", icent) );

		//---------------------------------------------------------------------------
		H2d_MvsPt_PairAcc[icent]    ->GetYaxis()->SetTitleSize(0.06);
		H2d_MvsPt_PairAcc[icent]    ->GetYaxis()->SetTitleOffset(0.80);
		H2d_MvsPt_PairAcc[icent]    ->GetXaxis()->SetTitleSize(0.05);
		H2d_MvsPt_PairAcc[icent]    ->GetXaxis()->SetTitleOffset(0.85);
		H2d_MvsPt_PairAcc[icent]    ->SetAxisRange(0., 0.5, "x");
		H2d_MvsPt_PairAcc[icent]    ->SetAxisRange(0., 0.5, "y");
		H2d_MvsPt_PairAcc[icent]    ->Draw("TEXT10");

		fRcut -> Draw("lsame");

		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_PairAcc_icent%d.png", icent) );
		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_PairAcc_icent%d.pdf", icent) );
		
		//---------------------------------------------------------------------------
		H2d_MvsPt_Reb_US[icent]    ->GetYaxis()->SetTitleSize(0.06);
		H2d_MvsPt_Reb_US[icent]    ->GetYaxis()->SetTitleOffset(0.80);
		H2d_MvsPt_Reb_US[icent]    ->GetXaxis()->SetTitleSize(0.05);
		H2d_MvsPt_Reb_US[icent]    ->GetXaxis()->SetTitleOffset(0.85);
		H2d_MvsPt_Reb_US[icent]    ->SetAxisRange(0., 0.5, "x");
		H2d_MvsPt_Reb_US[icent]    ->SetAxisRange(0., 0.5, "y");
		H2d_MvsPt_Reb_US[icent]    ->Draw("TEXT10");

		fRcut -> Draw("lsame");

		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_Reb_US_icent%d.png", icent) );
		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_Reb_US_icent%d.pdf", icent) );

		//---------------------------------------------------------------------------
		H2d_MvsPt_Reb_LSPos[icent]    ->GetYaxis()->SetTitleSize(0.06);
		H2d_MvsPt_Reb_LSPos[icent]    ->GetYaxis()->SetTitleOffset(0.80);
		H2d_MvsPt_Reb_LSPos[icent]    ->GetXaxis()->SetTitleSize(0.05);
		H2d_MvsPt_Reb_LSPos[icent]    ->GetXaxis()->SetTitleOffset(0.85);
		H2d_MvsPt_Reb_LSPos[icent]    ->SetAxisRange(0., 0.5, "x");
		H2d_MvsPt_Reb_LSPos[icent]    ->SetAxisRange(0., 0.5, "y");
		H2d_MvsPt_Reb_LSPos[icent]    ->Draw("TEXT10");

		fRcut -> Draw("lsame");

		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_Reb_LSPos_icent%d.png", icent) );
		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_Reb_LSPos_icent%d.pdf", icent) );

		//---------------------------------------------------------------------------
		H2d_MvsPt_Reb_LSNeg[icent]    ->GetYaxis()->SetTitleSize(0.06);
		H2d_MvsPt_Reb_LSNeg[icent]    ->GetYaxis()->SetTitleOffset(0.80);
		H2d_MvsPt_Reb_LSNeg[icent]    ->GetXaxis()->SetTitleSize(0.05);
		H2d_MvsPt_Reb_LSNeg[icent]    ->GetXaxis()->SetTitleOffset(0.85);
		H2d_MvsPt_Reb_LSNeg[icent]    ->SetAxisRange(0., 0.5, "x");
		H2d_MvsPt_Reb_LSNeg[icent]    ->SetAxisRange(0., 0.5, "y");
		H2d_MvsPt_Reb_LSNeg[icent]    ->Draw("TEXT10");

		fRcut -> Draw("lsame");

		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_Reb_LSNeg_icent%d.png", icent) );
		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_Reb_LSNeg_icent%d.pdf", icent) );

	}//icent

	delete c1;
}
//---------------------------------------------------------------------------
void draw_aftRcut()
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
	
	for(int icent=0; icent<nCentBins; icent++)
	{
		//---------------------------------------------------------------------------
		H2d_MvsPt_PairEff[icent]    ->GetYaxis()->SetTitleSize(0.06);
		H2d_MvsPt_PairEff[icent]    ->GetYaxis()->SetTitleOffset(0.80);
		H2d_MvsPt_PairEff[icent]    ->GetXaxis()->SetTitleSize(0.05);
		H2d_MvsPt_PairEff[icent]    ->GetXaxis()->SetTitleOffset(0.85);
		H2d_MvsPt_PairEff[icent]    ->SetAxisRange(0., 0.5, "x");
		H2d_MvsPt_PairEff[icent]    ->SetAxisRange(0., 0.5, "y");
		H2d_MvsPt_PairEff[icent]    ->Draw("TEXT10");

		fRcut -> Draw("lsame");

		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_PairEff_aftRcut_icent%d.png", icent) );
		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_PairEff_aftRcut_icent%d.pdf", icent) );

		//---------------------------------------------------------------------------
		H2d_MvsPt_PairAcc[icent]    ->GetYaxis()->SetTitleSize(0.06);
		H2d_MvsPt_PairAcc[icent]    ->GetYaxis()->SetTitleOffset(0.80);
		H2d_MvsPt_PairAcc[icent]    ->GetXaxis()->SetTitleSize(0.05);
		H2d_MvsPt_PairAcc[icent]    ->GetXaxis()->SetTitleOffset(0.85);
		H2d_MvsPt_PairAcc[icent]    ->SetAxisRange(0., 0.5, "x");
		H2d_MvsPt_PairAcc[icent]    ->SetAxisRange(0., 0.5, "y");
		H2d_MvsPt_PairAcc[icent]    ->Draw("TEXT10");

		fRcut -> Draw("lsame");

		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_PairAcc_aftRcut_icent%d.png", icent) );
		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_PairAcc_aftRcut_icent%d.pdf", icent) );
		
		//---------------------------------------------------------------------------
		H2d_MvsPt_Reb_US[icent]    ->GetYaxis()->SetTitleSize(0.06);
		H2d_MvsPt_Reb_US[icent]    ->GetYaxis()->SetTitleOffset(0.80);
		H2d_MvsPt_Reb_US[icent]    ->GetXaxis()->SetTitleSize(0.05);
		H2d_MvsPt_Reb_US[icent]    ->GetXaxis()->SetTitleOffset(0.85);
		H2d_MvsPt_Reb_US[icent]    ->SetAxisRange(0., 0.5, "x");
		H2d_MvsPt_Reb_US[icent]    ->SetAxisRange(0., 0.5, "y");
		H2d_MvsPt_Reb_US[icent]    ->Draw("TEXT10");

		fRcut -> Draw("lsame");

		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_Reb_US_aftRcut_icent%d.png", icent) );
		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_Reb_US_aftRcut_icent%d.pdf", icent) );

		//---------------------------------------------------------------------------
		H2d_MvsPt_Reb_LSPos[icent]    ->GetYaxis()->SetTitleSize(0.06);
		H2d_MvsPt_Reb_LSPos[icent]    ->GetYaxis()->SetTitleOffset(0.80);
		H2d_MvsPt_Reb_LSPos[icent]    ->GetXaxis()->SetTitleSize(0.05);
		H2d_MvsPt_Reb_LSPos[icent]    ->GetXaxis()->SetTitleOffset(0.85);
		H2d_MvsPt_Reb_LSPos[icent]    ->SetAxisRange(0., 0.5, "x");
		H2d_MvsPt_Reb_LSPos[icent]    ->SetAxisRange(0., 0.5, "y");
		H2d_MvsPt_Reb_LSPos[icent]    ->Draw("TEXT10");

		fRcut -> Draw("lsame");

		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_Reb_LSPos_aftRcut_icent%d.png", icent) );
		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_Reb_LSPos_aftRcut_icent%d.pdf", icent) );

		//---------------------------------------------------------------------------
		H2d_MvsPt_Reb_LSNeg[icent]    ->GetYaxis()->SetTitleSize(0.06);
		H2d_MvsPt_Reb_LSNeg[icent]    ->GetYaxis()->SetTitleOffset(0.80);
		H2d_MvsPt_Reb_LSNeg[icent]    ->GetXaxis()->SetTitleSize(0.05);
		H2d_MvsPt_Reb_LSNeg[icent]    ->GetXaxis()->SetTitleOffset(0.85);
		H2d_MvsPt_Reb_LSNeg[icent]    ->SetAxisRange(0., 0.5, "x");
		H2d_MvsPt_Reb_LSNeg[icent]    ->SetAxisRange(0., 0.5, "y");
		H2d_MvsPt_Reb_LSNeg[icent]    ->Draw("TEXT10");

		fRcut -> Draw("lsame");

		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_Reb_LSNeg_aftRcut_icent%d.png", icent) );
		c1->SaveAs( Form("./out4RcutPlots/H2d_MvsPt_Reb_LSNeg_aftRcut_icent%d.pdf", icent) );

	}//icent
	
	delete c1;
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

