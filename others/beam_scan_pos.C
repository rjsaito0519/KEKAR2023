double SAC_x = 145;
double SAC_y = 113;

double BAC_x = 115;
double BAC_y = 115;

double KVC_x = 26;
double KVC_y = 120;

double Trig_x = 10;
double Trig_y = 13;

double x_dis = 16;
double y_dis = 18;

int x_scan = 9;
int y_scan = 7;

void position_drawing(){

   // いろんな設定
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980392156862745, 0.054901960784313725);
    gROOT->GetColor(kBlue)->SetRGB(0.12156862745098039, 0.4666666666666667, 0.7058823529411765);
    gROOT->GetColor(kGreen)->SetRGB((Double_t) 44/255, (Double_t) 160/255, (Double_t) 44/255);
    gROOT->GetColor(kRed)->SetRGB((Double_t) 214/255, (Double_t) 39/255, (Double_t) 40/255);
    
    
    gStyle->SetOptStat(0);
    // gStyle->SetLabelSize(0.05, "XY");
    // gStyle->SetTitleSize(1., "XY");
    // gStyle->SetTitleFontSize(0.08);
    gROOT->GetColor(0)->SetAlpha(0.01);
    gStyle->SetPadBottomMargin(0.15);

  TMultiGraph *mg_geo = new TMultiGraph();

  TGraph *SAC_geo = new TGraph();

  SAC_geo->SetPoint(SAC_geo->GetN(),SAC_x/2,SAC_y/2);
  SAC_geo->SetPoint(SAC_geo->GetN(),SAC_x/2,-SAC_y/2);
  SAC_geo->SetPoint(SAC_geo->GetN(),-SAC_x/2,-SAC_y/2);
  SAC_geo->SetPoint(SAC_geo->GetN(),-SAC_x/2,SAC_y/2);
  SAC_geo->SetPoint(SAC_geo->GetN(),SAC_x/2,SAC_y/2);

  TGraph *BAC_geo = new TGraph();

  BAC_geo->SetPoint(BAC_geo->GetN(),BAC_x/2,BAC_y/2);
  BAC_geo->SetPoint(BAC_geo->GetN(),BAC_x/2,-BAC_y/2);
  BAC_geo->SetPoint(BAC_geo->GetN(),-BAC_x/2,-BAC_y/2);
  BAC_geo->SetPoint(BAC_geo->GetN(),-BAC_x/2,BAC_y/2);
  BAC_geo->SetPoint(BAC_geo->GetN(),BAC_x/2,BAC_y/2);

  TGraph *KVC_geo[4];
  for(int i=0;i<4;i++){
    KVC_geo[i] = new TGraph();


    KVC_geo[i]->SetPoint(KVC_geo[i]->GetN(),KVC_x/2-KVC_x*(i-1),KVC_y/2);
    KVC_geo[i]->SetPoint(KVC_geo[i]->GetN(),KVC_x/2-KVC_x*(i-1),-KVC_y/2);
    KVC_geo[i]->SetPoint(KVC_geo[i]->GetN(),-KVC_x/2-KVC_x*(i-1),-KVC_y/2);
    KVC_geo[i]->SetPoint(KVC_geo[i]->GetN(),-KVC_x/2-KVC_x*(i-1),KVC_y/2);
    KVC_geo[i]->SetPoint(KVC_geo[i]->GetN(),KVC_x/2-KVC_x*(i-1),KVC_y/2);

    KVC_geo[i]->SetLineColor(kGreen);
    KVC_geo[i]->SetLineWidth(4);

    mg_geo->Add(KVC_geo[i]);
    
  }

  TGraph *Beam_geo[x_scan][y_scan];


  for(int i=0;i<x_scan;i++){
    
    for(int j=0;j<y_scan;j++){
      Beam_geo[i][j] = new TGraph();


      Beam_geo[i][j]->SetPoint(Beam_geo[i][j]->GetN(),Trig_x/2+x_dis*(i-x_scan/2),Trig_y/2+y_dis*(j-y_scan/2));
      Beam_geo[i][j]->SetPoint(Beam_geo[i][j]->GetN(),Trig_x/2+x_dis*(i-x_scan/2),-Trig_y/2+y_dis*(j-y_scan/2));
      Beam_geo[i][j]->SetPoint(Beam_geo[i][j]->GetN(),-Trig_x/2+x_dis*(i-x_scan/2),-Trig_y/2+y_dis*(j-y_scan/2));
      Beam_geo[i][j]->SetPoint(Beam_geo[i][j]->GetN(),-Trig_x/2+x_dis*(i-x_scan/2),Trig_y/2+y_dis*(j-y_scan/2));
      Beam_geo[i][j]->SetPoint(Beam_geo[i][j]->GetN(),Trig_x/2+x_dis*(i-x_scan/2),Trig_y/2+y_dis*(j-y_scan/2));
      
      Beam_geo[i][j]->SetLineColor(kBlack);
      Beam_geo[i][j]->SetLineWidth(3);
      Beam_geo[i][j]->SetLineStyle(8);

      mg_geo->Add(Beam_geo[i][j]);
      
    }
  }
  SAC_geo->SetLineColor(kRed);
  SAC_geo->SetLineWidth(4);

  
  BAC_geo->SetLineColor(kBlue);
  BAC_geo->SetLineWidth(4);

  mg_geo->Add(SAC_geo);
  mg_geo->Add(BAC_geo);

  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
    c1->cd(1);
//   mg_geo->SetTitle("Detector Size & Beam position;X [mm];Y [mm]");
  mg_geo->Draw("a");
  c1->Print( Form("./beam_scan_pos.jpg") );
    
}
