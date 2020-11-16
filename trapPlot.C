#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TGraph.h"
#include "TCanvas.h"

#ifdef __MAKECINT__
#pragma link C++ class vector<double>+;
#endif

void trapPlot() {
    
    double e=1.6e-19;
    double r_0 = 0.03;
    double m= 2.1801714e-25;
    
    
    TFile *inputFile = new TFile("trapOutput_vacuum_e7.root");
    
    TTree *trapTree;
    trapTree = (TTree*)inputFile->Get("trap");
    
    std::vector<double> *vector_normR = 0;
    std::vector<double> *vector_normZ = 0;
    std::vector<double> *vector_Vac = 0;
    std::vector<double> *vector_Vdc = 0;
    std::vector<double> *vector_w = 0;
    
    trapTree->SetBranchAddress("norm_R", &vector_normR);
    trapTree->SetBranchAddress("norm_Z", &vector_normZ);
    trapTree->SetBranchAddress("Vdc", &vector_Vdc);
    trapTree->SetBranchAddress("Vac", &vector_Vac);
    trapTree->SetBranchAddress("w", &vector_w);
    
    
    trapTree->GetEntry(0);
    
    TGraph *gr = new TGraph();
    TGraph *gz = new TGraph();
    TGraph *gboth = new TGraph();
    TH1F *stable_Vac = new TH1F("stable_Vac","Vac",20 ,0 ,500);
    TH1F *stable_Vdc = new TH1F("stable_Vdc","Vdc",100 ,-150 ,150);
    TH1F *stable_w = new TH1F("stable_w","w",999 ,3000 ,3000000);
    
    TMultiGraph *multi = new TMultiGraph();
    
    for (int i = 0; i < vector_Vac->size(); i++) {
        
        double r_norm = vector_normR->at(i);
        double z_norm = vector_normZ->at(i);
        double thisVac = vector_Vac->at(i);
        double thisVdc = vector_Vdc->at(i);
        double thisW = vector_w->at(i);
        double a = (4*e*thisVdc)/(m*r_0 *r_0 * thisW*thisW);
        double q = (2*e*thisVac)/(m*r_0 *r_0 * thisW*thisW);
        
        if(r_norm < 3){
            
            gr->SetPoint(i, q, a);
        }
        
        if(z_norm < 3){
            
            gz->SetPoint(i, q, a);
            
        }
        
        if((r_norm < 3) && (z_norm < 3)){
            
            gboth->SetPoint(i, q, -a);
            stable_Vac->Fill(thisVac);
            stable_Vdc->Fill(thisVdc);
            stable_w->Fill(thisW);
        }
    }
    
    TCanvas *c1 = new TCanvas("c1", "c1", 0, 300);
    /*
    stable_w->Draw();
    
    TCanvas *c2 = new TCanvas("c2", "c2", 0, 500);
    stable_Vac->Draw();
    
    TCanvas *c3 = new TCanvas("c3", "c3", -150, 150);
    stable_Vdc->Draw();
    */
     //multi->Add(gz);
     //multi->Add(gr);
     multi->Add(gboth);
     multi->GetXaxis()->SetTitle("q");
     multi->GetYaxis()->SetTitle("a");
     multi->Draw("AP");
     c1->RedrawAxis();
     
     /*
     gz->SetMarkerColor(5);
     gz->SetMarkerStyle(3);
     
     gr->SetMarkerColor(6);
     gr->SetMarkerStyle(3);
     */
     gboth->SetMarkerColor(7);
     gboth->SetMarkerStyle(3);
     
     
    
    
}

