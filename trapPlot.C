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
    
    
    //Defining some constants
    double e= 1.6e-19;
    double r_0 = 3e-3;
    double m= 2.28e-25;

    double drag = 3.47e-16;
    
    //Creating my root file
    TFile *inputFile = new TFile("trapOutput_precious_withdrag.root");
    TTree *trapTree;
    trapTree = (TTree*)inputFile->Get("trap");
    
    //Setting my variables and creating the branchs for them
    std::vector<double> *vector_normR = 0;
    std::vector<double> *vector_normZ = 0;
    std::vector<double> *vector_Vac = 0;
    std::vector<double> *vector_Vdc = 0;
    std::vector<double> *vector_frequency = 0;
    std::vector<double> *vector_avgV = 0;
    
    trapTree->SetBranchAddress("norm_R", &vector_normR);
    trapTree->SetBranchAddress("norm_Z", &vector_normZ);
    trapTree->SetBranchAddress("Vdc", &vector_Vdc);
    trapTree->SetBranchAddress("Vac", &vector_Vac);
    trapTree->SetBranchAddress("w", &vector_frequency);
    trapTree->SetBranchAddress("vAvg", &vector_avgV);
    trapTree->GetEntry(0);
    
    
    //Creating the 3D histogram for the stability bands
    TH3* h3 = new TH3D("h3", "Stability plot", 200, -10.0, 10.0, 200, -10.0, 10.0,
                       100, 0.0, 1.0);
    
    //##############################################//
    //######### Old block of code , ignore #########//
    //##############################################//
    
    // TGraph *gr = new TGraph();
    //TGraph *gz = new TGraph();
    //TGraph *gboth = new TGraph();
    //TMultiGraph *multi = new TMultiGraph();
    
    //##############################################//
    //##############################################//
    //##############################################//
    
    //Looping over all the elements in my TTree and retrieving the information from them
    for (int i = 0; i < vector_Vac->size(); i++) {
        
        double r_norm = vector_normR->at(i);
        double z_norm = vector_normZ->at(i);
        double thisVac = vector_Vac->at(i);
        double thisVdc = vector_Vdc->at(i);
        double thisW = vector_frequency->at(i);
        double thisvAvg = vector_avgV->at(i);
        
        //Calculating my dimensionless variables
        double a = (-16*e*thisVdc)/(m*r_0 *r_0 * thisW*thisW);
        double q = (8*e*thisVac)/(m*r_0 *r_0 * thisW*thisW);
        double vTerm = m*g/drag;
        
        
        //This variable is 1 if the particle is in a unstable setting and 0 if it is perfectly stable
        double x = thisvAvg/vTerm;
        
        h3->Fill(a,q,x);
        
        
        //##############################################//
        //######### Old block of code , ignore #########//
        //##############################################//
        /*  if(r_norm < (r_0)){
         
         gr->SetPoint(i, q, a);
         
         }
         
         if(z_norm < (r_0)){
         
         gz->SetPoint(i, q, a);
         
         }
         
         if((r_norm < (r_0)) && (z_norm < (r_0))){
         
         gboth->SetPoint(i, q, -a);
         
         }*/
        //##############################################//
        //##############################################//
        //##############################################//
        
    } // -> closing the loop
    
    //Drawing my histogram
    h3->Draw();
    TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
    
    //##############################################//
    //######### Old block of code , ignore #########//
    //##############################################//
    // multi->Add(gz);
    //multi->Add(gr);
    //multi->Add(gboth);
    //multi->GetXaxis()->SetTitle("q");
    //multi->GetYaxis()->SetTitle("a");
    //multi->Draw("AP");
    //c1->RedrawAxis();
    /*
     gz->SetMarkerColor(5);
     gz->SetMarkerStyle(3);
     
     gr->SetMarkerColor(6);
     gr->SetMarkerStyle(3);
     
     gboth->SetMarkerColor(9);
     gboth->SetMarkerStyle(3);
     
     */
    //##############################################//
    //##############################################//
    //##############################################//
    
}

