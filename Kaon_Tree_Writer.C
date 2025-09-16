#include <cstdlib>
#include <iostream>
#include <sstream>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TChain.h>
#include <TBenchmark.h>
#include <vector>

// Macro name
void Kaon_Tree_Writer()
{

    //////////////////////////////////////////////////////////////////////////////
    //// Setting up input tree and variables    //////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////

    // Read file with information on vectors
    // gROOT->ProcessLine(".L /home/nics/work/York/Scaling/Loader.C+");
    gROOT->ProcessLine(".L /u/home/matthewn/Documents/Macros/Loader.C+");

    // Read input root file and assign it to 'f'
    TFile *f = new TFile("/volatile/clas12/matthewn/Scaling/eK+_sim_RGA_Fall18_In_Strangeness2_Flat0_Slope1-9650-0.root");
    // TFile *f = new TFile("/home/nics/work/York/Scaling/eK+_005196_Tree.root");
    TTree *t1 = (TTree *)f->Get("Tree");

    // Creating components to read from TTree
    // Particle information
    vector<TLorentzVector> *v_kp = new vector<TLorentzVector>; // 4-vectors for the detected particles
    TLorentzVector *missmassS1 = nullptr, *missmassS2 = nullptr, *missmassS3 = nullptr;
    vector <double> *M_KP;
    Double_t MMS1, MMS2, MMS3;

    // Setting the branch addresses to read from
    t1->SetBranchAddress("kp", &v_kp);
    t1->SetBranchAddress("M_kp", &M_KP);
    t1->SetBranchAddress("MissMassS1", &missmassS1);
    t1->SetBranchAddress("MissMassS2", &missmassS2);
    t1->SetBranchAddress("MissMassS3", &missmassS3);

    // K^+ variables
    Double_t M_kp_1;
    Double_t M_kp_2;
    Double_t M_kp_3;
    Double_t P_kp_1;
    Double_t P_kp_2;
    Double_t P_kp_3;
    Int_t kaonno;

    // Path and name for the output file to save
    TFile fileOutput1("~/eK+_sim_RGA_Fall18_In_Strangeness2_Flat0_Slope1-9650-0_Kaon_Tree.root", "recreate");
    // TFile fileOutput1("/home/nics/work/York/Scaling/eK+_005196_Kaon_Tree.root", "recreate");
    TTree Tree("T", "it's a tree!");
    Tree.Branch("kaonno", &kaonno);
    Tree.Branch("M_kp_1", &M_kp_1);
    Tree.Branch("M_kp_2", &M_kp_2);
    Tree.Branch("M_kp_3", &M_kp_3);
    Tree.Branch("P_kp_1", &P_kp_1);
    Tree.Branch("P_kp_2", &P_kp_2);
    Tree.Branch("P_kp_3", &P_kp_3);
    Tree.Branch("MMS1", &MMS1);
    Tree.Branch("MMS2", &MMS2);
    Tree.Branch("MMS3", &MMS3);

    //////////////////////////////////////////////////////////////////////////////
    //// Create histograms here    ///////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////

    auto *h_S1_Miss_Mass = new TH1F("h_S1_Miss_Mass", "Missing Mass", 300, 0, 3);
    auto *h_S2_Miss_Mass = new TH1F("h_S2_Miss_Mass", "Missing Mass", 300, 0, 3);
    auto *h_S3_Miss_Mass = new TH1F("h_S3_Miss_Mass", "Missing Mass", 300, 2, 5);

    // Histograms looking at kaon properties
    // auto *h_delta_beta_S1_kp_1 = new TH2F("h_delta_beta_S1_kp_1", "#Delta#Beta of K^{+} (1);P [GeV];#Delta#Beta", 480, 0, 12, 400, -0.4, 0.4);
    // auto *h_delta_beta_S2_kp_1 = new TH2F("h_delta_beta_S2_kp_1", "#Delta#Beta of K^{+} (1);P [GeV];#Delta#Beta", 480, 0, 12, 400, -0.4, 0.4);
    // auto *h_delta_beta_S2_kp_2 = new TH2F("h_delta_beta_S2_kp_2", "#Delta#Beta of K^{+} (2);P [GeV];#Delta#Beta", 480, 0, 12, 400, -0.4, 0.4);
    // auto *h_delta_beta_S3_kp_1 = new TH2F("h_delta_beta_S3_kp_1", "#Delta#Beta of K^{+} (1);P [GeV];#Delta#Beta", 480, 0, 12, 400, -0.4, 0.4);
    // auto *h_delta_beta_S3_kp_2 = new TH2F("h_delta_beta_S3_kp_2", "#Delta#Beta of K^{+} (2);P [GeV];#Delta#Beta", 480, 0, 12, 400, -0.4, 0.4);
    // auto *h_delta_beta_S3_kp_3 = new TH2F("h_delta_beta_S3_kp_3", "#Delta#Beta of K^{+} (3);P [GeV];#Delta#Beta", 480, 0, 12, 400, -0.4, 0.4);

    //////////////////////////////////////////////////////////////////////////////
    //// Looping over events in the tree    //////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////

    // Reads the total number of entries in the TTree
    Long64_t nentries = t1->GetEntries();

    // This is used to print out the percentage of events completed so far
    Int_t Percentage = nentries / 100;

    // This loops over all the entries in the TTree
    for (Long64_t i = 0; i < nentries; i++)
    {

        t1->GetEntry(i);

        // // This prints out the percentage of events completed so far
        // if (i % Percentage == 0)
        // {
            // cout << i / Percentage << endl;
        //     fprintf(stderr, "%lld\r", i / Percentage);
        //     fflush(stderr);
        // }

        ////////////////////////////////////////////////////////////////////////////
        // Looping over particles in the current event    //////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        // This reads the number of particles in the current entry/event
        kaonno = v_kp->size();

        M_kp_1 = M_KP->at(0);
        P_kp_1 = v_kp->at(0).Rho();

        if (kaonno == 1)
        {
            MMS1 = missmassS1->M();
            h_S1_Miss_Mass->Fill(MMS1);
        }

        else if (kaonno == 2)
        {
            MMS2 = missmassS2->M();
            h_S2_Miss_Mass->Fill(MMS2);
            M_kp_2 = M_KP->at(1);
            P_kp_2 = v_kp->at(1).Rho();
        }

        else if (kaonno == 3)
        {
            MMS3 = missmassS3->M();
            h_S3_Miss_Mass->Fill(MMS3);
            M_kp_2 = M_KP->at(1);
            P_kp_2 = v_kp->at(1).Rho();
            M_kp_3 = M_KP->at(2);
            P_kp_3 = v_kp->at(2).Rho();
        }
        Tree.Fill();
    }
    Tree.Write();
    fileOutput1.Write();
    f->Close();
}
