#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "HipoChain.h"
#include "clas12reader.h"
#include <vector>
#include "rcdb_reader.h"

using namespace clas12;
using namespace std;


void Tree_Maker_PureMC(TString loader_file, TString inFileName, TString outFileName, const std::string databaseF){
   auto start = std::chrono::high_resolution_clock::now();
   gBenchmark->Start("timer");
   int counter=0;

   /////////////////////////////////////////////////////////////////////////
   // Loading the various files required

   TString loadF = loader_file;
   gROOT->ProcessLine(loadF); // Uses Loader.C file, make sure Loader.C is in this file path

   TString outputFile = outFileName;
   TFile fileOutput1(outputFile,"recreate");

   TTree Tree("Tree","it's a tree!");

   // Access the PDG database to get information usind PID (e.g db->GetParticle(211)->Mass() for pi^{+} mass)
   auto db=TDatabasePDG::Instance();

   // Connecting to the RCDB database file for the input file
   const std::string DBfile = databaseF;
   clas12databases::SetRCDBRootConnection(DBfile);

   // Creating a chain for the data from different files
   TChain fake("hipo");

   // Data files to process
   TString inputFile = inFileName;

   /////////////////////////////////////////////////////////////////////////
   // Creating histograms

   auto* hkaonpno=new TH1F("hkaonpno","Number of K^{+};No. of K^{+};Counts",100,0,10);
   auto* hPID=new TH1F("hPID","PID ;PID;Counts",2000,-1000,1000);
   auto* hmass_n=new TH1F("hmass_n","Mass of negatives;Mass [GeV];Counts",1000,-10000,10000);
   auto* hmass_p=new TH1F("hmass_p","Mass of negatives;Mass [GeV];Counts",1000,-10000,10000);

   // Information to save to the tree
   // Any information specific to an event
   Int_t eventno; // Records the event number
   Int_t runno; // Records the run number
   Int_t triggerno; // Records the trigger number
   Double_t start_time; // Records the start time for the event
   Int_t Tree_Events=0; // Records the number of events that are put in the tree

   // Any information specific to an individual particle
   Double_t Photon_Energy; // Photon energy
   TLorentzVector Photon; // TLorentzVector for four vector of electron in FD (Px,Py,Pz,E)
   TLorentzVector el; // TLorentzVector for four vector of electron in FD (Px,Py,Pz,E)
   TLorentzVector p4; // TLorentzVector for four vector of each particle (Px,Py,Pz,E)
   TLorentzVector vertex;   // Vertex position and time (Vx,Vy,Vz,Vt)
   Double_t beta; // beta from time of flight (TOF)
   Double_t status; // Gives information on which detectors and how many picked it up
   Double_t energy; // Energy measured by detector
   Double_t charge; // Charge measured by drift chambers
   Double_t PID; // Records PID deterimined by TOF
   Double_t chi2PID; // Chi^2 for PID of TOF
   Double_t time; // Time recorded by by FTCAL,FTOF or CTOF
   Double_t path; // Path measured by FTCAL,FTOF or CTOF
   Double_t vertex_time; // Calculated vertex time from TOF information
   Double_t Region; // Records 0.0 for FT, 1.0 for FD and 2.0 for CD
   Double_t Mass; // Records the calculated mass from beta and momentum

   // Pure MC information
   TLorentzVector MC_el, MC_Kp1, MC_Kp2, MC_Kp3, MC_Omega, MC_Pim; // TLorentzVector for four vector of simulation particles (Px,Py,Pz,E)
   Double_t SimPx, SimPy, SimPz, SimE, SimMass; // Pure MC bank information

   // Vectors of particle measurables for when you have more than one of a type of particle (e.g 2 pi^{-})
   vector<TLorentzVector> v_p4;
   vector<TLorentzVector> v_vertex;
   vector<double> v_beta;
   vector<double> v_status;
   vector<double> v_energy;
   vector<double> v_charge;
   vector<double> v_PID;
   vector<double> v_chi2PID;
   vector<double> v_time;
   vector<double> v_path;
   vector<double> v_region;

   // Record the number of particles in an event
   Int_t elno; // electrons
   Int_t positronno; // positrons
   Int_t protonno; // protons
   Int_t antiprotonno; // antiprotons
   Int_t neutronno; // neutrons
   Int_t photonno; // photons
   Int_t pipno; // pi^{+}
   Int_t pimno; // pi^{-}
   Int_t pi0no; // pi^{0}
   Int_t kaonpno; // K^{+}
   Int_t Kaonp_FD_Chi3; // K^{+} in FD
   Int_t electronFD; // e' in FD
   Int_t kaonmno; // K^{-}
   Int_t positive_charge_tracks; // Count the number of positive charge tracks
   Int_t negative_charge_tracks; // Count the number of negative charge tracks

   // Setting TLorentzVectors for beam and target
   TLorentzVector beam(0,0,0,0); // Set 4 vector four the beam, all momentum in z-direction
   Double_t beam_E; // Used for the beam energy obtained from the RCDB
   TLorentzVector target(0,0,0,0.93827); // Set 4 vector for target, stationary so no momentum
   // TLorentzVector target(0,0,0,1.8756); // Set 4 vector for target, stationary so no momentum

   Double_t c = 30; // Speed of light, used to calculate vertex time

   // Assign a branch to each measurable and name it
   Tree.Branch("eventno",&eventno);
   Tree.Branch("runno",&runno,"runno/I");
   Tree.Branch("triggerno",&triggerno,"triggerno/I");
   Tree.Branch("start_time",&start_time);
   Tree.Branch("p4",&v_p4);
   Tree.Branch("el",&el);
   Tree.Branch("Photon_Energy",&Photon_Energy,"Photon_Energy/D");
   Tree.Branch("vertex",&v_vertex);
   Tree.Branch("beta",&v_beta);
   Tree.Branch("status",&v_status);
   Tree.Branch("energy",&v_energy);
   Tree.Branch("charge",&v_charge);
   Tree.Branch("PID",&v_PID);
   Tree.Branch("chi2PID",&v_chi2PID);
   Tree.Branch("region",&v_region);
   Tree.Branch("time",&v_time);
   Tree.Branch("path",&v_path);
   Tree.Branch("beam",&beam);
   Tree.Branch("target",&target);
   Tree.Branch("elno",&elno);
   Tree.Branch("negative_charge_tracks",&negative_charge_tracks, "neg");
   Tree.Branch("positive_charge_tracks",&positive_charge_tracks);
   Tree.Branch("Kaonp_FD_Chi3",&Kaonp_FD_Chi3);
   Tree.Branch("electronFD",&electronFD);


   //Creating a chain for the data from different filesc12->queryRcdb();
   clas12root::HipoChain chain;
   chain.Add(inputFile);
   chain.SetReaderTags({0});  //create clas12reader with just tag 0 events
   auto config_c12=chain.GetC12Reader();

   chain.db()->turnOffQADB();

   auto& rcdbData= config_c12->rcdb()->current();

   auto& c12=chain.C12ref();


   // This loop goes over the events within each file
   while (chain.Next()){
      auto c12=chain.GetC12Reader();

      counter++;
      // Clear the vectors from the previous event
      v_p4.clear();
      v_vertex.clear();
      v_beta.clear();
      v_status.clear();
      v_energy.clear();
      v_charge.clear();
      v_PID.clear();
      v_chi2PID.clear();
      v_time.clear();
      v_path.clear();
      v_region.clear();

      elno = 0;
      positronno = 0;
      protonno = 0;
      antiprotonno = 0;
      neutronno = 0;
      photonno = 0;
      pipno = 0;
      pimno = 0;
      pi0no = 0;
      kaonpno = 0;
      Kaonp_FD_Chi3 = 0;
      electronFD = 0;
      kaonmno = 0;
      positive_charge_tracks = 0;
      negative_charge_tracks = 0;

      // Define how to access the information for each event
      auto mceve = c12->mcevent(); // Getting the event number

      // Setting beam energy depending on which run it is
      // beam_E = rcdbData.beam_energy*0.001;
      // beam_E = mceve->getEbeam();
      beam_E = 10.6;
      beam.SetXYZM(0, 0, beam_E, 0.000511);

      // cout<<" beam energy  "<<mceve->getEbeam()<<" type "<<mceve->getBtype()<<endl;

      auto mcpbank=c12->mcparts();
      const Int_t  Ngen=mcpbank->getRows();

      for(Int_t i=0;i<Ngen;i++){
         mcpbank->setEntry(i);

         auto px=mcpbank->getPx();
         auto py=mcpbank->getPy();
         auto pz=mcpbank->getPz();
         auto pm=mcpbank->getMass();
         p4.SetXYZM(px,py,pz,pm);

         PID = mcpbank->getPid();

         // cout<<" particle "<<i<<" "<<PID<<" p4 = "<<p4.X()<<","<<p4.Y()<<","<<p4.Z()<<","<<p4.T()<<" and mass "<<p4.M()<<endl;


         // Looking at pid for particles
         if(PID == 11) elno++; // Count number of electrons
         if(PID==321) kaonpno++; // Count the number of positive kaons
         if(PID==2212)protonno++; // Count the number of protons
         if(PID==-211) pimno++; // Count the number of negative pions
         if(PID==-321) kaonmno++; // Count the number of negative pions

         // Save the particle information in the vectors by pushing it back
         v_PID.push_back(PID);
         v_p4.push_back(p4); // Recording the 4 vector for each neutron

         hPID->Fill(PID);

      }
      // Checking number of kaons per event
      hkaonpno->Fill(kaonpno);


      // Here you can apply a basic skim for events you want to save in your tree
      if(kaonpno > 0 && elno == 1){
         Photon_Energy = 0;
         Photon = beam - el;
         Photon_Energy = Photon.E();

         Tree.Fill();
         Tree_Events++;
      }
   }


   Tree.Write(); // Write information to the TTree
   auto finish = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = finish - start;
   std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";
   cout<<"events in tree "<<Tree_Events<<endl;
   fileOutput1.Write(); // Write information to the root file
   fileOutput1.Close(); // Close the root file at the end
}
