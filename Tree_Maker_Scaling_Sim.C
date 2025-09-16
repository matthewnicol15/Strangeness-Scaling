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

using namespace clas12;
using namespace std;

void Tree_Maker_Scaling_Sim(TString loader_file, TString inFileName, TString outFileName)
{
    auto start = std::chrono::high_resolution_clock::now();
    gBenchmark->Start("timer");
    int counter = 0;

    /////////////////////////////////////////////////////////////////////////
    // Loading the various files required

    TString loadF = loader_file;
    gROOT->ProcessLine(loadF); // Uses Loader.C file, make sure Loader.C is in this file path

    TString outputFile = outFileName;
    TFile fileOutput1(outputFile, "recreate");

    TTree Tree_Sim("Tree", "it's a tree!");
    TTree Tree_MC("Tree", "it's a tree!");

    // Access the PDG database to get information usind PID (e.g db->GetParticle(211)->Mass() for pi^{+} mass)
    auto db = TDatabasePDG::Instance();

    // Creating a chain for the data from different files
    TChain fake("hipo");

    // Data files to process
    TString inputFile = inFileName;

    /////////////////////////////////////////////////////////////////////////
    // Creating histograms

    // Information to save to the tree
    // Any information specific to an event
    Int_t eventno;   // Records the event number
    Int_t triggerno; // Records the trigger number
    Int_t Tree_Events_MC = 0;
    Int_t Tree_Events_Sim = 0;

    // Any information specific to an individual particle
    Double_t Photon_Energy_MC, Photon_Energy_Sim; // Photon energy
    TLorentzVector el_MC, el_Sim;      // TLorentzVector for four vector of electron in FD (Px,Py,Pz,E)
    TLorentzVector p4_MC, p4_Sim;      // TLorentzVector for four vector of each particle (Px,Py,Pz,E)
    vector <TLorentzVector> v_kp_MC, v_kp_Sim;
    Double_t beta;          // beta from time of flight (TOF)
    Double_t status;        // Gives information on which detectors and how many picked it up
    Double_t energy;        // Energy measured by detector
    Double_t charge;        // Charge measured by drift chambers
    Double_t PID;           // Records PID deterimined by TOF
    Double_t chi2PID;       // Chi^2 for PID of TOF
    Double_t time;          // Time recorded by by FTCAL,FTOF or CTOF
    Double_t path;          // Path measured by FTCAL,FTOF or CTOF
    Double_t vertex_time;   // Calculated vertex time from TOF information
    Double_t Region;        // Records 0.0 for FT, 1.0 for FD and 2.0 for CD
    Double_t Mass;          // Records the calculated mass from beta and momentum

    // Pure MC information
    TLorentzVector MC_el, MC_Kp1, MC_Kp2, MC_Kp3; // TLorentzVector for four vector of simulation particles (Px,Py,Pz,E)
    Double_t SimPx, SimPy, SimPz, SimE, SimMass;  // Pure MC bank information

    // Vectors of particle measurables for when you have more than one of a type of particle (e.g 2 pi^{-})
    vector<TLorentzVector> v_p4;
    vector<TLorentzVector> v_kp;
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
    Int_t elno;       // electrons
    Int_t kaonpno;    // K^{+}
    Int_t Kaonp_FD;   // K^{+} in FD
    Int_t electronFD; // e' in FD

    // Setting TLorentzVectors for beam and target
    TLorentzVector beam(0, 0, 0, 0);         // Set 4 vector four the beam, all momentum in z-direction
    Double_t beam_E;                         // Used for the beam energy
    TLorentzVector target(0, 0, 0, 0.93827); // Set 4 vector for target, stationary so no momentum
    // TLorentzVector target(0,0,0,1.8756); // Set 4 vector for target, stationary so no momentum

    Double_t c = 30; // Speed of light, used to calculate vertex time

    // Assign a branch to each measurable and name it
    // Pure MC tree
    Tree_MC.Branch("p4", &v_p4_MC);
    Tree_MC.Branch("el", &el_MC);
    Tree_MC.Branch("kp", &v_kp_MC);
    Tree_MC.Branch("eventno", &eventno_MC);
    Tree_MC.Branch("beam", &beam_MC);
    Tree_MC.Branch("target", &target_MC);
    Tree_MC.Branch("elno", &elno_MC);
    Tree_MC.Branch("Photon_Energy", &Photon_Energy_MC, "Photon_Energy/D");
    Tree_MC.Branch("M_kp", &v_M_kp_MC);
    Tree_MC.Branch("MissMassS1", &MissMassS1_MC);
    Tree_MC.Branch("MissMassS2", &MissMassS2_MC);
    Tree_MC.Branch("MissMassS3", &MissMassS3_MC);
    Tree_MC.Branch("PID", &v_PID_MC);

    // Processed sim tree
    Tree_Sim.Branch("p4", &v_p4_Sim);
    Tree_Sim.Branch("el", &el_Sim);
    Tree_Sim.Branch("kp", &v_kp_Sim);
    Tree_Sim.Branch("eventno", &eventno_Sim);
    Tree_Sim.Branch("beam", &beam_Sim);
    Tree_Sim.Branch("target", &target_Sim);
    Tree_Sim.Branch("elno", &elno_Sim);
    Tree_Sim.Branch("Photon_Energy", &Photon_Energy_Sim, "Photon_Energy/D");
    Tree_Sim.Branch("M_kp", &v_M_kp_Sim);
    Tree_Sim.Branch("MissMassS1", &MissMassS1_Sim);
    Tree_Sim.Branch("MissMassS2", &MissMassS2_Sim);
    Tree_Sim.Branch("MissMassS3", &MissMassS3_Sim);
    Tree_Sim.Branch("PID", &v_PID_Sim);
    Tree_Sim.Branch("Kaonp_FD", &Kaonp_FD_Sim);
    Tree_Sim.Branch("electronFD", &electronFD_Sim);
    Tree_Sim.Branch("chi2PID", &v_chi2PID_Sim);
    Tree_Sim.Branch("region", &v_region_Sim);

    //    Tree.Branch("runno",&runno,"runno/I");
    //    Tree.Branch("triggerno",&triggerno,"triggerno/I");
    //    Tree.Branch("start_time",&start_time);
    // Tree.Branch("vertex", &v_vertex);
    // Tree.Branch("beta", &v_beta);
    // Tree.Branch("status", &v_status);
    // Tree.Branch("energy", &v_energy);
    // Tree.Branch("charge", &v_charge);
    // Tree.Branch("time", &v_time);
    // Tree.Branch("path", &v_path);

    clas12root::HipoChain chain;
    chain.Add(inputFile);
    chain.SetReaderTags({0}); // create clas12reader with just tag 0 events
    auto config_c12 = chain.GetC12Reader();

    auto &c12 = chain.C12ref();

    // This loop goes over the events within each file
    while (chain.Next())
    {
        auto c12 = chain.GetC12Reader();

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
        // beam_E = mceve->getEbeam();
        beam_E = 10.6;
        beam.SetXYZM(0, 0, beam_E, 0.000511);

        // cout<<" beam energy  "<<mceve->getEbeam()<<" type "<<mceve->getBtype()<<endl;

        auto mcpbank = c12->mcparts();
        const Int_t Ngen = mcpbank->getRows();

        for (Int_t i = 0; i < Ngen; i++)
        {
            mcpbank->setEntry(i);

            auto px = mcpbank->getPx();
            auto py = mcpbank->getPy();
            auto pz = mcpbank->getPz();
            auto pm = mcpbank->getMass();
            p4.SetXYZM(px, py, pz, pm);

            PID = mcpbank->getPid();

            // cout<<" particle "<<i<<" "<<PID<<" p4 = "<<p4.X()<<","<<p4.Y()<<","<<p4.Z()<<","<<p4.T()<<" and mass "<<p4.M()<<endl;

            // Looking at pid for particles
            if (PID == 11)
                elno++; // Count number of electrons
            if (PID == 321)
                kaonpno++; // Count the number of positive kaons
            if (PID == 2212)
                protonno++; // Count the number of protons
            if (PID == -211)
                pimno++; // Count the number of negative pions
            if (PID == -321)
                kaonmno++; // Count the number of negative pions

            // Save the particle information in the vectors by pushing it back
            v_PID.push_back(PID);
            v_p4.push_back(p4); // Recording the 4 vector for each neutron

            hPID->Fill(PID);
        }
        // Checking number of kaons per event
        hkaonpno->Fill(kaonpno);

        // Here you can apply a basic skim for events you want to save in your tree
        if (kaonpno > 0 && elno == 1)
        {
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
    std::cout << "Elapsed time: " << elapsed.count() << " events = " << counter << " s\n";
    cout << "events in tree " << Tree_Events << endl;
    fileOutput1.Write(); // Write information to the root file
    fileOutput1.Close(); // Close the root file at the end
}
// Define how to access the information for each particle
p4.SetXYZM(p->par()->getPx(), p->par()->getPy(), p->par()->getPz(), 0);               // Sets the four vector for each particle, mass will be set once PID determined
time = p->getTime();                                                                  // Gets the TOF
path = p->getPath();                                                                  // Gets the path
beta = p->par()->getBeta();                                                           // Gets the beta from TOF
vertex_time = time - path / (beta * c);                                               // Calculated vertex time
vertex.SetXYZT(p->par()->getVx(), p->par()->getVy(), p->par()->getVz(), vertex_time); // Sets the vertex vector for each particle
status = p->par()->getStatus();                                                       // Information on which detectors are involved for this particle
energy = p->getDetEnergy();                                                           // Energy measured by detector
charge = p->par()->getCharge();                                                       // Charge of the particle measured (FTB isn't recognised!!!)
PID = p->par()->getPid();                                                             // PID determined by TOF
chi2PID = p->par()->getChi2Pid();                                                     // Chi^2 for PID from TOF
Mass = sqrt((pow(p4.Rho(), 2) / (pow(beta, 2))) - pow(p4.Rho(), 2));

// Save the particle information in the vectors by pushing it back
v_vertex.push_back(vertex);
v_beta.push_back(beta);
v_status.push_back(status);
v_energy.push_back(energy);
v_charge.push_back(charge);
v_PID.push_back(PID);
v_chi2PID.push_back(chi2PID);
v_time.push_back(time);
v_path.push_back(path);
// cout<<p4.Rho()<<endl;
v_p4.push_back(p4); // Recording the 4 vector for each neutron

hPID->Fill(PID);

if (PID == 11)
    elno++; // Count number of electrons

// Looking at positive particles
if (charge > 0)
{
    if (PID == 321)
        kaonpno++; // Count the number of positive kaons
    if (PID == 2212)
        protonno++;           // Count the number of protons
    positive_charge_tracks++; // Count number of positive charges

    // Calculated mass for positive particles
    hmass_p->Fill(Mass);
}

// Looking at negative particles
else if (charge < 0)
{
    if (PID == -211)
        pimno++; // Count the number of negative pions
    if (PID == -321)
        kaonmno++; // Count the number of negative pions
    negative_charge_tracks++;

    // Calculated mass for negative particles
    hmass_n->Fill(Mass);
}

// Recording the region the particles hit
if (p->getRegion() == FT)
{
    Region = 0.0;
}
else if (p->getRegion() == FD)
{
    Region = 1.0;
    if (PID == 321)
    {
        kp.SetXYZM(p->par()->getPx(), p->par()->getPy(), p->par()->getPz(), db->GetParticle(321)->Mass());
        v_kp.push_back(kp);
        v_M_kp.push_back(Mass);
        Kaonp_FD++;
    }
    if (PID == 11)
    {
        electronFD++;
        if (electronFD == 1)
            el.SetXYZM(p->par()->getPx(), p->par()->getPy(), p->par()->getPz(), db->GetParticle(11)->Mass());
    }
}
else if (p->getRegion() == CD)
{
    Region = 2.0;
}
v_region.push_back(Region); // pushing back the region of all particles in the event
}

// Checking number of kaons per event
hkaonpno->Fill(kaonpno);

// Here you can apply a basic skim for events you want to save in your tree
if (electronFD > 0 && Kaonp_FD > 0)
{
    Photon_Energy = 0;
    Photon = beam - el;
    Photon_Energy = Photon.E();

    // put the kaons in order of momentum
    if (Kaonp_FD == 1)
    {
        MissMassS1 = beam + target - el - v_kp.at(0);
    }
    else if (Kaonp_FD == 2)
    {
        MissMassS2 = beam + target - el - v_kp.at(0) - v_kp.at(1);
    }
    else if (Kaonp_FD == 3)
    {
        MissMassS3 = beam + target - el - v_kp.at(0) - v_kp.at(1) - v_kp.at(2);
    }
    // else if (Kaonp_FD > 3)
    // {
    //     for (int k = 0;k < v_kp.size(); k++)
    //     {

    //     }
    // }

    Tree.Fill();
    Tree_Events++;
}
}

Tree.Write(); // Write information to the TTree
auto finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed = finish - start;
std::cout << "Elapsed time: " << elapsed.count() << " events = " << counter << " s\n";
cout << "events in tree " << Tree_Events << endl;

Golden = 100. * evCountOK / evCount;
// Accumulated_Charge = qa->GetAccumulatedCharge();
// print fraction of events which pass QA cuts
printf("\nrun = %d\n", runno);
printf("number of events analyzed = %d\n", evCount);
printf("number of events in golden DSTs = %d  (%f%%)\n",
       evCountOK, Golden);
// Filling and saving the data tree
DataTree.Fill();
DataTree.Write();

fileOutput1.Write(); // Write information to the root file
fileOutput1.Close(); // Close the root file at the end
