#include <string>
#include <utility>
#include <vector>

#include "Riostream.h"
#include <string>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>
#include <limits>
#include <getopt.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TRandom.h"
#include <TMath.h>
#include "TSpline.h"

using namespace std;
//histograms{
TH2F * hist_signal = new TH2F("hist_signal", "signal;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_out3DST = new TH2F("hist_bkg_out3DST", "out3DST background;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_NC = new TH2F("hist_bkg_NC", "NC background;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_1 = new TH2F("hist_bkg_1", "secondary background;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_1_pion = new TH2F("hist_bkg_1_pion", "secondary background comming from pion;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_1_neutron = new TH2F("hist_bkg_1_neutron", "secondary background comming from neutron;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_1_proton = new TH2F("hist_bkg_1_proton", "secondary background comming from proton;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_1_other = new TH2F("hist_bkg_1_other", "secondary background comming from other;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_out3DST_largeTime = new TH2F("hist_bkg_out3DST_lt", "out3DST background;Lever Arm [cm]; Time [ns]", 20, 0, 200, 100, 0, 10000);
TH2F * hist_bkg_NC_largeTime = new TH2F("hist_bkg_NC_lt", "NC background;Lever Arm [cm]; Time [ns]", 20, 0, 200, 100, 0, 10000);
TH2F * hist_bkg_1_largeTime = new TH2F("hist_bkg_1_lt", "secondary background;Lever Arm [cm]; Time [ns]", 20, 0, 200, 100, 0, 10000);

TH3F * hist_neutron_hit = new TH3F("asdf","asdf",240,-120,120,240,-120,120,200,-100,100);

TH1F * neutronParentPDG = new TH1F("PDG","PDG",3500,-500,3000);
TH1F * neutronParentPDG_case4 = new TH1F("PDG_case4","PDG_case4",6,0,6);

TH1D * KE_secondary = new TH1D("seconday","seconday",100,0,200);
TH1D * KE_primary = new TH1D("primary","primary",100,0,200);

TH2F * first_n_position_XY = new TH2F("XY", "XY", 240, -120, 120, 240,-120, 120);
TH2F * first_n_position_YZ = new TH2F("YZ", "YZ", 240, -120, 120, 200,-100, 100);
TH2F * first_n_position_XZ = new TH2F("XZ", "XZ", 240, -120, 120, 200,-100, 100);

TH1F * dist_sp_vtx = new TH1F("sp_vtx" ,"from starting point to vertex;[cm]", 150, 0, 150);
TH1F * dist_sp_nh = new TH1F("sp_nh" ,"from starting point to neutron hit;[cm]", 150, 0, 150);
TH1F * dist_sig_sp_vtx = new TH1F("sp_sig_vtx" ,"from starting point to vertex;[cm]", 150, 0, 150);
TH1F * dist_sig_sp_nh = new TH1F("sp_sig_nh" ,"from starting point to neutron hit;[cm]", 150, 0, 150);
TH1F * dist_sig_sp_vtx1 = new TH1F("sp_sig_vtx1" ,"from starting point to vertex, primary;[cm]", 150, 0, 150);
TH1F * dist_sig_sp_nh1 = new TH1F("sp_sig_nh1" ,"from starting point to neutron hit, primary;[cm]", 150, 0, 150);
TH1F * dist_sig_sp_vtx2 = new TH1F("sp_sig_vtx2" ,"from starting point to vertex, primary;[cm]", 150, 0, 150);
TH1F * dist_sig_sp_nh2 = new TH1F("sp_sig_nh2" ,"from starting point to neutron hit, primary;[cm]", 150, 0, 150);
TH1F * dist_sig_sp_vtx3 = new TH1F("sp_sig_vtx3" ,"from starting point to vertex, secondary;[cm]", 150, 0, 150);
TH1F * dist_sig_sp_nh3 = new TH1F("sp_sig_nh3" ,"from starting point to neutron hit, secondary;[cm]", 150, 0, 150);

TH1F * dist_sig_sp_vtx_pion = new TH1F("sp_sig_vtx_pion" ,"from starting point to vertex, from pion;[cm]", 150, 0, 150);
TH1F * dist_sig_sp_nh_pion = new TH1F("sp_sig_nh_pion" ,"from starting point to neutron hit, from pion;[cm]", 150, 0, 150);
TH1F * dist_sig_sp_vtx_proton = new TH1F("sp_sig_vtx_proton" ,"from starting point to vertex, from proton;[cm]", 150, 0, 150);
TH1F * dist_sig_sp_nh_proton = new TH1F("sp_sig_nh_proton" ,"from starting point to neutron hit, from proton;[cm]", 150, 0, 150);
TH1F * dist_sig_sp_vtx_neutron = new TH1F("sp_sig_vtx_neutron" ,"from starting point to vertex, from neutron;[cm]", 150, 0, 150);
TH1F * dist_sig_sp_nh_neutron = new TH1F("sp_sig_nh_neutron" ,"from starting point to neutron hit, from neutron;[cm]", 150, 0, 150);

TH1F * dist_vtx_to_nh_secondary = new TH1F("vtx_nh_secondary","from vertex to neutron hit, secondary;[cm]", 150, 0, 150);
TH1F * dist_vtx_to_nh_pion = new TH1F("vtx_nh_pion","from vertex to neutron hit, from pion;[cm]", 150, 0, 150);
TH1F * dist_vtx_to_nh_neutron = new TH1F("vtx_nh_neutron","from vertex to neutron hit, from neutron;[cm]", 150, 0, 150);
TH1F * dist_vtx_to_nh_proton = new TH1F("vtx_nh_proton","from vertex to neutron hit, from proton;[cm]", 150, 0, 150);

TH1F * angle_event = new TH1F("asdf","angle; pi",10,0,1);
TH1F * angle_event_accumulated = new TH1F("asd","test; pi",10,0,1);

TH1F * angle_vtx_signal = new TH1F("a","angle between vtx and neutron signal; pi",10,0,1);
TH1F * angle_vtx_secondary = new TH1F("b","angle between vtx and secondary neutron; pi",10,0,1);

TH1F * angle_piDeath_neutron_hit = new TH1F("c","angle between pi death point and neutron hit from it; pi",10,0,1);

bool is_inFV = false;       //check if vertex is in FV
bool is_in3DST = false;     //check if vertex is in 3DST

int number_of_CC = 0;

class Hit_t 
{
    public:
    float timeWindow,           // time windows of the hit
          timeSmear,        // smear time
          energyDeposit,        // energy deposited by the neutron
          trackLength,          // lever arm
          trueRec,      // true reconstructed energy
          smearRec,
          vtxSignal[3],     // neutrino vertex position of the neutron
          vtxTime,      // neutrino  vertex time
          neutronTrueE,    //neutron true energy
          neutronTrueT,    //neutron true time
          piDeath[3],      //pion death
          protonDeath[3];      //proton Death

    //neutron hit position
    float neutronHit[3];

    float neutronStartingPoint[3];

    int bkgLoc,         // neutrino vertex position
        neutronParentId,    // Where the neutron come from
        neutronParentPdg;   // PDG of neutron parent

    bool isTherePion50,     // Is there a pion with KE > 50 MeV in FS particles
         isThereProton300;      // Is there a proton with KE > 300 MeV in FS particles
    bool isEmpty;

    Hit_t()
    {
        for(int i = 0; i < 3; i++)
        {
            this->vtxSignal[i] = 0; 
            this->piDeath[i] = 0;   
            this->protonDeath[i] = 0; 
            this->neutronHit[i] = 0;
            this->neutronStartingPoint[i] = 0;
        }
        this->timeWindow = 0;
        this->timeSmear = 0;    
        this->energyDeposit = 0;
        this->trackLength = 0;  
        this->trueRec = 0;      
        this->smearRec = 0;
        this->vtxTime = 0;      
        this->neutronTrueE = 0; 
        this->neutronTrueT = 0; 
        this->bkgLoc = 125124123;        
        this->neutronParentId = 123124123;
        this->neutronParentPdg = 123123123;
        this->isTherePion50 = 0;  
        this->isThereProton300 = 0;
        this->isEmpty = 1;
    }
};

double kineticEnergy(float arm, float time)
{
    double mass = 1.67492729;     //1.674*10^-27 kg
    double velocity = arm/time;      //10^7 m/s
    //cout<<"v:"<<velocity<<endl;
    double KE_J = mass*pow(velocity,2)/2;    //10^-13 kg*m/s (J)
    //cout<<"energy(J):"<<kineticEnergy_J<<endl;
    double KE_MeV = KE_J/1.60217646;     //MeV
    //cout<<"energy(MeV):"<<kineticEnergy_MeV<<endl;

    return KE_MeV;
}

float energyHitCut = 0.5; //energy deposit threshold for cube

int num_file = 0;
int all_interaction = 0;

void num_interaction(string file)
{
    auto _file = new TFile(TString(file));
    auto tree = (TTree*)_file->Get("tree");

    if(tree == NULL)
    {
        _file->Close();
        return;
    }
    else
    {
        num_file++;
    }

    float t_vtx[3], t_vtxTime;
    int nevents = tree->GetEntries();
    tree->SetBranchAddress("vtx", &t_vtx);

    for(int event = 0; event < nevents; event++)
    {
        tree->GetEntry(event);

        if(abs(t_vtx[0]) < 120 && abs(t_vtx[1]) < 120 && abs(t_vtx[2]) < 100)
        {
            all_interaction++;
        }
    }
    _file->Close();
}

//a, b are vector
double GetAngle(float a[], float b[])
{
    float norm_a[3];
    float norm_b[3];
    //normalize
    norm_a[0] = a[0]/(pow(pow(a[0],2)+pow(a[1],2)+pow(a[2],2),0.5));
    norm_a[1] = a[1]/(pow(pow(a[0],2)+pow(a[1],2)+pow(a[2],2),0.5));
    norm_a[2] = a[2]/(pow(pow(a[0],2)+pow(a[1],2)+pow(a[2],2),0.5));

    norm_b[0] = b[0]/(pow(pow(b[0],2)+pow(b[1],2)+pow(b[2],2),0.5));
    norm_b[1] = b[1]/(pow(pow(b[0],2)+pow(b[1],2)+pow(b[2],2),0.5));
    norm_b[2] = b[2]/(pow(pow(b[0],2)+pow(b[1],2)+pow(b[2],2),0.5));

    //get angle
    return TMath::ACos(norm_a[0]*norm_b[0]+ norm_a[1]*norm_b[1]+norm_a[2]*norm_b[2])/TMath::Pi();
}
    
//a, b are vector
double GetDistance(float a[], float b[])
{
    return pow(pow(a[0]-b[0],2)+pow(a[1]-b[1],2)+pow(a[2]-b[2],2),0.5);
}


void test_test_analyze(string file)
{
    auto _file = new TFile(TString(file));
    auto tree = (TTree*)_file->Get("tree");

    if(tree == NULL)
    {
        _file->Close();
        return;
    }

    float t_neutronHitX[1000], t_neutronHitY[1000], t_neutronHitZ[1000];
    float t_neutronStartingPointX[1000], t_neutronStartingPointY[1000], t_neutronStartingPointZ[1000];
    float t_neutronHitT[1000], t_neutronParentId[1000], t_neutronParentPDG[1000];
    float t_neutronHitE[1000], t_neutronTrueE[1000];
    float t_vtx[3], t_vtxTime;
    float t_piDeath[3], t_protonDeath[3];

    float vec_vtx_to_secondary_vertex[3], vec_secondary_vertex_to_neutron_hit[3];
    float vec_vtx_to_sig[3], vec_vtx_to_secondary_neutron[3];
    float norm_vec_vtx_to_sig[3], norm_vec_vtx_to_secondary_neutron[3];
    float norm_vec_vtx_to_secondary_vertex[3], norm_vec_secondary_vertex_to_neutron_hit[3];

    int PDG = 0;
    int t_nFS, t_fsPdg[1000];

    tree->SetBranchAddress("neutronHitX", &t_neutronHitX);
    tree->SetBranchAddress("neutronHitY", &t_neutronHitY);
    tree->SetBranchAddress("neutronHitZ", &t_neutronHitZ);
    tree->SetBranchAddress("neutronStartingPointX", &t_neutronStartingPointX);
    tree->SetBranchAddress("neutronStartingPointY", &t_neutronStartingPointY);
    tree->SetBranchAddress("neutronStartingPointZ", &t_neutronStartingPointZ);
    tree->SetBranchAddress("neutronHitT", &t_neutronHitT);
    tree->SetBranchAddress("neutronParentId", &t_neutronParentId);
    tree->SetBranchAddress("neutronParentPDG", &t_neutronParentPDG);
    tree->SetBranchAddress("neutronHitE", &t_neutronHitE);
    tree->SetBranchAddress("neutronTrueE", &t_neutronTrueE);
    tree->SetBranchAddress("vtx", &t_vtx);
    tree->SetBranchAddress("vtxTime", &t_vtxTime);
    tree->SetBranchAddress("nFS", &t_nFS);
    tree->SetBranchAddress("fsPdg", &t_fsPdg);
    tree->SetBranchAddress("piDeath", &t_piDeath);
    tree->SetBranchAddress("protonDeath", &t_protonDeath);

    int nevents = tree->GetEntries();

    for(int event = 0; event < nevents; event++)
    {
        Hit_t earliest_neutron_hit;
        Hit_t earliest_neutron_hit_for_pion;
        for(int i = 0; i < 3; i++)
        {
            vec_vtx_to_secondary_vertex[i] = 0; 
            vec_secondary_vertex_to_neutron_hit[i] = 0;
            vec_vtx_to_sig[i] = 0; 
            vec_vtx_to_secondary_neutron[i] = 0;
            norm_vec_vtx_to_sig[i] = 0; 
            norm_vec_vtx_to_secondary_neutron[i] = 0;
            norm_vec_vtx_to_secondary_vertex[i] = 0; 
            norm_vec_secondary_vertex_to_neutron_hit[i] = 0;
        }

        tree->GetEntry(event);
        //if(abs(t_vtx[0]) < 50 && abs(t_vtx[1]) < 50 && abs(t_vtx[2]) < 50)
        if(abs(t_vtx[0]) < 50 && abs(t_vtx[1]) < 50 && t_vtx[2] < 100 && t_vtx[2] >0)
            //if(1)
        {
            float temp_earliest_time = 1000000;
            for(int n_neutronHit = 0; n_neutronHit < 1000; n_neutronHit++)
            {
                //if(t_neutronHitX[n_neutronHit] != 0 && t_neutronHitT[n_neutronHit] < temp_earliest_time && t_neutronHitE[n_neutronHit] > energyHitCut)
                if(t_neutronHitX[n_neutronHit] != 0 && t_neutronStartingPointX[n_neutronHit] == t_piDeath[0] && t_neutronHitT[n_neutronHit] < temp_earliest_time && t_neutronHitE[n_neutronHit] > energyHitCut)
                {
                    temp_earliest_time = t_neutronHitT[n_neutronHit];
                    //look for a neutron hit in 3DST
                    /*
                       if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                       abs(t_neutronHitY[n_neutronHit]) < 120 && abs(t_neutronHitZ[n_neutronHit]) < 100) */
                    if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                            abs(t_neutronHitY[n_neutronHit]) < 120 && 
                            t_neutronHitZ[n_neutronHit] < 150 &&
                            t_neutronHitZ[n_neutronHit] > 50)
                    {
                        //calculate lever arm
                        float trackLength = pow(
                                pow(t_neutronHitX[n_neutronHit] - t_vtx[0],2)+
                                pow(t_neutronHitY[n_neutronHit] - t_vtx[1],2)+
                                pow(t_neutronHitZ[n_neutronHit] - t_vtx[2],2),0.5);

                        //calculate signal window; time of flight
                        float signalWindow = t_neutronHitT[n_neutronHit] - t_vtxTime;

                        //Fix a bug from edep-sim
                        if(signalWindow == 1)
                            signalWindow = 0.5;

                        if(signalWindow > 0)
                        {
                            earliest_neutron_hit.timeWindow = signalWindow;
                            earliest_neutron_hit.trackLength = trackLength;
                            earliest_neutron_hit.energyDeposit = t_neutronHitE[n_neutronHit];

                            earliest_neutron_hit.vtxSignal[0] = t_vtx[0];
                            earliest_neutron_hit.vtxSignal[1] = t_vtx[1];
                            earliest_neutron_hit.vtxSignal[2] = t_vtx[2];

                            earliest_neutron_hit.piDeath[0] = t_piDeath[0];
                            earliest_neutron_hit.piDeath[1] = t_piDeath[1];
                            earliest_neutron_hit.piDeath[2] = t_piDeath[2];

                            earliest_neutron_hit.protonDeath[0] = t_protonDeath[0];
                            earliest_neutron_hit.protonDeath[1] = t_protonDeath[1];
                            earliest_neutron_hit.protonDeath[2] = t_protonDeath[2];

                            earliest_neutron_hit.neutronHit[0] = t_neutronHitX[n_neutronHit];
                            earliest_neutron_hit.neutronHit[1] = t_neutronHitY[n_neutronHit];
                            earliest_neutron_hit.neutronHit[2] = t_neutronHitZ[n_neutronHit];
                            earliest_neutron_hit.neutronTrueT = t_neutronHitT[n_neutronHit];

                            earliest_neutron_hit.neutronStartingPoint[0] = t_neutronStartingPointX[n_neutronHit];
                            earliest_neutron_hit.neutronStartingPoint[1] = t_neutronStartingPointY[n_neutronHit];
                            earliest_neutron_hit.neutronStartingPoint[2] = t_neutronStartingPointZ[n_neutronHit];

                            earliest_neutron_hit.neutronParentId = t_neutronParentId[n_neutronHit];
                            earliest_neutron_hit.neutronParentPdg = t_neutronParentPDG[n_neutronHit];

                            earliest_neutron_hit.vtxTime = t_vtxTime;
                            earliest_neutron_hit.isEmpty = 0;
                        }
                    }
                }
            }   //end of n_neutronhit iterate

            if(earliest_neutron_hit.isEmpty == false)
            {
                if(earliest_neutron_hit.neutronStartingPoint[0] != -1 
                        &&earliest_neutron_hit.neutronStartingPoint[1] != -1
                        &&earliest_neutron_hit.neutronStartingPoint[2] != -1)
                {
                    if(earliest_neutron_hit.neutronParentId == -1 ||earliest_neutron_hit.neutronParentId == 0)
                    {
                        hist_signal->Fill(earliest_neutron_hit.trackLength,earliest_neutron_hit.timeWindow);
                        dist_sig_sp_vtx1->Fill(GetDistance(earliest_neutron_hit.neutronStartingPoint,earliest_neutron_hit.vtxSignal));
                        dist_sig_sp_nh1->Fill(GetDistance(earliest_neutron_hit.neutronStartingPoint,earliest_neutron_hit.neutronHit));

                        //vector from vtx to neutron hit
                        vec_vtx_to_sig[0] = earliest_neutron_hit.neutronHit[0]-earliest_neutron_hit.vtxSignal[0];
                        vec_vtx_to_sig[1] = earliest_neutron_hit.neutronHit[1]-earliest_neutron_hit.vtxSignal[1];
                        vec_vtx_to_sig[2] = earliest_neutron_hit.neutronHit[2]-earliest_neutron_hit.vtxSignal[2];

                        //z-axis
                        float b[3];
                        b[0] = 0;
                        b[1] = 0;
                        b[2] = 1;

                        angle_vtx_signal->Fill(GetAngle(vec_vtx_to_sig,b));
                    }

                    if(earliest_neutron_hit.neutronParentId > 0)
                    {
                        dist_vtx_to_nh_secondary->Fill(GetDistance(earliest_neutron_hit.neutronHit,earliest_neutron_hit.vtxSignal));
                        //hist_bkg_1->Fill(earliest_neutron_hit.trackLength,earliest_neutron_hit.timeWindow);
                        if(abs(earliest_neutron_hit.neutronParentPdg) == 211 || earliest_neutron_hit.neutronParentPdg == 111) //pion
                        {
                            dist_vtx_to_nh_pion->Fill(GetDistance(earliest_neutron_hit.neutronHit,earliest_neutron_hit.vtxSignal));
                            dist_sig_sp_vtx_pion->Fill(GetDistance(earliest_neutron_hit.neutronStartingPoint,earliest_neutron_hit.vtxSignal));
                            dist_sig_sp_nh_pion->Fill(GetDistance(earliest_neutron_hit.neutronStartingPoint,earliest_neutron_hit.neutronHit));
                            //vector from pi death point to earliest neutron hit from it
                            float vec_piDeath_neutron_hit[3];
                            vec_piDeath_neutron_hit[0] = earliest_neutron_hit.neutronHit[0]-earliest_neutron_hit.piDeath[0];
                            vec_piDeath_neutron_hit[1] = earliest_neutron_hit.neutronHit[1]-earliest_neutron_hit.piDeath[1];
                            vec_piDeath_neutron_hit[2] = earliest_neutron_hit.neutronHit[2]-earliest_neutron_hit.piDeath[2];
                            //vector from FV vertex to pi death point
                            float vec_vtx_piDeath[3];
                            vec_vtx_piDeath[0] = earliest_neutron_hit.piDeath[0]-earliest_neutron_hit.vtxSignal[0];
                            vec_vtx_piDeath[1] = earliest_neutron_hit.piDeath[1]-earliest_neutron_hit.vtxSignal[1];
                            vec_vtx_piDeath[2] = earliest_neutron_hit.piDeath[2]-earliest_neutron_hit.vtxSignal[2];
                            angle_piDeath_neutron_hit->Fill(GetAngle(vec_piDeath_neutron_hit,vec_vtx_piDeath));
                        }

                        if(earliest_neutron_hit.neutronParentPdg == 2112)    //neutron
                        {
                            hist_bkg_1->Fill(earliest_neutron_hit.trackLength,earliest_neutron_hit.timeWindow);
                            dist_vtx_to_nh_neutron->Fill(GetDistance(earliest_neutron_hit.neutronHit,earliest_neutron_hit.vtxSignal));
                            dist_sig_sp_vtx_neutron->Fill(GetDistance(earliest_neutron_hit.neutronStartingPoint,earliest_neutron_hit.vtxSignal));
                            dist_sig_sp_nh_neutron->Fill(GetDistance(earliest_neutron_hit.neutronStartingPoint,earliest_neutron_hit.neutronHit));
                        }

                        if(earliest_neutron_hit.neutronParentPdg == 2212)    //proton
                        {
                            dist_vtx_to_nh_proton->Fill(GetDistance(earliest_neutron_hit.neutronHit,earliest_neutron_hit.vtxSignal));
                            dist_sig_sp_vtx_proton->Fill(GetDistance(earliest_neutron_hit.neutronStartingPoint,earliest_neutron_hit.vtxSignal));
                            dist_sig_sp_nh_proton->Fill(GetDistance(earliest_neutron_hit.neutronStartingPoint,earliest_neutron_hit.neutronHit));
                        }
                        
                        dist_sig_sp_vtx3->Fill(GetDistance(earliest_neutron_hit.neutronStartingPoint,earliest_neutron_hit.vtxSignal));
                        dist_sig_sp_nh3->Fill(GetDistance(earliest_neutron_hit.neutronStartingPoint,earliest_neutron_hit.neutronHit));

                        //do when if it's secondary
                        //vector from FV vertex to secondary vertex
                        vec_vtx_to_secondary_vertex[0] = earliest_neutron_hit.neutronStartingPoint[0]-earliest_neutron_hit.vtxSignal[0];
                        vec_vtx_to_secondary_vertex[1] = earliest_neutron_hit.neutronStartingPoint[1]-earliest_neutron_hit.vtxSignal[1];
                        vec_vtx_to_secondary_vertex[2] = earliest_neutron_hit.neutronStartingPoint[2]-earliest_neutron_hit.vtxSignal[2];

                        //vector from secondary vertex to neutron hit
                        vec_secondary_vertex_to_neutron_hit[0] = earliest_neutron_hit.neutronHit[0]-earliest_neutron_hit.neutronStartingPoint[0];
                        vec_secondary_vertex_to_neutron_hit[1] = earliest_neutron_hit.neutronHit[1]-earliest_neutron_hit.neutronStartingPoint[1];
                        vec_secondary_vertex_to_neutron_hit[2] = earliest_neutron_hit.neutronHit[2]-earliest_neutron_hit.neutronStartingPoint[2];

                        //vector from vtx to secondary neutron hit
                        vec_vtx_to_secondary_neutron[0] = earliest_neutron_hit.neutronHit[0]-earliest_neutron_hit.vtxSignal[0];
                        vec_vtx_to_secondary_neutron[1] = earliest_neutron_hit.neutronHit[1]-earliest_neutron_hit.vtxSignal[1];
                        vec_vtx_to_secondary_neutron[2] = earliest_neutron_hit.neutronHit[2]-earliest_neutron_hit.vtxSignal[2];

                        float z[3];
                        z[0] = 0;
                        z[1] = 0;
                        z[2] = 1;

                        angle_event->Fill(GetAngle(vec_vtx_to_secondary_vertex,vec_secondary_vertex_to_neutron_hit));
                        angle_vtx_secondary->Fill(GetAngle(vec_vtx_to_secondary_neutron,z));
                    }
                }
            }
        }
    }       //end of event iterate
    _file->Close();
}

void neutron()
{
    int endPROD, beginPROD, filenum;
    cout<<"PROD begin :"<<endl;
    cin>>beginPROD;
    cout<<"PROD end :"<<endl;
    cin>>endPROD;
    cout<<"filenum :"<<endl;
    cin>>filenum;
    cout<<"start"<<endl;
    for(int j = beginPROD; j <endPROD+1; j++)
    {
        for(int i = 2; i <filenum; i++) //test_1 is not
        {
            cout<<"\033[1APROD"<<j<<": "<<(double)(i*100/filenum)<<"%\033[1000D"<<endl;
            //analyze(Form("/pnfs/dune/persistent/users/gyang/3DST/dump/standardGeo12/PROD%d/FHC_%d.root",j,i));
            //analyze(Form("/pnfs/dune/persistent/users/gyang/3DST/dump/standardGeo12/PROD%d/RHC_%d.root",j,i));
            //        test_analyze(Form("/pnfs/dune/persistent/users/gyang/3DST/dump/standardGeo12/PROD%d/FHC_%d.root",j,i));
            //        test_analyze(Form("/pnfs/dune/persistent/users/gyang/3DST/dump/standardGeo12/PROD%d/RHC_%d.root",j,i));
            //analyze(Form("/pnfs/dune/persistent/users/gyang/3DST/dump/standardGeo12/PROD%d/FHC_%d_test.root",j,i));
            //analyze(Form("/pnfs/dune/persistent/users/gyang/3DST/dump/standardGeo12/PROD%d/RHC_%d_test.root",j,i));
            //num_interaction(Form("/Users/gwon/Geo12/PROD%d/FHC_%d.root",j,i));
            //num_interaction(Form("/Users/gwon/Geo12/PROD%d/RHC_%d.root",j,i));
            //test_test_analyze(Form("/pnfs/dune/persistent/users/gyang/3DST/dump/standardGeo12/PROD%d/RHC_%d_test.root",j,i));
            test_test_analyze(Form("/Users/gwon/Geo12/PROD%d/RHC_%d_test.root",j,i));
        }
        cout<<endl;
    }
    cout<<"end"<<endl;
    cout<<"nubmer_of_CC: "<<number_of_CC<<endl;
    cout<<"number_of_file: "<<num_file<<endl;
    cout<<"number_of_interaction: "<<all_interaction<<endl;


    TFile * fi1 = new TFile("background.root","RECREATE");
    hist_signal->Write();
    hist_bkg_out3DST->Write();
    hist_bkg_NC->Write();
    hist_bkg_1->Write();
    hist_bkg_1_pion->Write();
    hist_bkg_1_neutron->Write();
    hist_bkg_1_proton->Write();
    hist_bkg_1_other->Write();
    hist_bkg_out3DST_largeTime->Write();
    hist_bkg_NC_largeTime->Write();
    hist_bkg_1_largeTime->Write();
    KE_primary->Write();
    KE_secondary->Write();
    neutronParentPDG->Write();
    neutronParentPDG_case4->Write();
    first_n_position_XY->Write();
    first_n_position_YZ->Write();
    first_n_position_XZ->Write();
    dist_sp_vtx->Write();
    dist_sp_nh->Write();
    dist_sig_sp_vtx->Write();
    dist_sig_sp_nh->Write();
    dist_sig_sp_vtx1->Write();
    dist_sig_sp_nh1->Write();
    dist_sig_sp_vtx2->Write();
    dist_sig_sp_nh2->Write();
    dist_sig_sp_vtx3->Write();
    dist_sig_sp_nh3->Write();
    dist_sig_sp_vtx_pion->Write();
    dist_sig_sp_nh_pion->Write();
    dist_sig_sp_vtx_neutron->Write();
    dist_sig_sp_nh_neutron->Write();
    dist_sig_sp_vtx_proton->Write();
    dist_sig_sp_nh_proton->Write();
    dist_vtx_to_nh_secondary->Write();
    dist_vtx_to_nh_pion->Write();
    dist_vtx_to_nh_neutron->Write();
    dist_vtx_to_nh_proton->Write();
    angle_vtx_signal->Write();
    angle_vtx_secondary->Write();
    angle_piDeath_neutron_hit->Write();
    fi1->Close();

    TCanvas * can = new TCanvas;
    can->Divide(2,2);
    can->cd(1);
    hist_signal->Draw("colz");
    can->cd(2);
    hist_bkg_out3DST->Draw("colz");
    can->cd(3);
    hist_bkg_NC->Draw("colz");
    can->cd(4);
    hist_bkg_1->Draw("colz");
    can->SaveAs("4plots.pdf");

    TCanvas * can1 = new TCanvas;
    can1->Divide(2,2);
    can1->cd(1);
    hist_bkg_1->Draw("colz");
    can1->cd(2);
    hist_bkg_1_pion->Draw("colz");
    can1->cd(3);
    hist_bkg_1_neutron->Draw("colz");
    can1->cd(4);
    hist_bkg_1_proton->Draw("colz");
    can1->SaveAs("3_1.pdf");

    TCanvas * can4 = new TCanvas;
    hist_bkg_1_other->Draw("colz");
    can4->SaveAs("other.pdf");

    TCanvas * can2 = new  TCanvas;
    KE_primary->Draw();
    TCanvas * can3 = new  TCanvas;
    KE_secondary->Draw();

    TCanvas * can5 = new TCanvas("asdf","asdf",1500,600);
    can5->Divide(3,1);
    can5->cd(1);
    first_n_position_XY->Draw("colz");
    can5->cd(2);
    first_n_position_YZ->Draw("colz");
    can5->cd(3);
    first_n_position_XZ->Draw("colz");
    can5->SaveAs("neutron_position.pdf");

    TCanvas * can6 = new TCanvas;
    can6->Divide(3,2);
    can6->cd(1);
    dist_sig_sp_vtx1->Draw();
    can6->cd(2);
    dist_sig_sp_nh1->Draw();
    can6->cd(4);
    dist_sig_sp_vtx3->Draw();
    can6->cd(5);
    dist_sig_sp_nh3->Draw();
    can6->cd(6);
    dist_vtx_to_nh_secondary->Draw();
    /*
       can6->cd(5);
       dist_sig_sp_vtx3->Draw();
       can6->cd(6);
       dist_sig_sp_nh3->Draw();
     */
    can6->SaveAs("test.pdf");

    TCanvas * can7 = new TCanvas;
    can7->Divide(3,3);
    can7->cd(1);
    dist_sig_sp_vtx_pion->Draw();
    can7->cd(2);
    dist_sig_sp_nh_pion->Draw();
    can7->cd(3);
    dist_vtx_to_nh_pion->Draw();
    can7->cd(4);
    dist_sig_sp_vtx_neutron->Draw();
    can7->cd(5);
    dist_sig_sp_nh_neutron->Draw();
    can7->cd(6);
    dist_vtx_to_nh_neutron->Draw();
    can7->cd(7);
    dist_sig_sp_vtx_proton->Draw();
    can7->cd(8);
    dist_sig_sp_nh_proton->Draw();
    can7->cd(9);
    dist_vtx_to_nh_proton->Draw();
    can7->SaveAs("secondary distribution.pdf");

    angle_event->Scale(1/angle_event->GetEntries(),"nosw2");
    double x[10];
    double y[10];
    for(int i = 0; i < 11; i++)
    {
        x[i] = angle_event->GetBinContent(i);
        for(int j = 1; j < i+1; j++)
        {
            y[i] = y[i]+x[j];
        }
    }

    for(int i = 1; i < 11; i++)
    {
        angle_event_accumulated->SetBinContent(i,y[i]);
    }

    TCanvas * can8 = new TCanvas;
    angle_event->SetStats(false);
    angle_event->Draw();
    can8->SaveAs("angle.pdf");

    TCanvas * can9 = new TCanvas;
    angle_event_accumulated->SetStats(false);
    angle_event_accumulated->Draw();
    can9->SaveAs("angle_accumulated.pdf");

    TCanvas * can10 = new TCanvas;
    angle_vtx_signal->Scale(1/angle_vtx_signal->GetEntries(),"nosw2");
    angle_vtx_signal->SetStats(false);
    angle_vtx_signal->Draw();
    can10->SaveAs("angle_vtx_signal.pdf");

    TCanvas * can11 = new TCanvas;
    angle_vtx_secondary->Scale(1/angle_vtx_secondary->GetEntries(),"nosw2");
    angle_vtx_secondary->SetStats(false);
    angle_vtx_secondary->Draw();
    can11->SaveAs("angle_vtx_secondary.pdf");

    TCanvas * can12 = new TCanvas;
    angle_piDeath_neutron_hit->Scale(1/angle_piDeath_neutron_hit->GetEntries(), "nosw2 ");
    angle_piDeath_neutron_hit->SetStats(false);
    angle_piDeath_neutron_hit->Draw();
    can12->SaveAs("angle_piDeath_neutron_hit.pdf");
}
