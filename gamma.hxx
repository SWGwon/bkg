#include <string>
#include <utility>
#include <vector>
#include <thread>

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
#include "TSystem.h"
#include <TError.h>



using namespace std;
//histograms{
TH1F * energy_of_signal = new TH1F("energy_of_signal" , "energy of signal;energy [MeV]",100,0,10);
TH1F * energy_of_gamma = new TH1F("energy_of_gamma" , "energy of gamma;energy [MeV]",100,0,10);
TH1F * energy_of_secondary = new TH1F("energy_of_secondary" , "energy of secondary;energy [MeV]",100,0,10);

TH1F * gammaPDG = new TH1F("gammaPDG", "PDG", 5000, -2500, 2500);
TH1F * num_of_primary_gamma = new TH1F;
TH1F * num_of_other_gamma = new TH1F;

TH1F * hitPDG_signal_big = new TH1F("","",3000,-500,2500);
TH1F * hitPDG_signal_small = new TH1F("","",3000,-500,2500);

TH2F * hitPDG = new TH2F("","",3000,-500,2500,500,0,15);


TH1F * beta_of_signal = new TH1F("beta_of_signal", "velocity/c of signal;beta", 10, 0, 1);
TH1F * beta_of_gamma = new TH1F("beta_of_gamma", "velocity/c of gamma;beta", 10, 0, 1);
TH1F * beta_of_secondary = new TH1F("beta_of_secondary", "velocity/c of secondary;beta", 10, 0, 1);

TH2F * hist_bkg_out3DST = new TH2F("hist_bkg_out3DST", "out3DST background;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_NC = new TH2F("hist_bkg_NC", "NC background;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);

//lever arm vs time 
TH2F * hist_sig_arm_vs_time = new TH2F("hist_sig_arm_vs_time", "signal;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_1_arm_vs_time = new TH2F("hist_bkg_1_arm_vs_time", "secondary background;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_gamma_arm_vs_time = new TH2F("hist_bkg_gamma_arm_vs_time", "gamma background;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_1_gamma_arm_vs_time = new TH2F("hist_bkg_1_gamma_arm_vs_time", "secondary+gamma background;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);

//lever arm vs time with cut
TH2F * hist_sig_arm_vs_time_linear_cut = new TH2F("hist_sig_arm_vs_time_linear_cut", "signal;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_1_arm_vs_time_linear_cut = new TH2F("hist_bkg_1_arm_vs_time_linear_cut", "secondary background;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_gamma_arm_vs_time_linear_cut = new TH2F("hist_bkg_gamma_arm_vs_time_linear_cut", "gamma background;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_1_gamma_arm_vs_time_linear_cut = new TH2F("hist_bkg_1_gamma_arm_vs_time_linear_cut", "secondary+gamma background;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);

//angla vs lever arm 
TH2F * hist_sig_ang_vs_dis = new TH2F("hist_sig_ang_vs_dis", "signal;angle(pi);Lever arm(cm)",10,0,1,10,0,200);
TH2F * hist_bkg_1_ang_vs_dis = new TH2F("hist_bkg_1_ang_vs_dis", "secondary background;angle(pi);Lever arm(cm)",10,0,1,10,0,200);
TH2F * hist_bkg_gamma_ang_vs_dis = new TH2F("hist_bkg_gamma_ang_vs_dis", "gamma background;angle(pi);Lever arm(cm)",10,0,1,10,0,200);
TH2F * hist_bkg_1_gamma_ang_vs_dis = new TH2F("hist_bkg_1_gamma_ang_vs_dis", "secondary+gamma background;angle(pi);Lever arm(cm)",10,0,1,10,0,200);

//angle vs lever arm with cut
TH2F * hist_sig_ang_vs_dis_linear_cut = new TH2F("hist_sig_ang_vs_dis_linear_cut", "signal;angle(pi);Lever arm(cm)",10,0,1,10,0,200);
TH2F * hist_bkg_gamma_ang_vs_dis_linear_cut = new TH2F("hist_bkg_gamma_ang_vs_dis_linear_cut", "gamma background;angle(pi);Lever arm(cm)",10,0,1,10,0,200);
TH2F * hist_bkg_1_ang_vs_dis_linear_cut = new TH2F("hist_bkg_1_ang_vs_dis_linear_cut", "secondagry background;angle(pi);Lever arm(cm)",10,0,1,10,0,200);
TH2F * hist_bkg_1_gamma_ang_vs_dis_linear_cut = new TH2F("hist_bkg_1_gamma_ang_vs_dis_linear_cut", "secondary+gamma background;angle(pi);Lever arm(cm)",10,0,1,10,0,200);

TH1F * distance_vtx_to_deathpoint = new TH1F("distance_vtx_to_deathpoint","distance;distance",20,0,200);


bool is_inFV = false;       //check if vertex is in FV
bool is_in3DST = false;     //check if vertex is in 3DST

float energyHitCut = 0.1; //energy deposit threshold for cube

int number_of_CC = 0;
int number_of_secondary_pion = 0;
int number_of_secondary_proton = 0;
int number_of_secondary_neutron = 0;
int number_of_secondary_other = 0;

//change this part to do slope, intercept test
int cut_slope = 0;
int cut_y_intercept = 0;
//channel type
int num_fspi = 0;   //number of fs charged pion
int num_fsp = 1;    //number of fs proton

int num_primary_gamma = 0;
int num_secondary_gamma = 0;
int num_gamma_from_pion = 0;
int num_gamma_from_muon = 0;
int num_gamma_from_other = 0;

int iftest = 0;
double test_cut_angle = 0.95;
double test_cut_distance = 20;

const double c_velocity = 29.9792458;

class Hit_t 
{
    public:
        float timeWindow,           // time windows of the hit
              timeWindowSmear,      // smeared time window 
              timeSmear,        // smear time
              energyDeposit,        // energy deposited by the neutron
              trackLength,          // lever arm
              trueRec,      // true reconstructed energy
              smearRec,
              vtxSignal[3],     // neutrino vertex position of the neutron
              vtxTime,      // neutrino  vertex time
              trueE,    //neutron true energy
              CubeE,    //neutron cube energy
              trueT,    //neutron true time
              hitPDG,    //hit PDG
              piDeath[3],      //pion death
              protonDeath[3];      //proton Death

        //neutron hit position
        float hit[3];

        float startingPoint[3];

        int bkgLoc,         // neutrino vertex position
            parentId,    // Where the hit come from
            parentPdg;   // PDG of parent

        bool isTherePion50,     // Is there a pion with KE > 50 MeV in FS particles
             isThereProton300;      // Is there a proton with KE > 300 MeV in FS particles
        bool isEmpty;
        bool isNeutron;
        bool isGamma;
        bool isFromPion;
        bool isFromProton;

        Hit_t()
            :timeWindow(0),
            timeSmear(0),    
            energyDeposit(0),
            trackLength(0),  
            trueRec(0),      
            smearRec(0),
            vtxTime(0),      
            trueE(0), 
            CubeE(0),
            trueT(0), 
            hitPDG(0),
            bkgLoc(125124123),        
            parentId(123124123),
            parentPdg(123123123),
            isTherePion50(0),  
            isThereProton300(0),
            isEmpty(1),
            isFromPion(0),
            isFromProton(0),
            isNeutron(0),
            isGamma(0)
    {
        for(int i = 0; i < 3; i++)
        {
            this->vtxSignal[i] = 0; 
            this->piDeath[i] = 0;   
            this->protonDeath[i] = 0; 
            this->hit[i] = 0;
            this->startingPoint[i] = 0;
        }
    }

        ~Hit_t() {}
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

    template <class T>
T save(T x,TCanvas* can)
{
    x->Draw("colz");
    x->Write();
    can->SaveAs(Form("%s.pdf",TO_STRING(x)));
    can->Clear();
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


void analyze(string file)
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
    float t_neutronCubeE[1000];
    float t_neutronHitSmearT[1000];
    float t_neutronHitPDG[1000];

    float t_vtx[3], t_vtxTime;
    float t_piDeath[3], t_protonDeath[3];

    float vec_piDeath_to_hit[3];
    float vec_protonDeath_to_hit[3];
    float vec_vtx_to_piDeath[3];
    float vec_vtx_to_protonDeath[3];

    float t_gammaHitX[1000], t_gammaHitY[1000], t_gammaHitZ[1000];
    float t_gammaStartingPointX[1000], t_gammaStartingPointY[1000], t_gammaStartingPointZ[1000];
    float t_gammaHitT[1000], t_gammaParentId[1000], t_gammaParentPDG[1000];
    float t_gammaHitE[1000], t_gammaTrueE[1000];
    float t_gammaCubeE[1000];
    float t_gammaHitSmearT[1000];
    float t_gammaHitPDG[1000];

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
    tree->SetBranchAddress("neutronCubeE", &t_neutronCubeE);
    tree->SetBranchAddress("neutronHitSmearT", &t_neutronHitSmearT);
    tree->SetBranchAddress("neutronHitPDG", &t_neutronHitPDG);
    tree->SetBranchAddress("vtx", &t_vtx);
    tree->SetBranchAddress("vtxTime", &t_vtxTime);
    tree->SetBranchAddress("nFS", &t_nFS);
    tree->SetBranchAddress("fsPdg", &t_fsPdg);
    tree->SetBranchAddress("piDeath", &t_piDeath);
    tree->SetBranchAddress("protonDeath", &t_protonDeath);

    tree->SetBranchAddress("gammaHitX", &t_gammaHitX);
    tree->SetBranchAddress("gammaHitY", &t_gammaHitY);
    tree->SetBranchAddress("gammaHitZ", &t_gammaHitZ);
    tree->SetBranchAddress("gammaStartingPointX", &t_gammaStartingPointX);
    tree->SetBranchAddress("gammaStartingPointY", &t_gammaStartingPointY);
    tree->SetBranchAddress("gammaStartingPointZ", &t_gammaStartingPointZ);
    tree->SetBranchAddress("gammaHitT", &t_gammaHitT);
    tree->SetBranchAddress("gammaParentId", &t_gammaParentId);
    tree->SetBranchAddress("gammaParentPDG", &t_gammaParentPDG);
    tree->SetBranchAddress("gammaHitE", &t_gammaHitE);
    tree->SetBranchAddress("gammaTrueE", &t_gammaTrueE);
    tree->SetBranchAddress("gammaCubeE", &t_gammaCubeE);
    tree->SetBranchAddress("gammaHitT", &t_gammaHitT);
    tree->SetBranchAddress("gammaHitSmearT", &t_gammaHitSmearT);
    tree->SetBranchAddress("gammaHitPDG", &t_gammaHitPDG);

    int nevents = tree->GetEntries();

    for(int event = 0; event < nevents; event++)
    {
        Hit_t earliest_hit;
        Hit_t earliest_neutron_hit;
        Hit_t earliest_gamma_hit;

        int num_pi = 0;
        int num_proton = 0;
        tree->GetEntry(event);
        //if(abs(t_vtx[0]) < 50 && abs(t_vtx[1]) < 50 && abs(t_vtx[2]) < 50)
        if(abs(t_vtx[0]) < 100 && abs(t_vtx[1]) < 100 && t_vtx[2] < 140 && t_vtx[2] > -40)
            //if(1)
        {
            bool is_CC = false;
            bool is_pion = false;
            bool is_proton = false;

            for(int inFS = 0; inFS < t_nFS; inFS++)
            {
                if(abs(t_fsPdg[inFS]) == 11 || abs(t_fsPdg[inFS]) == 13)    //electronPDG=11,muonPDG=13
                {
                    is_CC = true;
                    break;
                }
            }

            for(int inFS = 0; inFS < t_nFS; inFS++)
            {
                if(abs(t_fsPdg[inFS]) == 211)    //pionPDG=+-211
                {
                    is_pion = true;
                    num_pi++;
                }
            }

            for(int inFS = 0; inFS < t_nFS; inFS++)
            {
                if(abs(t_fsPdg[inFS]) == 2212)    //protonPDG=+-211
                {
                    is_proton = true;
                    num_proton++;
                }
            }


            if(!is_CC)
                continue;
            if(num_pi != num_fspi)
                continue;
            if(num_proton != num_fsp)
                continue;

            float temp_earliest_time = 1000000;
            for(int n_neutronHit = 0; n_neutronHit < 1000; n_neutronHit++)
            {
                //if(t_neutronHitSmearT[n_neutronHit] != 0 && t_neutronHitX[n_neutronHit] != 0 && t_neutronHitSmearT[n_neutronHit] < temp_earliest_time && t_neutronCubeE[n_neutronHit] > energyHitCut)
                if(t_neutronHitT[n_neutronHit] != 0 && t_neutronHitX[n_neutronHit] != 0 && t_neutronHitT[n_neutronHit] < temp_earliest_time && t_neutronHitE[n_neutronHit] > energyHitCut)
                {
                    temp_earliest_time = t_neutronHitT[n_neutronHit];
                    //look for a neutron hit in 3DST
                    /*
                       if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                       abs(t_neutronHitY[n_neutronHit]) < 120 && abs(t_neutronHitZ[n_neutronHit]) < 100) */
                    if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                            abs(t_neutronHitY[n_neutronHit]) < 120 && 
                            t_neutronHitZ[n_neutronHit] < 150 &&
                            t_neutronHitZ[n_neutronHit] > -50)
                    {
                        //calculate lever arm
                        float trackLength = pow(
                                pow(t_neutronHitX[n_neutronHit] - t_vtx[0],2)+
                                pow(t_neutronHitY[n_neutronHit] - t_vtx[1],2)+
                                pow(t_neutronHitZ[n_neutronHit] - t_vtx[2],2),0.5);

                        //calculate signal window; time of flight
                        float signalWindow = t_neutronHitT[n_neutronHit] - t_vtxTime;
                        float signalWindowSmear = t_neutronHitSmearT[n_neutronHit] - t_vtxTime;

                        //Fix a bug from edep-sim
                        if(signalWindow == 1)
                            signalWindow = 0.5;

                        if(signalWindow > 0)
                        {
                            earliest_neutron_hit.timeWindow = signalWindow;
                            earliest_neutron_hit.timeWindowSmear = signalWindowSmear;
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

                            earliest_neutron_hit.hit[0] = t_neutronHitX[n_neutronHit];
                            earliest_neutron_hit.hit[1] = t_neutronHitY[n_neutronHit];
                            earliest_neutron_hit.hit[2] = t_neutronHitZ[n_neutronHit];
                            earliest_neutron_hit.trueT = t_neutronHitT[n_neutronHit];
                            earliest_neutron_hit.CubeE = t_neutronCubeE[n_neutronHit];

                            earliest_neutron_hit.startingPoint[0] = t_neutronStartingPointX[n_neutronHit];
                            earliest_neutron_hit.startingPoint[1] = t_neutronStartingPointY[n_neutronHit];
                            earliest_neutron_hit.startingPoint[2] = t_neutronStartingPointZ[n_neutronHit];

                            earliest_neutron_hit.parentId = t_neutronParentId[n_neutronHit];
                            earliest_neutron_hit.parentPdg = t_neutronParentPDG[n_neutronHit];
                            earliest_neutron_hit.hitPDG = t_neutronHitPDG[n_neutronHit];

                            earliest_neutron_hit.vtxTime = t_vtxTime;
                            earliest_neutron_hit.isEmpty = 0;
                            if(t_neutronStartingPointX[n_neutronHit] == t_piDeath[0])
                                earliest_neutron_hit.isFromPion = 1;
                            else
                                earliest_neutron_hit.isFromPion = 0;
                            if(t_neutronStartingPointX[n_neutronHit] == t_protonDeath[0])
                                earliest_neutron_hit.isFromProton = 1;
                            else
                                earliest_neutron_hit.isFromProton = 0;
                            earliest_neutron_hit.isNeutron = 1;
                        }
                    }
                }
            }

            float temp_earliest_time_for_gamma = 1000000;
            for(int n_gammaHit = 0; n_gammaHit < 1000; n_gammaHit++)
            {
                //if(t_gammaHitSmearT[n_gammaHit] != 0 && t_gammaHitX[n_gammaHit] != 0 && t_gammaHitSmearT[n_gammaHit] < temp_earliest_time_for_gamma && t_gammaCubeE[n_gammaHit] > energyHitCut)
                if(t_gammaHitT[n_gammaHit] != 0 && t_gammaHitX[n_gammaHit] != 0 && t_gammaHitT[n_gammaHit] < temp_earliest_time_for_gamma && t_gammaHitE[n_gammaHit] > energyHitCut)
                {
                    temp_earliest_time_for_gamma = t_gammaHitT[n_gammaHit];
                    if(abs(t_gammaHitX[n_gammaHit]) < 120 && 
                            abs(t_gammaHitY[n_gammaHit]) < 120 && 
                            t_gammaHitZ[n_gammaHit] < 150 &&
                            t_gammaHitZ[n_gammaHit] > -50)
                    {
                        //calculate lever arm
                        float trackLength = pow(
                                pow(t_gammaHitX[n_gammaHit] - t_vtx[0],2)+
                                pow(t_gammaHitY[n_gammaHit] - t_vtx[1],2)+
                                pow(t_gammaHitZ[n_gammaHit] - t_vtx[2],2),0.5);

                        //calculate signal window; time of flight
                        float signalWindow = t_gammaHitT[n_gammaHit] - t_vtxTime;
                        float signalWindowSmear = t_gammaHitSmearT[n_gammaHit] - t_vtxTime;

                        //Fix a bug from edep-sim
                        if(signalWindow == 1)
                            signalWindow = 0.5;

                        if(signalWindow > 0)
                        {
                            earliest_gamma_hit.timeWindow = signalWindow;
                            earliest_neutron_hit.timeWindowSmear = signalWindowSmear;
                            earliest_gamma_hit.trackLength = trackLength;
                            earliest_gamma_hit.energyDeposit = t_gammaHitE[n_gammaHit];

                            earliest_gamma_hit.vtxSignal[0] = t_vtx[0];
                            earliest_gamma_hit.vtxSignal[1] = t_vtx[1];
                            earliest_gamma_hit.vtxSignal[2] = t_vtx[2];

                            earliest_gamma_hit.piDeath[0] = t_piDeath[0];
                            earliest_gamma_hit.piDeath[1] = t_piDeath[1];
                            earliest_gamma_hit.piDeath[2] = t_piDeath[2];

                            earliest_gamma_hit.protonDeath[0] = t_protonDeath[0];
                            earliest_gamma_hit.protonDeath[1] = t_protonDeath[1];
                            earliest_gamma_hit.protonDeath[2] = t_protonDeath[2];

                            earliest_gamma_hit.hit[0] = t_gammaHitX[n_gammaHit];
                            earliest_gamma_hit.hit[1] = t_gammaHitY[n_gammaHit];
                            earliest_gamma_hit.hit[2] = t_gammaHitZ[n_gammaHit];
                            earliest_gamma_hit.trueT = t_gammaHitT[n_gammaHit];
                            earliest_gamma_hit.CubeE = t_gammaCubeE[n_gammaHit];

                            earliest_gamma_hit.startingPoint[0] = t_gammaStartingPointX[n_gammaHit];
                            earliest_gamma_hit.startingPoint[1] = t_gammaStartingPointY[n_gammaHit];
                            earliest_gamma_hit.startingPoint[2] = t_gammaStartingPointZ[n_gammaHit];

                            earliest_gamma_hit.parentId = t_gammaParentId[n_gammaHit];
                            earliest_gamma_hit.parentPdg = t_gammaParentPDG[n_gammaHit];
                            earliest_gamma_hit.hitPDG = t_gammaHitPDG[n_gammaHit];

                            earliest_gamma_hit.vtxTime = t_vtxTime;
                            earliest_gamma_hit.isEmpty = 0;
                            if(t_gammaStartingPointX[n_gammaHit] == t_piDeath[0])
                                earliest_gamma_hit.isFromPion = 1;
                            else
                                earliest_gamma_hit.isFromPion = 0;
                            if(t_gammaStartingPointX[n_gammaHit] == t_protonDeath[0])
                                earliest_gamma_hit.isFromProton = 1;
                            else
                                earliest_gamma_hit.isFromProton = 0;
                            earliest_gamma_hit.isGamma = 1;
                        }
                    }
                }
            }

            if(earliest_gamma_hit.isEmpty == false && temp_earliest_time_for_gamma < temp_earliest_time)
            {
                earliest_hit = earliest_gamma_hit;
                //energy_of_gamma->Fill(earliest_hit.CubeE);
                energy_of_gamma->Fill(earliest_hit.energyDeposit);
            }
            if(earliest_neutron_hit.isEmpty == false && temp_earliest_time_for_gamma > temp_earliest_time)
            {
                earliest_hit = earliest_neutron_hit;

                if(earliest_hit.parentId == -1 || earliest_hit.parentId == 0)
                {
                    //energy_of_signal->Fill(earliest_hit.CubeE);
                    energy_of_signal->Fill(earliest_hit.energyDeposit);
                    //if(earliest_hit.CubeE > 1.5 && earliest_hit.CubeE < 2.2)
                    //    hitPDG_signal_big->Fill(earliest_hit.hitPDG);
                    //if(earliest_hit.CubeE < 0.8)
                    //    hitPDG_signal_small->Fill(earliest_hit.hitPDG);
                    if(earliest_hit.energyDeposit > 1.5 && earliest_hit.energyDeposit < 2.2)
                        hitPDG_signal_big->Fill(earliest_hit.hitPDG);
                    if(earliest_hit.energyDeposit < 0.8)
                        hitPDG_signal_small->Fill(earliest_hit.hitPDG);
                    hitPDG->Fill(earliest_hit.hitPDG,earliest_hit.energyDeposit);
                }
                else
                    //energy_of_secondary->Fill(earliest_hit.CubeE);
                    energy_of_secondary->Fill(earliest_hit.energyDeposit);
            }

            if(earliest_hit.isGamma)
            {
                gammaPDG->Fill(earliest_hit.parentPdg);
                if(earliest_gamma_hit.parentId == 0 || earliest_gamma_hit.parentId == -1)
                {
                    num_primary_gamma += 1;
                    num_of_primary_gamma->Fill(0);
                }
                if(earliest_gamma_hit.parentId > 0)
                {
                    num_secondary_gamma += 1;
                    num_of_other_gamma->Fill(0);
                }

                if(abs(earliest_gamma_hit.parentPdg) == 211)
                    num_gamma_from_pion += 1;
                else if(earliest_gamma_hit.parentPdg == 13)
                    num_gamma_from_muon += 1;
                else
                    num_gamma_from_other += 1;
            }

            if(earliest_hit.isEmpty == false && earliest_hit.CubeE != 0)
            {
                if(earliest_hit.parentId > 0)
                {
                    if(abs(earliest_hit.parentPdg) == 211)
                        number_of_secondary_pion += 1;
                    else if(earliest_hit.parentPdg == 2212)
                        number_of_secondary_proton += 1;
                    else if(earliest_hit.parentPdg == 2112)
                        number_of_secondary_neutron += 1;
                    else
                        number_of_secondary_other += 1;
                }
                if(earliest_hit.startingPoint[0] != -1 
                        &&earliest_hit.startingPoint[1] != -1
                        &&earliest_hit.startingPoint[2] != -1)
                {
                    for(int i = 0; i < 3; i++)
                    {
                        vec_piDeath_to_hit[i] = earliest_hit.hit[i]-earliest_hit.piDeath[i];

                        vec_vtx_to_piDeath[i] = earliest_hit.piDeath[i]-earliest_hit.vtxSignal[i];

                        vec_protonDeath_to_hit[i] = earliest_hit.hit[i]-earliest_hit.protonDeath[i];

                        vec_vtx_to_protonDeath[i] = earliest_hit.protonDeath[i]-earliest_hit.vtxSignal[i];
                    }


                    ////pi case
                    if(num_fspi == 1 && num_fsp ==0)
                    {
                        ////distance, angle cut 2D plot
                        distance_vtx_to_deathpoint->Fill(GetDistance(earliest_hit.vtxSignal,earliest_hit.piDeath));

                        //signal
                        if(earliest_hit.isNeutron && (earliest_hit.parentId == -1 ||earliest_hit.parentId == 0))
                        {
                            hist_sig_arm_vs_time->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                            hist_sig_ang_vs_dis->Fill(GetAngle(vec_piDeath_to_hit,vec_vtx_to_piDeath),earliest_hit.trackLength);
                            //linear cut
                            if(iftest == 0)
                            {
                                if(earliest_hit.trackLength-cut_slope*GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit)-cut_y_intercept > 0)
                                {
                                    hist_sig_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                                    hist_sig_ang_vs_dis_linear_cut->Fill(GetAngle(vec_piDeath_to_hit,vec_vtx_to_piDeath),earliest_hit.trackLength);
                                    //energy_of_signal->Fill(earliest_hit.CubeE);
                                    beta_of_signal->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                }
                            }
                            if(iftest == 1)
                            {
                                if(earliest_hit.trackLength < test_cut_distance && GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit) > test_cut_angle)
                                {
                                    hist_sig_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                                    hist_sig_ang_vs_dis_linear_cut->Fill(GetAngle(vec_piDeath_to_hit,vec_vtx_to_piDeath),earliest_hit.trackLength);
                                    //energy_of_signal->Fill(earliest_hit.CubeE);
                                    beta_of_signal->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                }
                            }
                        }

                        //background
                        if(earliest_hit.isGamma)
                        {
                            hist_bkg_gamma_arm_vs_time->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                            hist_bkg_1_gamma_arm_vs_time->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);

                            hist_bkg_gamma_ang_vs_dis->Fill(GetAngle(vec_piDeath_to_hit,vec_vtx_to_piDeath),earliest_hit.trackLength);
                            hist_bkg_1_gamma_ang_vs_dis->Fill(GetAngle(vec_piDeath_to_hit,vec_vtx_to_piDeath),earliest_hit.trackLength);

                            if(iftest == 0)
                            {
                                if(earliest_hit.trackLength-cut_slope*GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit)-cut_y_intercept > 0)
                                {
                                    hist_bkg_gamma_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                                    hist_bkg_1_gamma_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);

                                    hist_bkg_gamma_ang_vs_dis_linear_cut->Fill(GetAngle(vec_piDeath_to_hit,vec_vtx_to_piDeath),earliest_hit.trackLength);
                                    hist_bkg_1_gamma_ang_vs_dis_linear_cut->Fill(GetAngle(vec_piDeath_to_hit,vec_vtx_to_piDeath),earliest_hit.trackLength);
                                    //energy_of_gamma->Fill(earliest_hit.CubeE);
                                    //energy_of_secondary->Fill(earliest_hit.CubeE);
                                    beta_of_gamma->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                    beta_of_secondary->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                }
                            }
                            if(iftest == 1)
                            {
                                if(earliest_hit.trackLength < test_cut_distance && GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit) > test_cut_angle)
                                {
                                    hist_bkg_gamma_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                                    hist_bkg_1_gamma_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);

                                    hist_bkg_gamma_ang_vs_dis_linear_cut->Fill(GetAngle(vec_piDeath_to_hit,vec_vtx_to_piDeath),earliest_hit.trackLength);
                                    hist_bkg_1_gamma_ang_vs_dis_linear_cut->Fill(GetAngle(vec_piDeath_to_hit,vec_vtx_to_piDeath),earliest_hit.trackLength);
                                    //energy_of_gamma->Fill(earliest_hit.CubeE);
                                    //energy_of_secondary->Fill(earliest_hit.CubeE);
                                    beta_of_gamma->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                    beta_of_secondary->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                }
                            }
                        }
                        if(earliest_hit.parentId > 0)
                        {
                            if(abs(earliest_hit.piDeath[0]) < 120 && abs(earliest_hit.piDeath[1]) < 120 && earliest_hit.piDeath[2] < 150 && earliest_hit.piDeath[2] > -50)
                            {
                                if(earliest_hit.isFromPion && earliest_hit.isNeutron)
                                {
                                    hist_bkg_1_arm_vs_time->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                                    hist_bkg_1_gamma_arm_vs_time->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);

                                    hist_bkg_1_ang_vs_dis->Fill(GetAngle(vec_piDeath_to_hit,vec_vtx_to_piDeath),earliest_hit.trackLength);
                                    hist_bkg_1_gamma_ang_vs_dis->Fill(GetAngle(vec_piDeath_to_hit,vec_vtx_to_piDeath),earliest_hit.trackLength);

                                    //linear cut
                                    if(iftest == 0)
                                    {
                                        if(earliest_hit.trackLength-cut_slope*GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit)-cut_y_intercept > 0)
                                        {
                                            hist_bkg_1_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                                            hist_bkg_1_gamma_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);

                                            hist_bkg_1_ang_vs_dis_linear_cut->Fill(GetAngle(vec_piDeath_to_hit,vec_vtx_to_piDeath),earliest_hit.trackLength);
                                            hist_bkg_1_gamma_ang_vs_dis_linear_cut->Fill(GetAngle(vec_piDeath_to_hit,vec_vtx_to_piDeath),earliest_hit.trackLength);
                                            //energy_of_secondary->Fill(earliest_hit.CubeE);
                                            beta_of_secondary->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                        }
                                    }
                                    if(iftest == 1)
                                    {
                                        if(earliest_hit.trackLength < test_cut_distance && GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit) > test_cut_angle)
                                        {
                                            hist_bkg_1_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                                            hist_bkg_1_gamma_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);

                                            hist_bkg_1_ang_vs_dis_linear_cut->Fill(GetAngle(vec_piDeath_to_hit,vec_vtx_to_piDeath),earliest_hit.trackLength);
                                            hist_bkg_1_gamma_ang_vs_dis_linear_cut->Fill(GetAngle(vec_piDeath_to_hit,vec_vtx_to_piDeath),earliest_hit.trackLength);
                                            //energy_of_secondary->Fill(earliest_hit.CubeE);
                                            beta_of_secondary->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                        }
                                    }
                                }
                            }
                        }
                    }


                    ////proton case
                    if(num_fspi == 0 && num_fsp ==1)
                    {
                        ////distance, angle cut 2D plot
                        distance_vtx_to_deathpoint->Fill(GetDistance(earliest_hit.vtxSignal,earliest_hit.protonDeath));
                        //signal
                        if(earliest_hit.isNeutron && (earliest_hit.parentId == -1 ||earliest_hit.parentId == 0))
                        {
                            hist_sig_arm_vs_time->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                            hist_sig_ang_vs_dis->Fill(GetAngle(vec_protonDeath_to_hit,vec_vtx_to_protonDeath),earliest_hit.trackLength);
                            //linear cut
                            if(iftest == 0)
                            {
                                if(earliest_hit.trackLength-cut_slope*GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit)-cut_y_intercept < 0)
                                {
                                    hist_sig_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                                    hist_sig_ang_vs_dis_linear_cut->Fill(GetAngle(vec_protonDeath_to_hit,vec_vtx_to_protonDeath),earliest_hit.trackLength);
                                    //energy_of_signal->Fill(earliest_hit.CubeE);
                                    beta_of_signal->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                }
                            }
                            if(iftest == 1)
                            {
                                if(earliest_hit.trackLength < test_cut_distance && GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit) > test_cut_angle)
                                {
                                    hist_sig_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                                    hist_sig_ang_vs_dis_linear_cut->Fill(GetAngle(vec_protonDeath_to_hit,vec_vtx_to_protonDeath),earliest_hit.trackLength);
                                    //energy_of_signal->Fill(earliest_hit.CubeE);
                                    beta_of_signal->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                }
                            }
                        }
                        //background
                        if(earliest_hit.isGamma)
                        {
                            hist_bkg_gamma_arm_vs_time->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                            hist_bkg_1_gamma_arm_vs_time->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);

                            hist_bkg_gamma_ang_vs_dis->Fill(GetAngle(vec_protonDeath_to_hit,vec_vtx_to_protonDeath),earliest_hit.trackLength);
                            hist_bkg_1_gamma_ang_vs_dis->Fill(GetAngle(vec_protonDeath_to_hit,vec_vtx_to_protonDeath),earliest_hit.trackLength);

                            if(iftest == 0)
                            {
                                if(earliest_hit.trackLength-cut_slope*GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit)-cut_y_intercept < 0)
                                {
                                    hist_bkg_gamma_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                                    hist_bkg_1_gamma_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);

                                    hist_bkg_gamma_ang_vs_dis_linear_cut->Fill(GetAngle(vec_protonDeath_to_hit,vec_vtx_to_protonDeath),earliest_hit.trackLength);
                                    hist_bkg_1_gamma_ang_vs_dis_linear_cut->Fill(GetAngle(vec_protonDeath_to_hit,vec_vtx_to_protonDeath),earliest_hit.trackLength);
                                    //energy_of_gamma->Fill(earliest_hit.CubeE);
                                    //energy_of_secondary->Fill(earliest_hit.CubeE);
                                    beta_of_gamma->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                    beta_of_secondary->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                }
                            }
                            if(iftest == 1)
                            {
                                if(earliest_hit.trackLength < test_cut_distance && GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit) > test_cut_angle)
                                {
                                    hist_bkg_gamma_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                                    hist_bkg_1_gamma_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);

                                    hist_bkg_gamma_ang_vs_dis_linear_cut->Fill(GetAngle(vec_protonDeath_to_hit,vec_vtx_to_protonDeath),earliest_hit.trackLength);
                                    hist_bkg_1_gamma_ang_vs_dis_linear_cut->Fill(GetAngle(vec_protonDeath_to_hit,vec_vtx_to_protonDeath),earliest_hit.trackLength);
                                    //energy_of_gamma->Fill(earliest_hit.CubeE);
                                    //energy_of_secondary->Fill(earliest_hit.CubeE);
                                    beta_of_gamma->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                    beta_of_secondary->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                }
                            }
                        }
                        if(earliest_hit.parentId > 0)
                        {
                            if(abs(earliest_hit.protonDeath[0]) < 120 && abs(earliest_hit.protonDeath[1]) < 120 && earliest_hit.protonDeath[2] < 150 && earliest_hit.protonDeath[2] > -50)
                            {
                                if(earliest_hit.isFromProton && earliest_hit.isNeutron)
                                {
                                    hist_bkg_1_arm_vs_time->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                                    hist_bkg_1_gamma_arm_vs_time->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);

                                    hist_bkg_1_ang_vs_dis->Fill(GetAngle(vec_protonDeath_to_hit,vec_vtx_to_protonDeath),earliest_hit.trackLength);
                                    hist_bkg_1_gamma_ang_vs_dis->Fill(GetAngle(vec_protonDeath_to_hit,vec_vtx_to_protonDeath),earliest_hit.trackLength);

                                    //linear cut
                                    if(iftest == 0)
                                    {
                                        if(earliest_hit.trackLength-cut_slope*GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit)-cut_y_intercept < 0)
                                        {
                                            hist_bkg_1_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                                            hist_bkg_1_gamma_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);

                                            hist_bkg_1_ang_vs_dis_linear_cut->Fill(GetAngle(vec_protonDeath_to_hit,vec_vtx_to_protonDeath),earliest_hit.trackLength);
                                            hist_bkg_1_gamma_ang_vs_dis_linear_cut->Fill(GetAngle(vec_protonDeath_to_hit,vec_vtx_to_protonDeath),earliest_hit.trackLength);
                                            //energy_of_secondary->Fill(earliest_hit.CubeE);
                                            beta_of_secondary->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                        }
                                    }
                                    if(iftest == 1)
                                    {
                                        if(earliest_hit.trackLength < test_cut_distance && GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit) > test_cut_angle)
                                        {
                                            hist_bkg_1_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);
                                            hist_bkg_1_gamma_arm_vs_time_linear_cut->Fill(earliest_hit.trackLength,earliest_hit.timeWindow);

                                            hist_bkg_1_ang_vs_dis_linear_cut->Fill(GetAngle(vec_protonDeath_to_hit,vec_vtx_to_protonDeath),earliest_hit.trackLength);
                                            hist_bkg_1_gamma_ang_vs_dis_linear_cut->Fill(GetAngle(vec_protonDeath_to_hit,vec_vtx_to_protonDeath),earliest_hit.trackLength);
                                            //energy_of_secondary->Fill(earliest_hit.CubeE);
                                            beta_of_secondary->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }       //end of event iterate
    _file->Close();
    delete _file;
}

void gamma()
{
    int endPROD, beginPROD, filenum;
    
    //cout<<"filenum :"<<endl;
    //cin>>filenum;
    filenum = 300;
    cout<<"start"<<endl;
    for(int i = 2; i <filenum; i++) //test_1 is not
    {
        cout<<"\033[1APROD"<<101<<": "<<(double)(i*100/filenum)<<"%\033[1000D"<<endl;
        analyze(Form("/Users/gwon/Geo12/PROD101/RHC_%d_wGamma_2ndVersion.root",i));
        //analyze(Form("/pnfs/dune/persistent/users/gyang/3DST/dump/standardGeo12/PROD101/RHC_%d_wGamma_2ndVersion.root",i));
    }

    cout<<endl;
    cout<<"end"<<endl;

    gErrorIgnoreLevel = kWarning;
    TString folder_name = TString::Format("cc%dpi%dp_slope_%d_yintercept_%d",num_fspi,num_fsp,cut_slope,cut_y_intercept);
    if(num_fspi == 1)
    {
        gSystem->mkdir("pion");
        gSystem->cd("pion");
    }
    if(num_fsp == 1)
    {
        gSystem->mkdir("proton");
        gSystem->cd("proton");
    }
    gSystem->mkdir(folder_name);
    gSystem->cd(folder_name);

    double efficiency = hist_sig_ang_vs_dis_linear_cut->GetEntries()/hist_sig_ang_vs_dis->GetEntries();
    double purity = hist_sig_ang_vs_dis_linear_cut->GetEntries()/(hist_sig_ang_vs_dis_linear_cut->GetEntries()+hist_bkg_1_gamma_ang_vs_dis_linear_cut->GetEntries());
    cout<<"purity: "<<purity<<endl;
    cout<<"efficiency :"<<efficiency<<endl;

    cout<<"num_primary_gamma: "<<num_primary_gamma<<endl;
    cout<<"num_secondary_gamma: "<<num_secondary_gamma<<endl;
    cout<<"num_gamma_from_muon: "<<num_gamma_from_muon<<endl;
    cout<<"num_gamma_from_pion: "<<num_gamma_from_pion<<endl;
    
    TFile a(TString::Format("purity_%f",purity),"RECREATE");
    a.Close();
    TFile b(TString::Format("efficiency_%f",efficiency),"RECREATE");
    b.Close();
    TFile c(TString::Format("p*e_%f",purity*efficiency),"RECREATE");
    c.Close();


    TFile * fi1 = new TFile("background.root","RECREATE");

    TCanvas * can = new TCanvas;
    can->Divide(2,2);
    can->cd(1);
    hist_sig_arm_vs_time->Draw("colz");
    can->cd(2);
    hist_bkg_out3DST->Draw("colz");
    can->cd(3);
    hist_bkg_NC->Draw("colz");
    can->cd(4);
    hist_bkg_1_arm_vs_time->Draw("colz");
    can->SaveAs("4plots.pdf");
    can->Clear();


    hitPDG_signal_small->Draw();
    hitPDG_signal_small->SetStats(0);
    hitPDG_signal_small->Scale(1/hitPDG_signal_small->GetEntries(),"nosw2");
    hitPDG_signal_small->Write();
    can->SaveAs("hitPDG_signal_small.pdf");
    can->Clear();

    cout<<"debug2"<<endl;
    hitPDG_signal_big->Draw();
    hitPDG_signal_big->SetStats(0);
    hitPDG_signal_big->Scale(1/hitPDG_signal_big->GetEntries(),"nosw2");
    hitPDG_signal_big->Write();
    can->SaveAs("hitPDG_signal_big.pdf");
    can->Clear();

    gammaPDG->Draw();
    gammaPDG->Write();
    can->SaveAs("gammaPDG.pdf");
    can->Clear();

    num_of_primary_gamma->Draw();
    num_of_primary_gamma->Write();
    can->SaveAs("num_of_primary_gamma.pdf");
    can->Clear();

    num_of_other_gamma->Draw();
    num_of_other_gamma->Write();
    can->SaveAs("num_of_other_gamma.pdf");
    can->Clear();

    //beta
    //{
    beta_of_signal->Draw();
    beta_of_signal->Write();
    can->SaveAs("beta_of_signal.pdf");
    can->Clear();

    beta_of_gamma->Draw();
    beta_of_gamma->Write();
    can->SaveAs("beta_of_gamma.pdf");
    can->Clear();

    beta_of_secondary->Draw();
    beta_of_secondary->Write();
    can->SaveAs("beta_of_secondary.pdf");
    can->Clear();

    beta_of_signal->SetStats(0);
    beta_of_gamma->SetStats(0);
    beta_of_secondary->SetStats(0);
    beta_of_signal->Scale(1/beta_of_signal->GetEntries(),"nosw2");
    beta_of_gamma->Scale(1/beta_of_gamma->GetEntries(),"nosw2");
    beta_of_secondary->Scale(1/beta_of_secondary->GetEntries(),"nosw2");
    beta_of_signal->SetLineColor(6);        //purple
    beta_of_gamma->SetLineColor(8);     //green
    beta_of_secondary->SetLineColor(4);     //blue
    beta_of_signal->Draw();
    beta_of_gamma->Draw("same");
    beta_of_secondary->Draw("same");
    can->SaveAs("all_beta.pdf");
    can->Clear();
    //}
    
    //1D energy
    //{
    energy_of_gamma->Draw();
    energy_of_gamma->Write();
    can->SaveAs("energy_of_gamma.pdf");
    can->Clear();

    energy_of_signal->Draw();
    energy_of_signal->Write();
    can->SaveAs("energy_of_signal.pdf");
    can->Clear();

    energy_of_secondary->Draw();
    energy_of_secondary->Write();
    can->SaveAs("energy_of_secondary.pdf");
    can->Clear();

    energy_of_gamma->Scale(1/energy_of_gamma->GetEntries(),"nosw2");
    energy_of_gamma->SetStats(0);
    energy_of_gamma->SetLineColor(8);
    energy_of_gamma->Draw();
    energy_of_secondary->Scale(1/energy_of_secondary->GetEntries(),"nosw2");
    energy_of_secondary->SetStats(0);
    energy_of_secondary->SetLineColor(4);
    energy_of_secondary->Draw("same");
    energy_of_signal->Scale(1/energy_of_signal->GetEntries(),"nosw2");
    energy_of_signal->SetStats(0);
    energy_of_signal->SetLineColor(6);
    energy_of_signal->Draw("same");
    can->SaveAs("all_energy.pdf");
    can->Clear();
    //}


    distance_vtx_to_deathpoint->Draw();
    distance_vtx_to_deathpoint->Write();
    can->SaveAs("distance_vtx_to_deathpoint.pdf");
    can->Clear();

    //lever arm vs time 
    hist_sig_arm_vs_time->Draw("colz");
    hist_sig_arm_vs_time->Write();
    can->SaveAs("hist_sig_arm_vs_time.pdf");
    can->Clear();

    hist_bkg_1_arm_vs_time->Draw("colz");
    hist_bkg_1_arm_vs_time->Write();
    can->SaveAs("hist_bkg_1_arm_vs_time.pdf");
    can->Clear();

    hist_bkg_gamma_arm_vs_time->Draw("colz");
    hist_bkg_gamma_arm_vs_time->Write();
    can->SaveAs("hist_bkg_gamma_arm_vs_time.pdf");
    can->Clear();

    hist_bkg_1_gamma_arm_vs_time->Draw("colz");
    hist_bkg_1_gamma_arm_vs_time->Write();
    can->SaveAs("hist_bkg_1_gamma_arm_vs_time.pdf");
    can->Clear();
    
    //lever arm vs time with cut
    hist_sig_arm_vs_time_linear_cut->Draw("colz");
    hist_sig_arm_vs_time_linear_cut->Write();
    can->SaveAs("hist_sig_arm_vs_time_linear_cut.pdf");
    can->Clear();

    hist_bkg_1_arm_vs_time_linear_cut->Draw("colz");
    hist_bkg_1_arm_vs_time_linear_cut->Write();
    can->SaveAs("hist_bkg_1_arm_vs_time_linear_cut.pdf");
    can->Clear();

    hist_bkg_gamma_arm_vs_time_linear_cut->Draw("colz");
    hist_bkg_gamma_arm_vs_time_linear_cut->Write();
    can->SaveAs("hist_bkg_gamma_arm_vs_time_linear_cut.pdf");
    can->Clear();

    hist_bkg_1_gamma_arm_vs_time_linear_cut->Draw("colz");
    hist_bkg_1_gamma_arm_vs_time_linear_cut->Write();
    can->SaveAs("hist_bkg_1_gamma_arm_vs_time_linear_cut.pdf");
    can->Clear();
    
    //angla vs lever arm 
    hist_sig_ang_vs_dis->Draw("colz");
    hist_sig_ang_vs_dis->Write();
    can->SaveAs("hist_sig_ang_vs_dis.pdf");
    can->Clear();

    hist_bkg_1_ang_vs_dis->Draw("colz");
    hist_bkg_1_ang_vs_dis->Write();
    can->SaveAs("hist_bkg_1_ang_vs_dis.pdf");
    can->Clear();

    hist_bkg_gamma_ang_vs_dis->Draw("colz");
    hist_bkg_gamma_ang_vs_dis->Write();
    can->SaveAs("hist_bkg_gamma_ang_vs_dis.pdf");
    can->Clear();

    hist_bkg_1_gamma_ang_vs_dis->Draw("colz");
    hist_bkg_1_gamma_ang_vs_dis->Write();
    can->SaveAs("hist_bkg_1_gamma_ang_vs_dis.pdf");
    can->Clear();
    
    //angle vs lever arm with cut
    hist_sig_ang_vs_dis_linear_cut->Draw("colz");
    hist_sig_ang_vs_dis_linear_cut->Write();
    can->SaveAs("hist_sig_ang_vs_dis_linear_cut.pdf");
    can->Clear();

    hist_bkg_gamma_ang_vs_dis_linear_cut->Draw("colz");
    hist_bkg_gamma_ang_vs_dis_linear_cut->Write();
    can->SaveAs("hist_bkg_gamma_ang_vs_dis_linear_cut.pdf");
    can->Clear();

    hist_bkg_1_ang_vs_dis_linear_cut->Draw("colz");
    hist_bkg_1_ang_vs_dis_linear_cut->Write();
    can->SaveAs("hist_bkg_1_ang_vs_dis_linear_cut.pdf");
    can->Clear();

    hist_bkg_1_gamma_ang_vs_dis_linear_cut->Draw("colz");
    hist_bkg_1_gamma_ang_vs_dis_linear_cut->Write();
    can->SaveAs("hist_bkg_1_gamma_ang_vs_dis_linear_cut.pdf");
    can->Clear();
    
    //gamma
    //{
    TH2F * purity_distribution_gamma_ang_vs_dis_linear_cut = (TH2F*)hist_sig_ang_vs_dis_linear_cut->Clone();
    hist_bkg_gamma_ang_vs_dis_linear_cut->Add(hist_sig_ang_vs_dis_linear_cut);
    purity_distribution_gamma_ang_vs_dis_linear_cut->SetStats(0);
    purity_distribution_gamma_ang_vs_dis_linear_cut->SetMinimum(0);
    purity_distribution_gamma_ang_vs_dis_linear_cut->Divide(hist_bkg_gamma_ang_vs_dis_linear_cut);
    purity_distribution_gamma_ang_vs_dis_linear_cut->SetTitle("purity");
    purity_distribution_gamma_ang_vs_dis_linear_cut->Draw("colz");
    can->SaveAs("purity_distribution_gamma_ang_vs_dis_linear_cut.pdf");
    can->Clear();

    TH2F * purity_distribution_gamma_ang_vs_dis = (TH2F*)hist_sig_ang_vs_dis->Clone();
    hist_bkg_gamma_ang_vs_dis->Add(hist_sig_ang_vs_dis);
    purity_distribution_gamma_ang_vs_dis->SetStats(0);
    purity_distribution_gamma_ang_vs_dis->SetMinimum(0);
    purity_distribution_gamma_ang_vs_dis->Divide(hist_bkg_gamma_ang_vs_dis);
    purity_distribution_gamma_ang_vs_dis->SetTitle("purity(gamma)");
    purity_distribution_gamma_ang_vs_dis->Draw("colz");
    can->SaveAs("purity_distribution_gamma_ang_vs_dis.pdf");
    can->Clear();

    TH2F * purity_gamma_linear_cut = (TH2F*)hist_sig_arm_vs_time_linear_cut->Clone();
    hist_bkg_gamma_arm_vs_time_linear_cut->Add(hist_sig_arm_vs_time_linear_cut);
    purity_gamma_linear_cut->Divide(hist_bkg_gamma_arm_vs_time_linear_cut);
    purity_gamma_linear_cut->SetStats(0);
    purity_gamma_linear_cut->SetTitle("purity with cut");
    purity_gamma_linear_cut->SetMaximum(1);
    purity_gamma_linear_cut->Draw("colz");
    purity_gamma_linear_cut->Write();
    can->SaveAs("purity_gamma_linear_cut.pdf");
    can->Clear();

    TH2F * purity_gamma = (TH2F*)hist_sig_arm_vs_time->Clone();
    hist_bkg_gamma_arm_vs_time->Add(hist_sig_arm_vs_time);
    purity_gamma->Divide(hist_bkg_gamma_arm_vs_time);
    purity_gamma->SetStats(0);
    purity_gamma->SetTitle("purity(gamma)");
    purity_gamma->SetMaximum(1);
    purity_gamma->Draw("colz");
    purity_gamma->Write();
    can->SaveAs("purity_gamma.pdf");
    can->Clear();
    //}

    //secondary
    //{
    TH2F * purity_distribution_1_ang_vs_dis_linear_cut = (TH2F*)hist_sig_ang_vs_dis_linear_cut->Clone();
    hist_bkg_1_ang_vs_dis_linear_cut->Add(hist_sig_ang_vs_dis_linear_cut);
    purity_distribution_1_ang_vs_dis_linear_cut->SetStats(0);
    purity_distribution_1_ang_vs_dis_linear_cut->SetMinimum(0);
    purity_distribution_1_ang_vs_dis_linear_cut->Divide(hist_bkg_1_ang_vs_dis_linear_cut);
    purity_distribution_1_ang_vs_dis_linear_cut->SetTitle("purity");
    purity_distribution_1_ang_vs_dis_linear_cut->Draw("colz");
    can->SaveAs("purity_distribution_1_ang_vs_dis_linear_cut.pdf");
    can->Clear();

    TH2F * purity_distribution_1_ang_vs_dis = (TH2F*)hist_sig_ang_vs_dis->Clone();
    hist_bkg_1_ang_vs_dis->Add(hist_sig_ang_vs_dis);
    purity_distribution_1_ang_vs_dis->SetStats(0);
    purity_distribution_1_ang_vs_dis->SetMinimum(0);
    purity_distribution_1_ang_vs_dis->Divide(hist_bkg_1_ang_vs_dis);
    purity_distribution_1_ang_vs_dis->SetTitle("purity(secondary)");
    purity_distribution_1_ang_vs_dis->Draw("colz");
    can->SaveAs("purity_distribution_1_ang_vs_dis.pdf");
    can->Clear();

    TH2F * purity_1_linear_cut = (TH2F*)hist_sig_arm_vs_time_linear_cut->Clone();
    hist_bkg_1_arm_vs_time_linear_cut->Add(hist_sig_arm_vs_time_linear_cut);
    purity_1_linear_cut->Divide(hist_bkg_1_arm_vs_time_linear_cut);
    purity_1_linear_cut->SetStats(0);
    purity_1_linear_cut->SetTitle("purity with cut");
    purity_1_linear_cut->SetMaximum(1);
    purity_1_linear_cut->Draw("colz");
    purity_1_linear_cut->Write();
    can->SaveAs("purity_1_linear_cut.pdf");
    can->Clear();

    TH2F * purity_1 = (TH2F*)hist_sig_arm_vs_time->Clone();
    hist_bkg_1_arm_vs_time->Add(hist_sig_arm_vs_time);
    purity_1->Divide(hist_bkg_1_arm_vs_time);
    purity_1->SetStats(0);
    purity_1->SetTitle("purity(secondary)");
    purity_1->SetMaximum(1);
    purity_1->Draw("colz");
    purity_1->Write();
    can->SaveAs("purity_1.pdf");
    can->Clear();
    //}

    //secondary+gamma
    //{
    TH2F * purity_distribution_1_gamma_ang_vs_dis_linear_cut = (TH2F*)hist_sig_ang_vs_dis_linear_cut->Clone();
    hist_bkg_1_gamma_ang_vs_dis_linear_cut->Add(hist_sig_ang_vs_dis_linear_cut);
    purity_distribution_1_gamma_ang_vs_dis_linear_cut->SetStats(0);
    purity_distribution_1_gamma_ang_vs_dis_linear_cut->SetMinimum(0);
    purity_distribution_1_gamma_ang_vs_dis_linear_cut->Divide(hist_bkg_1_gamma_ang_vs_dis_linear_cut);
    purity_distribution_1_gamma_ang_vs_dis_linear_cut->SetTitle("purity");
    purity_distribution_1_gamma_ang_vs_dis_linear_cut->Draw("colz");
    can->SaveAs("purity_distribution_1_gamma_ang_vs_dis_linear_cut.pdf");
    can->Clear();

    TH2F * purity_distribution_1_gamma_ang_vs_dis = (TH2F*)hist_sig_ang_vs_dis->Clone();
    hist_bkg_1_gamma_ang_vs_dis->Add(hist_sig_ang_vs_dis);
    purity_distribution_1_gamma_ang_vs_dis->SetStats(0);
    purity_distribution_1_gamma_ang_vs_dis->SetMinimum(0);
    purity_distribution_1_gamma_ang_vs_dis->Divide(hist_bkg_1_gamma_ang_vs_dis);
    purity_distribution_1_gamma_ang_vs_dis->SetTitle("purity(secondary+gamma)");
    purity_distribution_1_gamma_ang_vs_dis->Draw("colz");
    can->SaveAs("purity_distribution_1_gamma_ang_vs_dis.pdf");
    can->Clear();

    TH2F * purity_1_gamma_linear_cut = (TH2F*)hist_sig_arm_vs_time_linear_cut->Clone();
    hist_bkg_1_gamma_arm_vs_time_linear_cut->Add(hist_sig_arm_vs_time_linear_cut);
    purity_1_gamma_linear_cut->Divide(hist_bkg_1_gamma_arm_vs_time_linear_cut);
    purity_1_gamma_linear_cut->SetStats(0);
    purity_1_gamma_linear_cut->SetTitle("purity with cut");
    purity_1_gamma_linear_cut->SetMaximum(1);
    purity_1_gamma_linear_cut->Draw("colz");
    purity_1_gamma_linear_cut->Write();
    can->SaveAs("purity_1_gamma_linear_cut.pdf");
    can->Clear();

    cout<<"debug3"<<endl;
    TH2F * purity_1_gamma = (TH2F*)hist_sig_arm_vs_time->Clone();
    hist_bkg_1_gamma_arm_vs_time->Add(hist_sig_arm_vs_time);
    purity_1_gamma->Divide(hist_bkg_1_gamma_arm_vs_time);
    purity_1_gamma->SetStats(0);
    purity_1_gamma->SetTitle("purity(secondary+gamma)");
    purity_1_gamma->SetMaximum(1);
    purity_1_gamma->Draw("colz");
    purity_1_gamma->Write();
    can->SaveAs("purity_1_gamma.pdf");
    can->Clear();
    //}
    
    //gamma 
    can->Divide(2,1);
    can->SetCanvasSize(700,300);
    can->cd(1);
    purity_gamma->Draw("colz");
    can->cd(2);
    purity_gamma_linear_cut->Draw("colz");
    can->SaveAs("2purity(gamma).pdf");
    can->Clear();
    
    //secondary 
    can->Divide(2,1);
    can->SetCanvasSize(700,300);
    can->cd(1);
    purity_1->Draw("colz");
    can->cd(2);
    purity_1_linear_cut->Draw("colz");
    can->SaveAs("2purity(secondary).pdf");
    can->Clear();

    //secondary+gamma
    can->Divide(2,1);
    can->SetCanvasSize(700,300);
    can->cd(1);
    purity_1_gamma->Draw("colz");
    can->cd(2);
    purity_1_gamma_linear_cut->Draw("colz");
    can->SaveAs("2purity(secondary+gamma).pdf");
    can->Clear();

    fi1->Close();


    delete can;
    delete fi1;
    delete energy_of_signal; 
    delete energy_of_gamma; 
    delete energy_of_secondary; 
    delete gammaPDG; 
    delete num_of_primary_gamma; 
    delete num_of_other_gamma; 
    delete hitPDG_signal_big; 
    delete hitPDG_signal_small; 
    delete beta_of_signal; 
    delete beta_of_gamma; 
    delete beta_of_secondary; 
    delete hist_bkg_out3DST; 
    delete hist_bkg_NC; 
    delete hist_sig_arm_vs_time; 
    delete hist_bkg_1_arm_vs_time; 
    delete hist_bkg_gamma_arm_vs_time; 
    delete hist_bkg_1_gamma_arm_vs_time; 
    delete hist_sig_arm_vs_time_linear_cut; 
    delete hist_bkg_1_arm_vs_time_linear_cut; 
    delete hist_bkg_gamma_arm_vs_time_linear_cut; 
    delete hist_bkg_1_gamma_arm_vs_time_linear_cut; 
    delete hist_sig_ang_vs_dis; 
    delete hist_bkg_1_ang_vs_dis; 
    delete hist_bkg_gamma_ang_vs_dis; 
    delete hist_bkg_1_gamma_ang_vs_dis; 
    delete hist_sig_ang_vs_dis_linear_cut; 
    delete hist_bkg_gamma_ang_vs_dis_linear_cut; 
    delete hist_bkg_1_ang_vs_dis_linear_cut; 
    delete hist_bkg_1_gamma_ang_vs_dis_linear_cut; 
    delete distance_vtx_to_deathpoint; 
}
