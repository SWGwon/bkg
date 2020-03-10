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

float leverArm;
float angle;
float beta;
float distanceCHit;
float tof;
float cubeE;
float category;

float signalLeverArm;
float signalAngle;
float signalBeta;
float signalDistanceCHit;
float signalTOF;
float signalCubeE;
float secondaryNeutronLeverArm;
float secondaryNeutronAngle;
float secondaryNeutronBeta;
float secondaryNeutronDistanceCHit;
float secondaryNeutronTOF;
float secondaryNeutronCubeE;
float primaryGammaLeverArm;
float primaryGammaAngle;
float primaryGammaBeta;
float primaryGammaDistanceCHit;
float primaryGammaTOF;
float primaryGammaCubeE;
float secondaryGammaLeverArm;
float secondaryGammaAngle;
float secondaryGammaBeta;
float secondaryGammaDistanceCHit;
float secondaryGammaTOF;
float secondaryGammaCubeE;

int num_secondary_neutron = 0;
int num_secondary_neutron_C_out3dst  = 0;

using namespace std;
//histograms{

TH2F * beta_vs_leverarm_parentId_0 = new TH2F("beta_vs_leverarm_parentId_0", "",30,0,1.5,20,0,200);
TH2F * beta_vs_leverarm_parentId_1 = new TH2F("beta_vs_leverarm_parentId_1", "",30,0,1.5,20,0,200);

//signal{
TH1F * leverarm_signal = new TH1F("leverarm_signal","lever arm of signal",20,0,200);
TH1F * angle_signal = new TH1F("angle_signal","angle between C and signal hit",20,0,1);
TH1F * beta_signal = new TH1F("beta_signal","beta of signal",30,0,1.5);
TH1F * distance_signal = new TH1F("distance_signal","distance b/w C and signal hit",20,0,200);
TH1F * TOF_signal = new TH1F("TOF_signal","time of flight of signal", 25,0,25);
TH1F * CubeE_signal = new TH1F("CubeE_signal", "CubeE of signal", 30, 0, 15);
//}

//secondary neutron{
TH1F * leverarm_secondary_neutron = new TH1F("leverarm_secondary_neutron","lever arm of secondary_neutron",20,0,200);
TH1F * angle_secondary_neutron = new TH1F("angle_secondary_neutron","angle between C and secondary_neutron hit",20,0,1);
TH1F * beta_secondary_neutron = new TH1F("beta_secondary_neutron","beta of secondary_neutron",30,0,1.5);
TH1F * distance_secondary_neutron = new TH1F("distance_secondary_neutron","distance b/w C and secondary_neutron hit",20,0,200);
TH1F * TOF_secondary_neutron = new TH1F("TOF_secondary_neutron","time of flight of secondary_neutron", 25,0,25);
TH1F * CubeE_secondary_neutron = new TH1F("CubeE_secondary_neutron", "CubeE of secondary_neutron", 30, 0, 15);
//}

//primary gamma{
TH1F * leverarm_primary_gamma = new TH1F("leverarm_primary_gamma","lever arm of primary_gamma",20,0,200);
TH1F * angle_primary_gamma = new TH1F("angle_primary_gamma","angle between C and primary_gamma hit",20,0,1);
TH1F * beta_primary_gamma = new TH1F("beta_primary_gamma","beta of primary_gamma",30,0,1.5);
TH1F * distance_primary_gamma = new TH1F("distance_primary_gamma","distance b/w C and primary_gamma hit",20,0,200);
TH1F * TOF_primary_gamma = new TH1F("TOF_primary_gamma","time of flight of primary_gamma", 25,0,25);
TH1F * CubeE_primary_gamma = new TH1F("CubeE_primary_gamma", "CubeE of primary_gamma", 30, 0, 15);
//}

//secondary gamma{
TH1F * leverarm_secondary_gamma = new TH1F("leverarm_secondary_gamma","lever arm of secondary_gamma",20,0,200);
TH1F * angle_secondary_gamma = new TH1F("angle_secondary_gamma","angle between C and secondary_gamma hit",20,0,1);
TH1F * beta_secondary_gamma = new TH1F("beta_secondary_gamma","beta of secondary_gamma",30,0,1.5);
TH1F * distance_secondary_gamma = new TH1F("distance_secondary_gamma","distance b/w C and secondary_gamma hit",20,0,200);
TH1F * TOF_secondary_gamma = new TH1F("TOF_secondary_gamma","time of flight of secondary_gamma", 25,0,25);
TH1F * CubeE_secondary_gamma = new TH1F("CubeE_secondary_gamma", "CubeE of secondary_gamma", 30, 0, 15);
//}

//}

bool is_inFV = false;       //check if vertex is in FV
bool is_in3DST = false;     //check if vertex is in 3DST

float energyHitCut = 0; //energy deposit threshold for cube
bool isCubeE = 1;

int number_of_CC = 0;

//change this part to do slope, intercept test
int cut_slope = 0;
int cut_y_intercept = 100000;
//channel type
int num_fspi = 1;   //number of fs charged pion
int num_fsp = 0;    //number of fs proton


int iftest = 1;
double test_cut_angle = 0;
double test_cut_distance = 1000;

const double c_velocity = 29.9792458;

class Hit_t 
{
    public:
        float timeWindow,           // time windows of the hit
              timeStartT,       //time between hit and startT
              lengthStart,       //length between hit and start position
              timeWindowSmear,      // smeared time window 
              timeSmear,        // smear time
              energyDeposit,        // energy deposited by the neutron
              trackLength,          // lever arm
              trueRec,      // true reconstructed energy
              smearRec,
              gammaTrueE,
              vtxSignal[3],     // neutrino vertex position of the neutron
              vtxTime,      // neutrino  vertex time
              trueE,    //neutron true energy
              CubeE,    //neutron cube energy
              trueT,    //neutron true time
              hitPDG,    //hit PDG
              startingPointT,
              category,     //kind of hit(1: signal, 2: secondary neutron, 3: primary gamma, 4: secondary gamma)
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
            timeStartT(0),
            lengthStart(0),
            energyDeposit(0),
            trackLength(0),  
            trueRec(0),      
            smearRec(0),
            vtxTime(0),      
            trueE(0), 
            CubeE(0),
            trueT(0), 
            hitPDG(0),
            startingPointT(0),
            gammaTrueE(0),
            bkgLoc(125124123),        
            parentId(123124123),
            parentPdg(123123123),
            category(-1),
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

void gamma()
{
    int endPROD, beginPROD, filenum;
    
    //cout<<"filenum :"<<endl;
    //cin>>filenum;
    filenum = 1000;
    cout<<"start"<<endl;

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

    TFile * outfile = new TFile("variables.root","RECREATE");
    TTree * output_tree = new TTree("output_tree","output_tree");

    /*
    output_tree->Branch("leverArm",&leverArm, "lever arm");
    output_tree->Branch("angle",&angle, "angle between C and hit");
    output_tree->Branch("beta",&beta, "beta");
    output_tree->Branch("distanceCHit",&distanceCHit, "distance C and hit");
    output_tree->Branch("tof",&tof, "time of flight");
    output_tree->Branch("cubeE",&cubeE, "CubeE");
    */
    output_tree->Branch("category", &category, "category");
    output_tree->Branch("signalLeverArm", &signalLeverArm, "signalLeverArm");
    output_tree->Branch("signalAngle", &signalAngle, "signalAngle");
    output_tree->Branch("signalBeta", &signalBeta, "signalBeta");
    output_tree->Branch("signalDistanceCHit", &signalDistanceCHit, "signalDistanceCHit");
    output_tree->Branch("signalTOF", &signalTOF, "signalTOF");
    output_tree->Branch("signalCubeE", &signalCubeE, "signalCubeE");
    output_tree->Branch("secondaryNeutronLeverArm", &secondaryNeutronLeverArm, "secondaryNeutronLeverArm");
    output_tree->Branch("secondaryNeutronAngle", &secondaryNeutronAngle, "secondaryNeutronAngle");
    output_tree->Branch("secondaryNeutronBeta", &secondaryNeutronBeta, "secondaryNeutronBeta");
    output_tree->Branch("secondaryNeutronDistanceCHit", &secondaryNeutronDistanceCHit, "secondaryNeutronDistanceCHit");
    output_tree->Branch("secondaryNeutronTOF", &secondaryNeutronTOF, "secondaryNeutronTOF");
    output_tree->Branch("secondaryNeutronCubeE", &secondaryNeutronCubeE, "secondaryNeutronCubeE");
    output_tree->Branch("primaryGammaLeverArm", &primaryGammaLeverArm, "primaryGammaLeverArm");
    output_tree->Branch("primaryGammaAngle", &primaryGammaAngle, "primaryGammaAngle");
    output_tree->Branch("primaryGammaBeta", &primaryGammaBeta, "primaryGammaBeta");
    output_tree->Branch("primaryGammaDistanceCHit", &primaryGammaDistanceCHit, "primaryGammaDistanceCHit");
    output_tree->Branch("primaryGammaTOF", &primaryGammaTOF, "primaryGammaTOF");
    output_tree->Branch("primaryGammaCubeE", &primaryGammaCubeE, "primaryGammaCubeE");
    output_tree->Branch("secondaryGammaLeverArm", &secondaryGammaLeverArm, "secondaryGammaLeverArm");
    output_tree->Branch("secondaryGammaAngle", &secondaryGammaAngle, "secondaryGammaAngle");
    output_tree->Branch("secondaryGammaBeta", &secondaryGammaBeta, "secondaryGammaBeta");
    output_tree->Branch("secondaryGammaDistanceCHit", &secondaryGammaDistanceCHit, "secondaryGammaDistanceCHit");
    output_tree->Branch("secondaryGammaTOF", &secondaryGammaTOF, "secondaryGammaTOF");
    output_tree->Branch("secondaryGammaCubeE", &secondaryGammaCubeE, "secondaryGammaCubeE");
                                                                                
    for(int i = 1; i != filenum; i++)
    {
        //cout<<"\033[1APROD"<<101<<": "<<(double)(i*100/filenum)<<"%\033[1000D"<<endl;
        cout<<i<<"th file"<<endl;
        string file = Form("/Users/gwon/Geo12/PROD101/RHC_%d_wGamma_2ndVersion.root",i);
        //analyze(Form("/Users/gwon/Geo12/PROD101/RHC_%d_wGamma_2ndVersion.root",i));
        //string file = Form("/pnfs/dune/persistent/users/gyang/3DST/dump/standardGeo12/PROD101/RHC_%d_wGamma_2ndVersion.root",i);

        auto _file = new TFile(TString(file));
        auto tree = (TTree*)_file->Get("tree");

        if(tree == NULL)
        {
            _file->Close();
            continue;
        }

        float t_neutronHitX[1000], t_neutronHitY[1000], t_neutronHitZ[1000];
        float t_neutronStartingPointX[1000], t_neutronStartingPointY[1000], t_neutronStartingPointZ[1000], t_neutronStartingPointT[1000];
        float t_neutronHitT[1000], t_neutronParentId[1000], t_neutronParentPDG[1000];
        float t_neutronHitE[1000], t_neutronTrueE[1000];
        float t_neutronCubeE[1000];
        float t_neutronHitSmearT[999];
        float t_neutronHitPDG[1000];

        float t_vtx[3], t_vtxTime;
        float t_piDeath[3], t_protonDeath[3];

        float vec_piDeath_to_hit[3];
        float vec_protonDeath_to_hit[3];
        float vec_vtx_to_piDeath[3];
        float vec_vtx_to_protonDeath[3];

        float t_gammaHitX[1000], t_gammaHitY[1000], t_gammaHitZ[1000];
        float t_gammaStartingPointX[1000], t_gammaStartingPointY[1000], t_gammaStartingPointZ[1000], t_gammaStartingPointT[1000];
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
        tree->SetBranchAddress("neutronStartingPointT", &t_neutronStartingPointT);
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
        tree->SetBranchAddress("gammaStartingPointT", &t_gammaStartingPointT);
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
            leverArm = -10;
            angle = -1;
            beta = -1;
            distanceCHit = -10;
            tof = -1;
            cubeE = -10;
            category = -1;

            signalLeverArm = -1000;
            signalAngle = -1000;
            signalBeta = -1000;
            signalDistanceCHit = -1000;
            signalTOF = -1000;
            signalCubeE = -1000;
            secondaryNeutronLeverArm = -1000;
            secondaryNeutronAngle = -1000;
            secondaryNeutronBeta = -1000;
            secondaryNeutronDistanceCHit = -1000;
            secondaryNeutronTOF = -1000;
            secondaryNeutronCubeE = -1000;
            primaryGammaLeverArm = -1000;
            primaryGammaAngle = -1000;
            primaryGammaBeta = -1000;
            primaryGammaDistanceCHit = -1000;
            primaryGammaTOF = -1000;
            primaryGammaCubeE = -1000;
            secondaryGammaLeverArm = -1000;
            secondaryGammaAngle = -1000;
            secondaryGammaBeta = -1000;
            secondaryGammaDistanceCHit = -1000;
            secondaryGammaTOF = -1000;
            secondaryGammaCubeE = -1000;

            Hit_t earliest_hit;
            Hit_t earliest_neutron_hit;
            Hit_t earliest_gamma_hit;

            int num_pi = 0;
            int num_proton = 0;

            tree->GetEntry(event);

            //out of fiducial volume
            if(abs(t_vtx[0]) > 50 || abs(t_vtx[1]) > 50 || abs(t_vtx[2]) > 50)
                continue;

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
                //out of 3dst
                if(abs(t_neutronHitX[n_neutronHit]) > 120 || 
                        abs(t_neutronHitY[n_neutronHit]) > 120 || 
                        abs(t_neutronHitZ[n_neutronHit]) > 100)
                    continue;

                //starting point == -1 is default value so remove it
                if(t_neutronStartingPointX[n_neutronHit] == -1 || 
                        t_neutronStartingPointY[n_neutronHit] == -1 || 
                        t_neutronStartingPointZ[n_neutronHit] == -1 )
                    continue;

                //energy threshold
                if(t_neutronCubeE[n_neutronHit] < energyHitCut || t_neutronCubeE[n_neutronHit] == 0)
                    continue;

                if(t_neutronHitT[n_neutronHit] != 0 && t_neutronHitX[n_neutronHit] != 0 && t_neutronHitT[n_neutronHit] < temp_earliest_time)
                {
                    temp_earliest_time = t_neutronHitT[n_neutronHit];
                    //calculate lever arm
                    float trackLength = pow(
                            pow(t_neutronHitX[n_neutronHit] - t_vtx[0],2)+
                            pow(t_neutronHitY[n_neutronHit] - t_vtx[1],2)+
                            pow(t_neutronHitZ[n_neutronHit] - t_vtx[2],2),0.5);
                    //length from starting point
                    float trackLengthStart = pow(
                            pow(t_neutronHitX[n_neutronHit] - t_neutronStartingPointX[n_neutronHit],2)+
                            pow(t_neutronHitY[n_neutronHit] - t_neutronStartingPointY[n_neutronHit],2)+
                            pow(t_neutronHitZ[n_neutronHit] - t_neutronStartingPointZ[n_neutronHit],2),0.5);

                    //calculate signal window; time of flight
                    float signalWindow = t_neutronHitT[n_neutronHit] - t_vtxTime;
                    //time b/w starting point and hit
                    float signalWindowStart = t_neutronHitT[n_neutronHit] - t_vtxTime - t_neutronStartingPointT[n_neutronHit];
                    float signalWindowSmear = t_neutronHitSmearT[n_neutronHit] - t_vtxTime;

                    if(signalWindow > 0)
                    {
                        earliest_neutron_hit.timeWindow = signalWindow;
                        earliest_neutron_hit.timeStartT = signalWindowStart;
                        earliest_neutron_hit.lengthStart = trackLengthStart;
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

            float temp_earliest_time_for_gamma = 1000000;
            for(int n_gammaHit = 0; n_gammaHit < 1000; n_gammaHit++)
            {
                //out of 3dst
                if(abs(t_gammaHitX[n_gammaHit]) > 120 || 
                        abs(t_gammaHitY[n_gammaHit]) > 120 || 
                        abs(t_gammaHitZ[n_gammaHit]) > 100)
                    continue;

                //starting point == -1 is default value so remove it
                if(t_gammaStartingPointX[n_gammaHit] == -1 || 
                        t_gammaStartingPointY[n_gammaHit] == -1 || 
                        t_gammaStartingPointZ[n_gammaHit] == -1 )
                    continue;

                //energy threshold
                if(t_gammaCubeE[n_gammaHit] < energyHitCut || t_gammaCubeE[n_gammaHit] == 0)
                    continue;

                //remove parentid == 0
                if(t_gammaParentId[n_gammaHit] == 0)
                    continue;

                if(t_gammaHitT[n_gammaHit] != 0 && t_gammaHitX[n_gammaHit] != 0 && t_gammaHitT[n_gammaHit] < temp_earliest_time_for_gamma)
                {
                    temp_earliest_time_for_gamma = t_gammaHitT[n_gammaHit];

                    //calculate lever arm
                    float trackLength = pow(
                            pow(t_gammaHitX[n_gammaHit] - t_vtx[0],2)+
                            pow(t_gammaHitY[n_gammaHit] - t_vtx[1],2)+
                            pow(t_gammaHitZ[n_gammaHit] - t_vtx[2],2),0.5);
                    //length from starting point
                    float trackLengthStart = pow(
                            pow(t_gammaHitX[n_gammaHit] - t_gammaStartingPointX[n_gammaHit],2)+
                            pow(t_gammaHitY[n_gammaHit] - t_gammaStartingPointY[n_gammaHit],2)+
                            pow(t_gammaHitZ[n_gammaHit] - t_gammaStartingPointZ[n_gammaHit],2),0.5);

                    //calculate signal window; time of flight
                    float signalWindow = t_gammaHitT[n_gammaHit] - t_vtxTime;

                    //time b/w starting point and hit
                    float signalWindowStart = t_gammaHitT[n_gammaHit] - t_vtxTime - t_gammaStartingPointT[n_gammaHit];
                    float signalWindowSmear = t_gammaHitSmearT[n_gammaHit] - t_vtxTime;

                    if(signalWindow > 0)
                    {
                        earliest_gamma_hit.timeWindow = signalWindow;
                        earliest_gamma_hit.timeStartT = signalWindowStart;
                        earliest_gamma_hit.lengthStart = trackLengthStart;
                        earliest_gamma_hit.timeWindowSmear = signalWindowSmear;
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
                        earliest_gamma_hit.gammaTrueE = t_gammaTrueE[n_gammaHit]+939.565;
                        earliest_gamma_hit.CubeE = t_gammaCubeE[n_gammaHit];

                        earliest_gamma_hit.startingPoint[0] = t_gammaStartingPointX[n_gammaHit];
                        earliest_gamma_hit.startingPoint[1] = t_gammaStartingPointY[n_gammaHit];
                        earliest_gamma_hit.startingPoint[2] = t_gammaStartingPointZ[n_gammaHit];
                        earliest_gamma_hit.startingPointT = t_gammaStartingPointT[n_gammaHit];

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

            //select earliest hit 
            if(earliest_gamma_hit.isEmpty == false && earliest_gamma_hit.timeWindow < earliest_neutron_hit.timeWindow)
            {
                earliest_hit = earliest_gamma_hit;
            }
            if(earliest_neutron_hit.isEmpty == false && earliest_gamma_hit.timeWindow > earliest_neutron_hit.timeWindow)
            {
                earliest_hit = earliest_neutron_hit;
            }

            if(earliest_hit.isEmpty)
                continue;

            if(earliest_hit.isNeutron)
            {
                if(earliest_hit.parentId == -1)
                    earliest_hit.category = 1;
                if(earliest_hit.parentId >= 0)
                    earliest_hit.category = 2;
            }
            if(earliest_hit.isGamma)
            {
                if(earliest_hit.parentId == -1)
                    earliest_hit.category = 3;
                if(earliest_hit.parentId > 0)
                    earliest_hit.category = 4;
            }

            if(earliest_hit.category == -1)
                continue;

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
                //signal
                if(earliest_hit.category == 1)
                {
                    leverarm_signal->Fill(earliest_hit.trackLength);
                    angle_signal->Fill(GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit));
                    beta_signal->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                    distance_signal->Fill(GetDistance(earliest_hit.piDeath,earliest_hit.hit));
                    TOF_signal->Fill(earliest_hit.timeWindow);
                    CubeE_signal->Fill(earliest_hit.CubeE);
                    signalLeverArm = earliest_hit.trackLength;
                    signalAngle = GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit);
                    signalBeta = (earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity;
                    signalDistanceCHit = GetDistance(earliest_hit.piDeath,earliest_hit.hit);
                    signalTOF = earliest_hit.timeWindow;
                    signalCubeE = earliest_hit.CubeE;
                }

                //secondary neutron
                if(earliest_hit.category == 2)
                {
                    if(abs(earliest_hit.piDeath[0]) < 120 && abs(earliest_hit.piDeath[1]) < 120 && abs(earliest_hit.piDeath[2]) < 100)
                    {
                        leverarm_secondary_neutron->Fill(earliest_hit.trackLength);
                        angle_secondary_neutron->Fill(GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit));
                        beta_secondary_neutron->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                        distance_secondary_neutron->Fill(GetDistance(earliest_hit.piDeath,earliest_hit.hit));
                        TOF_secondary_neutron->Fill(earliest_hit.timeWindow);
                        CubeE_secondary_neutron->Fill(earliest_hit.CubeE);
                        secondaryNeutronLeverArm = earliest_hit.trackLength;
                        secondaryNeutronAngle = GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit);
                        secondaryNeutronBeta = (earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity;
                        secondaryNeutronDistanceCHit = GetDistance(earliest_hit.piDeath,earliest_hit.hit);
                        secondaryNeutronTOF = earliest_hit.timeWindow;
                        secondaryNeutronCubeE = earliest_hit.CubeE;
                        num_secondary_neutron += 1;
                    }
                  
                    if(abs(earliest_hit.piDeath[0]) > 120 || abs(earliest_hit.piDeath[1]) > 120 || abs(earliest_hit.piDeath[2]) > 100)
                        num_secondary_neutron_C_out3dst += 1;
                }

                //primary gamma 
                if(earliest_hit.category == 3)
                {
                    leverarm_primary_gamma->Fill(earliest_hit.trackLength);
                    angle_primary_gamma->Fill(GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit));
                    beta_primary_gamma->Fill((earliest_hit.trackLength/(earliest_hit.timeWindow-1))/c_velocity);
                    distance_primary_gamma->Fill(GetDistance(earliest_hit.piDeath,earliest_hit.hit));
                    TOF_primary_gamma->Fill(earliest_hit.timeWindow);
                    CubeE_primary_gamma->Fill(earliest_hit.CubeE);
                    primaryGammaLeverArm = earliest_hit.trackLength;
                    primaryGammaAngle = GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit);
                    primaryGammaBeta = (earliest_hit.trackLength/(earliest_hit.timeWindow-1))/c_velocity;
                    primaryGammaDistanceCHit = GetDistance(earliest_hit.piDeath,earliest_hit.hit);
                    primaryGammaTOF = earliest_hit.timeWindow;
                    primaryGammaCubeE = earliest_hit.CubeE;
                }

                //secondary gamma
                if(earliest_hit.category == 4)
                {
                    leverarm_secondary_gamma->Fill(earliest_hit.trackLength);
                    angle_secondary_gamma->Fill(GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit));
                    beta_secondary_gamma->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                    distance_secondary_gamma->Fill(GetDistance(earliest_hit.piDeath,earliest_hit.hit));
                    TOF_secondary_gamma->Fill(earliest_hit.timeWindow);
                    CubeE_secondary_gamma->Fill(earliest_hit.CubeE);
                    secondaryGammaLeverArm = earliest_hit.trackLength;
                    secondaryGammaAngle = GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit);
                    secondaryGammaBeta = (earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity;
                    secondaryGammaDistanceCHit = GetDistance(earliest_hit.piDeath,earliest_hit.hit);
                    secondaryGammaTOF = earliest_hit.timeWindow;
                    secondaryGammaCubeE = earliest_hit.CubeE;
                }
            }    //end of pion case


            ////proton case
            if(num_fspi == 0 && num_fsp ==1)
            {
                //signal
                if(earliest_hit.category == 1)
                {
                    leverarm_signal->Fill(earliest_hit.trackLength);
                    angle_signal->Fill(GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit));
                    beta_signal->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                    distance_signal->Fill(GetDistance(earliest_hit.protonDeath,earliest_hit.hit));
                    TOF_signal->Fill(earliest_hit.timeWindow);
                    CubeE_signal->Fill(earliest_hit.CubeE);
                    signalLeverArm = earliest_hit.trackLength;
                    signalAngle = GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit);
                    signalBeta = (earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity;
                    signalDistanceCHit = GetDistance(earliest_hit.protonDeath,earliest_hit.hit);
                    signalTOF = earliest_hit.timeWindow;
                    signalCubeE = earliest_hit.CubeE;
                }

                //secondary neutron
                if(earliest_hit.category == 2)
                {
                    if(abs(earliest_hit.protonDeath[0]) < 120 && abs(earliest_hit.protonDeath[1]) < 120 && earliest_hit.protonDeath[2] < 150 && earliest_hit.protonDeath[2] > -50)
                    {
                        if(earliest_hit.isFromProton)
                        {
                            leverarm_secondary_neutron->Fill(earliest_hit.trackLength);
                            angle_secondary_neutron->Fill(GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit));
                            beta_secondary_neutron->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                            distance_secondary_neutron->Fill(GetDistance(earliest_hit.protonDeath,earliest_hit.hit));
                            TOF_secondary_neutron->Fill(earliest_hit.timeWindow);
                            CubeE_secondary_neutron->Fill(earliest_hit.CubeE);
                            secondaryNeutronLeverArm = earliest_hit.trackLength;
                            secondaryNeutronAngle = GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit);
                            secondaryNeutronBeta = (earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity;
                            secondaryNeutronDistanceCHit = GetDistance(earliest_hit.protonDeath,earliest_hit.hit);
                            secondaryNeutronTOF = earliest_hit.timeWindow;
                            secondaryNeutronCubeE = earliest_hit.CubeE;
                        }
                    }
                }

                //primary gamma
                if(earliest_hit.category == 3)
                {
                    leverarm_primary_gamma->Fill(earliest_hit.trackLength);
                    angle_primary_gamma->Fill(GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit));
                    beta_primary_gamma->Fill((earliest_hit.trackLength/(earliest_hit.timeWindow-1))/c_velocity);
                    distance_primary_gamma->Fill(GetDistance(earliest_hit.protonDeath,earliest_hit.hit));
                    TOF_primary_gamma->Fill(earliest_hit.timeWindow);
                    CubeE_primary_gamma->Fill(earliest_hit.CubeE);
                    primaryGammaLeverArm = earliest_hit.trackLength;
                    primaryGammaAngle = GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit);
                    primaryGammaBeta = (earliest_hit.trackLength/(earliest_hit.timeWindow-1))/c_velocity;
                    primaryGammaDistanceCHit = GetDistance(earliest_hit.protonDeath,earliest_hit.hit);
                    primaryGammaTOF = earliest_hit.timeWindow;
                    primaryGammaCubeE = earliest_hit.CubeE;
                }

                //secondary gamma
                if(earliest_hit.category == 4)
                {
                    leverarm_secondary_gamma->Fill(earliest_hit.trackLength);
                    angle_secondary_gamma->Fill(GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit));
                    beta_secondary_gamma->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                    distance_secondary_gamma->Fill(GetDistance(earliest_hit.protonDeath,earliest_hit.hit));
                    TOF_secondary_gamma->Fill(earliest_hit.timeWindow);
                    CubeE_secondary_gamma->Fill(earliest_hit.CubeE);
                    secondaryGammaLeverArm = earliest_hit.trackLength;
                    secondaryGammaAngle = GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit);
                    secondaryGammaBeta = (earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity;
                    secondaryGammaDistanceCHit = GetDistance(earliest_hit.protonDeath,earliest_hit.hit);
                    secondaryGammaTOF = earliest_hit.timeWindow;
                    secondaryGammaCubeE = earliest_hit.CubeE;
                }
            }   //end of proton case


            /*
            ////pi case
            if(num_fspi == 1 && num_fsp ==0)
            {
                //signal
                if(earliest_hit.category == 1)
                {
                    leverarm_signal->Fill(earliest_hit.trackLength);
                    leverArm = earliest_hit.trackLength;
                    angle_signal->Fill(GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit));
                    angle = GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit);
                    beta_signal->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                    beta = (earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity;
                    distance_signal->Fill(GetDistance(earliest_hit.piDeath,earliest_hit.hit));
                    distanceCHit = GetDistance(earliest_hit.piDeath,earliest_hit.hit);
                    TOF_signal->Fill(earliest_hit.timeWindow);
                    tof = earliest_hit.timeWindow;
                    CubeE_signal->Fill(earliest_hit.CubeE);
                    cubeE = earliest_hit.CubeE;
                }

                //secondary neutron
                if(earliest_hit.category == 2)
                {
                    if(abs(earliest_hit.piDeath[0]) < 120 && abs(earliest_hit.piDeath[1]) < 120 && abs(earliest_hit.piDeath[2]) < 100)
                    {
                        leverarm_secondary_neutron->Fill(earliest_hit.trackLength);
                        leverArm = earliest_hit.trackLength;
                        angle_secondary_neutron->Fill(GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit));
                        angle = GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit);
                        beta_secondary_neutron->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                        beta = (earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity;
                        distance_secondary_neutron->Fill(GetDistance(earliest_hit.piDeath,earliest_hit.hit));
                        distanceCHit = GetDistance(earliest_hit.piDeath,earliest_hit.hit);
                        TOF_secondary_neutron->Fill(earliest_hit.timeWindow);
                        tof = earliest_hit.timeWindow;
                        CubeE_secondary_neutron->Fill(earliest_hit.CubeE);
                        cubeE = earliest_hit.CubeE;
                        num_secondary_neutron += 1;
                    }
                  
                    if(abs(earliest_hit.piDeath[0]) > 120 || abs(earliest_hit.piDeath[1]) > 120 || abs(earliest_hit.piDeath[2]) > 100)
                        num_secondary_neutron_C_out3dst += 1;
                }

                //primary gamma 
                if(earliest_hit.category == 3)
                {
                    leverarm_primary_gamma->Fill(earliest_hit.trackLength);
                    leverArm = earliest_hit.trackLength;
                    angle_primary_gamma->Fill(GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit));
                    angle = GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit);
                    beta_primary_gamma->Fill((earliest_hit.trackLength/(earliest_hit.timeWindow-1))/c_velocity);
                    beta = (earliest_hit.trackLength/(earliest_hit.timeWindow-1))/c_velocity;
                    distance_primary_gamma->Fill(GetDistance(earliest_hit.piDeath,earliest_hit.hit));
                    distanceCHit = GetDistance(earliest_hit.piDeath,earliest_hit.hit);
                    TOF_primary_gamma->Fill(earliest_hit.timeWindow);
                    tof = earliest_hit.timeWindow;
                    CubeE_primary_gamma->Fill(earliest_hit.CubeE);
                    cubeE = earliest_hit.CubeE;
                }

                //secondary gamma
                if(earliest_hit.category == 4)
                {
                    leverarm_secondary_gamma->Fill(earliest_hit.trackLength);
                    leverArm = earliest_hit.trackLength;
                    angle_secondary_gamma->Fill(GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit));
                    angle = GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit);
                    beta_secondary_gamma->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                    beta = (earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity;
                    distance_secondary_gamma->Fill(GetDistance(earliest_hit.piDeath,earliest_hit.hit));
                    distanceCHit = GetDistance(earliest_hit.piDeath,earliest_hit.hit);
                    TOF_secondary_gamma->Fill(earliest_hit.timeWindow);
                    tof = earliest_hit.timeWindow;
                    CubeE_secondary_gamma->Fill(earliest_hit.CubeE);
                    cubeE = earliest_hit.CubeE;
                }
            }    //end of pion case


            ////proton case
            if(num_fspi == 0 && num_fsp ==1)
            {
                //signal
                if(earliest_hit.category == 1)
                {
                    leverarm_signal->Fill(earliest_hit.trackLength);
                    leverArm = earliest_hit.trackLength;
                    angle_signal->Fill(GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit));
                    angle = GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit);
                    beta_signal->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                    beta = (earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity;
                    distance_signal->Fill(GetDistance(earliest_hit.protonDeath,earliest_hit.hit));
                    distanceCHit = GetDistance(earliest_hit.protonDeath,earliest_hit.hit);
                    TOF_signal->Fill(earliest_hit.timeWindow);
                    tof = earliest_hit.timeWindow;
                    CubeE_signal->Fill(earliest_hit.CubeE);
                    cubeE = earliest_hit.CubeE;
                }

                //secondary neutron
                if(earliest_hit.category == 2)
                {
                    if(abs(earliest_hit.protonDeath[0]) < 120 && abs(earliest_hit.protonDeath[1]) < 120 && earliest_hit.protonDeath[2] < 150 && earliest_hit.protonDeath[2] > -50)
                    {
                        if(earliest_hit.isFromProton)
                        {
                            leverarm_secondary_neutron->Fill(earliest_hit.trackLength);
                            leverArm = earliest_hit.trackLength;
                            angle_secondary_neutron->Fill(GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit));
                            angle = GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit);
                            beta_secondary_neutron->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                            beta = (earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity;
                            distance_secondary_neutron->Fill(GetDistance(earliest_hit.protonDeath,earliest_hit.hit));
                            distanceCHit = GetDistance(earliest_hit.protonDeath,earliest_hit.hit);
                            TOF_secondary_neutron->Fill(earliest_hit.timeWindow);
                            tof = earliest_hit.timeWindow;
                            CubeE_secondary_neutron->Fill(earliest_hit.CubeE);
                            cubeE = earliest_hit.CubeE;
                        }
                    }
                }

                //primary gamma
                if(earliest_hit.category == 3)
                {
                    leverarm_primary_gamma->Fill(earliest_hit.trackLength);
                    leverArm = earliest_hit.trackLength;
                    angle_primary_gamma->Fill(GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit));
                    angle = GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit);
                    beta_primary_gamma->Fill((earliest_hit.trackLength/(earliest_hit.timeWindow-1))/c_velocity);
                    beta = (earliest_hit.trackLength/(earliest_hit.timeWindow-1))/c_velocity;
                    distance_primary_gamma->Fill(GetDistance(earliest_hit.protonDeath,earliest_hit.hit));
                    distanceCHit = GetDistance(earliest_hit.protonDeath,earliest_hit.hit);
                    TOF_primary_gamma->Fill(earliest_hit.timeWindow);
                    tof = earliest_hit.timeWindow;
                    CubeE_primary_gamma->Fill(earliest_hit.CubeE);
                    cubeE = earliest_hit.CubeE;
                }

                //secondary gamma
                if(earliest_hit.category == 4)
                {
                    leverarm_secondary_gamma->Fill(earliest_hit.trackLength);
                    leverArm = earliest_hit.trackLength;
                    angle_secondary_gamma->Fill(GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit));
                    angle = GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit);
                    beta_secondary_gamma->Fill((earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity);
                    beta = (earliest_hit.trackLength/earliest_hit.timeWindow)/c_velocity;
                    distance_secondary_gamma->Fill(GetDistance(earliest_hit.protonDeath,earliest_hit.hit));
                    distanceCHit = GetDistance(earliest_hit.protonDeath,earliest_hit.hit);
                    TOF_secondary_gamma->Fill(earliest_hit.timeWindow);
                    tof = earliest_hit.timeWindow;
                    CubeE_secondary_gamma->Fill(earliest_hit.CubeE);
                    cubeE = earliest_hit.CubeE;
                }
            }   //end of proton case
            */
            category = earliest_hit.category;
            output_tree->Fill();
        }       //end of event iterate
        _file->Close();
        delete _file;
    }
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<"num_secondary_neutron: "<<num_secondary_neutron<<endl;
    cout<<"num_secondary_neutron_C_out3dst: "<<num_secondary_neutron_C_out3dst<<endl;

    outfile->Write();
    outfile->Close();

    cout<<endl;
    cout<<"end"<<endl;


    TFile * fi1 = new TFile("background.root","RECREATE");

    TCanvas * can = new TCanvas;

    beta_vs_leverarm_parentId_0->Write();
    beta_vs_leverarm_parentId_0->Draw();
    can->SaveAs("beta_vs_leverarm_parentId_0.pdf");
    can->Clear();

    beta_vs_leverarm_parentId_1->Write();
    beta_vs_leverarm_parentId_1->Draw();
    can->SaveAs("beta_vs_leverarm_parentId_1.pdf");
    can->Clear();

    TLegend * legend = new TLegend(0.7,0.8,0.9,0.9);
    legend->AddEntry(leverarm_signal,"signal","l");
    legend->AddEntry(leverarm_secondary_neutron,"secondary neutron","l");
    legend->AddEntry(leverarm_primary_gamma,"primary gamma","l");
    legend->AddEntry(leverarm_secondary_gamma,"secondary gamma","l");
    /*
     * 2: red, signal
     * 4: blue, secondary neutron
     * 6: purple, primary gamma
     * 8: green, secondary gamma
     */

//leaverarm{
    leverarm_signal->Write();
    leverarm_signal->SetLineColor(2);
    leverarm_signal->SetStats(0);
    leverarm_signal->Scale(1/leverarm_signal->GetEntries(),"nosw2");
    leverarm_signal->GetYaxis()->SetRangeUser(0,0.4);
    leverarm_signal->SetTitle("lever arm");
    leverarm_signal->Draw();

    leverarm_secondary_neutron->Write();
    leverarm_secondary_neutron->SetLineColor(4);
    leverarm_secondary_neutron->SetStats(0);
    leverarm_secondary_neutron->Scale(1/leverarm_secondary_neutron->GetEntries(),"nosw2");
    leverarm_secondary_neutron->Draw("same");

    leverarm_primary_gamma->Write();
    leverarm_primary_gamma->SetLineColor(6);
    leverarm_primary_gamma->SetStats(0);
    leverarm_primary_gamma->Scale(1/leverarm_primary_gamma->GetEntries(),"nosw2");
    leverarm_primary_gamma->Draw("same");

    leverarm_secondary_gamma->Write();
    leverarm_secondary_gamma->SetLineColor(8);
    leverarm_secondary_gamma->SetStats(0);
    leverarm_secondary_gamma->Scale(1/leverarm_secondary_gamma->GetEntries(),"nosw2");
    leverarm_secondary_gamma->Draw("same");

    legend->Draw();
    can->SaveAs("leverarm.pdf");
    can->Clear();
//}

//angle{
    angle_signal->Write();
    angle_signal->SetLineColor(2);
    angle_signal->SetStats(0);
    angle_signal->Scale(1/angle_signal->GetEntries(),"nosw2");
    angle_signal->GetYaxis()->SetRangeUser(0,0.4);
    angle_signal->SetTitle("angle");
    angle_signal->Draw();

    angle_secondary_neutron->Write();
    angle_secondary_neutron->SetLineColor(4);
    angle_secondary_neutron->SetStats(0);
    angle_secondary_neutron->Scale(1/angle_secondary_neutron->GetEntries(),"nosw2");
    angle_secondary_neutron->Draw("same");

    angle_primary_gamma->Write();
    angle_primary_gamma->SetLineColor(6);
    angle_primary_gamma->SetStats(0);
    angle_primary_gamma->Scale(1/angle_primary_gamma->GetEntries(),"nosw2");
    angle_primary_gamma->Draw("same");

    angle_secondary_gamma->Write();
    angle_secondary_gamma->SetLineColor(8);
    angle_secondary_gamma->SetStats(0);
    angle_secondary_gamma->Scale(1/angle_secondary_gamma->GetEntries(),"nosw2");
    angle_secondary_gamma->Draw("same");

    legend->Draw();
    can->SaveAs("angle.pdf");
    can->Clear();
//}

//beta{
    beta_signal->Write();
    beta_signal->SetLineColor(2);
    beta_signal->SetStats(0);
    beta_signal->Scale(1/beta_signal->GetEntries(),"nosw2");
    beta_signal->GetYaxis()->SetRangeUser(0,0.3);
    beta_signal->SetTitle("beta");
    beta_signal->Draw();

    beta_secondary_neutron->Write();
    beta_secondary_neutron->SetLineColor(4);
    beta_secondary_neutron->SetStats(0);
    beta_secondary_neutron->Scale(1/beta_secondary_neutron->GetEntries(),"nosw2");
    beta_secondary_neutron->Draw("same");

    beta_primary_gamma->Write();
    beta_primary_gamma->SetLineColor(6);
    beta_primary_gamma->SetStats(0);
    beta_primary_gamma->Scale(1/beta_primary_gamma->GetEntries(),"nosw2");
    beta_primary_gamma->Draw("same");

    beta_secondary_gamma->Write();
    beta_secondary_gamma->SetLineColor(8);
    beta_secondary_gamma->SetStats(0);
    beta_secondary_gamma->Scale(1/beta_secondary_gamma->GetEntries(),"nosw2");
    beta_secondary_gamma->Draw("same");

    legend->Draw();
    can->SaveAs("beta.pdf");
    can->Clear();
//}

//distance{
    distance_signal->Write();
    distance_signal->SetLineColor(2);
    distance_signal->SetStats(0);
    distance_signal->Scale(1/distance_signal->GetEntries(),"nosw2");
    distance_signal->GetYaxis()->SetRangeUser(0,0.6);
    distance_signal->SetTitle("distance b/w C and hit");
    distance_signal->Draw();

    distance_secondary_neutron->Write();
    distance_secondary_neutron->SetLineColor(4);
    distance_secondary_neutron->SetStats(0);
    distance_secondary_neutron->Scale(1/distance_secondary_neutron->GetEntries(),"nosw2");
    distance_secondary_neutron->Draw("same");

    distance_primary_gamma->Write();
    distance_primary_gamma->SetLineColor(6);
    distance_primary_gamma->SetStats(0);
    distance_primary_gamma->Scale(1/distance_primary_gamma->GetEntries(),"nosw2");
    distance_primary_gamma->Draw("same");

    distance_secondary_gamma->Write();
    distance_secondary_gamma->SetLineColor(8);
    distance_secondary_gamma->SetStats(0);
    distance_secondary_gamma->Scale(1/distance_secondary_gamma->GetEntries(),"nosw2");
    distance_secondary_gamma->Draw("same");

    legend->Draw();
    can->SaveAs("distance.pdf");
    can->Clear();
//}

//TOF{
    TOF_signal->Write();
    TOF_signal->SetLineColor(2);
    TOF_signal->SetStats(0);
    TOF_signal->Scale(1/TOF_signal->GetEntries(),"nosw2");
    TOF_signal->GetYaxis()->SetRangeUser(0,0.7);
    TOF_signal->SetTitle("Time of flight");
    TOF_signal->Draw();

    TOF_secondary_neutron->Write();
    TOF_secondary_neutron->SetLineColor(4);
    TOF_secondary_neutron->SetStats(0);
    TOF_secondary_neutron->Scale(1/TOF_secondary_neutron->GetEntries(),"nosw2");
    TOF_secondary_neutron->Draw("same");

    TOF_primary_gamma->Write();
    TOF_primary_gamma->SetLineColor(6);
    TOF_primary_gamma->SetStats(0);
    TOF_primary_gamma->Scale(1/TOF_primary_gamma->GetEntries(),"nosw2");
    TOF_primary_gamma->Draw("same");

    TOF_secondary_gamma->Write();
    TOF_secondary_gamma->SetLineColor(8);
    TOF_secondary_gamma->SetStats(0);
    TOF_secondary_gamma->Scale(1/TOF_secondary_gamma->GetEntries(),"nosw2");
    TOF_secondary_gamma->Draw("same");

    legend->Draw();
    can->SaveAs("TOF.pdf");
    can->Clear();
//}
//
//CubeE{
    CubeE_signal->Write();
    CubeE_signal->SetLineColor(2);
    CubeE_signal->SetStats(0);
    CubeE_signal->Scale(1/CubeE_signal->GetEntries(),"nosw2");
    CubeE_signal->GetYaxis()->SetRangeUser(0,0.7);
    CubeE_signal->SetTitle("CubeE");
    CubeE_signal->Draw();

    CubeE_secondary_neutron->Write();
    CubeE_secondary_neutron->SetLineColor(4);
    CubeE_secondary_neutron->SetStats(0);
    CubeE_secondary_neutron->Scale(1/CubeE_secondary_neutron->GetEntries(),"nosw2");
    CubeE_secondary_neutron->Draw("same");

    CubeE_primary_gamma->Write();
    CubeE_primary_gamma->SetLineColor(6);
    CubeE_primary_gamma->SetStats(0);
    CubeE_primary_gamma->Scale(1/CubeE_primary_gamma->GetEntries(),"nosw2");
    CubeE_primary_gamma->Draw("same");

    CubeE_secondary_gamma->Write();
    CubeE_secondary_gamma->SetLineColor(8);
    CubeE_secondary_gamma->SetStats(0);
    CubeE_secondary_gamma->Scale(1/CubeE_secondary_gamma->GetEntries(),"nosw2");
    CubeE_secondary_gamma->Draw("same");

    legend->Draw();
    can->SaveAs("CubeE.pdf");
    can->Clear();

//}
    fi1->Close();
}
