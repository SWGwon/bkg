#include <ctime>
#include <string>
#include <utility>
#include <vector>
#include <thread>

#include "TChain.h"
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
#include "THistPainter.h"
#include <TStyle.h>
#include <TROOT.h>

float leverArm;
float angle;
float beta;
float distanceCHit;
float tof;
float cubeE;
float nCube;
float category;
float neutrinoE;
float hitPDG;

float neutronE;
float neutronAngle;

using namespace std;
//histograms{

//signal{
TH1F * leverarm_signal = new TH1F("leverarm_signal","lever arm of signal",20,0,200);
TH1F * angle_signal = new TH1F("angle_signal","angle between C and signal hit",20,0,1);
TH1F * beta_signal = new TH1F("beta_signal","beta of signal",30,0,1.5);
TH1F * distance_signal = new TH1F("distance_signal","distance b/w C and signal hit",20,0,200);
TH1F * TOF_signal = new TH1F("TOF_signal","time of flight of signal", 25,0,25);
TH1F * CubeE_signal = new TH1F("CubeE_signal", "CubeE of signal", 30, 0, 15);
TH1F * nCubeDis_signal = new TH1F("nCubeDis_signal","number of cubes, signal", 50, 0, 100);
//}
//secondary neutron{ 
TH1F * leverarm_secondary_neutron = new TH1F("leverarm_secondary_neutron","lever arm of secondary_neutron",20,0,200);
TH1F * angle_secondary_neutron = new TH1F("angle_secondary_neutron","angle between C and secondary_neutron hit",20,0,1);
TH1F * beta_secondary_neutron = new TH1F("beta_secondary_neutron","beta of secondary_neutron",30,0,1.5);
TH1F * distance_secondary_neutron = new TH1F("distance_secondary_neutron","distance b/w C and secondary_neutron hit",20,0,200);
TH1F * TOF_secondary_neutron = new TH1F("TOF_secondary_neutron","time of flight of secondary_neutron", 25,0,25);
TH1F * CubeE_secondary_neutron = new TH1F("CubeE_secondary_neutron", "CubeE of secondary_neutron", 30, 0, 15);
TH1F * nCubeDis_secondary_neutron = new TH1F("nCubeDis_secondary_neutron","number of cubes, secondary neutron", 50, 0, 100);
//}

//primary gamma{
TH1F * leverarm_primary_gamma = new TH1F("leverarm_primary_gamma","lever arm of primary_gamma",20,0,200);
TH1F * angle_primary_gamma = new TH1F("angle_primary_gamma","angle between C and primary_gamma hit",20,0,1);
TH1F * beta_primary_gamma = new TH1F("beta_primary_gamma","beta of primary_gamma",30,0,1.5);
TH1F * distance_primary_gamma = new TH1F("distance_primary_gamma","distance b/w C and primary_gamma hit",20,0,200);
TH1F * TOF_primary_gamma = new TH1F("TOF_primary_gamma","time of flight of primary_gamma", 25,0,25);
TH1F * CubeE_primary_gamma = new TH1F("CubeE_primary_gamma", "CubeE of primary_gamma", 30, 0, 15);
TH1F * nCubeDis_primary_gamma = new TH1F("nCubeDis_primary_gamma","number of cubes, primary gamma", 50, 0, 100);
//}

//secondary gamma{
TH1F * leverarm_secondary_gamma = new TH1F("leverarm_secondary_gamma","lever arm of secondary_gamma",20,0,200);
TH1F * angle_secondary_gamma = new TH1F("angle_secondary_gamma","angle between C and secondary_gamma hit",20,0,1);
TH1F * beta_secondary_gamma = new TH1F("beta_secondary_gamma","beta of secondary_gamma",30,0,1.5);
TH1F * distance_secondary_gamma = new TH1F("distance_secondary_gamma","distance b/w C and secondary_gamma hit",20,0,200);
TH1F * TOF_secondary_gamma = new TH1F("TOF_secondary_gamma","time of flight of secondary_gamma", 25,0,25);
TH1F * CubeE_secondary_gamma = new TH1F("CubeE_secondary_gamma", "CubeE of secondary_gamma", 30, 0, 15);
TH1F * nCubeDis_secondary_gamma = new TH1F("nCubeDis_secondary_gamma","number of cubes, secondary gamma", 50, 0, 100);
//}

TH1F * timeWindow = new TH1F("timeWindow", "",250,0,25);
TH2F * XYPlane = new TH2F("XY","XY;X;Y",240,-120,120,240,-120,120);
TH2F * XZPlane = new TH2F("XZ","XZ;X;Z",240,-120,120,200,-100,100);
TH2F * YZPlane = new TH2F("YZ","YZ;Y;Z",240,-120,120,200,-100,100);
TH2F * cube_XYPlane = new TH2F("cube_XY","XY;X;Y",240,-120,120,240,-120,120);
TH2F * cube_XZPlane = new TH2F("cube_XZ","XZ;X;Z",240,-120,120,200,-100,100);
TH2F * cube_YZPlane = new TH2F("cube_YZ","YZ;Y;Z",240,-120,120,200,-100,100);
TH2F * XYPlane_allhits = new TH2F("XY_allhits","XY;X;Y",240,-120,120,240,-120,120);
TH2F * XZPlane_allhits = new TH2F("XZ_allhits","XZ;X;Z",240,-120,120,200,-100,100);
TH2F * YZPlane_allhits = new TH2F("YZ_allhits","YZ;Y;Z",240,-120,120,200,-100,100);
//}

float energyHitCut = 0.2; //energy deposit threshold for cube

//channel type
//int num_fspi = 0;   //number of fs charged pion
//int num_fsp = 0;    //number of fs proton

const double c_velocity = 29.9792458;

class Hit
{
    private:
        float timeWindow;           // time windows of the hit
        float timeWindowSmear;      // smeared time window 
        float X;                    // X position of hit
        float Y;                    // Y position of hit
        float Z;                    // Z position of hit
        float T;                    // true T of hit
        float vtxX;                 // X position of vertex
        float vtxY;                 // Y position of vertex
        float vtxZ;                 // Z position of vertex
        float vtxT;                 // ture T of vertex
        float leverArm;             // lever arm
        float parentId;             // parentId of hit
        float cubeE;    //neutron cube energy
        float trueE;    //neutron true energy
        float trueT;    //neutron true time
        float hitPDG;    //hit PDG
        float energyDeposit;        // energy deposited by the neutron
        bool isNeutron;
        bool isGamma;

    public:
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetTimeWindow(float timeWindow){this->timeWindow = timeWindow;};
        float GetTimeWindow(){return this->timeWindow;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetTimeWindowSmear(float timeWindowSmear){this->timeWindowSmear = timeWindowSmear;};
        float GetTimeWindowSmear(){return this->timeWindowSmear;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetX(float X){this->X = X;};
        float GetX(){return this->X;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetY(float Y){this->Y = Y;};
        float GetY(){return this->Y;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetZ(float Z){this->Z = Z;};
        float GetZ(){return this->Z;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetT(float T){this->T = T;};
        float GetT(){return this->T;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetVtxX(float vtxX){this->vtxX = vtxX;};
        float GetVtxX(){return this->vtxX;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetVtxY(float vtxY){this->vtxY = vtxY;};
        float GetVtxY(){return this->vtxY;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetVtxZ(float vtxZ){this->vtxZ = vtxZ;};
        float GetVtxZ(){return this->vtxZ;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetVtxT(float vtxT){this->vtxT = vtxT;};
        float GetVtxT(){return this->vtxT;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetLeverArm(float leverArm){this->leverArm = leverArm;};
        float GetLeverArm(){return this->leverArm;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetParentId(float parentId){this->parentId = parentId;};
        float GetParentId(){return this->parentId;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetCubeE(float cubeE){this->cubeE = cubeE;};
        float GetCubeE(){return this->cubeE;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetTrueE(float trueE){this->trueE = trueE;};
        float GetTrueE(){return this->trueE;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetTrueT(float trueT){this->trueT = trueT;};
        float GetTrueT(){return this->trueT;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetHitPDG(float hitPDG){this->hitPDG = hitPDG;};
        float GetHitPDG(){return this->hitPDG;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetEnergyDeposit(float energyDeposit){this->energyDeposit = energyDeposit;};
        float GetEnergyDeposit(){return this->energyDeposit;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetIsNeutron(bool isNeutron){this->isNeutron = isNeutron;};
        bool IsNeutron(){return this->isNeutron;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        void SetIsGamma(bool isGamma){this->isGamma = isGamma;};
        bool IsGamma(){return this->isGamma;};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        float startingPointT;
        float category;     //kind of hit(1: signal; 2: secondary neutron; 3: primary gamma; 4: secondary gamma)
        float piDeath[3];      //pion death
        float protonDeath[3];      //proton Death

        bool isFromPion;
        bool isFromProton;

        float startingPoint[3];

        int parentPdg;   // PDG of parent

        bool isEmpty;

        Hit():
            timeWindow(-1),
            timeWindowSmear(-1),    
            X(0),
            Y(0),
            Z(0),
            parentId(0),
            isNeutron(false),
            isGamma(false),
            T(0),
            energyDeposit(0),
            leverArm(0),  
            trueE(0), 
            cubeE(0),
            trueT(0), 
            hitPDG(0),
            startingPointT(0),
            parentPdg(123123123),
            category(-1),
            isEmpty(1),
            isFromPion(0),
            isFromProton(0)
    {
        for(int i = 0; i < 3; i++)
        {
            this->piDeath[i] = 0;   
            this->protonDeath[i] = 0; 
            this->startingPoint[i] = 0;
        }
    }
        ~Hit() {}
};

bool tSort(Hit Hit1, Hit Hit2)
{
    return(Hit1.GetT() < Hit2.GetT());
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
    double angle = TMath::ACos(norm_a[0]*norm_b[0]+ norm_a[1]*norm_b[1]+norm_a[2]*norm_b[2])/TMath::Pi();
    if(angle >= 0 && angle <= 1)
        return angle;
    else
        return -1000;
}

//a, b are vector
double GetDistance(float a[], Hit hit)
{
    return pow(pow(a[0]-hit.GetX(),2)+pow(a[1]-hit.GetY(),2)+pow(a[2]-hit.GetZ(),2),0.5);
}

int main()
{
    int filenum;

    cout<<"filenum :"<<endl;
    cin>>filenum;

    gErrorIgnoreLevel = kWarning;
    //if(num_fspi == 1)
    //{
    //    gSystem->mkdir("pion");
    //    gSystem->cd("pion");
    //}
    //if(num_fsp == 1)
    //{
    //    gSystem->mkdir("proton");
    //    gSystem->cd("proton");
    //}
    //if(num_fsp == 0 && num_fspi == 0)
    //{
    //    gSystem->mkdir("0pi0p");
    //    gSystem->cd("0pi0p");
    //}

    TFile * outfile = new TFile("variables.root","RECREATE");
    TTree * output_tree = new TTree("output_tree", "output_tree");

    output_tree->Branch("leverArm",&leverArm, "lever arm");
    output_tree->Branch("angle",&angle, "angle between C and hit");
    output_tree->Branch("beta",&beta, "beta");
    output_tree->Branch("distanceCHit",&distanceCHit, "distance C and hit");
    output_tree->Branch("tof",&tof, "time of flight");
    output_tree->Branch("cubeE",&cubeE, "CubeE");
    output_tree->Branch("category", &category, "category");
    output_tree->Branch("nCube", &nCube, "nCube");
    output_tree->Branch("neutrinoE", &neutrinoE, "neutrinoE");
    output_tree->Branch("hitPDG", &hitPDG, "hitPDG");

    output_tree->Branch("neutronE", &neutronE, "neutronE");
    output_tree->Branch("neutronAngle", &neutronAngle, "neutronAngle");

    TChain tree("tree");
    cout<<"---------------------------"<<endl;
    cout<<"file loading..."<<endl;

    for(int i = 1; i != filenum+1; i++)
    {
        //cout<<"\033[1APROD"<<101<<": "<<(double)(i*100/filenum)<<"%\033[1000D"<<endl;
        string file = Form("/Users/gwon/Geo12/PROD101/RHC_%d_wGamma_2ndVersion.root",i);
        //string file = Form("/pnfs/dune/persistent/users/gyang/3DST/dump/standardGeo12/PROD101/RHC_%d_wGamma_2ndVersion.root",i);
        tree.Add(TString(file));
    }

    //tree.Add("/Users/gwon/Geo12/PROD101/RHC_*_wGamma_2ndVersion.root");

    float t_neutronHitX[1000], t_neutronHitY[1000], t_neutronHitZ[1000];
    float t_neutronStartingPointX[1000], t_neutronStartingPointY[1000], t_neutronStartingPointZ[1000], t_neutronStartingPointT[1000];
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
    float t_gammaStartingPointX[1000], t_gammaStartingPointY[1000], t_gammaStartingPointZ[1000], t_gammaStartingPointT[1000];
    float t_gammaHitT[1000], t_gammaParentId[1000], t_gammaParentPDG[1000];
    float t_gammaHitE[1000]; 
    float t_gammaCubeE[1000];
    float t_gammaHitSmearT[1000];
    float t_gammaHitPDG[1000];

    int PDG = 0;
    int t_nFS, t_fsPdg[1000];

    tree.SetBranchAddress("neutronHitX", &t_neutronHitX);
    tree.SetBranchAddress("neutronHitY", &t_neutronHitY);
    tree.SetBranchAddress("neutronHitZ", &t_neutronHitZ);
    tree.SetBranchAddress("neutronStartingPointX", &t_neutronStartingPointX);
    tree.SetBranchAddress("neutronStartingPointY", &t_neutronStartingPointY);
    tree.SetBranchAddress("neutronStartingPointZ", &t_neutronStartingPointZ);
    tree.SetBranchAddress("neutronStartingPointT", &t_neutronStartingPointT);
    tree.SetBranchAddress("neutronHitT", &t_neutronHitT);
    tree.SetBranchAddress("neutronParentId", &t_neutronParentId);
    tree.SetBranchAddress("neutronParentPDG", &t_neutronParentPDG);
    tree.SetBranchAddress("neutronHitE", &t_neutronHitE);
    tree.SetBranchAddress("neutronTrueE", &t_neutronTrueE);
    tree.SetBranchAddress("neutronCubeE", &t_neutronCubeE);
    tree.SetBranchAddress("neutronHitSmearT", &t_neutronHitSmearT);
    tree.SetBranchAddress("neutronHitPDG", &t_neutronHitPDG);
    tree.SetBranchAddress("vtx", &t_vtx);
    tree.SetBranchAddress("vtxTime", &t_vtxTime);
    tree.SetBranchAddress("nFS", &t_nFS);
    tree.SetBranchAddress("fsPdg", &t_fsPdg);
    tree.SetBranchAddress("piDeath", &t_piDeath);
    tree.SetBranchAddress("protonDeath", &t_protonDeath);

    tree.SetBranchAddress("gammaHitX", &t_gammaHitX);
    tree.SetBranchAddress("gammaHitY", &t_gammaHitY);
    tree.SetBranchAddress("gammaHitZ", &t_gammaHitZ);
    tree.SetBranchAddress("gammaStartingPointX", &t_gammaStartingPointX);
    tree.SetBranchAddress("gammaStartingPointY", &t_gammaStartingPointY);
    tree.SetBranchAddress("gammaStartingPointZ", &t_gammaStartingPointZ);
    tree.SetBranchAddress("gammaStartingPointT", &t_gammaStartingPointT);
    tree.SetBranchAddress("gammaHitT", &t_gammaHitT);
    tree.SetBranchAddress("gammaParentId", &t_gammaParentId);
    tree.SetBranchAddress("gammaParentPDG", &t_gammaParentPDG);
    tree.SetBranchAddress("gammaHitE", &t_gammaHitE);
    tree.SetBranchAddress("gammaCubeE", &t_gammaCubeE);
    tree.SetBranchAddress("gammaHitT", &t_gammaHitT);
    tree.SetBranchAddress("gammaHitSmearT", &t_gammaHitSmearT);
    tree.SetBranchAddress("gammaHitPDG", &t_gammaHitPDG);

    int nevents = tree.GetEntries();

    cout<<"file loading is done"<<endl;
    cout<<"---------------------------"<<endl;
    cout<<"total entries: "<<nevents<<endl;
    cout<<"event loop starts"<<endl;
    cout<<endl;

    int temp = 9907;
    for(int event = 0; event < nevents; event++)
        //for(int event = temp; event < temp+1; event++)
    {
        tree.GetEntry(event);

        cout<<"\033[1Aevent: "<<(double)(event*100/nevents)<<"% ,"<<event<<"\033[1000D"<<endl;
        leverArm = -1000;
        angle = -1000;
        beta = -1000;
        distanceCHit = -1000;
        tof = -1000;
        cubeE = -1000;
        category = -1000;
        nCube = -1000;
        hitPDG = -1000;

        neutronE = -1000;
        neutronAngle = -1000;

        //out of fiducial volume
        if(abs(t_vtx[0]) > 50 || abs(t_vtx[1]) > 50 || abs(t_vtx[2]) > 50)
            continue;

        //check whether it's CC event or not
        bool is_CC = false;
        for(int inFS = 0; inFS < t_nFS; inFS++)
        {
            if(abs(t_fsPdg[inFS]) == 11 || abs(t_fsPdg[inFS]) == 13)    //electronPDG=11,muonPDG=13
            {
                is_CC = true;
                break;
            }
        }
        //if it's not CC skip this event
        if(!is_CC)
            continue;

        //count # of charged pion,proton in FS
        int num_pi = 0;
        int num_proton = 0;
        for(int inFS = 0; inFS < t_nFS; inFS++)
        {
            if(abs(t_fsPdg[inFS]) == 211)    //pionPDG=+-211
            {
                num_pi++;
            }

            if(abs(t_fsPdg[inFS]) == 2212)    //protonPDG=+-211
            {
                num_proton++;
            }
        }

        //flags for channels
        bool _1pi0p = false;
        bool _0pi1p = false;
        bool _0pi0p = false;

        if(num_pi == 1 && num_proton == 0)
            _1pi0p = true;
        if(num_pi == 0 && num_proton == 1)
            _0pi1p = true;
        if(num_pi == 0 && num_proton == 0)
            _0pi0p = true;

        if(!_1pi0p && !_0pi1p && !_0pi0p)
            continue;

//        if(_0pi1p || _0pi0p)    //only looking for 1pi0p channel
//            continue;

        Hit temp_neutron_Hit;   
        std::vector<Hit> vectorHit;    //vector of all neutron+gamma hits

        //push_back satisying hits to vectorHit
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

            //calculate lever arm
            float leverArm = pow(
                    pow(t_neutronHitX[n_neutronHit] - t_vtx[0],2)+
                    pow(t_neutronHitY[n_neutronHit] - t_vtx[1],2)+
                    pow(t_neutronHitZ[n_neutronHit] - t_vtx[2],2),0.5);

            //calculate signal window; time of flight
            float signalWindow = t_neutronHitT[n_neutronHit] - t_vtxTime - 1;
            float signalWindowSmear = t_neutronHitSmearT[n_neutronHit] - t_vtxTime;

            //Fix a bug from edep-sim
            if(signalWindow == 1)
                signalWindow = 0.5;


            if(signalWindow > 0)
            {
                temp_neutron_Hit.SetTimeWindow(signalWindow);
                temp_neutron_Hit.SetTimeWindowSmear(signalWindowSmear);
                temp_neutron_Hit.SetLeverArm(leverArm);
                temp_neutron_Hit.SetEnergyDeposit(t_neutronHitE[n_neutronHit]);

                temp_neutron_Hit.SetVtxX(t_vtx[0]);
                temp_neutron_Hit.SetVtxY(t_vtx[1]);
                temp_neutron_Hit.SetVtxZ(t_vtx[2]);

                temp_neutron_Hit.piDeath[0] = t_piDeath[0];
                temp_neutron_Hit.piDeath[1] = t_piDeath[1];
                temp_neutron_Hit.piDeath[2] = t_piDeath[2];

                temp_neutron_Hit.protonDeath[0] = t_protonDeath[0];
                temp_neutron_Hit.protonDeath[1] = t_protonDeath[1];
                temp_neutron_Hit.protonDeath[2] = t_protonDeath[2];

                temp_neutron_Hit.SetTrueT(t_neutronHitT[n_neutronHit]);
                temp_neutron_Hit.SetTrueE(t_neutronTrueE[n_neutronHit]);
                temp_neutron_Hit.SetCubeE(t_neutronCubeE[n_neutronHit]);

                temp_neutron_Hit.startingPoint[0] = t_neutronStartingPointX[n_neutronHit];
                temp_neutron_Hit.startingPoint[1] = t_neutronStartingPointY[n_neutronHit];
                temp_neutron_Hit.startingPoint[2] = t_neutronStartingPointZ[n_neutronHit];

                temp_neutron_Hit.SetParentId(t_neutronParentId[n_neutronHit]);
                temp_neutron_Hit.parentPdg = t_neutronParentPDG[n_neutronHit];
                temp_neutron_Hit.SetHitPDG(t_neutronHitPDG[n_neutronHit]);

                temp_neutron_Hit.SetVtxT(t_vtxTime);
                temp_neutron_Hit.isEmpty = 0;
                if(t_neutronStartingPointX[n_neutronHit] == t_piDeath[0])
                    temp_neutron_Hit.isFromPion = 1;
                else
                    temp_neutron_Hit.isFromPion = 0;
                if(t_neutronStartingPointX[n_neutronHit] == t_protonDeath[0])
                    temp_neutron_Hit.isFromProton = 1;
                else
                    temp_neutron_Hit.isFromProton = 0;
                temp_neutron_Hit.SetT(t_neutronHitT[n_neutronHit]);
                temp_neutron_Hit.SetX(t_neutronHitX[n_neutronHit]);
                temp_neutron_Hit.SetY(t_neutronHitY[n_neutronHit]);
                temp_neutron_Hit.SetZ(t_neutronHitZ[n_neutronHit]);
                temp_neutron_Hit.SetParentId(t_neutronParentId[n_neutronHit]);
                temp_neutron_Hit.SetIsNeutron(true);
                temp_neutron_Hit.SetCubeE(t_neutronCubeE[n_neutronHit]);

                vectorHit.push_back(temp_neutron_Hit);
            }
        }

        Hit temp_gamma_Hit;

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

            //calculate lever arm
            float leverArm = pow(
                    pow(t_gammaHitX[n_gammaHit] - t_vtx[0],2)+
                    pow(t_gammaHitY[n_gammaHit] - t_vtx[1],2)+
                    pow(t_gammaHitZ[n_gammaHit] - t_vtx[2],2),0.5);

            //calculate signal window; time of flight
            float signalWindow = t_gammaHitT[n_gammaHit] - t_vtxTime - 1;
            float signalWindowSmear = t_gammaHitSmearT[n_gammaHit] - t_vtxTime;

            //Fix a bug from edep-sim
            if(signalWindow == 1)
                signalWindow = 0.5;


            if(signalWindow > 0)
            {
                temp_gamma_Hit.SetTimeWindow(signalWindow);
                temp_gamma_Hit.SetTimeWindowSmear(signalWindowSmear);
                temp_gamma_Hit.SetLeverArm(leverArm);
                temp_gamma_Hit.SetEnergyDeposit(t_gammaHitE[n_gammaHit]);

                temp_gamma_Hit.SetVtxX(t_vtx[0]);
                temp_gamma_Hit.SetVtxY(t_vtx[1]);
                temp_gamma_Hit.SetVtxZ(t_vtx[2]);

                temp_gamma_Hit.piDeath[0] = t_piDeath[0];
                temp_gamma_Hit.piDeath[1] = t_piDeath[1];
                temp_gamma_Hit.piDeath[2] = t_piDeath[2];

                temp_gamma_Hit.protonDeath[0] = t_protonDeath[0];
                temp_gamma_Hit.protonDeath[1] = t_protonDeath[1];
                temp_gamma_Hit.protonDeath[2] = t_protonDeath[2];

                temp_gamma_Hit.SetTrueT(t_gammaHitT[n_gammaHit]);
                temp_gamma_Hit.SetCubeE(t_gammaCubeE[n_gammaHit]);

                temp_gamma_Hit.startingPoint[0] = t_gammaStartingPointX[n_gammaHit];
                temp_gamma_Hit.startingPoint[1] = t_gammaStartingPointY[n_gammaHit];
                temp_gamma_Hit.startingPoint[2] = t_gammaStartingPointZ[n_gammaHit];

                temp_gamma_Hit.SetParentId(t_gammaParentId[n_gammaHit]);
                temp_gamma_Hit.parentPdg = t_gammaParentPDG[n_gammaHit];
                temp_gamma_Hit.SetHitPDG(t_gammaHitPDG[n_gammaHit]);

                temp_gamma_Hit.SetVtxT(t_vtxTime);
                temp_gamma_Hit.isEmpty = 0;
                if(t_gammaStartingPointX[n_gammaHit] == t_piDeath[0])
                    temp_gamma_Hit.isFromPion = 1;
                else
                    temp_gamma_Hit.isFromPion = 0;
                if(t_gammaStartingPointX[n_gammaHit] == t_protonDeath[0])
                    temp_gamma_Hit.isFromProton = 1;
                else
                    temp_gamma_Hit.isFromProton = 0;
                temp_gamma_Hit.SetT(t_gammaHitT[n_gammaHit]);
                temp_gamma_Hit.SetX(t_gammaHitX[n_gammaHit]);
                temp_gamma_Hit.SetY(t_gammaHitY[n_gammaHit]);
                temp_gamma_Hit.SetZ(t_gammaHitZ[n_gammaHit]);
                temp_gamma_Hit.SetParentId(t_gammaParentId[n_gammaHit]);
                temp_gamma_Hit.SetIsGamma(true);
                temp_gamma_Hit.SetCubeE(t_gammaCubeE[n_gammaHit]);

                vectorHit.push_back(temp_gamma_Hit);
            }
        }

        if(vectorHit.size() == 0)
            continue;

        //this makes beta <<< 1, 1 spill is 10us, skip this event
        if(vectorHit.at(0).GetTrueT() > 10000)
            continue;

        //sort by time
        if(vectorHit.size() != 0)
            std::sort(vectorHit.begin(), vectorHit.end(), tSort);

        //map<position,pair<true,time>>
        std::map<std::tuple<float,float,float>,std::pair<int,float>> cube_fired;

        for(auto t:vectorHit)
        {
            //cout<<"x,y,z: "<<t.GetX()<<", "<<t.GetY()<<", "<<t.GetZ()<<endl;
            XYPlane_allhits->Fill(t.GetX(),t.GetY());
            YZPlane_allhits->Fill(t.GetY(),t.GetZ());
            XZPlane_allhits->Fill(t.GetX(),t.GetZ());
            std::tuple<float, float, float> temp_position = std::make_tuple(t.GetX(),t.GetY(),t.GetZ());
            cube_fired.insert(make_pair(temp_position,make_pair(1,t.GetT())));  //모든 힛들이 제대로 들어가는것 확인완료
        }

        //cube cluster vector
        std::vector<Hit> cube_cluster;
        
        //push back the first hit
        cube_cluster.push_back(vectorHit.at(0));
        //erase that hit from cube_fired vector
        cube_fired.erase(std::make_tuple(cube_cluster.at(0).GetX(),cube_cluster.at(0).GetY(),cube_cluster.at(0).GetZ()));
        //cout<<"first hit: "<<cube_cluster.at(0).GetX()<<" ,"<<cube_cluster.at(0).GetY()<<" ,"<<cube_cluster.at(0).GetZ()<<", "<<cube_cluster.at(0).GetT()<<endl;

        bool size_check = true;
        while(size_check)
        {
            int size_before = cube_cluster.size();
            for(auto t:cube_cluster)
            {
                std::tuple<float, float, float> direction[26];
                direction[0]  = std::make_tuple(t.GetX()-1,t.GetY(),t.GetZ());
                direction[1]  = std::make_tuple(t.GetX()+1,t.GetY(),t.GetZ());
                direction[2]  = std::make_tuple(t.GetX(),t.GetY()-1,t.GetZ());
                direction[3]  = std::make_tuple(t.GetX(),t.GetY()+1,t.GetZ());
                direction[4]  = std::make_tuple(t.GetX(),t.GetY(),t.GetZ()-1);
                direction[5]  = std::make_tuple(t.GetX(),t.GetY(),t.GetZ()+1);
                direction[6]  = std::make_tuple(t.GetX()+1,t.GetY()+1,t.GetZ());
                direction[7]  = std::make_tuple(t.GetX()+1,t.GetY(),t.GetZ()+1);
                direction[8]  = std::make_tuple(t.GetX()+1,t.GetY()-1,t.GetZ());
                direction[9]  = std::make_tuple(t.GetX()+1,t.GetY(),t.GetZ()-1);
                direction[10] = std::make_tuple(t.GetX()-1,t.GetY()+1,t.GetZ());
                direction[11] = std::make_tuple(t.GetX()-1,t.GetY(),t.GetZ()+1);
                direction[12] = std::make_tuple(t.GetX()-1,t.GetY()-1,t.GetZ());
                direction[13] = std::make_tuple(t.GetX()-1,t.GetY(),t.GetZ()-1);
                direction[14] = std::make_tuple(t.GetX(),t.GetY()+1,t.GetZ()+1);
                direction[15] = std::make_tuple(t.GetX(),t.GetY()-1,t.GetZ()+1);
                direction[16] = std::make_tuple(t.GetX(),t.GetY()+1,t.GetZ()-1);
                direction[17] = std::make_tuple(t.GetX(),t.GetY()-1,t.GetZ()-1);
                direction[18] = std::make_tuple(t.GetX()+1,t.GetY()+1,t.GetZ()+1);
                direction[19] = std::make_tuple(t.GetX()-1,t.GetY()-1,t.GetZ()-1);
                direction[20] = std::make_tuple(t.GetX()-1,t.GetY()+1,t.GetZ()+1);
                direction[21] = std::make_tuple(t.GetX()+1,t.GetY()-1,t.GetZ()+1);
                direction[22] = std::make_tuple(t.GetX()+1,t.GetY()+1,t.GetZ()-1);
                direction[23] = std::make_tuple(t.GetX()+1,t.GetY()-1,t.GetZ()-1);
                direction[24] = std::make_tuple(t.GetX()-1,t.GetY()+1,t.GetZ()-1);
                direction[25] = std::make_tuple(t.GetX()-1,t.GetY()-1,t.GetZ()+1);
                for(int i = 0; i < 26; i++)
                {
                    if(cube_fired.find(direction[i])->second.first == 1 && abs(t.GetT()-cube_fired.find(direction[i])->second.second) < 1)
                    {
                        Hit temp_hit;
                        temp_hit.SetX(std::get<0>(direction[i]));
                        temp_hit.SetY(std::get<1>(direction[i]));
                        temp_hit.SetZ(std::get<2>(direction[i]));
                        temp_hit.SetT(cube_fired.find(direction[i])->second.second);
                        //cout<<std::get<0>(direction[i])<<", "<<std::get<1>(direction[i])<<", "<<std::get<2>(direction[i])<<endl;
                        cube_fired.erase(direction[i]);
                        cube_cluster.push_back(temp_hit);
                    }
                }
            }
            int size_after = cube_cluster.size();
            if(size_after == size_before)
                size_check = false;
        }

        for(auto t:cube_cluster)
        {
            cube_XYPlane->Fill(t.GetX(),t.GetY());
            cube_XZPlane->Fill(t.GetX(),t.GetZ());
            cube_YZPlane->Fill(t.GetY(),t.GetZ());
        }

        //gStyle->SetFrameFillColor(1);
        //if(XYPlane_allhits->GetEntries() != 0)
        //{
        //    TCanvas * can1 = new TCanvas();
        //    can1->Divide(3,2);
        //    can1->cd(1);
        //    XYPlane_allhits->SetStats(0);
        //    XYPlane_allhits->Draw("colz");
        //    can1->cd(2);
        //    XZPlane_allhits->SetStats(0);
        //    XZPlane_allhits->Draw("colz");
        //    can1->cd(3);
        //    YZPlane_allhits->SetStats(0);
        //    YZPlane_allhits->Draw("colz");

        //    can1->cd(4);
        //    cube_XYPlane->SetStats(0);
        //    cube_XYPlane->Draw("colz");
        //    can1->cd(5);
        //    cube_XZPlane->SetStats(0);
        //    cube_XZPlane->Draw("colz");
        //    can1->cd(6);
        //    cube_YZPlane->SetStats(0);
        //    cube_YZPlane->Draw("colz");
        //    can1->SaveAs(Form("cube_eventview_%d.pdf",event));
        //    can1->Clear();
        //}
        //XYPlane_allhits->Reset();
        //XZPlane_allhits->Reset();
        //YZPlane_allhits->Reset();
        //cube_XYPlane->Reset();
        //cube_XZPlane->Reset();
        //cube_YZPlane->Reset();

        //select earliest hit 
        Hit earliest_hit;
        earliest_hit = vectorHit.at(0);

        if(earliest_hit.IsNeutron())
        {
            if(earliest_hit.GetParentId() == -1)
            {
                earliest_hit.category = 1;
            }
            if(earliest_hit.GetParentId() >= 0)
                earliest_hit.category = 2;
        }
        if(earliest_hit.IsGamma())
        {
            if(earliest_hit.GetParentId() == -1)
                earliest_hit.category = 3;
            if(earliest_hit.GetParentId() > 0)
                earliest_hit.category = 4;
        }

        if(earliest_hit.category == -1000)
            continue;
        if(earliest_hit.GetHitPDG() > 10000)
            continue;

        for(int i = 0; i < 3; i++)
        {
            if(i == 0)
            {
                vec_piDeath_to_hit[i] = earliest_hit.GetX()-earliest_hit.piDeath[i];
                vec_vtx_to_piDeath[i] = earliest_hit.piDeath[i]-earliest_hit.GetVtxX();
                vec_protonDeath_to_hit[i] = earliest_hit.GetX()-earliest_hit.protonDeath[i];
                vec_vtx_to_protonDeath[i] = earliest_hit.protonDeath[i]-earliest_hit.GetVtxX();
            }
            if(i == 1)
            {
                vec_piDeath_to_hit[i] = earliest_hit.GetY()-earliest_hit.piDeath[i];
                vec_vtx_to_piDeath[i] = earliest_hit.piDeath[i]-earliest_hit.GetVtxY();
                vec_protonDeath_to_hit[i] = earliest_hit.GetY()-earliest_hit.protonDeath[i];
                vec_vtx_to_protonDeath[i] = earliest_hit.protonDeath[i]-earliest_hit.GetVtxY();
            }
            if(i == 2)
            {
                vec_piDeath_to_hit[i] = earliest_hit.GetZ()-earliest_hit.piDeath[i];
                vec_vtx_to_piDeath[i] = earliest_hit.piDeath[i]-earliest_hit.GetVtxZ();
                vec_protonDeath_to_hit[i] = earliest_hit.GetZ()-earliest_hit.protonDeath[i];
                vec_vtx_to_protonDeath[i] = earliest_hit.protonDeath[i]-earliest_hit.GetVtxZ();
            }
        }

        float Z[3] = {0,0,1};
        float vec_vtx_to_hit[3];
        vec_vtx_to_hit[0] = earliest_hit.GetX() - earliest_hit.GetVtxX();
        vec_vtx_to_hit[1] = earliest_hit.GetY() - earliest_hit.GetVtxY();
        vec_vtx_to_hit[2] = earliest_hit.GetZ() - earliest_hit.GetVtxZ();

        //signal
        if(earliest_hit.category == 1)
        {
            leverarm_signal->Fill(earliest_hit.GetLeverArm());
            leverArm = earliest_hit.GetLeverArm();
            if(_1pi0p)
            {
                angle_signal->Fill(GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit));
                angle = GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit);
                distance_signal->Fill(GetDistance(earliest_hit.piDeath,earliest_hit));
                distanceCHit = GetDistance(earliest_hit.piDeath,earliest_hit);
            }
            if(_0pi1p)
            {
                angle_signal->Fill(GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit));
                angle = GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit);
                distance_signal->Fill(GetDistance(earliest_hit.protonDeath,earliest_hit));
                distanceCHit = GetDistance(earliest_hit.protonDeath,earliest_hit);
            }
            if(_0pi0p)
            {
                angle = -1000;
                distanceCHit = -1000;
            }
            beta_signal->Fill((earliest_hit.GetLeverArm()/earliest_hit.GetTimeWindow())/c_velocity);
            beta = (earliest_hit.GetLeverArm()/earliest_hit.GetTimeWindow())/c_velocity;
            TOF_signal->Fill(earliest_hit.GetTimeWindow());
            tof = earliest_hit.GetTimeWindow();
            CubeE_signal->Fill(earliest_hit.GetCubeE());
            cubeE = earliest_hit.GetCubeE();
            nCubeDis_signal->Fill(cube_cluster.size());
            nCube = cube_cluster.size();
            hitPDG = earliest_hit.GetHitPDG();

            neutronE = earliest_hit.GetTrueE();
            neutronAngle = GetAngle(Z,vec_vtx_to_hit);
        }

        //secondary neutron
        if(earliest_hit.category == 2)
        {
            leverarm_secondary_neutron->Fill(earliest_hit.GetLeverArm());
            leverArm = earliest_hit.GetLeverArm();
            if(_1pi0p)
            {
                angle_secondary_neutron->Fill(GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit));
                angle = GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit);
                distance_secondary_neutron->Fill(GetDistance(earliest_hit.piDeath,earliest_hit));
                distanceCHit = GetDistance(earliest_hit.piDeath,earliest_hit);
            }
            if(_0pi1p)
            {
                angle_secondary_neutron->Fill(GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit));
                angle = GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit);
                distance_secondary_neutron->Fill(GetDistance(earliest_hit.protonDeath,earliest_hit));
                distanceCHit = GetDistance(earliest_hit.protonDeath,earliest_hit);
            }
            if(_0pi0p)
            {
                angle = -1000;
                distanceCHit = -1000;
            }
            beta_secondary_neutron->Fill((earliest_hit.GetLeverArm()/earliest_hit.GetTimeWindow())/c_velocity);
            beta = (earliest_hit.GetLeverArm()/earliest_hit.GetTimeWindow())/c_velocity;
            TOF_secondary_neutron->Fill(earliest_hit.GetTimeWindow());
            tof = earliest_hit.GetTimeWindow();
            CubeE_secondary_neutron->Fill(earliest_hit.GetCubeE());
            cubeE = earliest_hit.GetCubeE();
            nCubeDis_secondary_neutron->Fill(cube_cluster.size());
            nCube = cube_cluster.size();
            hitPDG = earliest_hit.GetHitPDG();

            neutronE = earliest_hit.GetTrueE();
            neutronAngle = GetAngle(Z,vec_vtx_to_hit);
        }

        //primary gamma 
        if(earliest_hit.category == 3)
        {
            leverarm_primary_gamma->Fill(earliest_hit.GetLeverArm());
            leverArm = earliest_hit.GetLeverArm();
            if(_1pi0p)
            {
                angle_primary_gamma->Fill(GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit));
                angle = GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit);
                distance_primary_gamma->Fill(GetDistance(earliest_hit.piDeath,earliest_hit));
                distanceCHit = GetDistance(earliest_hit.piDeath,earliest_hit);
            }
            if(_0pi1p)
            {
                angle_primary_gamma->Fill(GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit));
                angle = GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit);
                distance_primary_gamma->Fill(GetDistance(earliest_hit.protonDeath,earliest_hit));
                distanceCHit = GetDistance(earliest_hit.protonDeath,earliest_hit);
            }
            if(_0pi0p)
            {
                angle = -1000;
                distanceCHit = -1000;
            }
            beta_primary_gamma->Fill((earliest_hit.GetLeverArm()/(earliest_hit.GetTimeWindow()))/c_velocity);
            beta = (earliest_hit.GetLeverArm()/(earliest_hit.GetTimeWindow()))/c_velocity;
            TOF_primary_gamma->Fill(earliest_hit.GetTimeWindow());
            tof = earliest_hit.GetTimeWindow();
            CubeE_primary_gamma->Fill(earliest_hit.GetCubeE());
            cubeE = earliest_hit.GetCubeE();
            nCubeDis_primary_gamma->Fill(cube_cluster.size());
            nCube = cube_cluster.size();
            hitPDG = earliest_hit.GetHitPDG();
        }

        //secondary gamma
        if(earliest_hit.category == 4)
        {
            leverarm_secondary_gamma->Fill(earliest_hit.GetLeverArm());
            leverArm = earliest_hit.GetLeverArm();
            if(_1pi0p)
            {
                angle_secondary_gamma->Fill(GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit));
                angle = GetAngle(vec_vtx_to_piDeath,vec_piDeath_to_hit);
                distance_secondary_gamma->Fill(GetDistance(earliest_hit.piDeath,earliest_hit));
                distanceCHit = GetDistance(earliest_hit.piDeath,earliest_hit);
            }
            if(_0pi1p)
            {
                angle_secondary_gamma->Fill(GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit));
                angle = GetAngle(vec_vtx_to_protonDeath,vec_protonDeath_to_hit);
                distance_secondary_gamma->Fill(GetDistance(earliest_hit.protonDeath,earliest_hit));
                distanceCHit = GetDistance(earliest_hit.protonDeath,earliest_hit);
            }
            if(_0pi0p)
            {
                angle = -1000;
                distanceCHit = -1000;
            }
            beta_secondary_gamma->Fill((earliest_hit.GetLeverArm()/earliest_hit.GetTimeWindow())/c_velocity);
            beta = (earliest_hit.GetLeverArm()/earliest_hit.GetTimeWindow())/c_velocity;
            TOF_secondary_gamma->Fill(earliest_hit.GetTimeWindow());
            tof = earliest_hit.GetTimeWindow();
            CubeE_secondary_gamma->Fill(earliest_hit.GetCubeE());
            cubeE = earliest_hit.GetCubeE();
            nCubeDis_secondary_gamma->Fill(cube_cluster.size());
            nCube = cube_cluster.size();
            hitPDG = earliest_hit.GetHitPDG();
        }
        category = earliest_hit.category;
        output_tree->Fill();
    }//end of event iterate

    outfile->Write();
    outfile->Close();

    cout<<"event loop is aone"<<endl;
    cout<<"---------------------------"<<endl;
    cout<<"making output files"<<endl;

    TFile * fi1 = new TFile("background.root","RECREATE");
    gStyle->SetFrameFillColor(0);

    TCanvas * can = new TCanvas;

    TLegend * legend = new TLegend(0.6,0.7,0.9,0.9);
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
    leverarm_signal->Scale(1/leverarm_signal->Integral(),"nosw2");
    leverarm_signal->GetYaxis()->SetRangeUser(0,0.4);
    leverarm_signal->SetTitle("lever arm");
    leverarm_signal->Draw();

    leverarm_secondary_neutron->Write();
    leverarm_secondary_neutron->SetLineColor(4);
    leverarm_secondary_neutron->SetStats(0);
    leverarm_secondary_neutron->Scale(1/leverarm_secondary_neutron->Integral(),"nosw2");
    leverarm_secondary_neutron->Draw("same");

    leverarm_primary_gamma->Write();
    leverarm_primary_gamma->SetLineColor(6);
    leverarm_primary_gamma->SetStats(0);
    leverarm_primary_gamma->Scale(1/leverarm_primary_gamma->Integral(),"nosw2");
    leverarm_primary_gamma->Draw("same");

    leverarm_secondary_gamma->Write();
    leverarm_secondary_gamma->SetLineColor(8);
    leverarm_secondary_gamma->SetStats(0);
    leverarm_secondary_gamma->Scale(1/leverarm_secondary_gamma->Integral(),"nosw2");
    leverarm_secondary_gamma->Draw("same");

    legend->Draw();
    can->SaveAs("leverarm.pdf");
    can->Clear();
    //}

    //angle{
    angle_signal->Write();
    angle_signal->SetLineColor(2);
    angle_signal->SetStats(0);
    angle_signal->Scale(1/angle_signal->Integral(),"nosw2");
    angle_signal->GetYaxis()->SetRangeUser(0,0.4);
    angle_signal->SetTitle("angle");
    angle_signal->Draw();

    angle_secondary_neutron->Write();
    angle_secondary_neutron->SetLineColor(4);
    angle_secondary_neutron->SetStats(0);
    angle_secondary_neutron->Scale(1/angle_secondary_neutron->Integral(),"nosw2");
    angle_secondary_neutron->Draw("same");

    angle_primary_gamma->Write();
    angle_primary_gamma->SetLineColor(6);
    angle_primary_gamma->SetStats(0);
    angle_primary_gamma->Scale(1/angle_primary_gamma->Integral(),"nosw2");
    angle_primary_gamma->Draw("same");

    angle_secondary_gamma->Write();
    angle_secondary_gamma->SetLineColor(8);
    angle_secondary_gamma->SetStats(0);
    angle_secondary_gamma->Scale(1/angle_secondary_gamma->Integral(),"nosw2");
    angle_secondary_gamma->Draw("same");

    legend->Draw();
    can->SaveAs("angle.pdf");
    can->Clear();
    //}

    //beta{
    beta_signal->Write();
    beta_signal->SetLineColor(2);
    beta_signal->SetStats(0);
    beta_signal->Scale(1/beta_signal->Integral(),"nosw2");
    beta_signal->GetYaxis()->SetRangeUser(0,0.3);
    beta_signal->SetTitle("beta");
    beta_signal->Draw();

    beta_secondary_neutron->Write();
    beta_secondary_neutron->SetLineColor(4);
    beta_secondary_neutron->SetStats(0);
    beta_secondary_neutron->Scale(1/beta_secondary_neutron->Integral(),"nosw2");
    beta_secondary_neutron->Draw("same");

    beta_primary_gamma->Write();
    beta_primary_gamma->SetLineColor(6);
    beta_primary_gamma->SetStats(0);
    beta_primary_gamma->Scale(1/beta_primary_gamma->Integral(),"nosw2");
    beta_primary_gamma->Draw("same");

    beta_secondary_gamma->Write();
    beta_secondary_gamma->SetLineColor(8);
    beta_secondary_gamma->SetStats(0);
    beta_secondary_gamma->Scale(1/beta_secondary_gamma->Integral(),"nosw2");
    beta_secondary_gamma->Draw("same");

    legend->Draw();
    can->SaveAs("beta.pdf");
    can->Clear();
    //}

    //distance{
    distance_signal->Write();
    distance_signal->SetLineColor(2);
    distance_signal->SetStats(0);
    distance_signal->Scale(1/distance_signal->Integral(),"nosw2");
    distance_signal->GetYaxis()->SetRangeUser(0,0.6);
    distance_signal->SetTitle("distance b/w C and hit");
    distance_signal->Draw();

    distance_secondary_neutron->Write();
    distance_secondary_neutron->SetLineColor(4);
    distance_secondary_neutron->SetStats(0);
    distance_secondary_neutron->Scale(1/distance_secondary_neutron->Integral(),"nosw2");
    distance_secondary_neutron->Draw("same");

    distance_primary_gamma->Write();
    distance_primary_gamma->SetLineColor(6);
    distance_primary_gamma->SetStats(0);
    distance_primary_gamma->Scale(1/distance_primary_gamma->Integral(),"nosw2");
    distance_primary_gamma->Draw("same");

    distance_secondary_gamma->Write();
    distance_secondary_gamma->SetLineColor(8);
    distance_secondary_gamma->SetStats(0);
    distance_secondary_gamma->Scale(1/distance_secondary_gamma->Integral(),"nosw2");
    distance_secondary_gamma->Draw("same");

    legend->Draw();
    can->SaveAs("distance.pdf");
    can->Clear();
    //}

    //TOF{
    TOF_signal->Write();
    TOF_signal->SetLineColor(2);
    TOF_signal->SetStats(0);
    TOF_signal->Scale(1/TOF_signal->Integral(),"nosw2");
    TOF_signal->GetYaxis()->SetRangeUser(0,0.7);
    TOF_signal->SetTitle("Time of flight");
    TOF_signal->Draw();

    TOF_secondary_neutron->Write();
    TOF_secondary_neutron->SetLineColor(4);
    TOF_secondary_neutron->SetStats(0);
    TOF_secondary_neutron->Scale(1/TOF_secondary_neutron->Integral(),"nosw2");
    TOF_secondary_neutron->Draw("same");

    TOF_primary_gamma->Write();
    TOF_primary_gamma->SetLineColor(6);
    TOF_primary_gamma->SetStats(0);
    TOF_primary_gamma->Scale(1/TOF_primary_gamma->Integral(),"nosw2");
    TOF_primary_gamma->Draw("same");

    TOF_secondary_gamma->Write();
    TOF_secondary_gamma->SetLineColor(8);
    TOF_secondary_gamma->SetStats(0);
    TOF_secondary_gamma->Scale(1/TOF_secondary_gamma->Integral(),"nosw2");
    TOF_secondary_gamma->Draw("same");

    legend->Draw();
    can->SaveAs("TOF.pdf");
    can->Clear();
    //}
    
    //CubeE{
    CubeE_signal->Write();
    CubeE_signal->SetLineColor(2);
    CubeE_signal->SetStats(0);
    CubeE_signal->Scale(1/CubeE_signal->Integral(),"nosw2");
    CubeE_signal->GetYaxis()->SetRangeUser(0,0.7);
    CubeE_signal->SetTitle("CubeE");
    CubeE_signal->Draw();

    CubeE_secondary_neutron->Write();
    CubeE_secondary_neutron->SetLineColor(4);
    CubeE_secondary_neutron->SetStats(0);
    CubeE_secondary_neutron->Scale(1/CubeE_secondary_neutron->Integral(),"nosw2");
    CubeE_secondary_neutron->Draw("same");

    CubeE_primary_gamma->Write();
    CubeE_primary_gamma->SetLineColor(6);
    CubeE_primary_gamma->SetStats(0);
    CubeE_primary_gamma->Scale(1/CubeE_primary_gamma->Integral(),"nosw2");
    CubeE_primary_gamma->Draw("same");

    CubeE_secondary_gamma->Write();
    CubeE_secondary_gamma->SetLineColor(8);
    CubeE_secondary_gamma->SetStats(0);
    CubeE_secondary_gamma->Scale(1/CubeE_secondary_gamma->Integral(),"nosw2");
    CubeE_secondary_gamma->Draw("same");

    legend->Draw();
    can->SaveAs("CubeE.pdf");
    can->Clear();
    //}

    //nCube{
    nCubeDis_signal->Write();
    nCubeDis_signal->GetYaxis()->SetRangeUser(0,0.5);
    nCubeDis_signal->SetLineColor(2);
    nCubeDis_signal->SetStats(0);
    nCubeDis_signal->Scale(1/nCubeDis_signal->Integral(),"nosw2");
    nCubeDis_signal->GetYaxis()->SetRangeUser(0,0.5);
    nCubeDis_signal->SetTitle("number of cube");
    nCubeDis_signal->Draw();

    nCubeDis_secondary_neutron->Write();
    nCubeDis_secondary_neutron->SetLineColor(4);
    nCubeDis_secondary_neutron->SetStats(0);
    nCubeDis_secondary_neutron->Scale(1/nCubeDis_secondary_neutron->Integral(),"nosw2");
    nCubeDis_secondary_neutron->Draw("same");

    nCubeDis_primary_gamma->Write();
    nCubeDis_primary_gamma->SetLineColor(6);
    nCubeDis_primary_gamma->SetStats(0);
    nCubeDis_primary_gamma->Scale(1/nCubeDis_primary_gamma->Integral(),"nosw2");
    nCubeDis_primary_gamma->Draw("same");

    nCubeDis_secondary_gamma->Write();
    nCubeDis_secondary_gamma->SetLineColor(8);
    nCubeDis_secondary_gamma->SetStats(0);
    nCubeDis_secondary_gamma->Scale(1/nCubeDis_secondary_gamma->Integral(),"nosw2");
    nCubeDis_secondary_gamma->Draw("same");

    legend->Draw();
    can->SaveAs("nCubeDis.pdf");
    can->Clear();
    //}

    fi1->Close();

    cout<<"making output files is done"<<endl;
    cout<<"---------------------------"<<endl;
    cout<<"all done"<<endl;

    return 0;
}
