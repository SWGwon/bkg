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
TH2F * hist_signal = new TH2F("hist_signal", "hist_signal;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_out3DST = new TH2F("hist_bkg_out3DST", "hist_bkg_out3DST;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_NC = new TH2F("hist_bkg_NC", "hist_bkg_NC;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_1 = new TH2F("hist_bkg_1", "hist_bkg_1;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_out3DST_largeTime = new TH2F("hist_bkg_out3DST_lt", "hist_bkg_out3DST;Lever Arm [cm]; Time [ns]", 20, 0, 200, 100, 0, 10000);
TH2F * hist_bkg_NC_largeTime = new TH2F("hist_bkg_NC_lt", "hist_bkg_NC;Lever Arm [cm]; Time [ns]", 20, 0, 200, 100, 0, 10000);
TH2F * hist_bkg_1_largeTime = new TH2F("hist_bkg_1_lt", "hist_bkg_1;Lever Arm [cm]; Time [ns]", 20, 0, 200, 100, 0, 10000);

TH3F * hist_neutron_hit = new TH3F("asdf","asdf",240,-120,120,240,-120,120,200,-100,100);

TH1F * neutronParentPDG = new TH1F("PDG","PDG",3500,-500,3000);
TH1F * neutronParentPDG_case4 = new TH1F("PDG_case4","PDG_case4",6,0,6);

TH1D * KE_secondary = new TH1D("seconday","seconday",100,0,200);
TH1D * KE_primary = new TH1D("primary","primary",100,0,200);

bool is_inFV = false;       //check if vertex is in FV
bool is_in3DST = false;     //check if vertex is in 3DST


struct Hit_t 
{
    float timeWindow,           // time windows of the hit
          timeSmear,        // smear time
          energyDeposit,        // energy deposited by the neutron
          trackLength,          // lever arm
          trueRec,      // true reconstructed energy
          smearRec,
          vtxSignal[3],     // neutrino vertex position of the neutron
          vtxTime,      // neutrino  vertex time
          neutronTrueE,    //neutron true energy
          neutronTrueT;    //neutron true time

    int bkgLoc,         // neutrino vertex position
        neutronParentId,    // Where the neutron come from
        neutronParentPdg;   // PDG of neutron parent

    bool isTherePion50,     // Is there a pion with KE > 50 MeV in FS particles
         isThereProton300;      // Is there a proton with KE > 300 MeV in FS particles
};

struct neutronPostion
{
    float X,
          Y,
          Z;
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

void analyze(string file)
{
    //cout<<file<<endl;     //cout file name
    auto _file = new TFile(TString(file));
    auto tree = (TTree*)_file->Get("tree");

    if(tree == NULL)
    {
        _file->Close();
        return;
    }

    float t_neutronHitX[1000], t_neutronHitY[1000], t_neutronHitZ[1000];
    float t_neutronHitT[1000], t_neutronParentId[1000], t_neutronParentPDG[1000];
    float t_neutronHitE[1000], t_neutronTrueE[1000];
    float t_vtx[3], t_vtxTime;

    int PDG = 0;
    int t_nFS, t_fsPdg[1000];

    tree->SetBranchAddress("neutronHitX", &t_neutronHitX);
    tree->SetBranchAddress("neutronHitY", &t_neutronHitY);
    tree->SetBranchAddress("neutronHitZ", &t_neutronHitZ);
    tree->SetBranchAddress("neutronHitT", &t_neutronHitT);
    tree->SetBranchAddress("neutronParentId", &t_neutronParentId);
    tree->SetBranchAddress("neutronParentPDG", &t_neutronParentPDG);
    tree->SetBranchAddress("neutronHitE", &t_neutronHitE);
    tree->SetBranchAddress("neutronTrueE", &t_neutronTrueE);
    tree->SetBranchAddress("vtx", &t_vtx);
    tree->SetBranchAddress("vtxTime", &t_vtxTime);
    tree->SetBranchAddress("nFS", &t_nFS);
    tree->SetBranchAddress("fsPdg", &t_fsPdg);


    //flag 
    bool is_Sig = false,
         is_Bkg_out3DST = false,
         is_Bkg_NC = false,
         is_Bkg_1 = false,
         is_Bkg_outFV_in3DST_primary = false,
         is_Bkg_outFV_in3DST_secondary = false;

    int nevents = tree->GetEntries();

    //global git variables
    Hit_t signal,
          bkg_out3DST,
          bkg_NC,
          bkg_1,
          bkg_outFV_in3DST_primary,
          bkg_outFV_in3DST_secondary;


    //cout<<"number of event: "<<nevents<<endl;

    for(int event = 0; event < nevents; event++)
    {
        tree->GetEntry(event);

        if(abs(t_vtx[0]) < 50 && abs(t_vtx[1]) < 50 && abs(t_vtx[2]) < 50)
        {
            is_Sig = true;

            //flag to CC
            bool is_CC = false;
            int number_of_neutron = 0;

            //search for a muon or electron/positron, t_nFS = number of FS particle, inFS is iterator
            for(int inFS = 0; inFS < t_nFS; inFS++)
            {
                if(abs(t_fsPdg[inFS]) == 11 || abs(t_fsPdg[inFS]) == 13)    //electronPDG=11,muonPDG=13
                {
                    is_CC = true;
                    break;
                }
            }

            //if there's no CC event, skip to the next event
            if(!is_CC)
                continue;

            map<string,Hit_t> hitPerCube;

            /*SIGNAL : Neutron Information
              Look for neutron hit induced by a CC event in the FV across all the cube in 3DST
              Then select the earliest activated cube as a signal
             */

            //n_neutronHit = iterator
            for(int n_neutronHit = 0; n_neutronHit < 1000; n_neutronHit++)
            {
                if(t_neutronHitX[n_neutronHit] != 0)
                {
                    //look for a neutron hit in 3DST
                    if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                            abs(t_neutronHitY[n_neutronHit]) < 120 && 
                            abs(t_neutronHitZ[n_neutronHit]) < 100)
                    {
                        hist_neutron_hit->Fill(t_neutronHitX[n_neutronHit],t_neutronHitY[n_neutronHit],t_neutronHitZ[n_neutronHit]);
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

                        number_of_neutron++;

                        if(signalWindow > 0)
                        {
                            Hit_t temp;

                            temp.timeWindow = signalWindow;
                            temp.trackLength = trackLength;
                            temp.energyDeposit = t_neutronHitE[n_neutronHit];

                            temp.vtxSignal[0] = t_vtx[0];
                            temp.vtxSignal[1] = t_vtx[1];
                            temp.vtxSignal[2] = t_vtx[2];

                            temp.neutronParentId = t_neutronParentId[n_neutronHit];
                            temp.neutronParentPdg = t_neutronParentPDG[n_neutronHit];

                            temp.vtxTime = t_vtxTime;

                            //Positon(cm)
                            string key = string(Form("%d_%d_%d",
                                        (int)t_neutronHitX[n_neutronHit],
                                        (int)t_neutronHitY[n_neutronHit],
                                        (int)t_neutronHitZ[n_neutronHit]));

                            /*
                               +If the cube has been already activated by a neutron
                               -See which neutron hit the cube first
                               -Affect the first neutron hit to the cube
                               -Sum up the energy Deposit in the cube
                               +Affect the neutron hit to the cube otherwise
                             */
                            auto findKey_hitCubeEvent = hitPerCube.find(key);
                            if(findKey_hitCubeEvent != hitPerCube.end())
                            {
                                if(hitPerCube.at(key).timeWindow < temp.timeWindow)
                                {
                                    hitPerCube.at(key).energyDeposit += temp.energyDeposit;
                                }
                                else
                                {
                                    auto tempEnergy = hitPerCube.at(key).energyDeposit;
                                    hitPerCube.at(key) = temp;
                                    hitPerCube.at(key).energyDeposit += tempEnergy;
                                }
                            }
                            else
                            {
                                hitPerCube[key] = temp;     
                            }
                        }
                    }
                }
            }   //end of n_neutronhit iterate

            //cout<<"event:"<<event+1<<", number of neutron hit :"<<number_of_neutron<<endl;

            Hit_t temp_sig;
            temp_sig.timeWindow = 100000000;

            for(auto hit : hitPerCube)
            {
                if(temp_sig.timeWindow > hit.second.timeWindow && hit.second.energyDeposit > energyHitCut && (hit.second.neutronParentId == -1 || hit.second.neutronParentId == 0))
                {
                    temp_sig = hit.second;
                }
            }


            if(is_Sig  && temp_sig.timeWindow != 100000000)
            {
                //cout<<"there is signal"<<endl;
                signal = temp_sig;
                hist_signal->Fill(signal.trackLength,signal.timeWindow);
                //cout<<"arm :"<<signal.trackLength<<", time: "<<signal.timeWindow<<endl;
                //cout<<"vertex: "<<signal.vtxSignal[0]<<","<<signal.vtxSignal[1]<<","<<signal.vtxSignal[2]<<endl;
            }
        }
    }       //end of event iterate

    if(is_Sig)
    {
    }

    if(!is_Sig)
    {
        //cout<<"there is no signal"<<endl;
        _file->Close();
        return;
    }


    /*
BACKGROUND : Neutron information
+Search for neutron in 3DST 
1.vertex is out of 3DST
2.vertex is inFV && parentid != -1
3.vertex is outFV_in3DST && NC
*/
    for(int event = 0; event < nevents; event++)
    {
        tree->GetEntry(event);

        bool out3DST = false,
             outFV_in3DST = false;

        if(abs(t_vtx[0]) > 120 || abs(t_vtx[1]) > 120 || abs(t_vtx[2]) > 100)
        {
            out3DST = true;
        }

        if(abs(t_vtx[0]) < 120 && abs(t_vtx[1]) < 120 && abs(t_vtx[2]) < 100)
        {
            if(abs(t_vtx[0]) > 50 || abs(t_vtx[1]) > 50 || abs(t_vtx[2]) > 50)
            {
                outFV_in3DST = true;
            }
        }

        bool NC = false,
             CC = false;

        //look for NC event
        if(outFV_in3DST)
        {
            for(int inFS = 0; inFS < t_nFS; inFS++)
            {
                //look for muon or electron
                if(abs(t_fsPdg[inFS]) == 11 || abs(t_fsPdg[inFS]) == 13)
                {
                    CC = true;
                    break;
                }
            }

            //no muon or electron found
            if(!CC)
            {
                //look for pion
                for(int inFS = 0; inFS < t_nFS; inFS++)
                {
                    if(abs(t_fsPdg[inFS]) == 211 || abs(t_fsPdg[inFS]) == 111)
                    {
                        NC = true;
                        break;
                    }
                }
            }
        }   //end of NC found if

        //1.background from out of 3DST
        if(out3DST)
        {
            map<string,Hit_t> hitPerCube;

            for(int n_neutronHit = 0; n_neutronHit < 1000; n_neutronHit++)
            {
                if(t_neutronHitX[n_neutronHit] != 0)
                {
                    //look for a neutron hit in 3DST
                    if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                            abs(t_neutronHitY[n_neutronHit]) < 120 && 
                            abs(t_neutronHitZ[n_neutronHit]) < 100)
                    {
                        //calculate distance from FV vertex
                        float trackLength = pow(
                                pow(t_neutronHitX[n_neutronHit] - signal.vtxSignal[0],2)+
                                pow(t_neutronHitY[n_neutronHit] - signal.vtxSignal[1],2)+
                                pow(t_neutronHitZ[n_neutronHit] - signal.vtxSignal[2],2),0.5);

                        //calculate signal window; 
                        float backgroundWindow = t_neutronHitT[n_neutronHit] - signal.vtxTime;

                        //Fix a bug from edep-sim
                        if(backgroundWindow == 1)
                            backgroundWindow = 0.5;


                        if(backgroundWindow > 0)
                        {
                            Hit_t temp;

                            temp.timeWindow = backgroundWindow;
                            temp.trackLength = trackLength;
                            temp.energyDeposit = t_neutronHitE[n_neutronHit];

                            temp.vtxSignal[0] = t_vtx[0];
                            temp.vtxSignal[1] = t_vtx[1];
                            temp.vtxSignal[2] = t_vtx[2];

                            //temp.bkgLoc = 

                            temp.neutronParentId = t_neutronParentId[n_neutronHit];
                            temp.neutronParentPdg = t_neutronParentPDG[n_neutronHit];

                            temp.vtxTime = t_vtxTime;

                            //Positon(cm)
                            string key = string(Form("%d_%d_%d",
                                        (int)t_neutronHitX[n_neutronHit],
                                        (int)t_neutronHitY[n_neutronHit],
                                        (int)t_neutronHitZ[n_neutronHit]));

                            /*
                               +If the cube has been already activated by a neutron
                               -See which neutron hit the cube first
                               -Affect the first neutron hit to the cube
                               -Sum up the energy Deposit in the cube
                               +Affect the neutron hit to the cube otherwise
                             */
                            auto findKey_hitCubeEvent = hitPerCube.find(key);
                            if(findKey_hitCubeEvent != hitPerCube.end())
                            {
                                if(hitPerCube.at(key).timeWindow < temp.timeWindow)
                                {
                                    hitPerCube.at(key).energyDeposit += temp.energyDeposit;
                                }
                                else
                                {
                                    auto tempEnergy = hitPerCube.at(key).energyDeposit;
                                    hitPerCube.at(key) = temp;
                                    hitPerCube.at(key).energyDeposit += tempEnergy;
                                }
                            }
                            else
                            {
                                hitPerCube[key] = temp;
                            }
                        }
                    }
                }
            }   //end of n_neutronhit iterate

            Hit_t bkg_temp_out3DST;
            bkg_temp_out3DST.timeWindow = 10000000;

            //look for the earliest hit
            for(auto hit : hitPerCube)
            {
                if(bkg_temp_out3DST.timeWindow > hit.second.timeWindow && hit.second.energyDeposit > energyHitCut)
                {
                    bkg_temp_out3DST = hit.second;
                    is_Bkg_out3DST = true;
                }
            }

            if(is_Bkg_out3DST && bkg_temp_out3DST.timeWindow != 10000000)
            {
                bkg_out3DST = bkg_temp_out3DST;
                hist_bkg_out3DST->Fill(bkg_out3DST.trackLength,bkg_out3DST.timeWindow);
                hist_bkg_out3DST_largeTime->Fill(bkg_out3DST.trackLength,bkg_out3DST.timeWindow);
            }
        }       //end of if(out3DST)

        //2.parentid != -1
        if(abs(t_vtx[0]) < 50 && abs(t_vtx[1]) < 50 && abs(t_vtx[2]) < 50) 
        {
            map<string,Hit_t> hitPerCube;

            for(int n_neutronHit = 0; n_neutronHit < 1000; n_neutronHit++)
            {
                if(t_neutronHitX[n_neutronHit] != 0)
                {
                    //look for a neutron hit in 3DST
                    if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                            abs(t_neutronHitY[n_neutronHit]) < 120 && 
                            abs(t_neutronHitZ[n_neutronHit]) < 100)
                    {
                        //calculate distance from FV vertex
                        float trackLength = pow(
                                pow(t_neutronHitX[n_neutronHit] - signal.vtxSignal[0],2)+
                                pow(t_neutronHitY[n_neutronHit] - signal.vtxSignal[1],2)+
                                pow(t_neutronHitZ[n_neutronHit] - signal.vtxSignal[2],2),0.5);

                        //calculate signal window; time of flight
                        float backgroundWindow = t_neutronHitT[n_neutronHit] - signal.vtxTime;

                        //Fix a bug from edep-sim
                        if(backgroundWindow == 1)
                            backgroundWindow = 0.5;


                        if(backgroundWindow > 0)
                        {
                            Hit_t temp;

                            temp.timeWindow = backgroundWindow;
                            temp.trackLength = trackLength;
                            temp.energyDeposit = t_neutronHitE[n_neutronHit];

                            temp.vtxSignal[0] = t_vtx[0];
                            temp.vtxSignal[1] = t_vtx[1];
                            temp.vtxSignal[2] = t_vtx[2];

                            //temp.bkgLoc = 

                            temp.neutronParentId = t_neutronParentId[n_neutronHit];
                            temp.neutronParentPdg = t_neutronParentPDG[n_neutronHit];

                            temp.vtxTime = t_vtxTime;

                            //Positon(cm)
                            string key = string(Form("%d_%d_%d",
                                        (int)t_neutronHitX[n_neutronHit],
                                        (int)t_neutronHitY[n_neutronHit],
                                        (int)t_neutronHitZ[n_neutronHit]));

                            /*
                               +If the cube has been already activated by a neutron
                               -See which neutron hit the cube first
                               -Affect the first neutron hit to the cube
                               -Sum up the energy Deposit in the cube
                               +Affect the neutron hit to the cube otherwise
                             */
                            auto findKey_hitCubeEvent = hitPerCube.find(key);
                            if(findKey_hitCubeEvent != hitPerCube.end())
                            {
                                if(hitPerCube.at(key).timeWindow < temp.timeWindow)
                                {
                                    hitPerCube.at(key).energyDeposit += temp.energyDeposit;
                                }
                                else
                                {
                                    auto tempEnergy = hitPerCube.at(key).energyDeposit;
                                    hitPerCube.at(key) = temp;
                                    hitPerCube.at(key).energyDeposit += tempEnergy;
                                }
                            }
                            else
                            {
                                hitPerCube[key] = temp;
                            }
                        }
                    }
                }
            }   //end of n_neutronhit iterate

            Hit_t bkg_temp_1;
            bkg_temp_1.timeWindow = 100000000;

            //look for the earliest hit and parentid != -1
            for(auto hit : hitPerCube)
            {
                if(bkg_temp_1.timeWindow > hit.second.timeWindow && hit.second.energyDeposit > energyHitCut && hit.second.neutronParentId > 0)
                {
                    bkg_temp_1 = hit.second;
                    is_Bkg_1 = true;
                }
            }
            
            if(is_Bkg_1 && bkg_temp_1.timeWindow != 100000000)
            {
                bkg_1 = bkg_temp_1;
                hist_bkg_1->Fill(bkg_1.trackLength,bkg_1.timeWindow);
            }
        }       //end of if(abs(t_vtx[0]) < 50 && abs(t_vtx[1]) < 50 && abs(t_vtx[2]) < 50) 

        //3.for outFV_in3DST, NC
        if(outFV_in3DST && !CC)
        {
            map<string,Hit_t> hitPerCube;

            for(int n_neutronHit = 0; n_neutronHit < 1000; n_neutronHit++)
            {
                if(t_neutronHitX[n_neutronHit] != 0)
                {
                    //look for a neutron hit in 3DST
                    if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                            abs(t_neutronHitY[n_neutronHit]) < 120 && 
                            abs(t_neutronHitZ[n_neutronHit]) < 100)
                    {
                        //calculate distance from FV vertex
                        float trackLength = pow(
                                pow(t_neutronHitX[n_neutronHit] - signal.vtxSignal[0],2)+
                                pow(t_neutronHitY[n_neutronHit] - signal.vtxSignal[1],2)+
                                pow(t_neutronHitZ[n_neutronHit] - signal.vtxSignal[2],2),0.5);

                        //calculate signal window; time of flight
                        float backgroundWindow = t_neutronHitT[n_neutronHit] - signal.vtxTime;

                        //Fix a bug from edep-sim
                        if(backgroundWindow == 1)
                            backgroundWindow = 0.5;


                        if(backgroundWindow > 0)
                        {
                            Hit_t temp;

                            temp.timeWindow = backgroundWindow;
                            temp.trackLength = trackLength;
                            temp.energyDeposit = t_neutronHitE[n_neutronHit];

                            temp.vtxSignal[0] = t_vtx[0];
                            temp.vtxSignal[1] = t_vtx[1];
                            temp.vtxSignal[2] = t_vtx[2];

                            //temp.bkgLoc = 

                            temp.neutronParentId = t_neutronParentId[n_neutronHit];
                            temp.neutronParentPdg = t_neutronParentPDG[n_neutronHit];

                            temp.vtxTime = t_vtxTime;

                            //Positon(cm)
                            string key = string(Form("%d_%d_%d",
                                        (int)t_neutronHitX[n_neutronHit],
                                        (int)t_neutronHitY[n_neutronHit],
                                        (int)t_neutronHitZ[n_neutronHit]));

                            /*
                               +If the cube has been already activated by a neutron
                               -See which neutron hit the cube first
                               -Affect the first neutron hit to the cube
                               -Sum up the energy Deposit in the cube
                               +Affect the neutron hit to the cube otherwise
                             */
                            auto findKey_hitCubeEvent = hitPerCube.find(key);
                            if(findKey_hitCubeEvent != hitPerCube.end())
                            {
                                if(hitPerCube.at(key).timeWindow < temp.timeWindow)
                                {
                                    hitPerCube.at(key).energyDeposit += temp.energyDeposit;
                                }
                                else
                                {
                                    auto tempEnergy = hitPerCube.at(key).energyDeposit;
                                    hitPerCube.at(key) = temp;
                                    hitPerCube.at(key).energyDeposit += tempEnergy;
                                }
                            }
                            else
                            {
                                hitPerCube[key] = temp;
                            }
                        }
                    }
                }
            }   //end of n_neutronhit iterate

            Hit_t bkg_temp_NC;
            bkg_temp_NC.timeWindow = 100000000;

            //look for the earliest hit
            for(auto hit : hitPerCube)
            {
                if(bkg_temp_NC.timeWindow > hit.second.timeWindow && hit.second.energyDeposit > energyHitCut)
                {
                    bkg_temp_NC = hit.second;
                    is_Bkg_NC = true;
                }
            }

            if(is_Bkg_NC && bkg_temp_NC.timeWindow != 100000000)
            {
                bkg_NC = bkg_temp_NC;
                hist_bkg_NC->Fill(bkg_NC.trackLength,bkg_NC.timeWindow);
                hist_bkg_NC_largeTime->Fill(bkg_NC.trackLength,bkg_NC.timeWindow);
                //cout<<"arm: "<<bkg_NC.trackLength<<",time :"<<bkg_NC.timeWindow<<endl;
            }
        }       //end of if(NC)

        if(outFV_in3DST)
        {
            map<string,Hit_t> hitPerCube;

            for(int n_neutronHit = 0; n_neutronHit < 1000; n_neutronHit++)
            {
                if(t_neutronHitX[n_neutronHit] != 0)
                {
                    //look for a neutron hit in 3DST
                    if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                            abs(t_neutronHitY[n_neutronHit]) < 120 && 
                            abs(t_neutronHitZ[n_neutronHit]) < 100)
                    {
                        //calculate distance from FV vertex
                        float trackLength = pow(
                                pow(t_neutronHitX[n_neutronHit] - signal.vtxSignal[0],2)+
                                pow(t_neutronHitY[n_neutronHit] - signal.vtxSignal[1],2)+
                                pow(t_neutronHitZ[n_neutronHit] - signal.vtxSignal[2],2),0.5);

                        //calculate signal window; time of flight
                        float backgroundWindow = t_neutronHitT[n_neutronHit] - signal.vtxTime;

                        //Fix a bug from edep-sim
                        if(backgroundWindow == 1)
                            backgroundWindow = 0.5;


                        if(backgroundWindow > 0)
                        {
                            Hit_t temp;

                            temp.timeWindow = backgroundWindow;
                            temp.trackLength = trackLength;
                            temp.energyDeposit = t_neutronHitE[n_neutronHit];

                            temp.vtxSignal[0] = t_vtx[0];
                            temp.vtxSignal[1] = t_vtx[1];
                            temp.vtxSignal[2] = t_vtx[2];

                            //temp.bkgLoc = 

                            temp.neutronParentId = t_neutronParentId[n_neutronHit];
                            temp.neutronParentPdg = t_neutronParentPDG[n_neutronHit];

                            temp.vtxTime = t_vtxTime;

                            //Positon(cm)
                            string key = string(Form("%d_%d_%d",
                                        (int)t_neutronHitX[n_neutronHit],
                                        (int)t_neutronHitY[n_neutronHit],
                                        (int)t_neutronHitZ[n_neutronHit]));

                            /*
                               +If the cube has been already activated by a neutron
                               -See which neutron hit the cube first
                               -Affect the first neutron hit to the cube
                               -Sum up the energy Deposit in the cube
                               +Affect the neutron hit to the cube otherwise
                             */
                            auto findKey_hitCubeEvent = hitPerCube.find(key);
                            if(findKey_hitCubeEvent != hitPerCube.end())
                            {
                                if(hitPerCube.at(key).timeWindow < temp.timeWindow)
                                {
                                    hitPerCube.at(key).energyDeposit += temp.energyDeposit;
                                }
                                else
                                {
                                    auto tempEnergy = hitPerCube.at(key).energyDeposit;
                                    hitPerCube.at(key) = temp;
                                    hitPerCube.at(key).energyDeposit += tempEnergy;
                                }
                            }
                            else
                            {
                                hitPerCube[key] = temp;
                            }
                        }
                    }
                }
            }   //end of n_neutronhit iterate

            Hit_t bkg_temp_outFV_in3DST_secondary;
            bkg_temp_outFV_in3DST_secondary.timeWindow = 100000000;

            //look for the earliest hit
            for(auto hit : hitPerCube)
            {
                if(bkg_temp_outFV_in3DST_secondary.timeWindow > hit.second.timeWindow && hit.second.energyDeposit > energyHitCut && hit.second.neutronParentId > -1)
                {
                    bkg_temp_outFV_in3DST_secondary = hit.second;
                    is_Bkg_outFV_in3DST_secondary = true;
                }
            }

            if(is_Bkg_outFV_in3DST_secondary && bkg_temp_outFV_in3DST_secondary.timeWindow != 100000000)
            {
                bkg_outFV_in3DST_secondary = bkg_temp_outFV_in3DST_secondary;
                //cout<<"secondary KE:"<<kineticEnergy(bkg_earliestHit_outFV_in3DST_secondary.trackLength,bkg_earliestHit_outFV_in3DST_secondary.timeWindow)<<endl;
                KE_secondary->Fill(kineticEnergy(bkg_outFV_in3DST_secondary.trackLength,bkg_outFV_in3DST_secondary.timeWindow));
                neutronParentPDG->Fill(bkg_outFV_in3DST_secondary.neutronParentPdg);
                if(abs(bkg_outFV_in3DST_secondary.neutronParentPdg) == 211 || bkg_outFV_in3DST_secondary.neutronParentPdg == 111)
                    neutronParentPDG_case4->Fill(1);        //PDG = +-211,111
                if(bkg_outFV_in3DST_secondary.neutronParentPdg ==2112)
                    neutronParentPDG_case4->Fill(2);        //PDG = 2112
                if(bkg_outFV_in3DST_secondary.neutronParentPdg == 2212)
                    neutronParentPDG_case4->Fill(3);        //PDG = 2212
                if(abs(bkg_outFV_in3DST_secondary.neutronParentPdg) != 211 &&
                        bkg_outFV_in3DST_secondary.neutronParentPdg != 111 &&
                        bkg_outFV_in3DST_secondary.neutronParentPdg != 2112 &&
                        bkg_outFV_in3DST_secondary.neutronParentPdg != 2212 &&
                        bkg_outFV_in3DST_secondary.neutronParentPdg != 0)
                    neutronParentPDG_case4->Fill(4);        //PDG = others axcept 0
            }
        }   //end of if(outFV_in3DST)

        if(outFV_in3DST)
        {
            map<string,Hit_t> hitPerCube;

            for(int n_neutronHit = 0; n_neutronHit < 1000; n_neutronHit++)
            {
                if(t_neutronHitX[n_neutronHit] != 0)
                {
                    //look for a neutron hit in 3DST
                    if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                            abs(t_neutronHitY[n_neutronHit]) < 120 && 
                            abs(t_neutronHitZ[n_neutronHit]) < 100)
                    {
                        //calculate distance from FV vertex
                        float trackLength = pow(
                                pow(t_neutronHitX[n_neutronHit] - signal.vtxSignal[0],2)+
                                pow(t_neutronHitY[n_neutronHit] - signal.vtxSignal[1],2)+
                                pow(t_neutronHitZ[n_neutronHit] - signal.vtxSignal[2],2),0.5);

                        //calculate signal window; time of flight
                        float backgroundWindow = t_neutronHitT[n_neutronHit] - signal.vtxTime;

                        //Fix a bug from edep-sim
                        if(backgroundWindow == 1)
                            backgroundWindow = 0.5;


                        if(backgroundWindow > 0)
                        {
                            Hit_t temp;

                            temp.timeWindow = backgroundWindow;
                            temp.trackLength = trackLength;
                            temp.energyDeposit = t_neutronHitE[n_neutronHit];

                            temp.vtxSignal[0] = t_vtx[0];
                            temp.vtxSignal[1] = t_vtx[1];
                            temp.vtxSignal[2] = t_vtx[2];

                            //temp.bkgLoc = 

                            temp.neutronParentId = t_neutronParentId[n_neutronHit];
                            temp.neutronParentPdg = t_neutronParentPDG[n_neutronHit];

                            temp.vtxTime = t_vtxTime;

                            //Positon(cm)
                            string key = string(Form("%d_%d_%d",
                                        (int)t_neutronHitX[n_neutronHit],
                                        (int)t_neutronHitY[n_neutronHit],
                                        (int)t_neutronHitZ[n_neutronHit]));

                            /*
                               +If the cube has been already activated by a neutron
                               -See which neutron hit the cube first
                               -Affect the first neutron hit to the cube
                               -Sum up the energy Deposit in the cube
                               +Affect the neutron hit to the cube otherwise
                             */
                            auto findKey_hitCubeEvent = hitPerCube.find(key);
                            if(findKey_hitCubeEvent != hitPerCube.end())
                            {
                                if(hitPerCube.at(key).timeWindow < temp.timeWindow)
                                {
                                    hitPerCube.at(key).energyDeposit += temp.energyDeposit;
                                }
                                else
                                {
                                    auto tempEnergy = hitPerCube.at(key).energyDeposit;
                                    hitPerCube.at(key) = temp;
                                    hitPerCube.at(key).energyDeposit += tempEnergy;
                                }
                            }
                            else
                            {
                                hitPerCube[key] = temp;
                            }
                        }
                    }
                }
            }   //end of n_neutronhit iterate

            Hit_t      bkg_temp_outFV_in3DST_primary;
                       bkg_temp_outFV_in3DST_primary.timeWindow = 100000000;

            for(auto hit : hitPerCube)
            {
                if(bkg_temp_outFV_in3DST_primary.timeWindow > hit.second.timeWindow && hit.second.energyDeposit > energyHitCut && hit.second.neutronParentId == -1)
                {
                    bkg_temp_outFV_in3DST_primary = hit.second;
                    is_Bkg_outFV_in3DST_primary = true;
                }
            }

            if(is_Bkg_outFV_in3DST_primary && bkg_temp_outFV_in3DST_primary.timeWindow != 100000000)
            {
                //cout<<"primary KE:"<<kineticEnergy(bkg_earliestHit_outFV_in3DST_primary.trackLength,bkg_earliestHit_outFV_in3DST_primary.timeWindow)<<endl;
                bkg_outFV_in3DST_primary = bkg_temp_outFV_in3DST_primary;
                KE_primary->Fill(kineticEnergy(bkg_outFV_in3DST_primary.trackLength,bkg_outFV_in3DST_primary.timeWindow));
            }

        }       //end of if(outFV_in3DST)

        /*
        if(1)
        {
            map<string,Hit_t> hitPerCube;

            int number_of_neutron = 0;

            for(int n_neutronHit = 0; n_neutronHit < 1000; n_neutronHit++)
            {
                //look for a neutron hit in 3DST
                if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                        abs(t_neutronHitY[n_neutronHit]) < 120 && 
                        abs(t_neutronHitZ[n_neutronHit]) < 100)
                {
                    //calculate distance from FV vertex
                    float trackLength = pow(
                            pow(t_neutronHitX[n_neutronHit] - signal.vtxSignal[0],2)+
                            pow(t_neutronHitY[n_neutronHit] - signal.vtxSignal[1],2)+
                            pow(t_neutronHitZ[n_neutronHit] - signal.vtxSignal[2],2),0.5);

                    //calculate signal window; 
                    float backgroundWindow = t_neutronHitT[n_neutronHit] - signal.vtxTime;

                    //Fix a bug from edep-sim
                    if(backgroundWindow == 1)
                        backgroundWindow = 0.5;

                    if(t_neutronHitX[n_neutronHit] != 0)
                        number_of_neutron++;
                }
            }
            cout<<"event: "<<event+1<<", number of neutron :"<<number_of_neutron<<endl;
        }       //end of if(out3DST)
        */
    }       //end of for(int event = 0; event < nevents; event++)
    _file->Close();
}

void neutron(string filename)
{
    ifstream input(filename);

    auto outName = filename.substr(0,4);

    string file;

    while(input >> file)
    {
        analyze(file);
    }

    TFile * fi1 = new TFile("background.root","RECREATE");
    hist_signal->Write();
    hist_bkg_out3DST->Write();
    hist_bkg_NC->Write();
    hist_bkg_1->Write();
    KE_primary->Write();
    KE_secondary->Write();
    neutronParentPDG->Write();
    neutronParentPDG_case4->Write();
    fi1->Close();

    /*
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

    TCanvas * can2 = new  TCanvas;
    KE_primary->Draw();
    TCanvas * can3 = new  TCanvas;
    KE_secondary->Draw();
    */
}
