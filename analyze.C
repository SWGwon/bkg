//histograms{
TH2F * hist_signal = new TH2F("hist_signal", "hist_signal;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_out3DST = new TH2F("hist_bkg_out3DST", "hist_bkg_out3DST;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
TH2F * hist_bkg_NC = new TH2F("hist_bkg_NC", "hist_bkg_NC;Lever Arm [cm]; Time [ns]", 20, 0, 2000, 25, 0, 250);
TH2F * hist_bkg_1 = new TH2F("hist_bkg_1", "hist_bkg_1;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);

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
    //TFile * fi = new TFile("/Users/gwon/FHC_1.root");
    //cout<<file<<endl;     //cout file name
    auto fi = new TFile(TString(file));
    if(!fi->GetListOfKeys()->Contains("tree"))
        return;
    auto tree = (TTree*)fi->Get("tree");

    if(tree == NULL)
    {
        fi->Close();
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
         is_Bkg = false;

    int nevents = tree->GetEntries();

    //cout<<"number of event: "<<nevents<<endl;

    for(int event = 0; event < nevents; event++)
    {
        tree->GetEntry(event);

        if(abs(t_vtx[0]) < 50 && abs(t_vtx[1]) < 50 && abs(t_vtx[2]) < 50)
        {
            //flag to CC
            bool is_CC = false;

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
                //look for a neutron hit in 3DST
                if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                        abs(t_neutronHitY[n_neutronHit]) < 120 && 
                        abs(t_neutronHitZ[n_neutronHit]) < 100)
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
                            hitPerCube[key] = temp;      //problem
                        }


                        is_Sig = true;
                    }
                }
            }   //end of n_neutronhit iterate
            Hit_t sig_earliestHit;
            sig_earliestHit.timeWindow = 1000;

            for(auto hit : hitPerCube)
            {
                if(sig_earliestHit.timeWindow > hit.second.timeWindow && hit.second.energyDeposit > energyHitCut)
                {
                    sig_earliestHit = hit.second;
                }
            }

            if(is_Sig && sig_earliestHit.timeWindow != 1000)
            {
                //cout<<"there is signal"<<endl;
                //doing signal stuff using sig_earliestHit
                hist_signal->Fill(sig_earliestHit.trackLength,sig_earliestHit.timeWindow);
            }
        }
    }       //end of event iterate

    if(!is_Sig)
    {
        //cout<<"there is no signal"<<endl;
        fi->Close();
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

            for(int n_neutronHit = 0; n_neutronHit < 100; n_neutronHit++)
            {
                //look for a neutron hit in 3DST
                if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                        abs(t_neutronHitY[n_neutronHit]) < 120 && 
                        abs(t_neutronHitZ[n_neutronHit]) < 100)
                {
                    //calculate distance from FV vertex
                    float trackLength = pow(
                            pow(t_neutronHitX[n_neutronHit] - t_vtx[0],2)+
                            pow(t_neutronHitY[n_neutronHit] - t_vtx[1],2)+
                            pow(t_neutronHitZ[n_neutronHit] - t_vtx[2],2),0.5);

                    //calculate signal window; time of flight
                    float backgroundWindow = t_neutronHitT[n_neutronHit] - t_vtxTime;

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
            }   //end of n_neutronhit iterate

            Hit_t bkg_earliestHit_out3DST;
            bkg_earliestHit_out3DST.timeWindow = 1000;

            //look for the earliest hit
            for(auto hit : hitPerCube)
            {
                if(bkg_earliestHit_out3DST.timeWindow > hit.second.timeWindow && hit.second.energyDeposit > energyHitCut)
                {
                    bkg_earliestHit_out3DST = hit.second;
                }
            }

            if(bkg_earliestHit_out3DST.timeWindow != 1000)
            {
                hist_bkg_out3DST->Fill(bkg_earliestHit_out3DST.trackLength,bkg_earliestHit_out3DST.timeWindow);
            }
        }       //end of if(out3DST)

        //2.parentid != -1
        if(abs(t_vtx[0]) < 50 && abs(t_vtx[1]) < 50 && abs(t_vtx[2]) < 50) 
        {
            map<string,Hit_t> hitPerCube;

            for(int n_neutronHit = 0; n_neutronHit < 100; n_neutronHit++)
            {
                //look for a neutron hit in 3DST
                if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                        abs(t_neutronHitY[n_neutronHit]) < 120 && 
                        abs(t_neutronHitZ[n_neutronHit]) < 100)
                {
                    //calculate distance from FV vertex
                    float trackLength = pow(
                            pow(t_neutronHitX[n_neutronHit] - t_vtx[0],2)+
                            pow(t_neutronHitY[n_neutronHit] - t_vtx[1],2)+
                            pow(t_neutronHitZ[n_neutronHit] - t_vtx[2],2),0.5);

                    //calculate signal window; time of flight
                    float backgroundWindow = t_neutronHitT[n_neutronHit] - t_vtxTime;

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
            }   //end of n_neutronhit iterate

            Hit_t bkg_earliestHit_1;
            bkg_earliestHit_1.timeWindow = 1000;

            //look for the earliest hit and parentid != -1
            for(auto hit : hitPerCube)
            {
                if(bkg_earliestHit_1.timeWindow > hit.second.timeWindow && hit.second.energyDeposit > energyHitCut && hit.second.neutronParentId != -1)
                {
                    bkg_earliestHit_1 = hit.second;
                }
            }
            
            if(bkg_earliestHit_1.timeWindow != 1000)
            {
                hist_bkg_1->Fill(bkg_earliestHit_1.trackLength,bkg_earliestHit_1.timeWindow);
            }
        }       //end of if(abs(t_vtx[0]) < 50 && abs(t_vtx[1]) < 50 && abs(t_vtx[2]) < 50) 

        //3.for outFV_in3DST, NC
        if(NC)
        {
            map<string,Hit_t> hitPerCube;

            for(int n_neutronHit = 0; n_neutronHit < 100; n_neutronHit++)
            {
                //look for a neutron hit in 3DST
                if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                        abs(t_neutronHitY[n_neutronHit]) < 120 && 
                        abs(t_neutronHitZ[n_neutronHit]) < 100)
                {
                    //calculate distance from FV vertex
                    float trackLength = pow(
                            pow(t_neutronHitX[n_neutronHit] - t_vtx[0],2)+
                            pow(t_neutronHitY[n_neutronHit] - t_vtx[1],2)+
                            pow(t_neutronHitZ[n_neutronHit] - t_vtx[2],2),0.5);

                    //calculate signal window; time of flight
                    float backgroundWindow = t_neutronHitT[n_neutronHit] - t_vtxTime;

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
            }   //end of n_neutronhit iterate

            Hit_t bkg_earliestHit_NC;
            bkg_earliestHit_NC.timeWindow = 1000;

            //look for the earliest hit
            for(auto hit : hitPerCube)
            {
                if(bkg_earliestHit_NC.timeWindow > hit.second.timeWindow && hit.second.energyDeposit > energyHitCut)
                {
                    bkg_earliestHit_NC = hit.second;
                }
            }

            if(bkg_earliestHit_NC.timeWindow != 1000)
            {
                hist_bkg_NC->Fill(bkg_earliestHit_NC.trackLength,bkg_earliestHit_NC.timeWindow);
            }
        }       //end of if(NC)

        if(outFV_in3DST)
        {
            map<string,Hit_t> hitPerCube;

            for(int n_neutronHit = 0; n_neutronHit < 100; n_neutronHit++)
            {
                //look for a neutron hit in 3DST
                if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                        abs(t_neutronHitY[n_neutronHit]) < 120 && 
                        abs(t_neutronHitZ[n_neutronHit]) < 100)
                {
                    //calculate distance from FV vertex
                    float trackLength = pow(
                            pow(t_neutronHitX[n_neutronHit] - t_vtx[0],2)+
                            pow(t_neutronHitY[n_neutronHit] - t_vtx[1],2)+
                            pow(t_neutronHitZ[n_neutronHit] - t_vtx[2],2),0.5);

                    //calculate signal window; time of flight
                    float backgroundWindow = t_neutronHitT[n_neutronHit] - t_vtxTime;

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
            }   //end of n_neutronhit iterate

            Hit_t bkg_earliestHit_outFV_in3DST_secondary;
            bkg_earliestHit_outFV_in3DST_secondary.timeWindow = 1000;

            //look for the earliest hit
            for(auto hit : hitPerCube)
            {
                if(bkg_earliestHit_outFV_in3DST_secondary.timeWindow > hit.second.timeWindow && hit.second.energyDeposit > energyHitCut && hit.second.neutronParentId != -1)
                {
                    bkg_earliestHit_outFV_in3DST_secondary = hit.second;
                }
            }

            if(bkg_earliestHit_outFV_in3DST_secondary.timeWindow != 1000)
            {
                cout<<"secondary KE:"<<kineticEnergy(bkg_earliestHit_outFV_in3DST_secondary.trackLength,bkg_earliestHit_outFV_in3DST_secondary.timeWindow)<<endl;
                KE_secondary->Fill(kineticEnergy(bkg_earliestHit_outFV_in3DST_secondary.trackLength,bkg_earliestHit_outFV_in3DST_secondary.timeWindow));
            }
        }

        if(outFV_in3DST)
        {
            map<string,Hit_t> hitPerCube;

            for(int n_neutronHit = 0; n_neutronHit < 100; n_neutronHit++)
            {
                //look for a neutron hit in 3DST
                if(abs(t_neutronHitX[n_neutronHit]) < 120 && 
                        abs(t_neutronHitY[n_neutronHit]) < 120 && 
                        abs(t_neutronHitZ[n_neutronHit]) < 100)
                {
                    //calculate distance from FV vertex
                    float trackLength = pow(
                            pow(t_neutronHitX[n_neutronHit] - t_vtx[0],2)+
                            pow(t_neutronHitY[n_neutronHit] - t_vtx[1],2)+
                            pow(t_neutronHitZ[n_neutronHit] - t_vtx[2],2),0.5);

                    //calculate signal window; time of flight
                    float backgroundWindow = t_neutronHitT[n_neutronHit] - t_vtxTime;

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
            }   //end of n_neutronhit iterate

            Hit_t      bkg_earliestHit_outFV_in3DST_primary;
                       bkg_earliestHit_outFV_in3DST_primary.timeWindow = 1000;

            for(auto hit : hitPerCube)
            {
                if(bkg_earliestHit_outFV_in3DST_primary.timeWindow > hit.second.timeWindow && hit.second.energyDeposit > energyHitCut && hit.second.neutronParentId == -1)
                {
                    bkg_earliestHit_outFV_in3DST_primary = hit.second;
                }
            }

            if(bkg_earliestHit_outFV_in3DST_primary.timeWindow != 1000)
            {
                cout<<"primary KE:"<<kineticEnergy(bkg_earliestHit_outFV_in3DST_primary.trackLength,bkg_earliestHit_outFV_in3DST_primary.timeWindow)<<endl;
                KE_primary->Fill(kineticEnergy(bkg_earliestHit_outFV_in3DST_primary.trackLength,bkg_earliestHit_outFV_in3DST_primary.timeWindow));
            }

        }       //end of if(outFV_in3DST)
    }       //end of for(int event = 0; event < nevents; event++)
}
