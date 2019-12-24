# created by Guang Yang Jan 2019
#!/usr/bin/env python

import sys
import os.path
import os
import ROOT
import random
import math
from optparse import OptionParser
from array import array

lar_active_vols = [ "vol3DST" , "volCube" , "volcube" ]
neutron_mass = 939.5654133
speedOfLight = 30 #cm per ns


def loop( events, tgeo, tout, nfiles, okruns, GeoNow, doDetSim ):

    if len(okruns) == 0:
        print "There are no runs in this TTree...skipping!"
        return

    print "Inside event loop with %d files and first run %d" % (nfiles, okruns[0])

    #GeoNow=5;

    # updated geometry with less steel
    # cm
    # geometry 1 and 3 
    if GeoNow=="Geo1" or GeoNow=="Geo3":
    	offset = [ 0., 0., -50. ]
    	fvLo = [ -90., -90., -90. ]
    	fvHi = [ 90., 90., 90. ]
    	collarLo = [ -90., -90., -90. ]
    	collarHi = [ 90., 90., 90. ]
	size3DST = [200., 200., 200.]
        minEndInCM = [-100, -100, -100 ]
        maxEndInCM = [100, 100, 100]

    # geometry 2
    if GeoNow=="Geo2":
    	offset = [ 0., 0., -50. ]
    	fvLo = [ -140., -90., -90. ]
    	fvHi = [ 140., 90., 90. ]
    	collarLo = [ -140., -140., -90. ]
    	collarHi = [ 140., 140., 90. ]
        size3DST = [300., 200., 200.]
        minEndInCM = [-150, -100, -100 ]
        maxEndInCM = [150, 100, 100]

    # geometry 4 and 5
    if GeoNow=="Geo4":
    	offset = [ 0., 0., -50. ]
    	fvLo = [ -2000., -2000., -2000. ]
    	fvHi = [ 2000., 2000., 2000. ]
    	collarLo = [ -2000., -2000., -2000. ]
   	collarHi = [ 2000., 2000., 2000. ]
        size3DST = [200., 200., 200.]
        minEndInCM = [-100, -100, -100 ]
        maxEndInCM = [100, 100, 100]

    if GeoNow=="Geo5" or GeoNow=="Geo6" or GeoNow=="Geo7":
        offset = [ 0., 0., 550.-50. ]
        fvLo = [ -2000., -2000., -2000. ]
        fvHi = [ 2000., 2000., 2000. ]
        collarLo = [ -2000., -2000., -2000. ]
        collarHi = [ 2000., 2000., 2000. ]
        size3DST = [200., 200., 200.]
	minEndInCM = [-100, -100, -100 ]
	maxEndInCM = [100, 100, 100]

    if GeoNow=="Geo8":
        offset = [ 0., 0., -50. ]
        fvLo = [ -50., -50., -50. ]
        fvHi = [ 50., 50., 50. ]
        collarLo = [ -50., -50., -50. ]
        collarHi = [ 50., 50., 50. ]
        size3DST = [120., 120., 120.]
        minEndInCM = [-100, -100, -100 ]
        maxEndInCM = [100, 100, 100]

    if GeoNow=="Geo9":
        offset = [ 0., 0., 550.-50. ]
        fvLo = [ -2000., -2000., -2000. ]
        fvHi = [ 2000., 2000., 2000. ]
        collarLo = [ -2000., -2000., -2000. ]
        collarHi = [ 2000., 2000., 2000. ]
        size3DST = [120., 120., 120.]
        minEndInCM = [-100, -100, -100 ]
        maxEndInCM = [100, 100, 100]

    if GeoNow=="Geo10" or GeoNow=="Geo12":
        offset = [ 0., 0., 550.-50. ]
        fvLo = [ -2000., -2000., -2000. ]
        fvHi = [ 2000., 2000., 2000. ]
        collarLo = [ -2000., -2000., -2000. ]
        collarHi = [ 2000., 2000., 2000. ]
        size3DST = [240., 240., 200.]
        #minEndInCM = [-100, -100, -100 ]
        #maxEndInCM = [100, 100, 100]

    if GeoNow=="Geo11":
	offset = [ 0., 0., 550.-50. ]
        fvLo = [ -2000., -2000., -2000. ]
        fvHi = [ 2000., 2000., 2000. ]
        collarLo = [ -2000., -2000., -2000. ]
        collarHi = [ 2000., 2000., 2000. ]
        size3DST = [120., 120., 120.]

    if GeoNow == "Geo2":
    	sideExitParameter=150.
    elif GeoNow=="Geo8" or GeoNow=="Geo9" or GeoNow=="Geo11":
	sideExitParameter=60.
    elif GeoNow=="Geo10" or GeoNow=="Geo12":
        sideExitParameter=120.
    else:
    	sideExitParameter=100.

    print 'GeoNow is ', GeoNow

    event = ROOT.TG4Event()
    events.SetBranchAddress("Event",ROOT.AddressOf(event))
    print "Set branch address"

    N = events.GetEntries()
    evt_per_file = N/nfiles
    if N % nfiles:
        print "Files don't all have the same number of events!!!"
        print "\n\n\n\n\n\n\n\n\n\n\n"

    print "Starting loop over %d entries" % N
    for ient in range(N):
        if ient % 100 == 0:
            print "Event %d of %d..." % (ient,N)
        events.GetEntry(ient)
        for ivtx,vertex in enumerate(event.Primaries):

            ## initialize output variables
            fileidx = ient/evt_per_file
            t_ifileNo[0] = okruns[fileidx]
            t_ievt[0] = ient%evt_per_file;
            t_vtx[0]=0.0; t_vtx[1]=0.0; t_vtx[2]=0.0;
            t_p3lep[0]=0.0; t_p3lep[1]=0.0; t_p3lep[2]=0.0;
            t_lepDeath[0]=0.0; t_lepDeath[1]=0.0; t_lepDeath[2]=0.0;
            t_lepPdg[0] = 0
            t_lepKE[0] = 0.

            t_p3pi[0]=0.0; t_p3pi[1]=0.0; t_p3pi[2]=0.0;
            t_piDeath[0]=0.0; t_piDeath[1]=0.0; t_piDeath[2]=0.0;
            t_piPdg[0] = 0
            t_piKE[0] = 0.

            t_p3proton[0]=0.0; t_p3proton[1]=0.0; t_p3proton[2]=0.0;
            t_protonDeath[0]=0.0; t_protonDeath[1]=0.0; t_protonDeath[2]=0.0;
            t_protonPdg[0] = 0
            t_protonKE[0] = 0.

            t_muonExitPt[0] = 0.0; t_muonExitPt[1] = 0.0; t_muonExitPt[2] = 0.0; 
            t_muonExitMom[0] = 0.0; t_muonExitMom[1] = 0.0; t_muonExitMom[2] = 0.0; 
            t_muonReco[0] = -1;
            t_muGArLen[0]=0.0;

            t_hadTot[0] = 0.; t_hadTot_ECAL[0] = 0.; t_hadTot_3DST[0] = 0.; t_hadTot_TPC[0] = 0.;
	    t_hadTot_allECAL[0] = 0.; t_hadTot_leak[0] = 0.;
            t_hadCollar[0] = 0.; t_vtxTime[0] = 0.;
            t_nFS[0] = 0
	    for inii in range(1000):
	        t_neutronHitX[inii]=0.;
                t_neutronHitY[inii]=0.;
                t_neutronHitZ[inii]=0.;
                t_neutronHitT[inii]=0.;
                t_neutronHitSmearT[inii]=0.;
                t_neutronHitE[inii]=0.;
                t_neutronRecoE[inii]=0.;
                t_neutronHitS[inii]=0.;
		t_neutronTrueE[inii]=0.;
		t_neutronParentId[inii]=0.;
		t_neutronParentPDG[inii]=0.;
            ## done

            # now ID numucc
            reaction=vertex.Reaction        

            # set the vertex location for output
            for i in range(3): 
                t_vtx[i] = vertex.Position[i] / 10. - offset[i] # cm
		#print t_vtx[i]

            # fiducial vertex cut
            fvCut = False
            for i in range(3):
                if t_vtx[i] < fvLo[i] or t_vtx[i] > fvHi[i]:
                    fvCut = True
            vtxv = ROOT.TVector3( t_vtx[0], t_vtx[1], t_vtx[2] )
            if fvCut:
                continue

            ileptraj = -1
	    ipi0traj = -1
	    ipitraj = -1

            nfsp = 0
            # get the lepton kinematics from the edepsim file
            fsParticleIdx = {}
            for ipart,particle in enumerate(vertex.Particles):
                e = particle.Momentum[3]
                p = (particle.Momentum[0]**2 + particle.Momentum[1]**2 + particle.Momentum[2]**2)**0.5
                m = (e**2 - p**2)**0.5
	        print "fs pdg list ",particle.PDGCode
                t_fsPdg[nfsp] = particle.PDGCode
                t_fsPx[nfsp] = particle.Momentum[0]
                t_fsPy[nfsp] = particle.Momentum[1]
                t_fsPz[nfsp] = particle.Momentum[2]
                t_fsE[nfsp] = e
                fsParticleIdx[particle.TrackId] = nfsp
                nfsp += 1
                if abs(particle.PDGCode) in [11,12,13,14]:
                    ileptraj = particle.TrackId
                    t_lepPdg[0] = particle.PDGCode
                    # set the muon momentum for output
                    for i in range(3): t_p3lep[i] = particle.Momentum[i]
                    t_lepKE[0] = e - m

                if abs(particle.PDGCode) in [211]:
                    ipitraj = particle.TrackId
                    t_piPdg[0] = particle.PDGCode
                    # set the pion momentum for output
                    for i in range(3): t_p3pi[i] = particle.Momentum[i]
                    t_piKE[0] = e - m

                if abs(particle.PDGCode) in [111]:
                    ipi0traj = particle.TrackId

                if abs(particle.PDGCode) in [2212]:
                    iprotontraj = particle.TrackId
                    t_protonPdg[0] = particle.PDGCode
                    # set the proton momentum for output
                    for i in range(3): t_p3proton[i] = particle.Momentum[i]
                    t_protonKE[0] = e - m

            assert ileptraj != -1, "There isn't a lepton??"
            t_nFS[0] = nfsp


            # If there is a muon, determine how to reconstruct its momentum and charge
            t_muexit[0] = 0
            exitKE = 0.
            exitP = None
            endVolIdx = -1 # where does the muon die
            if abs(t_lepPdg[0]) == 13:
                leptraj = event.Trajectories[ileptraj]
                for p in leptraj.Points:
                    pt = p.Position
                    node = tgeo.FindNode( pt.X(), pt.Y(), pt.Z() )
                    volName = node.GetName()
                    active = False
                    for v in lar_active_vols:
                        if v in volName:
                            active = True
                            break
                    if active:
                        t_muonExitPt[0] = pt.X() / 10. - offset[0]
                        t_muonExitPt[1] = pt.Y() / 10. - offset[1]
                        t_muonExitPt[2] = pt.Z() / 10. - offset[2]
                        continue
                    # first hit outside 3DST -- determine exit
                    t_muonExitPt[0] = pt.X() / 10. - offset[0]
                    t_muonExitPt[1] = pt.Y() / 10. - offset[1]
                    t_muonExitPt[2] = pt.Z() / 10. - offset[2]
                    t_muonExitMom[0] = p.Momentum.x()
                    t_muonExitMom[1] = p.Momentum.y()
                    t_muonExitMom[2] = p.Momentum.z()
                    if abs(pt.X() / 10. - offset[0]) > sideExitParameter: 
			t_muexit[0] = 1 # side exit
                    elif abs(pt.Y() / 10. - offset[1]) > 100.: 
			t_muexit[0] = 2 # top/bottom exit
                    elif pt.Z() / 10. - offset[2] < -100.: 
			t_muexit[0] = 3 # upstream exit
                    elif pt.Z() / 10. - offset[2] > 100.: 
			t_muexit[0] = 4 # downstream exit
                    else:
                        print "Hit in %s at position (%1.1f, %1.1f, %1.1f) unknown exit!" % (volName, pt.X()/10.-offset[0], pt.Y()/10.-offset[1], pt.Z()/10.-offset[2])
                    exitP = p.Momentum
                    exitKE = (exitP.x()**2 + exitP.y()**2 + exitP.z()**2 + 105.3**2)**0.5 - 105.3
                    break

                endpt = leptraj.Points[-1].Position

                node = tgeo.FindNode( endpt.X(), endpt.Y(), endpt.Z() )

                t_lepDeath[0] = endpt.X()/10. - offset[0]
                t_lepDeath[1] = endpt.Y()/10. - offset[1]
                t_lepDeath[2] = endpt.Z()/10. - offset[2]

                endVolName = node.GetName()
		#print endVolName
		#print t_muexit[0]

                # dipole+HPGTPC
                if "volWorld" in endVolName or "volDetEnclosure" in endVolName: endVolIdx = 0 # outside detector components
                elif "volCube" in endVolName or "volcube" in endVolName: endVolIdx = 1 # 3DST
                elif "volTPC" in endVolName or "volTPCTop" in endVolName: endVolIdx = 2 # TPC
                elif "volECALStripx" in endVolName or "volECALStripy" in endVolName or "volRadiatorPlate" in endVolName: endVolIdx = 3 # ECAL
                elif "volMagnet" in endVolName: endVolIdx = 4 # Magnet

                # look for muon hits in the gas TPC
                hits = []
                for key in event.SegmentDetectors:
                    if key.first in ["volTPC", "volTPCTop"]:
                        hits += key.second

                tot_length = 0.0
                for hit in hits:
                    if hit.PrimaryId == ileptraj: # hit is due to the muon
                        # TG4HitSegment::TrackLength includes all delta-rays, which spiral in gas TPC and give ridiculously long tracks
                        hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                        hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                        tot_length += (hStop-hStart).Mag()

                t_muGArLen[0] = tot_length

                # muon reconstruction method
                # 1 = contained
                if endVolIdx == 1:
                    t_muonReco[0] = 1
                # 2 = gas TPC match
                #elif tot_length > 0.:
		elif endVolIdx == 2:
                    t_muonReco[0] = 2
                # 3 = ECAL stopping
                elif endVolIdx == 3: 
                    t_muonReco[0] = 3
                # 4 = magnet/coil stopper
                elif endVolIdx == 4: 
                    t_muonReco[0] = 4


            # If there is a pion, determine how to reconstruct its momentum and charge
            t_piexit[0] = 0
            exitKE = 0.
            exitP = None
            endVolIdx = -1 # where does the pion die
            if abs(t_piPdg[0]) == 211:
                pitraj = event.Trajectories[ipitraj]
                for p in pitraj.Points:
                    pt = p.Position
                    node = tgeo.FindNode( pt.X(), pt.Y(), pt.Z() )
                    volName = node.GetName()
                    active = False
                    for v in lar_active_vols:
                        if v in volName:
                            active = True
                            break
                    if active:
                        t_pionExitPt[0] = pt.X() / 10. - offset[0]
                        t_pionExitPt[1] = pt.Y() / 10. - offset[1]
                        t_pionExitPt[2] = pt.Z() / 10. - offset[2]
                        continue
                    # first hit outside 3DST -- determine exit
                    t_pionExitPt[0] = pt.X() / 10. - offset[0]
                    t_pionExitPt[1] = pt.Y() / 10. - offset[1]
                    t_pionExitPt[2] = pt.Z() / 10. - offset[2]
                    t_pionExitMom[0] = p.Momentum.x()
                    t_pionExitMom[1] = p.Momentum.y()
                    t_pionExitMom[2] = p.Momentum.z()
                    if abs(pt.X() / 10. - offset[0]) > 150.:
                        t_piexit[0] = 1 # side exit
                    elif abs(pt.Y() / 10. - offset[1]) > 100.:
                        t_piexit[0] = 2 # top/bottom exit
                    elif pt.Z() / 10. - offset[2] < -100.:
                        t_piexit[0] = 3 # upstream exit
                    elif pt.Z() / 10. - offset[2] > 100.:
                        t_piexit[0] = 4 # downstream exit
                    else:
                        print "Hit in %s at position (%1.1f, %1.1f, %1.1f) unknown exit!" % (volName, pt.X()/10.-offset[0], pt.Y()/10.-offset[1], pt.Z()/10.-offset[2])
                    exitP = p.Momentum
                    exitKE = (exitP.x()**2 + exitP.y()**2 + exitP.z()**2 + 105.3**2)**0.5 - 105.3
                    break

                endpt = pitraj.Points[-1].Position

                node = tgeo.FindNode( endpt.X(), endpt.Y(), endpt.Z() )

                t_piDeath[0] = endpt.X()/10. - offset[0]
                t_piDeath[1] = endpt.Y()/10. - offset[1]
                t_piDeath[2] = endpt.Z()/10. - offset[2]

                endVolName = node.GetName()
                #print endVolName
                #print t_piexit[0]

                # dipole+HPGTPC
                if "volWorld" in endVolName or "volDetEnclosure" in endVolName: endVolIdx = 0 # outside detector components
                elif "volCube" in endVolName or "volcube" in endVolName: endVolIdx = 1 # 3DST
                elif "volTPC" in endVolName or "volTPCTop" in endVolName: endVolIdx = 2 # TPC
                elif "volECALStripx" in endVolName or "volECALStripy" in endVolName or "volRadiatorPlate" in endVolName: endVolIdx = 3 # ECAL
                elif "volMagnet" in endVolName: endVolIdx = 4 # Magnet

                # look for pion hits in the gas TPC
                hits = []
                for key in event.SegmentDetectors:
                    if key.first in ["volTPC", "volTPCTop"]:
                        hits += key.second

                tot_length = 0.0
                for hit in hits:
                    if hit.PrimaryId == ipitraj: # hit is due to the pion
                        # TG4HitSegment::TrackLength includes all delta-rays, which spiral in gas TPC and give ridiculously long tracks
                        hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                        hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                        tot_length += (hStop-hStart).Mag()

                t_piGArLen[0] = tot_length

                # pion reconstruction method
                # 1 = contained
                if endVolIdx == 1:
                    t_pionReco[0] = 1
                # 2 = gas TPC match
                #elif tot_length > 0.:
                elif endVolIdx == 2:
                    t_pionReco[0] = 2
                # 3 = ECAL stopping
                elif endVolIdx == 3:
                    t_pionReco[0] = 3
                # 4 = magnet/coil stopper
                elif endVolIdx == 4:
                    t_pionReco[0] = 4


            # If there is a proton, determine how to reconstruct its momentum and charge
            t_protonexit[0] = 0
            exitKE = 0.
            exitP = None
            endVolIdx = -1 # where does the proton die
            if abs(t_protonPdg[0]) == 2212:
                protontraj = event.Trajectories[iprotontraj]
                for p in protontraj.Points:
                    pt = p.Position
                    node = tgeo.FindNode( pt.X(), pt.Y(), pt.Z() )
                    volName = node.GetName()
                    active = False
                    for v in lar_active_vols:
                        if v in volName:
                            active = True
                            break
                    if active:
                        t_protonExitPt[0] = pt.X() / 10. - offset[0]
                        t_protonExitPt[1] = pt.Y() / 10. - offset[1]
                        t_protonExitPt[2] = pt.Z() / 10. - offset[2]
                        continue
                    # first hit outside 3DST -- determine exit
                    t_protonExitPt[0] = pt.X() / 10. - offset[0]
                    t_protonExitPt[1] = pt.Y() / 10. - offset[1]
                    t_protonExitPt[2] = pt.Z() / 10. - offset[2]
                    t_protonExitMom[0] = p.Momentum.x()
                    t_protonExitMom[1] = p.Momentum.y()
                    t_protonExitMom[2] = p.Momentum.z()
                    if abs(pt.X() / 10. - offset[0]) > 150.:
                        t_protonexit[0] = 1 # side exit
                    elif abs(pt.Y() / 10. - offset[1]) > 100.:
                        t_protonexit[0] = 2 # top/bottom exit
                    elif pt.Z() / 10. - offset[2] < -100.:
                        t_protonexit[0] = 3 # upstream exit
                    elif pt.Z() / 10. - offset[2] > 100.:
                        t_protonexit[0] = 4 # downstream exit
                    else:
                        print "Hit in %s at position (%1.1f, %1.1f, %1.1f) unknown exit!" % (volName, pt.X()/10.-offset[0], pt.Y()/10.-offset[1], pt.Z()/10.-offset[2])
                    exitP = p.Momentum
                    exitKE = (exitP.x()**2 + exitP.y()**2 + exitP.z()**2 + 105.3**2)**0.5 - 105.3
                    break

                endpt = protontraj.Points[-1].Position

                node = tgeo.FindNode( endpt.X(), endpt.Y(), endpt.Z() )

                t_protonDeath[0] = endpt.X()/10. - offset[0]
                t_protonDeath[1] = endpt.Y()/10. - offset[1]
                t_protonDeath[2] = endpt.Z()/10. - offset[2]

                endVolName = node.GetName()
                #print endVolName
                #print t_protonexit[0]

                # dipole+HPGTPC
                if "volWorld" in endVolName or "volDetEnclosure" in endVolName: endVolIdx = 0 # outside detector components
                elif "volCube" in endVolName or "volcube" in endVolName: endVolIdx = 1 # 3DST
                elif "volTPC" in endVolName or "volTPCTop" in endVolName: endVolIdx = 2 # TPC
                elif "volECALStripx" in endVolName or "volECALStripy" in endVolName or "volRadiatorPlate" in endVolName: endVolIdx = 3 # ECAL
                elif "volMagnet" in endVolName: endVolIdx = 4 # Magnet

                # look for proton hits in the gas TPC
                hits = []
                for key in event.SegmentDetectors:
                    if key.first in ["volTPC", "volTPCTop"]:
                        hits += key.second

                tot_length = 0.0
                for hit in hits:
                    if hit.PrimaryId == iprotontraj: # hit is due to the proton
                        # TG4HitSegment::TrackLength includes all delta-rays, which sprotonral in gas TPC and give ridiculously long tracks
                        hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                        hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                        tot_length += (hStop-hStart).Mag()

                t_protonGArLen[0] = tot_length

                # proton reconstruction method
                # 1 = contained
                if endVolIdx == 1:
                    t_protonReco[0] = 1
                # 2 = gas TPC match
                #elif tot_length > 0.:
                elif endVolIdx == 2:
                    t_protonReco[0] = 2
                # 3 = ECAL stopprotonng
                elif endVolIdx == 3:
                    t_protonReco[0] = 3
                # 4 = magnet/coil stopper
                elif endVolIdx == 4:
                    t_protonReco[0] = 4



            # hadronic containment -- find hits in TPC
            hitsTPC = []
            for key in event.SegmentDetectors:
		if "volTPC" in key.first:
                    hitsTPC += key.second

            total_energy_TPC = 0.
            track_length_TPC = [0. for i in range(nfsp)]
            for hit in hitsTPC:
                hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                if event.Trajectories[hit.PrimaryId].ParentId == -1:
                    track_length_TPC[fsParticleIdx[hit.PrimaryId]] += (hStop-hStart).Mag()

                if hit.PrimaryId != ileptraj:
                    hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                    total_energy_TPC += hit.EnergyDeposit

            t_hadTot_TPC[0] = total_energy_TPC
            #print total_energy_TPC

	    # hits in ECAL without radiator
            hitsECAL = []
            for key2 in event.SegmentDetectors:
                if "volECAL" in key2.first:
                    hitsECAL += key2.second

            total_energy_ECAL = 0.
            track_length_ECAL = [0. for i in range(nfsp)]
            for hit2 in hitsECAL:
                hStart = ROOT.TVector3( hit2.Start[0]/10.-offset[0], hit2.Start[1]/10.-offset[1], hit2.Start[2]/10.-offset[2] )
                hStop = ROOT.TVector3( hit2.Stop[0]/10.-offset[0], hit2.Stop[1]/10.-offset[1], hit2.Stop[2]/10.-offset[2] )
                if event.Trajectories[hit2.PrimaryId].ParentId == -1:
                    track_length_ECAL[fsParticleIdx[hit2.PrimaryId]] += (hStop-hStart).Mag()

                if hit2.PrimaryId != ileptraj:
                    hStart = ROOT.TVector3( hit2.Start[0]/10.-offset[0], hit2.Start[1]/10.-offset[1], hit2.Start[2]/10.-offset[2] )
                    total_energy_ECAL += hit2.EnergyDeposit
	
            t_hadTot_ECAL[0] = total_energy_ECAL
	    #print total_energy_ECAL

            # hits in ECAL with radiator
            hitsECAL = []
            for key2 in event.SegmentDetectors:
                if "volECAL" in key2.first or "volRadiator" in key2.first:
                    hitsECAL += key2.second

            total_energy_ECAL = 0.
            track_length_ECAL = [0. for i in range(nfsp)]
            for hit2 in hitsECAL:
                hStart = ROOT.TVector3( hit2.Start[0]/10.-offset[0], hit2.Start[1]/10.-offset[1], hit2.Start[2]/10.-offset[2] )
                hStop = ROOT.TVector3( hit2.Stop[0]/10.-offset[0], hit2.Stop[1]/10.-offset[1], hit2.Stop[2]/10.-offset[2] )
                if event.Trajectories[hit2.PrimaryId].ParentId == -1:
                    track_length_ECAL[fsParticleIdx[hit2.PrimaryId]] += (hStop-hStart).Mag()

                if hit2.PrimaryId != ileptraj:
                    hStart = ROOT.TVector3( hit2.Start[0]/10.-offset[0], hit2.Start[1]/10.-offset[1], hit2.Start[2]/10.-offset[2] )
                    total_energy_ECAL += hit2.EnergyDeposit

            t_hadTot_allECAL[0] = total_energy_ECAL

            # hits outside sensitive detectors
            hitsLeak = []
            for key2 in event.SegmentDetectors:
                if "volECAL" in key2.first or "volRadiator" in key2.first or "volTPC" in key2.first or "volCube" in key2.first or "volcube" in key2.first:
		    print "contained in sensitive detector"
		else:
                    hitsLeak += key2.second

            total_energy_Leak = 0.
            track_length_Leak = [0. for i in range(nfsp)]
            for hit2 in hitsLeak:
                hStart = ROOT.TVector3( hit2.Start[0]/10.-offset[0], hit2.Start[1]/10.-offset[1], hit2.Start[2]/10.-offset[2] )
                hStop = ROOT.TVector3( hit2.Stop[0]/10.-offset[0], hit2.Stop[1]/10.-offset[1], hit2.Stop[2]/10.-offset[2] )
                if event.Trajectories[hit2.PrimaryId].ParentId == -1:
                    track_length_Leak[fsParticleIdx[hit2.PrimaryId]] += (hStop-hStart).Mag()

                if hit2.PrimaryId != ileptraj:
                    hStart = ROOT.TVector3( hit2.Start[0]/10.-offset[0], hit2.Start[1]/10.-offset[1], hit2.Start[2]/10.-offset[2] )
                    total_energy_Leak += hit2.EnergyDeposit

            t_hadTot_leak[0] = total_energy_Leak


            # event vtx time
            t_vtxTime[0] = random.uniform(0.0,10000.)

	    # Get neutron trajectory information PDG neutron 2112 proton 2212
	    nTraj = []

	    icounter = 0
            for ctraj in event.Trajectories:
                if event.Trajectories[ctraj.TrackId].PDGCode == 2112:

                    #print "a FS neutron with kinetic energy ", event.Trajectories[ctraj.TrackId].InitialMomentum.Energy()-939.565

		    #t_neutronTrueE[icounter]= event.Trajectories[ctraj.TrackId].InitialMomentum.Energy()-939.565

		    #print "Parent PDG ", event.Trajectories[ctraj.ParentId].PDGCode
		    #print event.Trajectories[ctraj.TrackId].InitialMomentum.Energy()-939.565
		    #print event.Trajectories[ctraj.TrackId].InitialMomentum.T()
		    icounter = icounter+1
                    #print ctraj.ParentId
                    #print ctraj.TrackId
                    #print ctraj.Name

	    parentDaughterMap = dict()
	    neutronMap = dict()

            for ctraj in event.Trajectories:
		if event.Trajectories[ctraj.ParentId].PDGCode == 2112:
		    #print "a neutron hit"
		    #print "parent is neutron, then neutron's parent is : ", event.Trajectories[event.Trajectories[ctraj.ParentId].TrackId].PDGCode
		    #print ctraj.ParentId
		    #print ctraj.TrackId
		    #print ctraj.Name
		    nTraj.append(ctraj.TrackId)
		    parentDaughterMap[ctraj.TrackId] = ctraj.ParentId

            	    for ctraj4 in event.Trajectories:
                        if ctraj4.TrackId == ctraj.ParentId:
			    neutronMap[ctraj.TrackId] = ctraj4.ParentId
			    if ctraj4.ParentId == -1:
			    	print "neutron induced ID ",ctraj.TrackId, " neutron parent Id ", ctraj4.ParentId
			    	print "neutron induced PDG ",event.Trajectories[ctraj.TrackId].PDGCode, " parent of parent PDG ", event.Trajectories[ctraj4.ParentId].PDGCode


                #if ctraj.PDGCode == 2112:
		    #print ctraj.Name
                    #nTraj.append(ctraj.TrackId)

	    print "----------------- ",ient

	    """
	    nTraj_1daughter = []
            for cctraj in event.Trajectories:
		#print cctraj.ParentId
                for iloop in nTraj:
                    if cctraj.ParentId == iloop:
                        #print cctraj.Name
			nTraj_1daughter.append(cctraj.TrackId)

            nTraj_2daughter = []
            for cctraj in event.Trajectories:
                for iloop in nTraj_1daughter:
                    if cctraj.ParentId == iloop:
                        #print cctraj.Name
			nTraj_2daughter.append(cctraj.TrackId)

            nTraj_3daughter = []
            for cctraj in event.Trajectories:
                for iloop in nTraj_2daughter:
                    if cctraj.ParentId == iloop:
                        #print cctraj.Name
                        nTraj_3daughter.append(cctraj.TrackId)
	    """

	    # neutron hit number looper that will be used soon
	    iN = 0;     
       
	    # hits in 3DST 

            hits = []
	    cubeInven = []
            for key3 in event.SegmentDetectors:
                if "volCube" in key3.first or "volcube" in key3.first or "vol" in key3.first:
                    hits += key3.second

            collar_energy = 0.
            total_energy_3DST = 0.
            track_length = [0. for i in range(nfsp)]
            for hit in hits:
                hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                if event.Trajectories[hit.PrimaryId].ParentId == -1:
                    track_length[fsParticleIdx[hit.PrimaryId]] += (hStop-hStart).Mag()

                if hit.PrimaryId != ileptraj:
                    hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                    total_energy_3DST += hit.EnergyDeposit
                    # check if hit is in collar region
                    if hStart.x() < collarLo[0] or hStart.x() > collarHi[0] or hStart.y() < collarLo[1] or hStart.y() > collarHi[1] or hStart.z() < collarLo[2] or hStart.z() > collarHi[2]:
                        collar_energy += hit.EnergyDeposit

		if hit.PrimaryId == ipi0traj:
		    print '---------------!!!!!!!!!!!!!!!!!' 

		#print 'random time stamp ',random.uniform(0.0,10000.)
		#print 'hit time is ',hit.Start.T()
		#print 'generated time is ',hit.Start.T()+random.uniform(0.0,10000.)

		#print hit.Contrib[0]
		#Neutron PDG code 2112, neutron should not give energy deposit itself
	 	for iloop in nTraj:
		    #print "checking"
		    #print iloop
#                    if hit.Contrib[0] == iloop:
#		    	thisMatch = False
#		    	cubeInfo = getCube(minEndInCM, maxEndInCM, (hStart.x()+hStop.x())/2., (hStart.y()+hStop.y())/2., (hStart.z()+hStop.z())/2.)
#	   		for cubeloop in len(cubeInven):
#			    if cubeInfo == cubeInven[cubeloop]:
#				thisMatch = True
#			if thisMatch == False:	
#			    cubeInven.append(cubeInfo) 
#			    iN += 1;

		    if hit.Contrib[0] == iloop and iN < 999:
			#print hStart.x()
			#print 'energy deposit by ',event.Trajectories[hit.PrimaryId].Name

                        #print iN, "--- neutron induced ID ",iloop, " neutron parent Id ", neutronMap.get(iloop)
                        #print "--- neutron induced PDG ",event.Trajectories[iloop].PDGCode, " parent of parent PDG ", event.Trajectories[neutronMap.get(iloop)].PDGCode

			if doDetSim == False: 
			    t_neutronHitX[iN] = hStart.x()
                            t_neutronHitY[iN] = hStart.y()
                            t_neutronHitZ[iN] = hStart.z()
			else: 
			    #detResSim
			    dist = math.sqrt( math.sqrt( math.pow(t_vtx[0]-hStart.x(),2) + math.pow(t_vtx[1]-hStart.y(),2) + math.pow(t_vtx[2]-hStart.z(),2) ))
			    calVec=detResSim( hStart.x(), hStart.y() , hStart.z(), hit.Start.T(), hit.EnergyDeposit, size3DST[0], size3DST[1], size3DST[2], dist )
			    #print 'checking location x y z and time and recoE and distance',' ',calVec[0],' ',calVec[1],' ',calVec[2],' ',calVec[3],' ',calVec[4],' ',dist
			    t_neutronHitX[iN] = calVec[0]
                            t_neutronHitY[iN] = calVec[1]
                            t_neutronHitZ[iN] = calVec[2]
                            t_neutronHitSmearT[iN] = calVec[3]+t_vtxTime[0]
			    t_neutronRecoE[iN] = calVec[4]

			t_neutronTrueE[iN]= event.Trajectories[iloop].InitialMomentum.Energy()-939.565
			t_neutronHitT[iN] = hit.Start.T()+t_vtxTime[0]
			t_neutronHitE[iN] = hit.EnergyDeposit
                        t_neutronHitS[iN] = event.Trajectories[hit.PrimaryId].PDGCode
			if (neutronMap.get(iloop) > -1000000):
			    t_neutronParentId[iN] = neutronMap.get(iloop)
			    t_neutronParentPDG[iN] = event.Trajectories[neutronMap.get(iloop)].PDGCode
			parID = parentDaughterMap.get(hit.PrimaryId)
			#print type(parID),' ',hit.PrimaryId,' ',parentDaughterMap.get(hit.PrimaryId)
                        #print 'checking',event.Trajectories[parID].PDGCode
			#print event.Trajectories[hit.PrimaryId].InitialMomentum.Energy()
			iN += 1;

		"""
                for iloop in nTraj_1daughter:
                    if hit.Contrib[0] == iloop:
                        #print hStart.x()
                        t_neutronHitX[iN] = hStart.x()
                        t_neutronHitY[iN] = hStart.y()
                        t_neutronHitZ[iN] = hStart.z()
                        t_neutronHitE[iN] = hit.EnergyDeposit
                        t_neutronHitS[iN] = event.Trajectories[hit.PrimaryId].PDGCode
                        iN += 1;

                for iloop in nTraj_2daughter:
                    if hit.Contrib[0] == iloop:
                        #print hStart.x()
                        t_neutronHitX[iN] = hStart.x()
                        t_neutronHitY[iN] = hStart.y()
                        t_neutronHitZ[iN] = hStart.z()
                        t_neutronHitE[iN] = hit.EnergyDeposit
                        t_neutronHitS[iN] = event.Trajectories[hit.PrimaryId].PDGCode
                        iN += 1;

                for iloop in nTraj_3daughter:
                    if hit.Contrib[0] == iloop:
                        #print hStart.x()
                        t_neutronHitX[iN] = hStart.x()
                        t_neutronHitY[iN] = hStart.y()
                        t_neutronHitZ[iN] = hStart.z()
                        t_neutronHitE[iN] = hit.EnergyDeposit
                        t_neutronHitS[iN] = event.Trajectories[hit.PrimaryId].PDGCode
                        iN += 1;
		"""

	    parentDaughterMap.clear()
            t_hadTot_3DST[0] = total_energy_3DST

            #print total_energy_3DST

            t_hadTot[0] = total_energy_3DST + total_energy_ECAL + total_energy_TPC
            t_hadCollar[0] = collar_energy


            for i in range(nfsp):
                t_fsTrkLen[i] = track_length[i]

            tout.Fill()


def getCube( start, end, posx, posy, posz ):

    pz= int(math.modf(posz)[1])
    py= int(math.modf(posy)[1])
    px= int(math.modf(posx)[1])

    res = []
    res.append(pz - start[2])
    res.append(py - start[1])
    res.append(px - start[0])
    res.append(res[0]*1000000 + res[1]*1000 + res[2])

    return res

def detResSim( threeVecX, threeVecY, threeVecZ, time, ene, size3dstX, size3dstY, size3dstZ, dist ):
    aveX = int(threeVecX+size3dstX/2)
    aveY = int(threeVecY+size3dstY/2)
    aveZ = int(threeVecZ+size3dstZ/2)
    xlocation = aveX + 0.5 - size3dstX/2;
    ylocation = aveY + 0.5 - size3dstY/2;
    zlocation = aveZ + 0.5 - size3dstZ/2;
    tRes = (1/math.sqrt(3))* math.sqrt((ene * 70.8 * 0.38)/50.)
    aTime = random.gauss(time, tRes)
    if dist/(speedOfLight * aTime)<0 or dist/(speedOfLight * aTime)>1 :
	recoE = 0;
    else:
	#print 'neutron_mass dist speedOfLight aTime time resolution ',neutron_mass,' ',dist,' ',speedOfLight,' ',aTime,' ',time,' ',tRes,' ',dist/(speedOfLight * aTime)
	recoE = neutron_mass/(math.sqrt(1-math.pow(dist/(speedOfLight * aTime),2))) - neutron_mass
        #print math.pow((72.3 * dist/100.)/aTime,2)
    vec = []
    vec.append(xlocation)
    vec.append(ylocation)
    vec.append(zlocation)
    vec.append(aTime)
    vec.append(recoE)
    #print 'doing the det. sim., the location, original, resolution and smeared time and reco. E are: ', xlocation,' ',ylocation, ' ',zlocation, ' ',time,' ',tRes,' ',aTime, ' ', recoE
    return vec

if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--topdir', help='Input file top directory', default="")
    parser.add_option('--first_run', type=int, help='First run number', default=0)
    parser.add_option('--last_run', type=int, help='Last run number', default=0)
    parser.add_option('--rhc', action='store_true', help='Reverse horn current', default=False)
    parser.add_option('--grid', action='store_true', help='grid mode', default=False)
    parser.add_option('--geo', help='geometry number', default="Geo4")
    parser.add_option('--doDetSim', help='if do detSim', default=False)

    (args, dummy) = parser.parse_args()

    # make an output ntuple
    fout = ROOT.TFile( args.outfile, "RECREATE" )
    tout = ROOT.TTree( "tree","tree" )
    t_ifileNo = array('i',[0])
    tout.Branch('ifileNo',t_ifileNo,'ifileNo/I')
    t_ievt = array('i',[0])
    tout.Branch('ievt',t_ievt,'ievt/I')
    t_p3lep = array('f',3*[0.0])
    tout.Branch('p3lep',t_p3lep,'p3lep[3]/F')
    t_p3pi = array('f',3*[0.0])
    tout.Branch('p3pi',t_p3pi,'p3pi[3]/F')
    t_p3proton = array('f',3*[0.0])
    tout.Branch('p3proton',t_p3proton,'p3proton[3]/F')

    t_vtx = array('f',3*[0.0])
    tout.Branch('vtx',t_vtx,'vtx[3]/F')
    t_lepDeath = array('f',3*[0.0])
    tout.Branch('lepDeath',t_lepDeath,'lepDeath[3]/F')
    t_lepPdg = array('i',[0])
    tout.Branch('lepPdg',t_lepPdg,'lepPdg/I')
    t_lepKE = array('f',[0])
    tout.Branch('lepKE',t_lepKE,'lepKE/F')

    t_piPdg = array('i',[0])
    tout.Branch('piPdg',t_piPdg,'piPdg/I')
    t_piKE = array('f',[0])
    tout.Branch('piKE',t_piKE,'piKE/F')
    t_piDeath = array('f',3*[0.0])
    tout.Branch('piDeath',t_piDeath,'piDeath[3]/F')

    t_protonPdg = array('i',[0])
    tout.Branch('protonPdg',t_protonPdg,'protonPdg/I')
    t_protonKE = array('f',[0])
    tout.Branch('protonKE',t_protonKE,'protonKE/F')
    t_protonDeath = array('f',3*[0.0])
    tout.Branch('protonDeath',t_protonDeath,'protonDeath[3]/F')

    t_muexit = array('i',[0])
    tout.Branch('muexit',t_muexit,'muexit/I')
    t_muonExitPt = array('f',3*[0.0])
    tout.Branch('muonExitPt',t_muonExitPt,'muonExitPt[3]/F')
    t_muonExitMom = array('f',3*[0.0])
    tout.Branch('muonExitMom',t_muonExitMom,'muonExitMom[3]/F')
    t_muonReco = array('i',[0])
    tout.Branch('muonReco',t_muonReco,'muonReco/I')
    t_muGArLen = array('f',[0])
    tout.Branch('muGArLen',t_muGArLen,'muGArLen/F')

    t_piexit = array('i',[0])
    tout.Branch('piexit',t_piexit,'piexit/I')
    t_pionExitPt = array('f',3*[0.0])
    tout.Branch('pionExitPt',t_pionExitPt,'pionExitPt[3]/F')
    t_pionExitMom = array('f',3*[0.0])
    tout.Branch('pionExitMom',t_pionExitMom,'pionExitMom[3]/F')
    t_pionReco = array('i',[0])
    tout.Branch('pionReco',t_pionReco,'pionReco/I')
    t_piGArLen = array('f',[0])
    tout.Branch('piGArLen',t_piGArLen,'piGArLen/F')

    t_protonexit = array('i',[0])
    tout.Branch('protonexit',t_protonexit,'protonexit/I')
    t_protonExitPt = array('f',3*[0.0])
    tout.Branch('protonExitPt',t_protonExitPt,'protonExitPt[3]/F')
    t_protonExitMom = array('f',3*[0.0])
    tout.Branch('protonExitMom',t_protonExitMom,'protonExitMom[3]/F')
    t_protonReco = array('i',[0])
    tout.Branch('protonReco',t_protonReco,'protonReco/I')
    t_protonGArLen = array('f',[0])
    tout.Branch('protonGArLen',t_protonGArLen,'protonGArLen/F')

    t_hadTot = array('f', [0.] )
    tout.Branch('hadTot', t_hadTot, 'hadTot/F' )
    t_hadTot_TPC = array('f', [0.] )
    tout.Branch('hadTot_TPC', t_hadTot_TPC, 'hadTot_TPC/F' )
    t_hadTot_3DST = array('f', [0.] )
    tout.Branch('hadTot_3DST', t_hadTot_3DST, 'hadTot_3DST/F' )
    t_hadTot_ECAL = array('f', [0.] )
    tout.Branch('hadTot_ECAL', t_hadTot_ECAL, 'hadTot_ECAL/F' )
    t_hadTot_allECAL = array('f', [0.] )
    tout.Branch('hadTot_allECAL', t_hadTot_allECAL, 'hadTot_allECAL/F' )
    t_hadTot_leak = array('f', [0.] )
    tout.Branch('hadTot_leak', t_hadTot_leak, 'hadTot_leak/F' )

    t_hadCollar = array('f', [0.] )
    tout.Branch('hadCollar', t_hadCollar, 'hadCollar/F' )
    t_nFS = array('i',[0])
    tout.Branch('nFS',t_nFS,'nFS/I')
    t_fsPdg = array('i',100*[0])
    tout.Branch('fsPdg',t_fsPdg,'fsPdg[nFS]/I')
    t_fsPx = array('f',100*[0.])
    tout.Branch('fsPx',t_fsPx,'fsPx[nFS]/F')
    t_fsPy = array('f',100*[0.])
    tout.Branch('fsPy',t_fsPy,'fsPy[nFS]/F')
    t_fsPz = array('f',100*[0.])
    tout.Branch('fsPz',t_fsPz,'fsPz[nFS]/F')
    t_fsE = array('f',100*[0.])
    tout.Branch('fsE',t_fsE,'fsE[nFS]/F')
    t_fsTrkLen = array('f',100*[0.])
    tout.Branch('fsTrkLen',t_fsTrkLen,'fsTrkLen[nFS]/F')

    t_neutronHitX = array('f' ,1000*[0.])
    tout.Branch('neutronHitX',t_neutronHitX,'neutronHitX[1000]/F')
    t_neutronHitY = array('f' ,1000*[0.])
    tout.Branch('neutronHitY',t_neutronHitY,'neutronHitY[1000]/F')
    t_neutronHitZ = array('f' ,1000*[0.])
    tout.Branch('neutronHitZ',t_neutronHitZ,'neutronHitZ[1000]/F')
    t_neutronHitT = array('f' ,1000*[0.])
    tout.Branch('neutronHitT',t_neutronHitT,'neutronHitT[1000]/F')
    t_neutronHitSmearT = array('f' ,1000*[0.])
    tout.Branch('neutronHitSmearT',t_neutronHitSmearT,'neutronHitSmearT[1000]/F')
    t_neutronHitE = array('f' ,1000*[0.])
    tout.Branch('neutronHitE',t_neutronHitE,'neutronHitE[1000]/F')
    t_neutronRecoE = array('f' ,1000*[0.])
    tout.Branch('neutronRecoE',t_neutronRecoE,'neutronRecoE[1000]/F')
    t_neutronHitS = array('f' ,1000*[0.])
    tout.Branch('neutronHitPDG',t_neutronHitS,'neutronHitPDG[1000]/F')
    t_neutronTrueE = array('f' ,1000*[0.])
    tout.Branch('neutronTrueE',t_neutronTrueE,'neutronTrueE[1000]/F')
    t_neutronParentId = array('f' ,1000*[0.])
    tout.Branch('neutronParentId',t_neutronParentId,'neutronParentId[1000]/F')
    t_neutronParentPDG = array('f' ,1000*[0.])
    tout.Branch('neutronParentPDG',t_neutronParentPDG,'neutronParentPDG[1000]/F')

    t_vtxTime = array('f', [0.] )
    tout.Branch('vtxTime', t_vtxTime, 'vtxTime/F' )

    loaded = False
    #if os.path.exists("EDepSimEvents/EDepSimEvents.so"):
    #    print "Found EDepSimEvents.so"
    #    ROOT.gSystem.Load("EDepSimEvents/EDepSimEvents.so")
    #    loaded = True

    tgeo = None

    events = ROOT.TChain( "EDepSimEvents", "main event tree" )

    neutrino = "neutrino"
    if args.rhc:
        neutrino = "antineutrino"

    print "Building TChains for runs %d-%d..." % (args.first_run, args.last_run)
    nfiles = 0
    okruns = []
    for run in range( args.first_run, args.last_run+1 ):
        if args.grid:
            fname = "%s/edep.%d.root" % (args.topdir,run)
        else:
            fname = "%s/%02d/full3DST.%s.%d.edepsim.root" % (args.topdir, run/1000, neutrino, run)
        print fname

        # see if it is an OK file
        if not os.access( fname, os.R_OK ):
            print "Can't access file: %s" % fname
            continue
        tf = ROOT.TFile( fname )
        if tf.TestBit(ROOT.TFile.kRecovered): # problem with file
            print "File is crap: %s" % fname
            continue
        nfiles += 1
        okruns.append( run )

        if not loaded:
            loaded = True
            tf.MakeProject("EDepSimEvents","*","RECREATE++")

        # add it to the tchain
        events.Add( fname )

        if tgeo is None: # first OK file, get geometry
            tgeo = tf.Get("EDepSimGeometry")
        tf.Close() # done with this one

    print "OK runs: ", sorted(okruns)
    print "got %d events in %d files = %1.1f events per file" % (events.GetEntries(), nfiles, 1.0*events.GetEntries()/nfiles)
    loop( events, tgeo, tout, nfiles, sorted(okruns), args.geo, args.doDetSim )

    fout.cd()
    tout.Write()




