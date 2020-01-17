#include "/Users/gwon/bkg/neutron.hxx"

int main()
{
    neutron();
    /*
    int endPROD, beginPROD;
    cout<<"PROD begin :"<<endl;
    cin>>beginPROD;
    cout<<"PROD end :"<<endl;
    cin>>endPROD;
    cout<<"start"<<endl;
    for(int j = beginPROD; j <endPROD+1; j++)
    {
        for(int i = 1; i <1001; i++)
        {
            cout<<"\033[1APROD"<<j<<": "<<i<<"\033[1000D"<<endl;
            //analyze(Form("/Users/gwon/Geo12/PROD%d/FHC_%d.root",j,i));
            //analyze(Form("/Users/gwon/Geo12/PROD%d/RHC_%d.root",j,i));
            //num_interaction(Form("/Users/gwon/Geo12/PROD%d/FHC_%d.root",j,i));
            //num_interaction(Form("/Users/gwon/Geo12/PROD%d/RHC_%d.root",j,i));
            //test_test_analyze(Form("/Users/gwon/Geo12/PROD%d/RHC_%d_test.root",j,i));
        }
        cout<<endl;
    }
    test_test_analyze(Form("/Users/gwon/Geo12/PROD%d/RHC_%d_test.root",101,1));
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
    can6->Divide(2,3);
    can6->cd(1);
    dist_sig_sp_vtx1->Draw();
    can6->cd(2);
    dist_sig_sp_nh1->Draw();
    can6->cd(3);
    dist_sig_sp_vtx2->Draw();
    can6->cd(4);
    dist_sig_sp_nh2->Draw();
    can6->cd(5);
    dist_sig_sp_vtx3->Draw();
    can6->cd(6);
    dist_sig_sp_nh3->Draw();
    can6->SaveAs("test.pdf");

    return 0;
    */
}
