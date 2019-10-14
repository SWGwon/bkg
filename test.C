#include "/Users/gwon/analyze.C"
//#include "/Users/gwon/neutron.hxx"

int main()
{
    cout<<"start"<<endl;
    for(int j = 1; j <11; j++)
    {
        for(int i = 1; i <1001; i++)
        {
            cout<<"\033[1APROD"<<j<<": "<<i<<"\033[1000D"<<endl;
            analyze(Form("/Users/gwon/Geo12/PROD%d/FHC_%d.root",j,i));
            analyze(Form("/Users/gwon/Geo12/PROD%d/RHC_%d.root",j,i));
        }
        cout<<"                    "<<endl;
    }
    cout<<"end"<<endl;


    TFile * fi1 = new TFile("background.root","RECREATE");
    hist_signal->Write();
    hist_bkg_out3DST->Write();
    hist_bkg_NC->Write();
    hist_bkg_1->Write();
    hist_bkg_out3DST_largeTime->Write();
    hist_bkg_NC_largeTime->Write();
    hist_bkg_1_largeTime->Write();
    KE_primary->Write();
    KE_secondary->Write();
    neutronParentPDG->Write();
    neutronParentPDG_case4->Write();
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

    TCanvas * can2 = new  TCanvas;
    KE_primary->Draw();
    TCanvas * can3 = new  TCanvas;
    KE_secondary->Draw();

    return 0;
}
