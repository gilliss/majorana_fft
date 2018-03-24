void waveform_plotter() {
    TFile* myfile = new TFile("/global/u1/s/stew/OR_run23521.root"); //("/global/project/projectdirs/majorana/data/mjd/surfmjd/data/built/P3KJR/OR_run17454.root");
	TTree* myttree = (TTree*)myfile->Get("MGTree");
	MGTEvent* myevent = 0;
	MGTWaveform* mywaveform = 0;
	myttree->SetBranchAddress("event",&myevent);
	TCanvas* c = new TCanvas("c", "c");
	c->cd();
    int NWFLimit = 100;
    TH1D* myhist[NWFLimit];
	int totalwaveforms = 0;
    int colorindex = 1;
	int totalentries = 0;
    for(int i = 0; i < myttree->GetEntries(); i++)
    {
        if (totalwaveforms < NWFLimit - 1)
        {
            myttree->GetEntry(i);
            totalentries += 1;
            for(int j = 0; j < myevent->GetNWaveforms(); j++)
            {
                
                cout << totalwaveforms << endl;
                
                mywaveform = myevent->GetWaveform(j);
        
                MGWFBaselineRemover* BLR = new MGWFBaselineRemover();
                BLR->SetStartSample(50);
                BLR->SetBaselineSamples(750);
                BLR->TransformInPlace(*mywaveform);

                MGWFExtremumFinder* WFEF = new MGWFExtremumFinder();
                WFEF->SetFindMaximum(true);
                WFEF->SetLocalMinimumTime(0);
                WFEF->SetLocalMaximumTime(20000);
                WFEF->Transform(mywaveform);
                double WFmax = WFEF->GetTheExtremumValue();
                int WFmaxpoint = WFEF->GetTheExtremumPoint();

                if(WFmax > 40 && totalwaveforms < NWFLimit)
                {
                    //cout << "Max #" << i << " = " << WFmax << " " << WFmaxpoint <<  endl;
        
                    MGWaveform mywaveformnorm = ((MGWaveform) *mywaveform) /= (WFmax);
                    mywaveform->SetData(mywaveformnorm.GetVectorData());

                    myhist[totalwaveforms] =  mywaveform->GimmeHist();
                    myhist[totalwaveforms]->SetName("h");
                    myhist[totalwaveforms]->SetLineColor(colorindex);
                    
                    if(totalwaveforms == 0)
                    {
                        myhist[totalwaveforms]->Draw();
                        cout << "Entry = " << i << endl;
                        cout << "WF = " << j << endl;
                    }
                    else
                    {
                        myhist[totalwaveforms]->Draw("same");
                    }
                    
                    totalwaveforms += 1;
                    colorindex += 1;
                }
            }
        }
        else
            break;
    }
cout << "Entries = " << totalentries << endl;
cout << "Waveforms = " << totalwaveforms << endl;
}
