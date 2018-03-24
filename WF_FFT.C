void WF_FFT() {
    
    TFile* myfile = new TFile("/global/project/projectdirs/majorana/data/mjd/surfmjd/data/built/P3KJR/OR_run17454.root"); //"/global/u1/s/stew/OR_run23521.root");
    TTree* myttree = (TTree*)myfile->Get("MGTree");
    MGTEvent* myevent = 0;
    MGTWaveform* mywaveform = 0;
    myttree->SetBranchAddress("event", &myevent);
    //int totalWavefroms = 0;
    //int totalEntries = 0;
    int startSample = 0;
    int endSample = 900;
    int sampleSize = endSample - startSample;
    //int NFWLimit = 1;
    vector<double> f_n;
    vector<double> power;
    
    //Get the waveform from myevent
    for(int i = 0; i < myttree->GetEntries(); i++)
    {
        myttree->GetEntry(i);
        
        for(int j = 0; j < myevent->GetNWaveforms(); j++)
        {
            mywaveform = myevent->GetWaveform(j);
            
            //See if waveform is a good one to use
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
            
            if(WFmax < 10)
            {
                //Start of FFT section . . .
                //Invoke FFT Transform
                MGWFFastFourierTransformDefault* FT = new MGWFFastFourierTransformDefault();
                //Invoke object to hold FFT result
                MGWaveformFT* wfFT = new MGWaveformFT();
                //Change WF type so we can use it
                MGWaveform* rawWF_MG = (MGWaveform*)mywaveform;
                cout<<"   rawWF_MG length (samples):"<<rawWF_MG->GetLength()<<endl;
                //Set WF region to study
                MGTWaveformRegion* wfRegion_MGT = new MGTWaveformRegion(startSample,endSample);
                MGWaveformRegion* wfRegion_MG = (MGWaveformRegion*)wfRegion_MGT;
                //Perform FFT
                FT->PerformFFT(rawWF_MG, wfFT, wfRegion_MG);
                //Pull out data and store it in vector
                vector<complex<double> > wfFT_vec = wfFT->GetVectorData();
                //Get sampling period for frequency calculation later
                double samplingPeriod = 10; //wfFT->GetSamplingPeriod();
                
                //cout << samplingPeriod << endl;
                cout << wfFT_vec.size() << endl;
                
                //Loop through FT result data to calculate power vs freq
                for(int k = 0; k < wfFT_vec.size(); k++)
                {
                    complex<double> z = complex<double>(wfFT_vec[k].real(), wfFT_vec[k].imag());
                    //Create power vector (y axis)
                    power.push_back((z * conj(z)).real());
                    //Create frequency vector (x axis)
                    f_n.push_back(k/(sampleSize*samplingPeriod)*1000);
                    //cout << f_n[k] << endl;
                }
                
                TCanvas* c = new TCanvas("waveform", "waveform");
                c->Divide(1,2);
                TGraph* wfFTgraph;
                TH1D* waveformHist = new TH1D;
                
                //Plot power and f_n on TGraph on canvas pad 1
                wfFTgraph = new TGraph(f_n.size(), &(f_n[0]), &(power[0]));
                c->cd(1);
                wfFTgraph->Draw();
                wfFTgraph->GetXaxis()->SetLimits(0, f_n[wfFT_vec.size() - 1]);
                //wfFTgraph->GetHistogram()->SetMaximum(500000);
                wfFTgraph->GetHistogram()->SetMinimum(0);
                wfFTgraph->GetHistogram()->SetMaximum(50000);
                //Plot raw waveform on canvas pad 2
                c->cd(2);
                waveformHist = mywaveform->GimmeHist();
                waveformHist->Draw();
                c->Update();
                
                cout << "entry: " << i << "     wf: " << j << endl;
                cout << "Next waveform? (N = no) or save canvas and continue? (S = save): ";
                string decision;
                getline(cin, decision);
                if(decision == "N")
                {
                    return;
                }
                else if(decision == "S")
                {
                    //Save canvas
                    string entryNumber = to_string(i);
                    string waveformNumber = to_string(j);
                    //string fileName = "run17525" + "entry" + entryNumber + "waveform" + waveformNumber + ".pdf";
                    //const char* fileName_p = fileName.c_str();
                    //c->SaveAs(fileName_p);
                    
                    cout << endl;
                    delete c, waveformHist, wfFTgraph;
                    power.clear();
                    f_n.clear();
                }
                else
                {
                    cout << endl;
                    delete c, waveformHist, wfFTgraph;
                    power.clear();
                    f_n.clear();
                }
                
            }
        }
    }
}
