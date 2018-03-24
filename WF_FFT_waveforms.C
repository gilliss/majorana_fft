void WF_FFT_waveforms() {
    
    //GATDataSet ds(20747);
    //TChain* mytchain = ds.GetChains();
    TChain* mytchain = new TChain("MGTree");
    mytchain->Add("/global/u1/s/stew/OR_run23522.root"); //"/global/project/projectdirs/majorana/data/mjd/surfmjd/data/built/P3KJR/OR_run17454.root");
    MGTEvent* myevent = 0;
    MGTWaveform* mywaveform = 0;
    //vector<double>* trapENFCal = 0;
    mytchain->SetBranchAddress("event", &myevent);
    //mytchain->SetBranchAddress("trapENFCal",&trapENFCal);
    int totalWaveforms = 0;
    int totalEntries = 0;
    int startSample = 0;
    int endSample = 900;
    int sampleSize = endSample - startSample;
    int NFWLimit = 100;
    vector<double> f_n;
    vector<double> power;
    TCanvas* c2 = new TCanvas("waterfall", "waterfall");
    TH2D* wfWaterfall = new TH2D("wfWaterfall", "wfWaterfall", 30000, 0, 30000, 50, 0, 56.8889);
    
    //Get the waveform from myevent
    for(int i = 0; i < mytchain->GetEntries(); i++)
    {
        
        //if(totalWaveforms < NFWLimit)
        //{
            mytchain->GetEntry(i);
            for(int j = 0; j < myevent->GetNWaveforms() ; j++)
            {
                //cout << j << endl;
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
                
                if(WFmax > 5.0)
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
                    
                    if(totalWaveforms == 0)
                    {
                        //Loop through FT result data to calculate power vs freq
                        for(int k = 0; k < wfFT_vec.size(); k++)
                        {
                            complex<double> z = complex<double>(wfFT_vec[k].real(), wfFT_vec[k].imag());
                            //Create power vector (y axis)
                            double temp_power = (z * conj(z)).real();
                            power.push_back(temp_power);
                            //Create frequency vector (x axis)
                            double temp_frequency = (k/(sampleSize*samplingPeriod)*1000*0.8789);
                            f_n.push_back(temp_frequency);
                            
                            //Create TH2D histogram of frequncy, event, and power (as color)
                            wfWaterfall->Fill(totalWaveforms, temp_frequency, temp_power);
                            c2->Update();
                        }
                    }
                    else
                    {
                        //Loop through FT result data to calculate power vs freq (without rederiving entire frequency values)
                        for(int k = 0; k < wfFT_vec.size(); k++)
                        {
                            complex<double> z = complex<double>(wfFT_vec[k].real(), wfFT_vec[k].imag());
                            //Create power vector (y axis)
                            double temp_power = (z * conj(z)).real();
                            power.push_back(temp_power);
                            
                            //Create TH2D histogram of frequncy, event, and power (as color)
                            wfWaterfall->Fill(totalWaveforms, f_n[k], temp_power);
                            c2->Update();
                        }
                    }
                    
                    totalWaveforms += 1;
                    cout << totalWaveforms << endl;
                    
                } //Energy lower limit if statement
                
            } //Get waveform for loop
        //} //if WFNLimit
        totalEntries += 1;
        cout << mytchain->GetEntries() << "       " << totalEntries << endl;
        
        if(totalEntries == mytchain->GetEntries())//totalWaveforms == (NFWLimit - 1))
        {
            c2->cd();
            wfWaterfall->SetStats(kFALSE);
            wfWaterfall->Draw("COLZ");
            return;
        }
        else
        {
            power.clear();
            f_n.clear();
        }
    } //Get entry for loop
} //void
