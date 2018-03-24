void WF_FFT_runs() {
    
    int StartRun = 20745;
    int EndRun = 20755;
    int runRange = EndRun - StartRun;
    GATDataSet ds(StartRun, EndRun);
    TChain* mytchain = ds.GetChains();
    MGTEvent* myevent = 0;
    MGTWaveform* mywaveform = 0;
    vector<double>* trapENFCal = 0;
    mytchain->SetBranchAddress("event", &myevent);
    mytchain->SetBranchAddress("trapENFCal",&trapENFCal);
    int totalWaveforms = 0;
    int totalEntries = 0;
    int startSample = 0;
    int endSample = 900;
    int sampleSize = endSample - startSample;
    //int NFWLimit = 1;
    vector<double> f_n;
    vector<double> powerSum;
    int runValue = 1;
    TCanvas* c2 = new TCanvas("waterfall", "waterfall");
    TH2D* wfWaterfall = new TH2D("wfWaterfall", "wfWaterfall", runRange + 1, StartRun, EndRun + 1, 512, 0, 56.8889);
    
    //Get the waveform from myevent
    for(int i = 0; i < mytchain->GetEntries(); i++)
    {
        mytchain->GetEntry(i);
        for(int j = 0; j < trapENFCal->size(); j++)
        {
            mywaveform = myevent->GetWaveform(j);
            //See if waveform is a good one to use
            MGWFBaselineRemover* BLR = new MGWFBaselineRemover();
            BLR->SetStartSample(50);
            BLR->SetBaselineSamples(750);
            BLR->TransformInPlace(*mywaveform);
            
            if(trapENFCal->at(j) > 5.0)
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
                double samplingPeriod = wfFT->GetSamplingPeriod();
                
                //Loop through FT result data to calculate power vs freq
                if(totalWaveforms == 0)
                {
                    for(int k = 0; k < wfFT_vec.size(); k++)
                    {
                        complex<double> z = complex<double>(wfFT_vec[k].real(), wfFT_vec[k].imag());
                        //Initialize power vector (y axis)
                        double temp_power = (z * conj(z)).real();
                        powerSum.push_back(temp_power);
                        //Create frequency vector (x axis)
                        double temp_frequency = (k/(sampleSize*samplingPeriod)*1000);
                        f_n.push_back(temp_frequency);
                    }
                }
                else
                {
                    for(int k = 0; k < wfFT_vec.size(); k++)
                    {
                        complex<double> z = complex<double>(wfFT_vec[k].real(), wfFT_vec[k].imag());
                        //Create summed power vector (y axis)
                        double temp_power = (z * conj(z)).real();
                        powerSum[k] += temp_power;
                        //Create frequency vector (x axis)
                        double temp_frequency = (k/(sampleSize*samplingPeriod)*1000);
                        f_n[k] = temp_frequency;
                    }
                }
                
                totalWaveforms += 1;
                
            } //if statement for energy
        } //for loop for waveforms
        
        totalEntries += 1;
        
        if(totalEntries == ((mytchain->GetTree()->GetEntries())))
        {
            cout << totalEntries << endl;
            c2->cd();
            //Create TH2D histogram of frequncy, event, and power (as color)
            for(int w = 0; w < f_n.size(); w++)
            {
                //normalize powerSum based on how many waveforms we looked at
                powerSum[w] = (powerSum[w] / totalWaveforms);
                wfWaterfall->SetBinContent(runValue, w + 1, powerSum[w]);
            }
            f_n.clear();
            powerSum.clear();
            totalWaveforms = 0;
            totalEntries = 0;
            runValue += 1;
        }
    } //for loop for entries
    
    //Draw TH2D
    wfWaterfall->SetStats(kFALSE);
    wfWaterfall->Draw("COLZ");
    c2->Update();

    return;
} //void
