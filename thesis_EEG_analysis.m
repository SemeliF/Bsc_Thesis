close all
clear
clc

%% SEMELI FRANGOPOULOU - SNR 2020452

%% SETUP Environment Variables
%Set working directory
w_path = "C:\Users\Semeli\Desktop\Thesis_Matlab";
ad_path = "C:\Users\Semeli\Desktop\Thesis_Matlab\AD";
hc_path = "C:\Users\Semeli\Desktop\Thesis_Matlab\HC";

%enter dir
cd(w_path);

%% CLEAN DATASET
ad_dataset = ingest(ad_path);
% hc_dataset = ingest(hc_path);

%% GET PSD
%UNCOMMENT FOR PSD
%ad_dataset_psd = get_ds_psd(ad_dataset);
%hc_dataset_psd = get_ds_psd(hc_dataset);
%figure
%plot(hc_psd_dataset(8).PSD_Freqs,hc_psd_dataset(8).PSD_Spectra)
% 
%plotting the PSD for each band
% dfig = figure('Name','Delta');
% 
% plot(freqs(deltaIdx),spectra(:,deltaIdx))
% tfig = figure('Name','Theta');
% plot(freqs(thetaIdx),spectra(:,thetaIdx))
% afig = figure('Name','Alpha');
% plot(freqs(alphaIdx),spectra(:,alphaIdx))
% bfig = figure('Name','Beta');
% plot(freqs(betaIdx),spectra(:,betaIdx))
% gfig = figure('Name','Gamma');
% plot(freqs(gammaIdx),spectra(:,gammaIdx))


%% GET Phase Sync.
%Loop for all files in dataset

for i = 1:size(ad_dataset,1)
    ad_dataset_ph(i,1) = get_ph_sync(ad_dataset(i));
end

% for i = 1:size(hc_dataset,1)
%     hc_dataset_ph(i,1) = get_ph_sync(hc_dataset(i));
% end

%% END

return

%% Functions

%This function reads all files and prepares the data
function processed_dataset = ingest(dataset_folder)

    %get a list of files in the directory
    files = dir(dataset_folder);

    %get size of list
    [r,c] = size(files);

    %remove . and .. dirs
    names = strings(r-2,1);

    %get full path for each file
    for i = 3:r
        names(i-2,1) = strcat(dataset_folder,'\',files(i).name);
    end

    %Load Data
    for i = length(names):-1:1
    raw_dataset(i,1) = importdata(flip(names(i,1)));    
    end

    %Guess sampling rate based on length + Downsample to 500 Hz
    raw_dataset_srate = zeros(length(raw_dataset),1);

    %Transfer header info
    downsampled_dataset = raw_dataset;
    %downsampled_dataset = raw_dataset;
    for i = 1:length(raw_dataset)
        if length(raw_dataset(i,1).data) > 15000
            raw_dataset_srate(i,1)=2000;
            downsampled_dataset(i,1).data=downsample(raw_dataset(i,1).data,4);
        else
            raw_dataset_srate(i,1)=1000;
            downsampled_dataset(i,1).data=downsample(raw_dataset(i,1).data,2);
        end
    end

    %List of pairs that are of interest and have to be extarcted from each file
    pair_list = [ "F8-F4" ; "F7-F3"; "F4-C4"; "F3-C3"; "F4-FZ"; "FZ-CZ"; "F3-FZ"; "T4-C4"; "T3-C3"; "C4-CZ"; "C3-CZ"; "CZ-PZ"; "C4-P4"; "C3-P3"; "T4-T6"; "T3-T5"; "P4-PZ"; "P3-PZ"; "T6-O2"; "T5-O1"; "P4-O2"; "P3-O1"; "O1-O2"];

    %for each file in dataset
    for i = 1:length(downsampled_dataset)

        data = [];
        colheaders = strings(1,1);

        %for each column label
        for j = 1:length(pair_list)
            found_counter=0;
            string_buffer = strings(50,1);
            col_idx = zeros(50,1);
            %for each column in each file
            for k =1:length(downsampled_dataset(i,1).colheaders)
                if contains(regexprep(downsampled_dataset(i,1).colheaders(1,k),'0','O','all'),pair_list(j,1),'IgnoreCase',true)
                    found_counter=found_counter+1;
                    string_buffer(found_counter,1)=downsampled_dataset(i,1).colheaders(1,k);
                    %string_buffer(found_counter,1)=regexprep(downsampled_dataset(i,1).colheaders(1,k),'0','O','all');
                    col_idx(found_counter,1)=k;
                end
            end
            %if no matches delete (by doing nothing)

            if found_counter == 1 %if there is one match for column type

                data(:,end+1) = downsampled_dataset(i,1).data(:,col_idx(1,1));
                colheaders(1,end+1) = pair_list(j,1);

            elseif found_counter > 1 %if there are multiple select the smallest

                %arbirary large value
                current_min = 4000000;
                smallest_idx = 0;

                for b = 1:found_counter
                    tmp = split(string_buffer(b,1),' ');
                    tmp = str2num(regexprep(tmp(1,1),'[^0-9]','','all'));
                    if current_min > tmp
                        current_min = tmp;
                        smallest_idx=col_idx(b,1);
                    end  
                end
                data(:,end+1) = downsampled_dataset(i,1).data(:,smallest_idx);            
                colheaders(1,end+1) = pair_list(j,1);
            end
        end
        %delete first columns header which is empty (little string array hack
        %fix)
        colheaders=colheaders(1,2:end);
        
        %Create a struct instance 
        struct_instance = struct('path',names(i),'data',data','time',downsampled_dataset(i,1).data(:,1)','labels',colheaders');
        %Push to return 
        processed_dataset(i,1)=struct_instance;

    end

end

%This function reads all files and derives PSD data
function dataset = get_ds_psd(cle_dataset)
    for k = 1:size(cle_dataset,1)

        imported = pop_importdata('setname', cle_dataset(k).path, 'data',cle_dataset(k).data,'srate',500);
        band_filtered = pop_eegfiltnew(imported, 'locutoff', 0.1, 'hicutoff', 100);
        %notch_filtered = pop_eegfiltnew(band_filtered, 'locutoff', 47.85,'hicutoff', 50.5  , 'revfilt', 1);
        notch_filtered = pop_eegfiltnew(band_filtered, 'locutoff', 48,'hicutoff', 52  , 'revfilt', 1);
        
        plot = 0;
        if plot == 1
            fig = figure('Name','PSD Full Spectrum');
            [spectra,freqs] = spectopo(notch_filtered.data, 0, 500); %psd for the entire signal
            %close(fig);
        else %no plotting
            [spectra,freqs] = spectopo(notch_filtered.data, 0, 500,'plot','off'); %psd for the entire signal
        end
        
        %Get indexes for values within the bands
        deltaIdx = find(freqs>=1 & freqs<=4);
        thetaIdx = find(freqs>=4 & freqs<=8);
        alphaIdx = find(freqs>=8 & freqs<=13);
        betaIdx  = find(freqs>=13 & freqs<=30);
        gammaIdx = find(freqs>=36 & freqs<=44);

        %Calculate mean power for each band
        dP = mean(10.^(spectra(deltaIdx) / 10));
        tP = mean(10.^(spectra(thetaIdx) / 10));
        aP = mean(10.^(spectra(alphaIdx) / 10));
        bP = mean(10.^(spectra(betaIdx) / 10));
        gP = mean(10.^(spectra(gammaIdx) / 10));

        %create storage struct

        struct_instance = struct('Datafile_Path',cle_dataset(k).path,'PSD_Spectra',spectra,'PSD_Freqs',freqs','Channel_Labels',cle_dataset(k).labels,'DeltaPower',dP,'ThetaPower',tP,'AlphaPower',aP,'BetaPower',bP,'GammaPower',gP);

        %push to struct array
        psd_dataset(k,1)=struct_instance;

    end
    
    dataset = psd_dataset;
    
end

%This function takes 1 file and calculates phase syncronization between all
%electrode pairs for all bands (adaptation of Sue's code)
function ph_sync_data = get_ph_sync(one_file_data)

    data = one_file_data.data;
    timevec = one_file_data.time;
    phase_data = zeros(2, length(timevec)); 
    real_data  = zeros(2, length(timevec));

    %(23^2 -23 )/2 possible unique pairs
    band = zeros(5, 253);

    %create a string array to store labels
    label_tracker = strings(1,253);
    
    %band is frex +-fwhm
    frex_list = [2.5 6 10.5 21.5 40];
    fwhm_list = [0.8 0.5 0.3 0.13 0.27];

    for s = 1:5
        srate = 500;
        time = -1.5: 1/srate : 1.5;
        frex = frex_list(s);
        fwhm = fwhm_list(s);

        sine_wave = exp( 1i*2*pi*frex*time );
        gaus_win = exp( (-4 * log(2)*time.^2) / fwhm^2 );

        wavelet = sine_wave .* gaus_win;
        half_wavN = (length(time)-1)/2;

        nWave = length(time);
        nData = size(data,2);
        nConv = nWave + nData - 1;

        waveletX = fft(wavelet,nConv);
        waveletX = waveletX ./ max(waveletX);

        band_idx = 1;

        for a = 1:22
            chan1idx = a;
            for j = a+1 : 23
                chan2idx = j;

                dataX = fft(data(chan1idx, :), nConv);
                as = ifft(waveletX.*dataX,nConv);
                as = as(half_wavN+1:end-half_wavN);

                phase_data(1,:) = angle(as); % extract phase angles
                real_data(1,:)  = real(as);  % extract the real part (projection onto real axis)


                dataX = fft(data(chan2idx, : ), nConv);
                as = ifft(waveletX.*dataX,nConv);
                as = as(half_wavN+1:end-half_wavN);

                phase_data(2,:) = angle(as);
                real_data(2,:)  = real(as);

                %pahse syncronization value for the current 2 channels for the
                %current band
                band(s, band_idx) = abs(mean(exp(1i * diff(phase_data))));
                
                %store label
                label_tracker(1,band_idx) = strcat('(',one_file_data.labels(a,1),')&(',one_file_data.labels(j,1),')');

                %increment band
                band_idx = band_idx + 1;

            end
        end
    end
    
    ph_sync_data = struct('path',one_file_data.path,'data',band,'combos',label_tracker);
end
