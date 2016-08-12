%codar_DP_imageFOLtest_07122016.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%codar_DP_imageFOLtest_07122016.m
%
%
%
%
%    Anthony Kirincich
%    WHOI-PO
%    akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set the current directory to that of this mfile
cd ~/matlab/mfiles/codar_DP/ImageFOLs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the test datasets as used in the paper
fnames={'CSS_HATY_16_04_28_0030.cs4';
    'CSS_LPWR_11_04_06_0315.cs4';
    'CSS_LPWR_11_07_19_2230.cs4';
    'CSS_LPWR_11_08_02_0645.cs4';
    'CSS_LPWR_14_08_06_1730.cs4';
    'CSS_NWTP_16_05_20_2300.cs4';
    'CSS_STV1_05_08_24_1620.cs'};

%user_param=[vel_scale max_vel snr_min];
user_param_list=[40 300 5;  %used for HATY
    20 200 5;               %used for LPWR
    20 200 5;
    20 200 5;
    20 200 5;               %note that 
    20 200 5;               %used for NWTP
    70 200 5];              %used for STV1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%plot the [results details]
goplot=[1 0];  %just the results
%goplot=[1 1];  %everything

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% set up the data structure and constants.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set constants and parameters for spectra processing
CONST=[];
%%% other constants that won't change
% speed of light
CONST.c=299792458;  %m/s   %value used by COS
%    c=299705000;   % in air from wiki
%    c=299644310;   % in air based on getting v=0 at ibragg
%    (rad vel difference of 2 mm/s from first value....
%number of antenna elements... not really needed here, but
CONST.N=3;

%%%%%% set a generic HEAD structure for this file %%%%
% normally all the info from the header.txt would go here
%  but this is not included/necessary for looking at the FOLs only
HEAD=[];
HEAD.Site_name='XXXX';
HEAD.Site_loc=[0 0];
patt=[];
%start processing steps list now that HEAD has been established.
HEAD.ProcessingSteps=[];
 


%%
%%%%%% loop over all files in the dataset
for jjj=1:length(fnames)
    

% set parameters for imageFOL
user_param=user_param_list(jjj,:);


%%%%%%%%%%%%%  load the CSS file %%%%%%%%%%%%%%%%%%%%

disp(char(fnames(jjj)))
%%%%%load CSS file
[CSS_Head,Data]=HFR_spect_load_v2(['data/' char(fnames(jjj))]);

%%%% calculated constants from the CSS file %%%
SpecHead=[];
SpecHead.rad_res=CSS_Head.fRangeCellDistKm;
SpecHead.Fc=CSS_Head.fStartFreqMHz.*1e6+CSS_Head.fBandwidthKHz.*1e3/2;
%establish frequency, period
SpecHead.Tr=1/CSS_Head.fRepFreqHz; %s
SpecHead.fr=CSS_Head.fRepFreqHz;  %Hz
SpecHead.fDmax=.5*SpecHead.fr;

%array size
[n m]=size(Data.a3);
SpecHead.rangecell=1:n;  SpecHead.rang=SpecHead.rangecell.*SpecHead.rad_res;
SpecHead.spectra_length=m;

% more frequencies
SpecHead.delta_f=2*SpecHead.fDmax/(SpecHead.spectra_length);     % correct .....
SpecHead.doppler_freq=(-SpecHead.fDmax+SpecHead.delta_f:SpecHead.delta_f:SpecHead.fDmax);  %correct......

%%% find what the bragg frequency is
SpecHead.L=(CONST.c/SpecHead.Fc)/2;
SpecHead.v_p=sqrt(9.81*SpecHead.L/2/pi);
SpecHead.FBragg=2*SpecHead.v_p*SpecHead.Fc/CONST.c;
[s,i]=sort(abs(abs(SpecHead.doppler_freq)-SpecHead.FBragg));
SpecHead.iFBragg=i(1:2);

%get doppler_vel
SpecHead.doppler_vel=SpecHead.doppler_freq*CONST.c/SpecHead.Fc/2;
SpecHead.c_vel=SpecHead.doppler_vel-[-SpecHead.v_p*ones(1,(SpecHead.spectra_length)/2) SpecHead.v_p*ones(1,(SpecHead.spectra_length)/2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%occasionally, there will be an Inf in the data, not sure why.
% fix by changing all indices with data.a3==Inf to 1e-17
% works out to be a amp of -204 dB
i=find(Data.a3==Inf);
%  set all of bad range bin to background for this bin
%  and the bin before it
Data.a3(i)=nan; Data.a1(i)=nan; Data.a2(i)=nan;
Data.a12(i)=nan; Data.a13(i)=nan; Data.a23(i)=nan;

j=find(isnan(mean(Data.a3,2))==1);
i=find(abs(SpecHead.c_vel)<0.75);
for jj=1:length(j)
    a=j(jj)-1:j(jj)+1;
    if a(end)>n; a(end)=n; end
    Data.a1(j(jj),:)=(nanmean(Data.a1(a,:)));
    Data.a2(j(jj),:)=(nanmean(Data.a2(a,:)));
    Data.a3(j(jj),:)=(nanmean(Data.a3(a,:)));
    Data.a23(j(jj),:)=(nanmean(Data.a23(a,:)));
    Data.a13(j(jj),:)=(nanmean(Data.a13(a,:)));
    Data.a12(j(jj),:)=(nanmean(Data.a12(a,:)));
    %now change all data from bad rows inside bragg region to min
    %values
    Data.a3(j(jj),i)=1e-17; Data.a1(j(jj),i)=1e-17; Data.a2(j(jj),i)=1e-17;
    Data.a12(j(jj),i)=1e-17; Data.a13(j(jj),i)=1e-17; Data.a23(j(jj),i)=1e-17;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the FOLs using the imageFOL method
iFBragg=SpecHead.iFBragg;
v_incr=median(diff(SpecHead.c_vel));

[FOreg, FOregi, Alims, HEAD]=imageFOLs_v07122016(Data.a3,iFBragg,v_incr,user_param,goplot,HEAD);

%%
 
pause

end






