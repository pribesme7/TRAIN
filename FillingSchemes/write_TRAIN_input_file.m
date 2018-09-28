clear;

fill = 5137;
nominalIntensity = 1.2e11;
scan = 1;

data=ExtractEmittanceScans(ConfigForIP5(), fill, ...
    {'LHC.BCTFR.A6R4.B1:BUNCH_INTENSITY','LHC.BCTFR.A6R4.B2:BUNCH_INTENSITY'});%, ...
%    'LHC.BQM.B1:BUNCH_INTENSITIES','LHC.BQM.B2:BUNCH_INTENSITIES', ...
%    'LHC.BCTDC.A6R4.B1:BEAM_INTENSITY','LHC.BCTDC.A6R4.B2:BEAM_INTENSITY'});


%bunchIntB1_bqm = BunchIntensityFromBQM(data{scan}.additional{3}, data{scan}.additional{5}, data{scan}.meta.filledSlots.b1);
%bunchIntB2_bqm = BunchIntensityFromBQM(data{scan}.additional{4}, data{scan}.additional{6}, data{scan}.meta.filledSlots.b2);

bunchIntensityToUse_B1 = data{scan}.additional{1};
bunchIntensityToUse_B2 = data{scan}.additional{2};

trainMatrix = zeros(3564,5);
trainMatrix(:,1) = 1:3564;
trainMatrix(data{scan}.meta.filledSlots.b1,2) = 1;
trainMatrix(data{scan}.meta.filledSlots.b2,3) = 1;
trainMatrix(:,4) = mean(bunchIntensityToUse_B1(:,2:end),1) ./ nominalIntensity;
trainMatrix(:,5) = mean(bunchIntensityToUse_B2(:,2:end),1) ./ nominalIntensity;

fileid = fopen(sprintf('train_%d_scan%d.in',fill,scan), 'w');
fprintf(fileid,'%d %d %d %.5f %.5f\n', trainMatrix');
fclose(fileid);