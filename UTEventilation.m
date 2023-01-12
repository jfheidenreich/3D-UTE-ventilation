% skript to analyse the UTE lung MRI data ...

%% import all data as a 3d-matrix
% choosing the primary path containing the images

directory = uigetdir('', 'Please choose directory containing DICOMs.');
cd(directory);

for i = 1:numel(directory)
    path = directory(i).folder; 
    cd(path)
    rawDat(:,:,:,i) = double(read3dDicom(1));
    cd('..');
end
    
%% perform segmentation
%for further information see DOI: 10.1186/s12880-021-00608-1
load('path to pretrained network','net','info')

for cnt = size(rawDat,4);
    SegMasks(:,:,:,cnt) = applySemSeg(rawDat(:,:,:,cnt),net);
end
      
%% calculate volume data
% volumes are calculated, for any other resolution the input-voxSize has to
% be changed

voxSize = 'input voxel size'
volumes = squeeze(sum(sum(sum(SegMasks, 1),2),3));
volumes = volumes.*voxSize;
[B I] = sort(volumes, 'ascend');
volumes = volumes(I,:);

%% sort data 
% most important thing! data is being sorted ascendingly! 
recoSort = rawDat(:,:,:,I); % Morphologie ist nach volumen geordnet
maskSort = SegMasks(:,:,:,I); % Masken sind nach Volumen geordnet
mask = squeeze(maskSort(:,:,:,3));
%% image registration 
%for further information see https://elastix.lumc.nl/
%find further info for matlab wrapper on https://github.com/raacampbell/matlab_elastix

fixed = recoSort(:,:,:,3);
[x y z] = size(fixed);
[~] = copyfile('path to elastix parameter file, [dir, '\regParams.txt']);
maskedAffine = zeros(x,y,z,stad);

for cnt = 1:stad
moving = maskedRaw(:,:,:,cnt);
[maskedAffine(:,:,:,cnt),~] = elastix(moving, fixed, 'C:\temp' , {'regParams.txt'});
end
affine = double(maskedAffine);

fixed = affine(:,:,:,3);
[~] = copyfile('path to BSpline elastix parameter file', [dir, '\regParams.txt']);

for cnt = 1:stad
moving = affine(:,:,:,cnt);
[Bspline(:,:,:,cnt),~] = elastix(moving, fixed, 'C:\temp' , {'regParams.txt'});
end
regDat = double(Bspline);
%% perform a rational fit to model data. 
[rows columns slices stads] = size(regDat);
% reshape regDat
regDatReshape = reshape(normedData,[rows*columns*slices,stads]);

szchk    = size(regDatReshape);
modeldata = regDatReshape.*0;
fitted   = zeros(szchk(1),2); 
     
for cnt = 1:szchk(1)
    if maskmid(cnt) ~=0
       fitted(cnt,:) = coeffvalues(fit(B, squeeze(regDatReshape(cnt,:))', 'rat01','StartPoint',[600 3]));            
    end
end
for cnt = 1:szchk(1)
    if maskmid(cnt) ~=0
       modeldata(cnt,:) = fitted(cnt,1) ./ (B + fitted(cnt,2));
    end
end
     
% build files for further analysis
modeldata    = reshape(modeldata, [ex,columns,slices,stads]);    
%% Evaluation of Ventilation ventresults
%3D ventilation map 
% calculate ventilation per voxel following Maren Zapke et al., 2007 doi:10.1186/1465-9921-7-106
VentNorm = (squeeze(modeldata(:,:,:,2))-squeeze(modeldata(:,:,:,stads-1))) ./ (squeeze(modeldata(:,:,:,2)));
%ventilation parameters
VentNormMean  = mean(VentNorm(find(maskmid)))
VentNormSTD   = std(VentNorm(find(maskmid)))
VentNormIqr   = iqr(VentNorm(find(maskmid))) 

%normalization to mean as denominator 
normedVent = VentNorm./VentNormMean;
normedMean = nanmean(normedVent(find(maskmid)))
normedSTD = nanstd(normedVent(find(maskmid)))
normedIQR = iqr(normedVent(find(maskmid)))
% end of skript!