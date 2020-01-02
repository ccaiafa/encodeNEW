function [fh, fe] = test_NANKAIdata()


%% (0) Check matlab dependencies and path settings.
if ~exist('vistaRootPath.m','file')
    disp('Vistasoft package either not installed or not on matlab path.')
    error('Please, download it from https://github.com/vistalab/vistasoft');
end

if ~exist('mbaComputeFibersOutliers','file')
    disp('ERROR: mba package either not installed or not on matlab path.')
    error('Please, download it from https://github.com/francopestilli/mba')
end

%% Select b=2000 from multishel data and save as a single shell data
split_shells('../../../../data/segmentation/dwi_data.nii.gz','../../../../data/segmentation/bvecs','../../../../data/segmentation/bvals',200,100,2000)

%% Build the file names for the diffusion data, the anatomical MRI.
t1File        = '../../../../data/segmentation/T1w.nii.gz';

%% Load brain connectome from disk. 
fgFileName = '../../../../data/segmentation/sift_1mio.tck';

%% Identify a DWI file from disk.
dwiFile       = '../../../../data/segmentation/dwi_2000.nii.gz';

% The final connectome and data astructure will be saved with this name:
feFileName    = 'fe_NANKAI';

%% Identify voxel segmentation from disk.
%segFile       = '../../../../data/segmentation/dwi_data_segmentation.nii.gz';
segFile       = '../../../../data/segmentation/change_hcp.nii.gz';


%% Encode connectome and data in a multidimensional tensor. 

L = 360; 

fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],[],t1File,segFile,L,[1,0]);

%% The encoding model is comprised by a large, sparse Phi tensor containing the connectome.
%
% To extract the Phi tensor you can use feGet.m
Phi = feGet(fe, 'Phi');

% The Phi tensor encodes fascicles' nodes orientation in mode 1 (see Caiafa
% Pestilli Figure 1)
Number_of_Orientations = feGet(fe,'n atoms');

% The Phi tensor encodes spatial location of nodes (voxel indices) in mode
% 2 (see Caiafa Pestilli Figure 1).
Number_of_voxels = feGet(fe,'n voxels');

% The Phi tensor encodes fascicles identify in mode 3 (see Caiafa Pestilli
% Figure 1).
Number_of_Fascicles = feGet(fe,'nfibers');

disp(['The size of the sparse tensor Phi is (Na,Nv,Nf) = (',num2str(Number_of_Orientations), ...
      ',',num2str(Number_of_voxels),',',num2str(Number_of_Fascicles),')'])

% The precomputed (demeaned) diffusion signals are stored in a Dictionary
% matrix D. Each column (atom) in the Dictionary corresponds to one spatial
% orientation of a fascicle's node
%
% To extract the Dictionary matrix you can use feGet.m
D = feGet(fe,'Dictionary');

Number_of_gradient_directions = feGet(fe,'nbvecs');
disp(['The size of the dictionary D is (Ntheta,Na) = (', ...
          num2str(Number_of_gradient_directions), ...
      ',',num2str(Number_of_Orientations),')'])

%% Example of identification of fibers connecting two regions: Region A and Region B, each one composed by the union of a series of subregions as included in the segmentation
%save('variables.mat','-v7.3')

%load('variables.mat')

RegionA = [1]; % Subregions of segmentation for Region A
RegionB = [181]; % Subregions of segmentation for Region A

voxelsA = feGetRegionVox(fe,RegionA); % find voxel indices for Region A
Phi_subtensorA = Phi(:,voxelsA,:);
[indsA, ~] = find(Phi_subtensorA); % find nonzero entries of subtensor
fibers_A_indx = unique(indsA(:,3)); % find fibers touching RegionA

voxelsB = feGetRegionVox(fe,RegionB); % find voxel indices for Region B
Phi_subtensorB = Phi(:,voxelsB,:);
[indsB, ~] = find(Phi_subtensorB); % find nonzero entries of subtensor
fibers_B_indx = unique(indsB(:,3)); % find fibers touching RegionB

fibers_AB_indx = intersect(fibers_A_indx,fibers_B_indx); % find fibers touching Regions A and B simultanously
 
% Finally, we generate a visualization of the fascicles and voxels
Visualize_fascicles(fe,fibers_AB_indx,voxelsA, voxelsB, ...
                    'Subset of fascicles connecting A and B')

end

