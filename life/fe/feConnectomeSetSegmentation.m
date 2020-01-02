function fe = feConnectomeSetSegmentation(fe,SegFileName)
% Set all the fields necessary to store the voxel segmentation info.
%
% Inputs:
%   fe          - an fe structure
%   SegFileName - the full path to a Segmentation nifti file
%        
% Outputs:
%   fe - fe structure with data added in the appropriate fields.

% load the segmentation data data
nifti_Seg = niftiRead(SegFileName);

% Check segmentation data dimensions
dwi = feGet(fe, sprintf('dwi'));
dwisize = size(dwi.nifti.data);
if any((size(nifti_Seg.data) ~= dwisize(1:end-1)))
    disp('segmentation data size do not match with diffusion data size')
    keyboard;
end

% Store voxel segmentation.
fe  = feSet(fe, sprintf('segmentation'), nifti_Seg, feGet(fe,'roi coords'));

% Extract the dwi signal at the coordinates of the connectome
%fe  = feSet(fe, sprintf('diffusion signal image %s',tag), ...
%            dwiGet(dwi, 'diffusion signal image',feGet(fe,'roi coords')) );

end
