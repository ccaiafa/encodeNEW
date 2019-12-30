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
nifti_Seg = niftiRead(SegFileName)

% Store voxel segmentation.
fe  = feSet(fe, sprintf('segmentaton'),   ...
            dwiGet(dwi, 'diffusionimagenums'));
        
dim = dwi.nifti.dim;
coords = fe.roi.coords;
coords(coords<1)=1;
coords(coords(:,1)>dim(1), 1)=dim(1);
coords(coords(:,2)>dim(2), 2)=dim(2);
coords(coords(:,3)>dim(3), 3)=dim(3);
fe.roi.coords = coords;

% Extract the dwi signal at the coordinates of the connectome
fe  = feSet(fe, sprintf('diffusion signal image %s',tag), ...
            dwiGet(dwi, 'diffusion signal image',feGet(fe,'roi coords')) );
          
% Extract the non-diffusion direction signal at the coordinates of the
% conenctome
fe  = feSet(fe, sprintf('S0 image %s',tag), ...
            dwiGet(dwi, 'S0 image',feGet(fe,'roi coords')));

if length(isrepeat)>1, keyboard;end
          
% Here I set the dimensions of the dwi file so that I have that
% available everytime when creating a map of parameters.
fe = feSet(fe,sprintf('img size %s',tag),dwiGet(dwi,'size'));

end
