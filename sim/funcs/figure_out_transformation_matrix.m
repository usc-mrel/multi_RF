%% figure_out_traj_info.m
% 
%  Ziwei Zhao    10092023

%% Get info from Siemens files
if 0
    twix = obj{end};
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dSag')
        dNormalSag = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dSag;
    else
        dNormalSag = 0;
    end
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dCor')
        dNormalCor = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dCor;
    else
        dNormalCor = 0;
    end
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dTra')
        dNormalTra = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dTra;
    else
        dNormalTra = 0;
    end
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}, 'dInPlaneRot')
        dRotAngle = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dInPlaneRot; % [rad]
    else
        dRotAngle = 0; % [rad]
    end
end

% assume axial plane % double check!! 
dNormalSag = 0;
dNormalCor = 0;
dNormalTra = 1;
dRotAngle  = 0;

%% Get a slice offset of a stack in the PCS from Siemens TWIX format
[R_gcs2pcs, phase_sign, read_sign, main_orientation] = siemens_calculate_transform_gcs_to_pcs(dNormalSag, dNormalCor, dNormalTra, dRotAngle);
patient_position = 'HFS';
R_pcs2dcs = siemens_calculate_transform_pcs_to_dcs(patient_position);
R_gcs2dcs = R_pcs2dcs * R_gcs2pcs;
