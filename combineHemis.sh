# Combine LR gifti files
wb_command -cifti-create-dense-scalar fullBrain.dscalar.nii \
  -left-metric L.atlasroi.32k_fs_LR.shape.gii \
  -right-metric R.atlasroi.32k_fs_LR.shape.gii
