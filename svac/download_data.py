!pip -q install stsynphot synphot astropy scipy matplotlib

!mkdir -p /content/trds

# Castelli-Kurucz 2004 atmosphere models: 約41 MB
!wget -q -nc -O /content/ck04.tar \
"https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_castelli-kurucz-2004-atlas_multi_v2_synphot3.tar"

!tar -xf /content/ck04.tar -C /content/trds

# PHOENIX models: 約1.8 GB
!wget -q -nc -O /content/phoenix.tar \
"https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_pheonix-models_multi_v3_synphot5.tar"

!tar -xf /content/phoenix.tar -C /content/trds
