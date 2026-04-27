#! /bin/bash
# Any modification to a parameter file must be done in the one corresponding to p65
reds=(10.700 10.440 10.190 9.941 9.707 9.471 9.235 9.010 8.794 8.579 8.372 8.166 7.961 7.764 7.576 7.389 7.203 7.026 6.849 8.675 6.508 6.342 6.184 6.022 5.873 5.720 5.575 5.431 5.289 5.150 5.017 4.882 4.757 4.627 4.507 4.385 4.266 4.152 4.038 3.929 3.819 3.715 3.61\
0 3.511 3.411 3.3414 3.221 3.217 3.037 2.949 2.862 2.778 2.695 2.614 2.525 2.458 2.382 2.308 2.235 2.165 2.095 2.028 1.961 1.896 1.833 1.771 1.710 1.650 1.593 1.535 1.480 1.425 1.372 1.321 1.270 1.220 1.172 1.124 1.077 1.032 0.987 0.944 0.901 0.859 0.819 0.779 0.74\
0 0.702 0.664 0.628 0.592 0.557 0.523 0.490 0.487 0.425 0.394 0.364 0.334 0.305 0.276 0.248 0.221 0.194 0.168 0.142 0.117 0.093 0.069 0.045 0.023 0.000)

for i in {0..111}
do
    new_red=${reds[i]}
    let new_index=$i+17
    echo $new_index, $new_red
    cp INIT_PARS_UNITSIM_ASSIGN/cosmicatlas_new_unitsim_p65.ini INIT_PARS_UNITSIM_ASSIGN/cosmicatlas_new_unitsim_p"$new_index".ini
    sed -i "s/Redshift = 3.037/Redshift = ${new_red}/g" INIT_PARS_UNITSIM_ASSIGN/cosmicatlas_new_unitsim_p"$new_index".ini
    sed -i "s/Unitsim_plabel = 65/Unitsim_plabel = ${new_index}/g" INIT_PARS_UNITSIM_ASSIGN/cosmicatlas_new_unitsim_p"$new_index".ini
    sed -i "s/SNAPSHOT_p65/SNAPSHOT_p${new_index}/g" INIT_PARS_UNITSIM_ASSIGN/cosmicatlas_new_unitsim_p"$new_index".ini
    sed -i "s/Name_Catalog_X = dm_field_p65/Name_Catalog_X = dm_field_p${new_index}/g" INIT_PARS_UNITSIM_ASSIGN/cosmicatlas_new_unitsim_p"$new_index".ini
    sed -i "s/Name_Catalog_X_NEW = dm_field_p65_Nres256_MAS1_downsampled_average.dat/Name_Catalog_X_NEW = dm_field_p${new_index}_Nres256_MAS1_downsampled_average.dat/g" INIT_PARS_UNITSIM_ASSIGN/cosmicatlas_new_unitsim_p"$new_index".ini
    sed -i "s/Name_Catalog_X = dm_field_p65_Nres256_MAS1_downsampled_average.dat/Name_Catalog_X = dm_field_p${new_index}_Nres256_MAS1_downsampled_average.dat/g" INIT_PARS_UNITSIM_ASSIGN/cosmicatlas_new_unitsim_p"$new_index".ini 
    sed -i "s/Name_Catalog_Y = UNITSIM_HALOS_NCounts_Nft256_MAS0/Name_Catalog_Y = UNITSIM_HALOS300dm_p${new_index}/g" INIT_PARS_UNITSIM_ASSIGN/cosmicatlas_new_unitsim_p"$new_index".ini
    sed -i "s/dm_p65/dm_p${new_index}/g" INIT_PARS_UNITSIM_ASSIGN/cosmicatlas_new_unitsim_p"$new_index".ini
    sed -i "s/_p65/_p${new_index}/g" INIT_PARS_UNITSIM_ASSIGN/cosmicatlas_new_unitsim_p"$new_index".ini
    sed -i "s/out_65p/out_${new_index}p/g" INIT_PARS_UNITSIM_ASSIGN/cosmicatlas_new_unitsim_p"$new_index".ini  
    sed -i "s/UNITSIM_p65/UNITSIM_p${new_index}/g" INIT_PARS_UNITSIM_ASSIGN/cosmicatlas_new_unitsim_p"$new_index".ini
    sed -i "s/CAL_SNAPSHOT_p65/CAL_SNAPSHOT_p${new_index}/g" INIT_PARS_UNITSIM_ASSIGN/cosmicatlas_new_unitsim_p"$new_index".ini
done;
