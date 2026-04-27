
#! /bin/bash
for i in {904..1003};
do
    cp INIT_PARS_SLICS/cosmicatlas_slics_voids_real903.ini INIT_PARS_SLICS/cosmicatlas_slics_voids_real"$i".ini
    sed -i "s/densDMALPTrS20.0TETCICz1.041G192V505.0S903.dat/densDMALPTrS20.0TETCICz1.041G192V505.0S"$i".dat/g" INIT_PARS_SLICS/cosmicatlas_slics_voids_real"$i".ini
    sed -i "s/l_LOS_new_dm = 903/l_LOS_new_dm = $i/g" INIT_PARS_SLICS/cosmicatlas_slics_voids_real"$i".ini
    sed -i "s/IC_index = 903/IC_index = $i/g" INIT_PARS_SLICS/cosmicatlas_slics_voids_real"$i".ini
    sed -i "s/Realization = 903/Realization = $i/g" INIT_PARS_SLICS/cosmicatlas_slics_voids_real"$i".ini
    sed -i "s/test_slics_903/test_slics_$i/g" INIT_PARS_SLICS/cosmicatlas_slics_voids_real"$i".ini
    sed -i "s/test_slics_offibinning_903/test_slics_offibinning_$i/g" INIT_PARS_SLICS/cosmicatlas_slics_voids_real"$i".ini

    sed -i "s/seed = 903/seed = $i/g" INIT_PARS_SLICS/cosmicatlas_slics_voids_real"$i".ini
    sed -i "s/LOS903/LOS$i/g" INIT_PARS_SLICS/cosmicatlas_slics_voids_real"$i".ini
    sed -i "s/los903/los$i/g" INIT_PARS_SLICS/cosmicatlas_slics_voids_real"$i".ini
    sed -i "s/R903/R$i/g" INIT_PARS_SLICS/cosmicatlas_slics_voids_real"$i".ini
    sed -i "s/CAT_R903/CAT_R$i/g" INIT_PARS_SLICS/cosmicatlas_slics_voids_real"$i".ini
    sed -i "s/IC_LOS903/IC_LOS$i/g" INIT_PARS_SLICS/cosmicatlas_slics_voids_real"$i".ini
done;

