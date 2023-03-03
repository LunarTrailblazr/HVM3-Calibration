# David R Thompson

# This script relies on the caller to have first dark-subtracted and pedestal-shift-corrected all 
# the SRF data from the monochromator sweeps
# It outputs text files, one per base file, corresponding to the SRFs

wl=/home/drt/src/hvm3/l1/data/HVM3_Wavelengths_20221016.txt

# Make srfs
for field in 173 319 465; do
for fin in `ls /beegfs/scratch/drt/20221031_HVM3_SRF/FieldPixel\#${field}/*_darksub`; do
   python /home/drt/src/emit-sds-l1b/utils/makesrf.py --wavelengths ${wl} --target_index ${field} $fin > ${fin}.txt
done
done


