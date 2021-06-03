PREF="M2"
cp ${PREF}.crd min/
cd min/
sander -O -i min.in -p ${PREF}.prmtop -c ${PREF}.crd -r ${PREF}_min.ncrst  -o min.out -x min.nc 
cp ${PREF}_min.ncrst ../heat/
cd ../heat/
pmemd.cuda -O -i heat.in -p ${PREF}.prmtop -c ${PREF}_min.ncrst -r ${PREF}_heat.ncrst -ref ${PREF}_min.ncrst -o heat.out -l heat.log -x heat.nc 
cp ${PREF}_heat.ncrst ../nvt/
cd ../nvt/
pmemd.cuda -O -i nvt.in -p ${PREF}.prmtop -c ${PREF}_heat.ncrst -r ${PREF}_nvt.ncrst -ref ${PREF}_heat.ncrst -o nvt.out -l nvt.log -x nvt.nc 
cp ${PREF}_nvt.ncrst ../npt/
cd ../npt/
pmemd.cuda -O -i npt.in -p ${PREF}.prmtop -c ${PREF}_nvt.ncrst -r ${PREF}_npt.ncrst -ref ${PREF}_nvt.ncrst -o npt.out -l npt.log -x npt.nc 
cd ../
echo "Done"
