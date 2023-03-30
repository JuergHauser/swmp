for dir in random  checkerboard  blobs
do
cd $dir
if test -f "run.py"; then
echo '=============================================================================='
pwd
echo '=============================================================================='
export my_path=`pwd`
export my_demo=`basename $my_path`
cp data/observed.dat ../../../demos/$my_demo/data/
cp output/start/mod.vel ../../../demos/$my_demo/output/start/
cp output/true/mod.vel ../../../demos/$my_demo/output/true/
cp output/current/mod.vel ../../../demos/$my_demo/output/current/
fi
cd ..
done
