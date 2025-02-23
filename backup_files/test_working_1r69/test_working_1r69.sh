mkdir first
cd first
awsem_create 1r69 --frag
cd 1r69
cp ../../forces_setup.py .
awsem_run 1r69 --platform CUDA --steps 6e6 --tempStart 800 --tempEnd 200 -f forces_setup.py --to run1
cd ../..

mkdir second
cd second
awsem_create 1r69 --frag
cd 1r69
cp ../../forces_setup.py .
awsem_run 1r69 --platform CUDA --steps 6e6 --tempStart 800 --tempEnd 200 -f forces_setup.py --to run1
cd ../..

mkdir third
cd third
awsem_create 1r69 --frag
cd 1r69
cp ../../forces_setup.py .
awsem_run 1r69 --platform CUDA --steps 6e6 --tempStart 800 --tempEnd 200 -f forces_setup.py --to run1
cd ../..

mkdir fourth
cd fourth
awsem_create 1r69 --frag
cd 1r69
cp ../../forces_setup.py .
awsem_run 1r69 --platform CUDA --steps 6e6 --tempStart 800 --tempEnd 200 -f forces_setup.py --to run1
cd ../..

mkdir fifth
cd fifth
awsem_create 1r69 --frag
cd 1r69
cp ../../forces_setup.py .
awsem_run 1r69 --platform CUDA --steps 6e6 --tempStart 800 --tempEnd 200 -f forces_setup.py --to run1
cd ../..

mkdir sixth
cd sixth
awsem_create 1r69 --frag
cd 1r69
cp ../../forces_setup.py .
awsem_run 1r69 --platform CUDA --steps 6e6 --tempStart 800 --tempEnd 200 -f forces_setup.py --to run1
cd ../..
