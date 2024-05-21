cp ../src/make.inc .
cp ../src/mod_bio_3D.F .
cp ../src/run/input_wqem/NUTRIENT_INI_* run/input_wqem/
cp ../src/run/input_wqem/GOMDOM.in run/input_wqem/
cp ../src/run/input_cgem/CGEM.in run/input_cgem/
cp ../src/run/input_tp/TP.in run/input_tp/
cd BIO_CGEM/
cp ../../src/BIO_CGEM/*.F .
cp ../../src/BIO_CGEM/*.F90 .
cp ../../src/BIO_CGEM/makefile .
cd ../BIO_WQEM/
cp ../../src/BIO_WQEM/*.F .
cp ../../src/BIO_WQEM/makefile .
cp ../../src/BIO_WQEM/GOMDOM_InputFile_LakeOntario .
cd ../BIO_TP/
cp ../../src/BIO_TP/*.F .
cp ../../src/BIO_TP/makefile .
cd ..
