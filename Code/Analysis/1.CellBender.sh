#python 3.7
#cellbender_0.3.0

#create virtual env for cellbender
module load python/3.7
virtualenv --no-download --clear ~/py_3.7_cellbender_0.3.0 && source ~/py_3.7_cellbender_0.3.0/bin/activate
pip install -U pip
pip install cellbender numexpr==2.7.0 torch==1.10.0


#activate virtual env
module load python/3.7
source ~/py_3.7_cellbender_0.3.0/bin/activate

#launch cellbender
cellbender remove-background \
     --cuda \
     --input </path/to/raw.h5> \
     --output </path/to/cellbender.h5> \

#convert to seurat format
ptrepack --complevel 5 </path/to/cellbender.h5>:/matrix </path/to/cellbender.seurat.h5>:/matrix