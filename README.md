# PrimerDesign
### IPython Notebook for Primer Design

##### Make a folder for the virtual environment
mkdir -p ~/venv/base
wget --no-check-certificate https://raw.github.com/pypa/virtualenv/master/virtualenv.py

##### Create a new virtual environment
python virtualenv.py --no-site-packages --no-setuptools ~/venv/base

##### Download get-pip.py for installing pip and setuptools in the new environment
wget --no-check-certificate https://bootstrap.pypa.io/get-pip.py

##### Activate the environment. This should be added to the .bashrc or run manually every time.
source ~/venv/base/bin/activate

##### Install pip and setuptools
python get-pip.py

##### Now everything is ready for installing the remaining dependencies
pip install ipython[all]
pip install regex
pip install pandas
pip install biopython
pip install runipy
pip install selenium

##### To install rpy2 it is necessary to load R 3.1.0
module load R/3.1.0
pip install rpy2

##### Installing the primer3 package
wget https://pypi.python.org/packages/source/p/primer3-py/primer3-py-0.3.1.tar.gz
tar -xvzf primer3-py-0.3.1.tar.gz
cd primer3-py-0.3.1
pip install primer3-py
cd ..

##### Download the latest firefox to your computer and scp the tar.gz file to your HPCC account and unzip it
##### I added the following lines to my ~/.bashrc file:
module load R/3.1.0
source ~/venv/base/bin/activate

##### Instal the R packages (after loading the appropriate version of R
##### R Code (type within R)
install.packages("ggplot2")
install.packages("plyr")
install.packages("reshape2")

#### To run the notebook interactively:
##### Run the ipython notebook without a browser:
ipython notebook --no-browser &

##### Run firefox in background:
./firefox/firefox &

##### In the address bar type: http://127.0.0.1:8888/tree
##### This will open the IPython notebook
##### I was able to run the notebook without any problem using this setup
##### If you don’t want to have to do this workaround everytime, you can create an IPython profile and set the browser appropriately.
##### Don’t use the firefox originally available at HPCC because it crashes very easily.
##### Any questions, feel free to contact me at teo.decarvalho@gmail.com

#### To run the notebook using qsub 
##### Install runipy
pip install runipy

##### See the file job_example.qsub for an example of running the Nitrosopumilus example.
##### You will have to comment out part of the notebook or use "if" statements to select
##### the section of the notebook that you want to run.
