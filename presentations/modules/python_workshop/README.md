## Pre workshop instructions

1. Install docker: https://www.docker.com/products/docker-desktop

2. Pull the workshop image
```
docker pull mpgagebioinformatics/python workshop
```

3. Run the image and enter the container:
```
mdkir -p ~/python_workshop
docker run -p 8787:8787 -p 8888:8888 -v ~/python_workshop/:/home/mpiage --name python_workshop -it mpgagebioinformatics/python_workshop:latest
```

4. Start jupyterhub inside the container
```
module load jupyterhub
jupyter notebook --ip=0.0.0.0
```

5. Access jupyter.
You will be shown a message like this:
```
   Copy/paste this URL into your browser when you connect for the first time,
    to login with a token:
        http://b4b294d48b2b:8888/?token=01d1532b8624becc5c7b288cf0c600b7805e802349e81478&token=01d1532b8624becc5c7b288cf0c600b7805e802349e81478
```
Replace the "b4b294d48b2b" string (ie. bettween the http:// and :8888/) with "localhost" (eg. http://localhost:8888/?token=01d1532b8624becc5c7b288cf0c600b7805e802349e81478&token=01d1532b8624becc5c7b288cf0c600b7805e802349e81478 ) and paste it into your web browser.


If you are not able to install docker you can still run the notebook by install Python3 - https://www.python.org - and jupyter - https://jupyter.org - asafterwards make sure you have installed all required packages:
```
Package            Version
------------------ -------
autograd           1.2    
backcall           0.1.0  
bleach             2.1.3  
Bottleneck         1.2.1  
cycler             0.10.0 
decorator          4.3.0  
entrypoints        0.2.3  
future             0.17.1 
html5lib           1.0.1  
ipykernel          4.8.2  
ipython            6.4.0  
ipython-genutils   0.2.0  
ipywidgets         7.2.1  
jedi               0.12.0 
Jinja2             2.10   
jsonschema         2.6.0  
jupyter            1.0.0  
jupyter-client     5.2.3  
jupyter-console    5.2.0  
jupyter-core       4.4.0  
kiwisolver         1.0.1  
lifelines          0.19.4 
MarkupSafe         1.0    
matplotlib         3.0.2  
matplotlib-venn    0.11.5 
mistune            0.8.3  
nbconvert          5.3.1  
nbformat           4.4.0  
notebook           5.5.0  
numpy              1.16.1 
pandas             0.24.1 
pandocfilters      1.4.2  
parso              0.2.1  
pexpect            4.6.0  
pickleshare        0.7.4  
pip                9.0.3  
prompt-toolkit     1.0.15 
ptyprocess         0.6.0  
Pygments           2.2.0  
pyparsing          2.3.1  
python-dateutil    2.7.3  
pytz               2018.9 
pyzmq              17.0.0 
qtconsole          4.3.1  
scikit-learn       0.20.2 
scipy              1.2.1  
seaborn            0.9.0  
Send2Trash         1.5.0  
setuptools         39.0.1 
simplegeneric      0.8.1  
six                1.11.0 
sklearn            0.0    
terminado          0.8.1  
testpath           0.3.1  
tornado            5.0.2  
traitlets          4.3.2  
wcwidth            0.1.7  
webencodings       0.5.1  
widgetsnbextension 3.2.1  
xlrd               1.2.0 
```

Check your packages with `pip3 list` and install with `pip3 install <package_name>==<version> --user`.

You can then download the notebook from: https://github.com/mpg-age-bioinformatics/presentations-tutorials/blob/gh-pages/presentations/modules/python_workshop/python_workshop.ipynb.

&nbsp;


---