## 기풍

```
conda create -n "py38" python=3.8
conda activate py38 

pip install --upgrade pip
python -m pip install --upgrade pip

pip install jupyterlab
pip install --upgrade jupyterlab

pip install jupyter_contrib_nbextensions
jupyter contrib nbextension install --user

jupyter labextension install @jupyterlab/toc
jupyter labextension install @jupyterlab/shortcutui
jupyter labextension install @lckr/jupyterlab_variableinspector
jupyter labextension install @krassowski/jupyterlab_go_to_definition
jupyter labextension install @jupyter-widgets/jupyterlab-manager
jupyter labextension install lineup_widget
jupyter labextension install @jupyterlab/github

pip install nbresuse
jupyter labextension install jupyterlab-topbar-extension jupyterlab-system-monitor
jupyter labextension install jupyterlab_filetree
conda install --yes jupyter-archive

jupyter lab build
jupyter labextension update --all



```



## 호영

### Update of PIP
```
pip install --upgrade pip
python -m pip install --upgrade pip
```


### Jupyter Nbextensions
```
pip install jupyter_contrib_nbextensions
jupyter contrib nbextension install --user
```

### Jupyter Lab
```
pip install jupyterlab
pip install --upgrade jupyterlab
```

### Jupyter Lab Extensions Package
```
pip install nodejs
conda install --yes nodejs
conda install -c conda-forge --yes nodejs
```

### Table of Contents
```
jupyter labextension install @jupyterlab/toc
```

### Shortcut UI
```
jupyter labextension install @jupyterlab/shortcutui
```

### Variable Inspector
```
jupyter labextension install @lckr/jupyterlab_variableinspector
```

### Go to Definition of Module
```
jupyter labextension install @krassowski/jupyterlab_go_to_definition
```

### Interactive Visualization
```
jupyter labextension install @jupyter-widgets/jupyterlab-manager
jupyter labextension install lineup_widget
```

### Connection to Github
```
jupyter labextension install @jupyterlab/github
```

### CPU+RAM Monitor
```
pip install nbresuse
jupyter labextension install jupyterlab-topbar-extension jupyterlab-system-monitor
```

### File Tree Viewer
```
jupyter labextension install jupyterlab_filetree
```

### Download Folder as Zip File
```
conda install --yes jupyter-archive
jupyter lab build
jupyter labextension update --all
```

