# blastocoel_review
Jupyter notebooks used for the review "Blastocoel morphogenesis: a biophysics perspective" by M. Le-Verge-Serandour and H. Turlier

## Requirements
This project is developped using Python3 and requires the following python libraries to run :
- numpy (np)
- scipy.integrate
- matplotlib.pyplot (plt)
- matplotlib.patches
- ipywidgets

## Content
- figures/                        : folder containing figures included in the notebooks.
- blastula_class.py               : main class used to plot the blastula with individual cells. 
- blastomeres_packing.ipynb       : notebook plotting a blastula, calculate the aspect ratio a single cells in an embryo (Fig. 4C).
- hemispherical_blastomeres.ipynb : notebook calculating the aspect ratio of a single cell spreading on a surface (Fig. 4B).
- pumping.ipynb                   : notebook that integrates the volume-concentration-pump model of Fig. 6A.

## About
Developped by M. Le Verge--Serandour and H. Turlier, 2021
