# landau-2022
Code accompanying our manuscript published in eLife 2022.

## Overview
This repository contains code for running neuron simulations and generating figures as presented in our eLife 2022 paper. The codebase provides tools for neural simulation and analysis, with a focus on [brief description of main research focus].

## Installation
To get started with this codebase:

1. Clone the repository:
   ```bash
   git clone https://github.com/landoskape/landau-2022.git
   cd landau-2022
   ```

2. Install dependencies in a new conda environment:
   ```bash
   conda env create -f environment.yml
   ```

3. Activate the environment:
   ```bash
   conda activate landau-2022
   ```

> Note: NEURON is installed via pip as instructed in the `environment.yml` file. This is a bit more finnicky than typical python packages. You might need to install NEURON manually and also might need to uninstall / reinstall to get it working. Note that some of the NEURON code requires version 8 or greater (maybe even 8.2 as specified in the `environment.yml` file). I haven't tested to see which versions it works with, but it definitely doesn't work with any version <8.

## How to use the code
The code base is built around a few components:
1. [collection_L23](./src/collection_L23.py) / [collection_uncageMapping](./src/collection_uncageMapping.py)
    - These files contain the definitions of the neuron models we used in our paper. The morphologies were identified through two-photon imaging of the Alex 594 fluorophore and documented in hoc files in the `src/morphologies_from_paper` folder.
    - The `collection_L23.py` file contains the definitions of additional neuron models from the literature (I used them to learn and develop NEURON in the beginning and kept them around). 
    - Both collections depend on hoc files in the `src/morphologies_from_*` folders. 
2. [morphologyFunctions](./src/morphologyFunctions.py)
    - This file contains some functions for analyzing neural morphologies in NEURON and a few other functions for manipulating activity in a site specific way. 
3. [neuronFunctions](./src/neuronFunctions.py)
    - This file contains a few very basic functions for simulating with NEURON and handling the hoc stuff in python. 
4. [NEURON Base Components](./src/mod.files/)
    - This folder contains the NEURON mechanisms that we used in our paper (and some other ones!). They are the base code that define channels and synapses.

### Notebooks that show how the code works

#### [KeyNeuronPlots](./KeyNeuronPlots.ipynb)
This notebook serves as a tutorial and demonstration of the core functionality. It includes:
- Examples of how to run neuron simulations using our codebase
- Visualization of neural activity and key parameters
- Sample plots demonstrating various analysis techniques
- These plots draw from the neuron models in the `src/collection_L23.py` file which are based on Neuron morphologies we found in the literature. The `src/collection_L23.py` file has information about where the cells came from in the comments of the ``returnSecListL23`` method. The citations should have been better but the code was written a long time ago and I've lost track. I can find it again if you ask (but they're all from modeldb and should be somewhat easy to find). 

#### [Optimization Notebooks](./optimizeKaDensity.ipynb)
The notebooks that start with the name "optimizeKaDensity" are the ones we used to perform analyses in our paper. They are a bit slow because they all start with an optimization procedure to select potassium densities. 
- These notebooks build on `src/collection_uncageMapping.py` which contains the definitions of the neuron models we used in our paper. The morphologies were identified through two-photon imaging of the Alex 594 fluorophore and documented in hoc files in the `src/morphologies_from_paper` folder.


### FYI
> Full disclosure: this code is old and was written in my early days of python programming so it could definitely be organized better. Apologies if it's opaque or hard to follow. Hopefully the notebooks show enough to get you started and moving from there. Otherwise send me an email and I can try to help out.

## Citation
If you use this code in your research, please cite our paper:
Andrew T Landau, Pojeong Park, J David Wong-Campos, He Tian, Adam E Cohen, Bernardo L Sabatini (2022) Dendritic branch structure compartmentalizes voltage-dependent calcium influx in cortical layer 2/3 pyramidal cells eLife 11:e76993

## License
This code is licensed under the GNU General Public License v3.0 (GPL-3.0). The GPL is a copyleft license that ensures the code remains free and open source. See the [LICENSE](./LICENSE) file for full details.

## Contact
If this is useful, not useful enough and you want it to be, or you have any other questions, please contact me at andrewtylerlandau at gmail dot com. I'm excited you're interested and trying out our code!

