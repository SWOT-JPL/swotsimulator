================================
SWOT Simulator for Ocean Science
================================

Description
-----------
This software simulates SWOT observations of sea surface height (SSH) that can be applied to an ocean general circulation model (OGCM), allowing the exploration of ideas and methods to optimize information retrieval from the SWOT Mission in the future. From OGCM SSH inputs, the software generates SWOT-like outputs on a swath along the orbit ground track and adds measurement error and noise, which are generated according to technical characteristics published by the SWOT project team. Not designed to directly simulate the payload instrument performance, this SWOT simulator aims at providing statistically realistic outputs for the science community with a simple software package released as an open source in Python. The software is scalable and designed to support future evolution of orbital parameters, error budget estimates from the project team and suggestions from the science community.

Data in the example directory
-------------------------------
WARNING: we had to remove the example data from the input_field directory as we can't (and don't want to) git data that are too large. This data are available under this link: http://swot.jpl.nasa.gov/science/resources/ 
Before running the example:
  - click on swotsimulator_largedata.zip to download the data,
  - unzip the downloaded zipped document
  - copy everything from the unzipped directory into example/input_field directory
  
Note that you only have to do that once when you clone the git swotsimulator repository. Updating the code by doing a 'git pull' won't modify the example data files. 

Manual
------------
Refer to [Full manual](https://github.com/SWOTsimulator/swotsimulator/blob/master/doc/source/science.rst) for more details.

Updates
------------

* August 22, 2016
 * updated Karin noise
