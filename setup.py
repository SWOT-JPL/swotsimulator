"""Build and install the SWOT Simulator for Ocean Science package."""
#from distutils.core import setup
from setuptools import setup, find_packages
import os

#- Read in the package version and author fields from the Python
#  source code package_version.py file:

with open(os.curdir+os.sep+'swotsimulator'+os.sep+'package_version.py') as f:
  code=compile(f.read(), os.curdir+os.sep+'swotsimulator'+os.sep+'package_version.py', 'exec')
  exec(code)
#execfile(os.curdir+os.sep+'swotsimulator'+os.sep+'package_version.py')

setup(name='swotsimulator',
          version=version,
          description=description,
	  author=author,
	  author_email=email,
	  url=url,
	  license = 'COPYING',
	  long_description=open('README').read(),
      packages=find_packages(),
	  install_requires=[
	  	"Numpy",
		"Scipy",
		],
      entry_points={                                                            
          'console_scripts': [                                                  
               'swotsim = swotsimulator.run:main',                        
               ]
      }
     )

