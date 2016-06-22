"""To use this file and amf.c do the following:

>python setup.py build
>ln -sf build/lib.linux-x86_64-2.6/amf.so amf.so

For more information, see:

http://docs.python.org/extending/building.html#building
from distutils.core import setup, Extension
http://docs.scipy.org/doc/numpy/reference/distutils.html
"""
import numpy
import numpy.numarray
import numpy.distutils
import numpy.distutils.core
import numpy.distutils.misc_util
from numpy.distutils.core import setup, Extension

module1 = Extension('amf',
                    sources = ['amf.c'])



def try_config():
    config = numpy.distutils.misc_util.Configuration(None, None, None)
    ## I should be able to use thsi get_numarray_include_dirs(), but
    ## it does not work
    config.add_include_dirs(numpy.numarray.get_numarray_include_dirs())
    #config.add_include_dirs('/usr/lib/python2.5/site-packages/numpy/numarray')
    return config

setup (name = 'PackageName',
       version = '1.0',
       description = 'Speed up for make_Radon',
       ext_modules = [module1],
       configuration=try_config)
