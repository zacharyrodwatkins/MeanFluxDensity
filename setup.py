from distutils.core import setup, Extension

setup(name = "FluxDen", maintainer = "Zachary Watkins" ,
maintainer_email = "watkinsz@nrc-cnrc.gc.ca", ext_modules=[Extension(
'FluxDen', sources = ['mean_flux_den.c'])])
