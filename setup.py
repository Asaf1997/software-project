from setuptools import setup, find_packages, Extension
setup(
    name='mykmeanssp',
    version='0.1.0',
    packages=find_packages(),
    author="Asaf & Ido Productions",
    ext_modules=[
        Extension(
            'mykmeanssp',
            ['spkmeansmodule.c', 'spkmeans.c']
        ),
    ]
)