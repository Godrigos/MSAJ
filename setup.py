import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='MSAJ',
    version='1.0.3',
    url='https://github.com/Godrigos/MSAJ',
    license='GPL-3.0',
    author='Rodrigo Aluizio',
    author_email='',
    description='Multi Sequence Alignment Joiner',
    long_description=long_description,
    entry_points={'gui_scripts': ["MSAJ = ui.py"]},
    package_data={
        '': ['icons/*.*'],
    },
    install_requires=['PyQt5', 'biopython'],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL-3.0",
        "Operating System :: OS Independent",
    ],
)

