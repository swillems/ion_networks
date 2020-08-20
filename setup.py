import setuptools
import ion_networks._version


with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()


setuptools.setup(
    name=ion_networks._version.__project__,
    version=ion_networks._version.__version__,
    license=ion_networks._version.__license__,
    description="Analysis of LC-[...]-MSMS data with ion-networks",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    author=ion_networks._version.__author__,
    author_email="sander.willems@ugent.be",
    url="https://github.com/swillems/ion_networks",
    project_urls={
        "Publication": "https://www.biorxiv.org/content/10.1101/726273v2",
    },
    keywords=[
        "MS",
        "Mass spectrometry",
        "Proteomics",
        "DIA",
        "Data-independent acquisition",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=setuptools.find_packages(),
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "ion_networks=ion_networks.ion_networks:main",
        ],
    },
    install_requires=[
        "numpy",
        "numexpr",
        "matplotlib",
        "PySimpleGUI",
        "pandas",
        "click",
        "scipy",
        "scikit-learn",
        "h5py",
        "pyteomics",
        "numba",
        "pyqt5-sip",
        "pyqtwebengine",
    ],
    extras_require={
        # MS2PIP only installable on Linux,  manual installation needed for OSX
        'ms2pip': [
            'ms2pip'
        ]
    },
    python_requires=">=3.6,<4",
)
