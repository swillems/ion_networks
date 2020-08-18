from setuptools import setup
import src._version


with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()


setup(
    name=ion_networks._version.__project__,
    version=ion_networks._version.__version__,
    license=ion_networks._version.__license__,
    description="Analysis of LC-[...]-MSMS data with ion-networks.",
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
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=["src", "lib"],
    include_package_data=True,
    package_data={'lib': ['*']},
    entry_points={
        "console_scripts": [
            "ion_networks=src.ion_networks:main",
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
        "ipykernel",
        "pyqt5-sip",
        "pyqtwebengine",
        "tables",
        "ms2pip",
    ],
    python_requires=">=3.6,<4",
)
