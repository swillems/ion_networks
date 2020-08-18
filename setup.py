from setuptools import setup
import ion_networks._version

VERSION = ion_networks._version.__version__

NAME = ion_networks._version.__project__
LICENSE = ion_networks._version.__license__
DESCRIPTION = "Analysis of LC-[...]-MSMS data with ion-networks."
AUTHOR = ion_networks._version.__author__
AUTHOR_EMAIL = "sander.willems@ugent.be"
URL = "https://github.com/swillems/ion_networks"
PROJECT_URLS = {
    "Source": "",
    "Publication": "https://www.biorxiv.org/content/10.1101/726273v2",
}
KEYWORDS = [
    "MS",
    "Mass spectrometry",
    "Proteomics",
    "DIA",
    "Data-independent acquisition",
]
CLASSIFIERS = [
    "Intended Audience :: Science/Research",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python :: 3 :: Only",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]
INSTALL_REQUIRES = [
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
    # "ms2pip",
]
PYTHON_REQUIRES = ">=3.6,<4"

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

#_ = [[os.remove(f) for f in glob(pat)] for pat in to_remove]

setup(
    name=NAME,
    version=VERSION,
    license=LICENSE,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    project_urls=PROJECT_URLS,
    keywords=KEYWORDS,
    classifiers=CLASSIFIERS,
    packages=["src", "lib"],
    include_package_data=True,
    package_data={'lib': ['*']},
    entry_points={
        "console_scripts": [
            "ion_networks=src.ion_networks:main",
        ],
    },
    install_requires=INSTALL_REQUIRES,
    python_requires=PYTHON_REQUIRES,
    # include_dirs=["lib"],
    # cmdclass={"build_ext": build_ext},
)
