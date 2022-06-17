import setuptools
import codecs
import os.path
from pybind11 import get_cmake_dir
# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_string(string, rel_path="actipy/__init__.py"):
    for line in read(rel_path).splitlines():
        if line.startswith(string):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError(f"Unable to find {string}.")

with open("README.md", "r") as fh:
    long_description = fh.read()

ext_modules = [
    Pybind11Extension("GENEActivReaderCPP",
        ["actipy/GENEActivReader.cxx"],
        # # Example: passing in the version to the compiled code
        # define_macros = [('VERSION_INFO', __version__)],
        ),
]

setuptools.setup(
    name="actipy",
    python_requires=">=3.7",
    version=get_string("__version__"),
    description="Python package to process wearable accelerometer data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/activityMonitoring/actipy",
    author=get_string("__author__"),
    author_email=get_string("__email__"),
    license=get_string("__license__"),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
    ],
    packages=setuptools.find_packages(exclude=("test",)),
    include_package_data=True,
    install_requires=[
        "numpy",
        "scipy",
        "pandas>=1.2.5",
        "statsmodels>=0.12.2",
        "Jpype1==1.3.0",
        "pybind11>=2.9.2"
    ],
    ext_modules=ext_modules,
    # Currently, build_ext only provides an optional "highest supported C++
    # level" feature, but in the future it may provide more features.
    cmdclass={"build_ext": build_ext},
)
