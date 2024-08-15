
from setuptools import setup, find_namespace_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="cgen2gmx",  
    version="1.0",  # Required
    long_description=long_description, 
    long_description_content_type="text/markdown", 
    url="https://github.com/chrispy67/cgen2gmx",
    author="Christian Phillips", 
    author_email="christian_phillips1@msn.com",
    classifiers=[
        "Development Status :: 1 - Planning",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    keywords="molecular dynamics, charmm, forcefield, computational chemistry",
    packages=find_namespace_packages(),

    python_requires="<4",
    extras_require={ 
        "test": ["coverage", 'pytest'],
    },

    package_data={},
    data_files=None, 
    entry_points={},
)