
from setuptools import setup, find_namespace_packages

# Load the long description from the README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="cgen2gmx",  
    version="1.1.0", # Needs to be updated EACH TIME A NEW RELEASE IS CREATED
    description="A small commandline tool for managing forcefield parameters used in molecular dynamics simulations", 
    long_description=long_description,  # Detailed description from README
    long_description_content_type="text/markdown",  # Type of content in long description
    url="https://github.com/chrispy67/cgen2gmx",  
    author="Christian Phillips ", 
    author_email="christian_phillips1@msn.com", 
    classifiers=[
        "Development Status :: 1 - Planning",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    keywords="molecular dynamics, charmm, forcefield, computational chemistry",  # Relevant keywords
    packages=find_namespace_packages(),  # Finds all packages under the namespace
    python_requires=">=3.7",  
    install_requires=[
        # List your package dependencies here, e.g.,
        # 'pandas', 'os', 'click'
    ],
    py_modules=['cgen2gmx'],
    entry_points={
        'console_scripts': [
            'cgen2gmx=cgen2gmx:main', 
        ],
    },
    include_package_data=True, 
    project_urls={ 
        "Bug Reports": "https://github.com/chrispy67/cgen2gmx/issues",
        "Source": "https://github.com/chrispy67/cgen2gmx",
    },
)





# setup(
#     name="cgen2gmx",  
#     version="1.0",  # Required
#     long_description=long_description, 
#     long_description_content_type="text/markdown", 
#     url="https://github.com/chrispy67/cgen2gmx",
#     author="Christian Phillips", 
#     author_email="christian_phillips1@msn.com",
#     classifiers=[
#         "Development Status :: 1 - Planning",
#         "License :: OSI Approved :: MIT License",
#         "Programming Language :: Python :: 3",
#     ],
#     keywords="molecular dynamics, charmm, forcefield, computational chemistry",
#     packages=find_namespace_packages(),

#     python_requires="<4",
#     extras_require={ 
#         "test": ["coverage", 'pytest'],
#     },

#     package_data={},
#     data_files=None, 
#     entry_points={},
# )