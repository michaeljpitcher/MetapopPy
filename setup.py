from setuptools import setup, find_packages

setup(
   name='MetapopPy',
   version='1.0',
   description='Metapopulation network framework',
   author='Michael Pitcher',
   author_email='mjp22@st-andrews.ac.uk',
   packages=find_packages(),
   install_requires=['epyc','matplotlib','networkx'], #external packages as dependencies
)