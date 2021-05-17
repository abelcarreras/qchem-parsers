from setuptools import setup
import os


def get_version_number():
    main_ns = {}
    for line in open('qcparsers/__init__.py', 'r').readlines():
        if not(line.find('__version__')):
            exec(line, main_ns)
            return main_ns['__version__']


def package_dirs(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        print(directories)
        if '__pycache__' in directories:
            paths.append(path[2:].replace('/', '.'))
    return paths


# Make python package
setup(name='qcparsers',
      version=get_version_number(),
      description='Python parsers for Q-Chem',
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      install_requires=['numpy'],
      author='Abel Carreras',
      author_email='abelcarreras83@gmail.com',
      packages=['qcparsers', 'qcparsers.tools', 'qcparsers.abstractions'] + package_dirs('./qcparsers/parsers'),
      url='https://github.com/abelcarreras/PyQchem',
      classifiers=[
          "Programming Language :: Python",
          "License :: OSI Approved :: MIT License"]
      )
