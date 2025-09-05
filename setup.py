from setuptools import setup

# Version
exec(open('pyhmmsearch/__init__.py').read())

# Read requirements from requirements.txt
with open("requirements.txt", "r") as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    
setup(
    name='pyhmmsearch',
    version=__version__,
    description='Fast implementation of HMMSEARCH optimized for high-memory systems using PyHmmer',
    url='https://github.com/jolespin/pyhmmsearch',
    author='Josh L. Espinoza',
    author_email='jol.espinoz@gmail.com',
    license='MIT License',
    packages=["pyhmmsearch"],
    install_requires=requirements,
    include_package_data=False,
    entry_points={
        'console_scripts': [
            'pyhmmsearch=pyhmmsearch.pyhmmsearch:main',   # Executes pyhmmsearch.main()
            'reformat_pyhmmsearch=pyhmmsearch.reformat_pyhmmsearch:main',  # Executes reformat_pyhmmsearch.main()
            'serialize_hmm_models=pyhmmsearch.serialize_hmm_models:main',  # Executes serialize_hmm_models.main()
        ],
    },
#      scripts=[
#          "pyhmmsearch/pyhmmsearch.py",
#          "pyhmmsearch/reformat_pyhmmsearch.py",
#          "pyhmmsearch/serialize_hmm_models.py",
#          ],
)

