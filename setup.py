from setuptools import setup

# Version
# version = None
# with open("./pykofamsearch/__init__.py", "r") as f:
#     for line in f.readlines():
#         line = line.strip()
#         if line.startswith("__version__"):
#             version = line.split("=")[-1].strip().strip('"')
# assert version is not None, "Check version in pykofamsearch/__init__.py"

exec(open('pyhmmsearch/__init__.py').read())

setup(
    name='pyhmmsearch',
    version=__version__,
    description='Fast implementation of HMMSEARCH optimized for high-memory systems using PyHmmer',
    url='https://github.com/jolespin/pyhmmsearch',
    author='Josh L. Espinoza',
    author_email='jolespin@newatlantis.io, jol.espinoz@gmail.com',
    license='MIT License',
    packages=["pyhmmsearch"],
    install_requires=[
        "pyhmmer >=0.10.12",
        "pandas",
        "tqdm",
        "requests",  # Assuming API interaction may still be useful
    ],
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

