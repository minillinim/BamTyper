from distutils.core import setup

setup(
    name='BamTyper',
    version='0.2.6',
    author='Michael Imelfort',
    author_email='mike@mikeimelfort.com',
    packages=['bamtyper', 'bamtyper.test'],
    scripts=['bin/bamtyper'],
    url='http://pypi.python.org/pypi/BamTyper/',
    license='GPLv3',
    description='Working with paired reads in BAM format',
    long_description=open('README.md').read(),
    install_requires=[
        "pysam >= 0.6",
    ],
)
