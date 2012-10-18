from distutils.core import setup

setup(
    name='BamTyper',
    version='0.1.0',
    author='Michael Imelfort',
    author_email='mike@mikeimelfort.com',
    packages=['bamtyper', 'bamtyper.test'],
    scripts=['bin/bamtyper'],
    url='http://pypi.python.org/pypi/BamTyper/',
    license='LICENSE.txt',
    description='Working with paired reads in BAM format',
    long_description=open('README.txt').read(),
    install_requires=[
        "pysam >= 0.6",
    ],
)
