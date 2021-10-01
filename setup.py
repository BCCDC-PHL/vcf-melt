from setuptools import setup, find_packages

from vcf_melt import __version__

setup(
    name='vcf-melt',
    python_requires='>=3, <3.6',
    version=__version__,
    packages=find_packages(),
    scripts=[
    ],
    package_data={
    },
    install_requires=[
        'pyvcf==0.6.*',
    ],
    description='Melt vcf files to resolve one-to-many relationships',
    url='https://github.com/BCCDC-PHL/vcf-melt',
    author='Dan Fornika',
    author_email='dan.fornika@bccdc.ca',
    entry_points="""
    [console_scripts]
    vcf-melt = vcf_melt.vcf_melt:main
    """,
    include_package_data=True,
    keywords=[],
    zip_safe=False,
)
