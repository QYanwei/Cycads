from setuptools import setup, find_packages
from cycads import __version__, __description__, __url__

with open('README.md') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    required = [x for x in f.read().splitlines() if x]

setup(
    name='cycads',
    version=__version__,
    packages=find_packages(),
    author="Yanwei Qi, Jia-Yuan Zhang",
    author_email="qiyanwei1@genomics.cn, zhangjiayuan@genomics.cn",
    description=__description__,
    long_description=long_description,
    url=__url__,
    entry_points={
        'console_scripts': [
            'cycads = cycads.cycads:main',
        ],
    },
    package_data={'cycads':  ['resources/*']},
    classifiers=[
        # Trove classifiers
        # Full list at https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
    ],
    install_requires=required,
)
