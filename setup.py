from setuptools import setup, find_packages

setup(
    name='biofkit',
    version='0.0.1',
    description=('With biofkit, sequence can be easily extract. And more functions will be added.'),
    long_description=open('README.md').read(),
    author='Zhang Yujian',
    author_email='Chouuken@outlook.com',
    maintainer='Zhang Yujian',
    maintainer_email='Chouuken@outlook.com',
    license='Apache License 2.0',
    url='https://github.com/Chou-Uken/biofkit',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Operating System :: OS Independent',
        'License :: OSI Approved :: Apache License 2.0',
    ],
    python_requires='>=3.9',
    install_requires=[
        'os'
    ]
)
