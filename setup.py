from setuptools import setup

setup(name='hsz',
        version='0.0.1.dev',
        description='Calculate the Stark and Zeeman maps in Rydberg helium using the Numerov method.',
        url='',
        author='Alex Morgan',
        author_email='alexandre.morgan.15@ucl.ac.uk',
        license='GPL-3.0',
        packages=['hsz'],
        install_requires=[
            'tqdm',
            'attrs'
        ],
        include_package_data=True,
        zip_safe=False)
