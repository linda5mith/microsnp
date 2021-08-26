import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="example-pkg-linda5mith",
    version="1.0.0",
    author="Linda Smith",
    author_email="linda.smith@ucc.ie",
    description="A package for finding SNPs in microbial genomic data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/linda5mith/microsnp",
    project_urls={
        "Bug Tracker": "https://github.com/linda5mith/microsnp/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src",exclude='data'),
    entry_points={
        'console_scripts': [
            'microsnp = microsnp.microsnp:main',
        ]
    },
    python_requires=">=3.6",
)

