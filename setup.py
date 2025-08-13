from setuptools import setup, find_packages
import pathlib

# Read the README file
here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="fasttrimmatic",
    version="1.0.0",  # Match your GitHub release tag
    author="Beckley Brown",
    author_email="brownbeckley94@gmail.com",
    description="Hybrid fastp + Trimmomatic-like FASTQ quality control and trimming tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bbeckley-hub/fasttrimmatic",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "biopython",
        "numpy",
        "tqdm",
    ],
    extras_require={
        "ml": ["scikit-learn"]
    },
    entry_points={
        "console_scripts": [
            "fasttrimmatic=fasttrimmatic.fasttrimmatic:main_cli",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # Pick the license you chose
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    project_urls={
        "Bug Tracker": "https://github.com/bbeckley-hub/fasttrimmatic/issues",
        "Documentation": "https://github.com/bbeckley-hub/fasttrimmatic#readme",
    },
)
