from setuptools import setup, find_packages

## install main application
desc = "utility script for queue MAGs against the MGnify db"
setup(
    name="mgnify_search",
    version='1.0.0',
    description=desc,
    long_description=desc + "\n See README for more information.",
    author="Santino Faack",
    author_email="santino_faack@gmx.de",
    license="GPL-3.0 license",
    packages=find_packages(),
    url="https://github.com/SantaMcCloud/MGnify-genome-search",
    entry_points={
        "console_scripts": [
            "mgnify_search = mgnify_search.cli:main",
        ]
    },
)
