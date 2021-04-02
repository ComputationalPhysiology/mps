#!/usr/bin/env python

"""The setup script."""

from setuptools import find_packages, setup

with open("README.md") as readme_file:
    readme = readme_file.read()

requirements = ["numpy", "scipy", "imageio", "tifffile", "typer"]

extras_require = {
    "zip": ["xmltodict"],
    "mp4": ["imageio-ffmpeg"],
    "plot": ["matplotlib"],
}
extras_require.update(
    {"all": [val for values in extras_require.values() for val in values]}
)


setup(
    author="Henrik Finsberg",
    author_email="henriknf@simula.no",
    python_requires=">=3.5",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    description="Tools for working with mps files",
    entry_points={"console_scripts": ["mps=mps.cli:main"]},
    install_requires=requirements,
    long_description=readme,
    include_package_data=True,
    keywords="mps",
    name="mps",
    extras_require=extras_require,
    packages=find_packages(include=["mps", "mps.*"]),
    url="https://github.com/ComputationalPhysiology/mps",
    version="0.1.0",
    zip_safe=False,
)
