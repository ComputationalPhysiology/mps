.PHONY: clean clean-test clean-pyc clean-build docs help
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python3 -c "$$BROWSER_PYSCRIPT"

help:
	@python3 -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .pytest_cache

lint: ## check style with flake8
	python3 -m flake8 mps tests

type: ## Run mypy
	python3 -m mypy mps tests

black: ## Run mypy
	python3 -m black mps tests

test: ## run tests quickly with the default Python
	python3 -m pytest

docs: ## generate Sphinx HTML documentation, including API docs
	rm -f docs/source/mps.rst
	rm -f docs/source/modules.rst
	sphinx-apidoc -o docs/source mps/ mps/tifffile.py mps/czifile.py mps/nd2file.py
	cp CONTRIBUTING.md docs/source/.
	$(MAKE) -C docs clean
	$(MAKE) -C docs html


show: ## show docs
	open docs/build/html/index.html


dist: clean ## builds source and wheel package
	python3 setup.py sdist
	python3 setup.py bdist_wheel
	ls -l dist

install: clean upgrade-pip ## install on unix
	python3 -m pip install "."

install-windows: clean ## install on windows usig pipwin
	python3 -m pip install --upgrade pip
	python3 -m pip install pipwin
	pipwin install -r requirements.txt
	python3 -m pip install "."

dev: clean upgrade-pip ## Developement install
	python3 -m pip install git+https://github.com/ComputationalPhysiology/ap_features.git@master
	python3 -m pip install -e ".[dev]"
	pre-commit install

upgrade-pip:
	python3 -m pip install pip --upgrade

dev-windows: clean upgrade-pip ## Developement install - windows
	python3 -m pip install pipwin
	python3 -m pip install ".[dev]"
	pre-commit install

installer: clean upgrade-pip ## make installer for unix
	python3 -m pip install ".[dev,motion]"
	python3 -m pip install pyinstaller
	python -m pip install "opencv-contrib-python==4.5.3.56"
	python -m pip install "opencv-python==4.5.3.56"
	pyinstaller -F mps/__main__.py -n mps --hidden-import=imageio_ffmpeg --hidden-import=matplotlib --hidden-import=scipy.special.cython_special --collect-submodules imageio --additional-hooks-dir=pyinstaller_hooks

installer-windows: clean upgrade-pip ## make installer for windows
	python3 -m pip install pipwin
	pipwin -m pip install ".[dev,motion]"
	pipwin install pyinstaller
	pyinstaller -F mps/__main__.py -n mps --hidden-import=imageio_ffmpeg --hidden-import=matplotlib --hidden-import=scipy.special.cython_special --collect-submodules imageio --additional-hooks-dir=pyinstaller_hooks


release: dist ## package and upload a release
	python3 -m twine upload dist/*
