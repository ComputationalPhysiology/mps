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

BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

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
	python -m flake8 mps tests

type: ## Run mypy
	python -m mypy mps tests

black: ## Run mypy
	python -m black mps tests

test: ## run tests quickly with the default Python
	python -m pytest

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
	python setup.py sdist
	python setup.py bdist_wheel
	ls -l dist

install: clean ## install on unix
	python -m pip install "."

install-windows: clean ## install on windows usig pipwin
	python -m pip install --upgrade pip
	python -m pip install pipwin
	pipwin install -r requirements.txt
	python -m pip install "."

dev: clean ## Developement install
	python -m pip install --upgrade pip
	python -m pip install git+https://github.com/ComputationalPhysiology/ap_features.git@master
	python -m pip install -e ".[dev]"
	pre-commit install

dev-windows: clean ## Developement install - windows
	python -m pip install --upgrade pip
	python -m pip install pipwin
	python -m pip install ".[dev]"
	pre-commit install

installer: clean  ## make installer for unix
	python -m pip install -r requirements.txt
	python -m pip install pyinstaller
	pyinstaller -F mps/__main__.py -n mps --hidden-import=imageio_ffmpeg --hidden-import=matplotlib --hidden-import=scipy.special.cython_special --additional-hooks-dir=pyinstaller_hooks

installer-windows: clean  ## make installer for windows
	python -m pip install pipwin
	pipwin install -r requirements.txt
	pipwin install pyinstaller
	pyinstaller -F mps/__main__.py -n mps --hidden-import=imageio_ffmpeg --hidden-import=matplotlib --hidden-import=scipy.special.cython_special --additional-hooks-dir=pyinstaller_hooks
